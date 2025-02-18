#include "dawn_2_obj.h"
#include <numeric>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

size_t model::add_unique_vertex(std::vector<std::array<double, 3>> &v, const std::array<double, 3> &p) {
	/* Adds a new vertex applying the origin offset and the basis transfomation. 
	 * Returns the index of the vertex for use later in the obj file. 
	 * 
	 * We check to make sure the vertex is unique before adding it to prevent the obj file from getting too big. 
	 * This is important when creating tracks as every child track will have shared verticies. 
	 * It also helps blender as it makes all fances part of the same mesh. 
	 */
	std::array<double,3> xp;
	for (size_t i=0; i<3; i++) {
		xp[i]  = std::inner_product(xform[i].begin(), xform[i].end(), p.begin(), static_cast<double>(0.) );
		xp[i] += origin[i];
	}
	
	auto it = std::find_if( v.begin(), v.end()
	                       , [&xp](const std::array<double, 3> &x)->bool{
		                       return xp[0]==x[0] && xp[1]==x[1] && xp[2]==x[2];
	                       }
	);
	
	if (it == v.end()) {
		v.push_back(xp);
		return v.size();
	} else {
		return std::distance(v.begin(), it)+1;
	}
}

void model::add_basis(const std::array<double, 3> &e1, const std::array<double, 3> &e2) {
	/* Defines the current rotation that the next object created will follow.
	 * \begin{verbatim}
	 * /BaseVector 1 0 0 0 1 0 
	 * \end{verbatim}
	 * 
	 * the obj file does not have an equivelent so all the transformations must be applied.
	 * 
	 * DAWN defines it's 'default' basis vectors as $e_1=(1, 0, 0)$, $e_2=(0, 1, 0)$ and $e_3 = e_1 \wedge e_2$. 
	 * any change to this basis defines a new transformed set of coordinates. 
	 * As such we need to transform any object in the new coordinates ($e'_i$) back to the standard ones before writing to the obj file. 
	 * 
	 * To do this we construct a transformation matix $\mathbf{X}$ such that $X e_i = e'_i$
	 * As such a transformation is a special orthogonal transformation, we cn use $X X^{-1} = X X^\intercal = \mathbb{1}$
	 */ 
	std::array<double, 3> e3 = { e1[1]*e2[2]-e1[2]*e2[1], -e1[0]*e2[2]+e1[2]*e2[0], e1[0]*e2[1]-e1[1]*e2[0] }; 	
	for (size_t i=0; i<3; i++){
		xform[i] = { e1[i], e2[i], e3[i] } ;
	}
}

void model::add_box(const std::array<double, 3> &s) {
	/* Creates a box object. 
	 * DAWN definition to create a 1$\times$1$\times$1 box: 
	 * \begin{verbatim}
	 * /Box 0.5 0.5 0.5
	 * \end{verbatim}
	 * 
	 * the obj file on the otherhand does not have a box primative and so we must add 6 square faces, one for each corner. 
	 * a minimal example:
	 * 
	 * \begin{verbatim}
	 * v -0.5 -0.5 -0.5 
	 * v  0.5 -0.5 -0.5
	 * v -0.5  0.5 -0.5
	 * v  0.5  0.5 -0.5
	 * v -0.5 -0.5  0.5
	 * v  0.5 -0.5  0.5
	 * v -0.5  0.5  0.5
	 * v  0.5  0.5  0.5
	 * g box 
	 * f 1 2 3 4
	 * f 1 2 5 6 
	 * f 5 6 7 8
	 * f 3 4 7 8
	 * f 2 4 6 8
	 * f 1 3 5 7
	 * \end{verbatim}
	 * 
	 * As you can see this gets a little complex and verbose.
	 */
	std::stringstream ss;
	ss << "g group_" << group_no++ << std::endl;
	ss << "usemtl colour_" << add_unique_vertex(colour_db, rgb) << std::endl;
	
	auto add_face = [this, &s, &ss](int f_){
		int f = std::abs(f_)-1;
		
		std::array<std::array<int, 4>, 3> A;
		A[f]       = {f_/std::abs(f_), f_/std::abs(f_), f_/std::abs(f_), f_/std::abs(f_)};
		A[(f+1)%3] = {-1, +1, +1, -1};
		A[(f+2)%3] = {-1, -1, +1, +1};
		
		ss << "f ";
		for (int i=0; i<4; i++) {
			ss << add_unique_vertex(vertex_db, {A[0][i]*s[0], A[1][i]*s[1], A[2][i]*s[2]}) << " ";
		}
		ss << std::endl;
	};
	add_face(-3); add_face(+3);
	add_face(-2); add_face(+2);
	add_face(-1); add_face(+1);
	
	obj_db.push_back(ss.str());
}

void model::add_line(const std::vector<std::array<double, 3>> &line) {
	/* Adds a line defined by a set of verticies. 
	 * \begin{verbatim}
	 * /Polyline 
	 * /PLVertex 0 0 0
	 * /PLVertex 0 0 1
	 * ...
	 * /EndPolyline
	 * \end{verbatim}
	 * 
	 * This is Implemented in the obj using the \verb|l| command.
	 */
	std::stringstream ss;
	ss << "g group_" << group_no++ << std::endl;
	ss << "usemtl colour_" << add_unique_vertex(colour_db, rgb) << std::endl;
	ss << "l";
	for (const auto &p : line ) {
		ss << " " << add_unique_vertex(vertex_db, p);
	}
	ss << std::endl;
	obj_db.push_back(ss.str());
}

void model::add_polyhedron(const polyhedron_t &p) {
	/* adds a generic polyhedron to the obj file. 
	 * The structure of the DAWN file is similar to the object file, but with a few differences. 
	 * 
	 * \begin{verbatim}
	 * /Polyhedron
	 * /Vertex -0.5 -0.5 -0.5
	 * /Vertex  0.5  0.5  0.5
	 * ...
	 * /Facet 1 2 3 4
	 * /EndPolyhedron
	 * \end{verbatim}
	 * 
	 * The major diffrences are the vertex nunber is unique to the polyhedron.
	 * Therefore if there are 2 polyhedron in a file they will both have a vertex 1. 
	 * 
	 * The other diffrences are the transformations such as origin offsets and basis vector rotations. 
	 * 
	 */
	std::vector<size_t> m;
	for ( const auto &v : p.v ) {
		m.push_back(add_unique_vertex(vertex_db, v));
	}
	
	std::stringstream ss;
	ss << "g group_" << group_no++ << std::endl;
	ss << "usemtl colour_" << add_unique_vertex(colour_db, rgb) << std::endl;
	for ( const auto &f : p.f ) {
		ss << "f " << m[f[0]-1] << " " << m[f[1]-1] << " " << m[f[2]-1] << std::endl;
	}
	obj_db.push_back(ss.str());
}

void model::write(const std::string &filename) {
	/*
	 * Writes the obj file to an actual file. 
	 * The precision is set to 8dp, this can be changed. 
	 */
	std::ofstream fp_obj(filename + ".obj");
	fp_obj << std::setprecision(8) << std::fixed;
	for (const auto &v : vertex_db) {
		fp_obj << "v " << v[0] << " " << v[1] << " " << v[2] << " " << std::endl;
	}
	for (const auto &s : obj_db ) { 
		fp_obj << s << std::endl;
	}
	fp_obj.close();
	
	std::ofstream fp_mat(filename + ".mtl");
	fp_mat << std::setprecision(8) << std::fixed;
	for ( size_t i=1; const auto &v : colour_db ) { 
		fp_mat << "newmtl color_" << i++ << std::endl;
		fp_mat << "Kd " << v[0] << " " << v[1] << " " << v[2] << " " << std::endl;
	}
	fp_mat.close();
}

template <typename T>
T *safe_cast(void *self) {
	if (self != NULL) return static_cast<T *>(self);
	std::fprintf(stderr, "Unable to open %s. You probs have a syntax error in your prim file\n", typeid(T).name());
	exit(1);
}

extern "C" void *construct_model() {
	return static_cast<void *>(new model);
}

extern "C" void set_rgb(void *_self, double r, double g, double b){
	auto self = safe_cast<model>(_self); 
	self->rgb = {r, g, b};
}

extern "C" void set_origin(void *_self, double x, double y, double z){
	auto self = safe_cast<model>(_self); 
	self->origin = {x, y, z};
}

extern "C" void set_basis(void *_self, double e1x, double e1y, double e1z, double e2x, double e2y, double e2z) { 
	auto self = safe_cast<model>(_self); 
	self->add_basis({e1x, e1y, e1z}, {e2x, e2y, e2z});
}

extern "C" void *polyline() {
	return static_cast<void *>(new std::vector<std::array<double, 3>>);
}

extern "C" void line_add_vertex(void * _line, double x, double y, double z){
	auto line = safe_cast<std::vector<std::array<double, 3>>>(_line); 
	line->push_back({x, y, z});
}

extern "C" void add_line(void *_self, void *_line) {
	auto self = safe_cast<model>(_self); 
	auto line = safe_cast<std::vector<std::array<double, 3>>>(_line); 
	self->add_line(*line);
	delete line;
}

extern "C" void add_box(void *_self, double dx, double dy, double dz) {
	auto self = safe_cast<model>(_self); 
	self->add_box({dz, dy, dz});
}

void *polyhedron() {
	return static_cast<void *>(new polyhedron_t());
}

void polyhedron_add_vertex(void *_polyhedron, double x, double y, double z) {
	auto polyhedron = safe_cast<polyhedron_t>(_polyhedron); 
	polyhedron->v.push_back({x, y, z});
}

void polyhedron_add_face(void *_polyhedron, long int i, long int j, long int k) {
	auto polyhedron = safe_cast<polyhedron_t>(_polyhedron); 
	polyhedron->f.push_back({i, j, k});
}

void add_polyhedron(void *_self, void *_polyhedron) {
	auto self = safe_cast<model>(_self); 
	auto polyhedron = safe_cast<polyhedron_t>(_polyhedron); 
	self->add_polyhedron(*polyhedron);
	delete polyhedron;
}

extern "C" void write_obj(void *_self, int argc, char **argv) {
	auto self = safe_cast<model>(_self); 
	std::string filename = (argc == 3) ? argv[2] : argv[1];
	filename = filename.substr(0, filename.find_last_of('.'));
	self->write( filename );
}
