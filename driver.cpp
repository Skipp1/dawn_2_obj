#include "driver.h"

#include <stdarg.h>

#include <numeric>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <vector>

namespace dawn {
size_t driver::add_unique_vertex(std::vector<std::array<double, 3>> &v, const std::array<double, 3> &p) {
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

void driver::add_basis(const std::array<double, 3> &e1, const std::array<double, 3> &e2) {
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

void driver::add_box(const std::array<double, 3> &oset, const std::array<double, 3> &s) {
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
	
	size_t colour_idx = add_unique_vertex(colour_db, rgb);
	ss << "o box_" << colour_idx << "." << group_no++ << std::endl;
	ss << "usemtl colour_" << colour_idx << std::endl;
	
	auto add_face = [this, &s, &ss, &oset](int f_){
		int f = std::abs(f_)-1;
		
		std::array<std::array<int, 4>, 3> A;
		A[f]       = {f_/std::abs(f_), f_/std::abs(f_), f_/std::abs(f_), f_/std::abs(f_)};
		A[(f+1)%3] = {-1, +1, +1, -1};
		A[(f+2)%3] = {-1, -1, +1, +1};
		
		ss << "f ";
		for (int i=0; i<4; i++) {
			ss << add_unique_vertex(vertex_db, {A[0][i]*s[0]+oset[0], A[1][i]*s[1]+oset[1], A[2][i]*s[2]+oset[2]}) << " ";
		}
		ss << std::endl;
	};
	add_face(-3); add_face(+3);
	add_face(-2); add_face(+2);
	add_face(-1); add_face(+1);
	
	obj_db.push_back(ss.str());
}

void driver::add_sphere( const std::array<double, 3> &o, double r, double res) {
	/* Add a sphere of radius $rR at the point $\vec o$. 
	 * 
	 * We currently construct a uv-sphere as that is the easiest for me to code up. 
	 * This is currently only used as a marker shape, but if I ever am bothered enough to make a geant simulation that has a sphere in it, I'll update the parsing code. 
	 * 
	 * The defintion for $(x, y, z)$ is taken from the conversion from spherical coordinates to cartesians. 
	 */
	
	std::stringstream ss;
	size_t colour_idx = add_unique_vertex(colour_db, rgb);
	ss << "o sphere_" << colour_idx << "." << group_no++ << std::endl;
	ss << "usemtl colour_" << colour_idx << std::endl;
	
	std::vector<size_t> cur_ring(std::lrint(2*res)+1);
	std::vector<size_t> old_ring(cur_ring.size());
	
	for (long i=0; i<=std::lrint(res); i++){
		for (size_t j=0; j<cur_ring.size(); j++){
			double z = r * std::cos(M_PI*(double)i/res);
			double y = r * std::sin(M_PI*(double)i/res) * std::sin(2*M_PI*(double)j/(2*res));
			double x = r * std::sin(M_PI*(double)i/res) * std::cos(2*M_PI*(double)j/(2*res));
			cur_ring[j] = add_unique_vertex(vertex_db, {x+o[0], y+o[1], z+o[2]});
		}
		if (i != 0) {
			for (size_t j=0; j<cur_ring.size(); j++){
				ss << "f " << cur_ring[j] 
				   <<  " " << cur_ring[(j+1)%old_ring.size()] 
				   <<  " " << old_ring[(j+1)%old_ring.size()]
				   <<  " " << old_ring[j] << std::endl;
			}
		}
		old_ring = cur_ring;
	}
	obj_db.push_back(ss.str());
}

void driver::add_line(const std::vector<std::array<double, 3>> &line) {
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
	size_t colour_idx = add_unique_vertex(colour_db, rgb);
	ss << "o line_" << colour_idx << "." << group_no++ << std::endl;
	ss << "usemtl colour_" << colour_idx << std::endl;
	ss << "l";
	for (const auto &p : line ) {
		ss << " " << add_unique_vertex(vertex_db, p);
	}
	ss << std::endl;
	obj_db.push_back(ss.str());
}

void driver::add_polyhedron(const std::vector<std::array<double, 3>> &v, const std::vector<std::vector<long int>> &f) {
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
	for ( const auto &v : v ) {
		m.push_back(add_unique_vertex(vertex_db, v));
	}
	
	std::stringstream ss;
	size_t colour_idx = add_unique_vertex(colour_db, rgb);
	ss << "o mesh_" << colour_idx << "." << group_no++ << std::endl;
	ss << "usemtl colour_" << colour_idx << std::endl;
	for ( const auto &idx : f ) {
		ss << "f";
		for ( const auto &i : idx ) { 
			ss << " " << m[i-1];
		}
		ss << std::endl;
	}
	obj_db.push_back(ss.str());
}

void driver::write() {
	/*
	 * Writes the obj file to an actual file. 
	 * The precision is set to 8dp, this can be changed. 
	 */
	std::ofstream fp_obj(filename_out + ".obj");
	fp_obj << std::setprecision(8) << std::fixed;
	for (const auto &v : vertex_db) {
		fp_obj << "v " << v[0] << " " << v[1] << " " << v[2] << " " << std::endl;
	}
	for (const auto &s : obj_db ) { 
		fp_obj << s << std::endl;
	}
	fp_obj.close();
	
	std::ofstream fp_mat(filename_out + ".mtl");
	fp_mat << std::setprecision(8) << std::fixed;
	for ( size_t i=1; const auto &v : colour_db ) { 
		fp_mat << "newmtl color_" << i++ << std::endl;
		fp_mat << "Kd " << v[0] << " " << v[1] << " " << v[2] << " " << std::endl;
	}
	fp_mat.close();
}

std::string driver::help() {
	/*
	 * Prints help and exits. 
	 * Typename is just for compatibility 
	 */
	std::fprintf(stderr, "./dawn_2_obj <input file> <output file>\n");
	exit(1);

}

}// NAMEPSACE dawn 
