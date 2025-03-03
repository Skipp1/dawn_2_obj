#include "driver.h"

#include <stdarg.h>

#include <numeric>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <utility>
#include <vector>
#include <map>
#include <fmt/core.h>

namespace dawn {
size_t driver::add_unique_vertex(const std::array<double, 3> &p) {
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
	
	auto it = std::find_if( vertex_db.begin(), vertex_db.end()
	                       , [&xp](const std::array<double, 3> &x)->bool{
		                       return xp[0]==x[0] && xp[1]==x[1] && xp[2]==x[2];
	                       }
	);
	
	if (it == vertex_db.end()) {
		vertex_db.push_back(xp);
		return vertex_db.size() + vertex_idx_offset;
	} else {
		return std::distance(vertex_db.begin(), it)+1 + vertex_idx_offset;
	}
}


size_t driver::get_rgb(const std::array<double, 3> &v) {
	auto it = std::find_if( colour_db.begin(), colour_db.end()
	                       , [&v](const std::array<double, 3> &x)->bool{
		                       return v[0]==x[0] && v[1]==x[1] && v[2]==x[2];
	                       }
	);
	if (it == colour_db.end()) {
		colour_db.push_back(v);
		return colour_db.size();
	} else {
		return std::distance(colour_db.begin(), it)+1;
	}
}

std::pair<std::string, bool> driver::o_name(const std::string &dflt, size_t colour_idx) {
	
	if (pv_name.size() == 0) {
		pv_name = fmt::format("{}_{}_", dflt, colour_idx);
	}
	
	std::string simplified_pv_name = pv_name.substr(0, pv_name.find_last_not_of("0123456789."));
	bool same_obj = simplified_pv_name == old_pv_name;
	old_pv_name   = simplified_pv_name;
	
	if ( o_name_db.find(simplified_pv_name) == o_name_db.end() )  {
		o_name_db[simplified_pv_name] = 0;
	}
	
	std::string retval = fmt::format("{}.{}", simplified_pv_name, o_name_db[simplified_pv_name]++);
	
	pv_name.clear();
	return std::make_pair(retval, same_obj);
}

std::stringstream driver::add_object(const std::string &dflt) {
	std::stringstream ss;
	size_t colour_idx = get_rgb(rgb);
	
	auto name = o_name(dflt, colour_idx);
	if (!name.second) {
		ss << "o " << name.first << std::endl;
		ss << "usemtl colour_" << colour_idx << std::endl;
	}
	return ss;
}

void driver::set_basis(const std::array<double, 3> &e1, const std::array<double, 3> &e2) {
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


void driver::add_torus(double rmin, double rmax, double rtor, double phi_0, double phi) {
	std::printf("skipped torus %s\n", pv_name.c_str());
}


void driver::add_cons(std::array<double, 2> rmin, std::array<double, 2> rmax, double dz, double phi_0, double phi) {

	std::stringstream ss = add_object("cons");
	
	std::array<std::vector<size_t>, 4> rings;
	
	for (size_t ring_idx =0; auto &ring : rings) {
		for (size_t i=0; i<ndiv; i++) {
			double r = (ring_idx==0) ? rmax[0] 
			         : (ring_idx==1) ? rmin[0]
			         : (ring_idx==2) ? rmin[1] : rmax[1]; 
			double x = r * sin(phi_0 + i * phi/((double)ndiv) );
			double y = r * cos(phi_0 + i * phi/((double)ndiv) );
			double z = ((ring_idx/2)%2==0)? dz : -dz; 
			ring.push_back(add_unique_vertex({x, y, z}));
		}
		ring_idx++;
	}
	
	for (size_t j=0; j<rings.size(); j++)
		for (size_t i=0; i<ndiv; i++) {
			ss << "f " << rings[j][i]
			   <<  " " << rings[j][(i+1)%ndiv]
			   <<  " " << rings[(j+1)%4][(i+1)%ndiv]
			   <<  " " << rings[(j+1)%4][i] << std::endl;
		}
	obj_db = ss.str();
	write_tmp();
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

	
	std::stringstream ss = add_object("box");
	std::array<size_t, 8> vs;
	vs[0] = add_unique_vertex({ s[0]+oset[0],  s[1]+oset[1],  s[2]+oset[2]});
	vs[1] = add_unique_vertex({-s[0]+oset[0],  s[1]+oset[1],  s[2]+oset[2]});
	vs[2] = add_unique_vertex({ s[0]+oset[0], -s[1]+oset[1],  s[2]+oset[2]});
	vs[3] = add_unique_vertex({-s[0]+oset[0], -s[1]+oset[1],  s[2]+oset[2]});
	vs[4] = add_unique_vertex({ s[0]+oset[0],  s[1]+oset[1], -s[2]+oset[2]});
	vs[5] = add_unique_vertex({-s[0]+oset[0],  s[1]+oset[1], -s[2]+oset[2]});
	vs[6] = add_unique_vertex({ s[0]+oset[0], -s[1]+oset[1], -s[2]+oset[2]});
	vs[7] = add_unique_vertex({-s[0]+oset[0], -s[1]+oset[1], -s[2]+oset[2]});
	
	eight_faces(ss, vs);
	
	obj_db = ss.str();
	write_tmp();
}

void driver::eight_faces(std::stringstream &ss, const std::array<size_t, 8> &vs) const {
	ss << "f " << vs[0] << " " << vs[1] << " " << vs[3] << " " << vs[2] << std::endl;
	ss << "f " << vs[4] << " " << vs[5] << " " << vs[1] << " " << vs[0] << std::endl;
	ss << "f " << vs[5] << " " << vs[7] << " " << vs[3] << " " << vs[1] << std::endl;
	ss << "f " << vs[7] << " " << vs[6] << " " << vs[2] << " " << vs[3] << std::endl;
	ss << "f " << vs[6] << " " << vs[4] << " " << vs[0] << " " << vs[2] << std::endl;
	ss << "f " << vs[6] << " " << vs[7] << " " << vs[5] << " " << vs[4] << std::endl;
}

void driver::add_sphere( const std::array<double, 3> &o, double r, double res) {
	/* Add a sphere of radius $rR at the point $\vec o$. 
	 * 
	 * We currently construct a uv-sphere as that is the easiest for me to code up. 
	 * This is currently only used as a marker shape, but if I ever am bothered enough to make a geant simulation that has a sphere in it, I'll update the parsing code. 
	 * 
	 * The defintion for $(x, y, z)$ is taken from the conversion from spherical coordinates to cartesians. 
	 */
	

	
	std::stringstream ss = add_object("sphere");
	std::vector<size_t> cur_ring(std::lrint(2*res)+1);
	std::vector<size_t> old_ring(cur_ring.size());
	
	for (long i=0; i<=std::lrint(res); i++){
		for (size_t j=0; j<cur_ring.size(); j++){
			double z = r * std::cos(M_PI*(double)i/res);
			double y = r * std::sin(M_PI*(double)i/res) * std::sin(2*M_PI*(double)j/(2*res));
			double x = r * std::sin(M_PI*(double)i/res) * std::cos(2*M_PI*(double)j/(2*res));
			cur_ring[j] = add_unique_vertex({x+o[0], y+o[1], z+o[2]});
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
	obj_db = ss.str();
	write_tmp();
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

	std::stringstream ss = add_object("line");
	ss << "l";
	for (const auto &p : line ) {
		ss << " " << add_unique_vertex(p);
	}
	ss << std::endl;
	obj_db = ss.str();
	write_tmp();
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
		m.push_back(add_unique_vertex(v));
	}
	

	std::stringstream ss = add_object("mesh");
	for ( const auto &idx : f ) {
		ss << "f";
		for ( const auto &i : idx ) { 
			ss << " " << m[i-1];
		}
		ss << std::endl;
	}
	obj_db = ss.str();
	write_tmp();
}

void driver::add_tubs(double rmin, double rmax, double dz, double phi_0, double phi) {
	/*
	/Tubs describes a tube or its segment in azimuthal angle with constant minimum (inside) and maximum (outside) radii.
	Its height extends along the z axis.
	The top facet is on plane $z = +\dd z$, and the bottom facet on plane $z = -\dd z$.
	Centers of the top and bottom facets are $(0, 0, +\dd z)$ and (0, 0, -\dd z)$, respectively.
	/Tubs is equivalent to /Cons with rmin1 = rmin2 and rmax1 = rmax2.
	/Tubs corresponds to class G4Tubs of GEANT4.
	\begin{tabular}{ll}
	Format & /Tubs rmin rmax dz sphi dphi                     \\\hline
	rmin   & minimum (inside) radius                          \\
	rmax   & maximum (outside) radius                         \\
	dz     & half height along the $z$ axis                   \\
	sphi   & starting azimuthal angle, $[−2\pi, +2\pi]$       \\
	dphi   & extension of azimuthal angle, dphi = $[0, 2\pi]$ \\
	\end{tabular}
	*/

	std::stringstream ss = add_object("tubs");
	std::array<std::vector<size_t>, 4> rings;
	
	for (size_t ring_idx =0; auto &ring : rings) {
		for (size_t i=0; i<ndiv; i++) {
			double r = (((ring_idx+1)/2)%2==0) ? rmin : rmax;
			double x = r * sin(phi_0 + i * phi/((double)ndiv) );
			double y = r * cos(phi_0 + i * phi/((double)ndiv) );
			double z = ((ring_idx/2)%2==0)? dz : -dz; 
			ring.push_back(add_unique_vertex({x, y, z}));
		}
		ring_idx++;
	}
	
	for (size_t j=0; j<rings.size(); j++)
		for (size_t i=0; i<ndiv; i++) {
			ss << "f " << rings[j][i]
			   <<  " " << rings[j][(i+1)%ndiv]
			   <<  " " << rings[(j+1)%4][(i+1)%ndiv]
			   <<  " " << rings[(j+1)%4][i] << std::endl;
		}
	obj_db = ss.str();
	write_tmp();
}

void driver::add_trap(double dz, double theta, double phi, double h1, double bl1, double tl1, double alpha1
                                                         , double h2, double bl2, double tl2, double alpha2 ) {
	/*
	/Trap is a skewed version of /Trd, i.e., asymmetric pyramid with its upper part cut away.
	Its top and bottom facets are asymmetric trapezoids.
	Its skew is expressed with direction of a line joining the centers of the top and bottom trapezoids.
	(Here we deﬁne a center of a trapezoid by an intersection point of two lines passing through middle points of opponent edges.)
	This line should pass through the origin.
	The top trapezoid is on plane z = +dz, and the bottom trapezoid on plane z = -dz.
	Note that there are 11 parameters, but only 9 are really independent.
	Meanings of some parameters are similar to those of /Parallelepiped.
	/Trap corresponds to class G4Trap of GEANT4.
	dz     & half height of this shape along the z axis
	theta  & polar angle of the line expressing skew
	phi    & azimuthal angle of the line expressing skew
	h1     & half height of the bottom trapezoid along the y axis
	tl1    & half length along the x direction of the side at minumum y of the bottom trapezoid
	alpha1 & angle formed the y axis and a line joining middle points of the two x-directional edges of the bottom trapezoid
	h2     & half height of the top trapezoid along the y axis
	bl2    & half length along the x direction of the side at minimum y of the top trapezoid
	tl2    & half length along the x direction of the side at maximum y of the top trapezoid
	alpha2 & angle formed the y axis and a line joining middle points of the two x-directional edges of the bottom trapezoid
	*/
	

	
	std::stringstream ss = add_object("trap");
	double cx = dz * tan(theta)*cos(phi);
	double cy = dz * tan(theta)*sin(phi);
	
	double dx1 = h1*tan(alpha1);
	double dx2 = h2*tan(alpha2);
	
	std::array<size_t, 8> vs;
	vs[0] = add_unique_vertex({ +tl2-dx2+cx, +h2+cy, +dz });
	vs[1] = add_unique_vertex({ -tl2-dx2+cx, +h2+cy, +dz });
	vs[2] = add_unique_vertex({ +bl2+dx2+cx, -h2+cy, +dz });
	vs[3] = add_unique_vertex({ -bl2+dx2+cx, -h2+cy, +dz });
	vs[4] = add_unique_vertex({ +tl1-dx1-cx, +h1-cy, -dz });
	vs[5] = add_unique_vertex({ -tl1-dx1-cx, +h1-cy, -dz });
	vs[6] = add_unique_vertex({ +bl1+dx1-cx, -h1-cy, -dz });
	vs[7] = add_unique_vertex({ -bl1+dx1-cx, -h1-cy, -dz });
	
	eight_faces(ss, vs);
	obj_db = ss.str();
	write_tmp();
}


void driver::add_trd(double dx1, double dx2, double dy1, double dy2, double dz) {

	
	std::stringstream ss = add_object("trd");
	std::array<size_t, 8> vs;
	vs[0] = add_unique_vertex({ +dx2, +dy2,  dz });
	vs[1] = add_unique_vertex({ -dx2, +dy2,  dz });
	vs[2] = add_unique_vertex({ +dx2, -dy2,  dz });
	vs[3] = add_unique_vertex({ -dx2, -dy2,  dz });
	vs[4] = add_unique_vertex({ +dx1, +dy1, -dz });
	vs[5] = add_unique_vertex({ -dx1, +dy1, -dz });
	vs[6] = add_unique_vertex({ +dx1, -dy1, -dz });
	vs[7] = add_unique_vertex({ -dx1, -dy1, -dz });
	
	eight_faces(ss, vs);
	obj_db = ss.str();
}

void driver::write_tmp() {
	/*
	 * Writes the obj file to an actual file. 
	 * The precision is set to 8dp, this can be changed. 
	 */
	if ( !filter_object() ) {
		for (const auto &v : vertex_db) {
			fp_out << "v " << v[0] << " " << v[1] << " " << v[2] << " " << std::endl;
		}
		vertex_idx_offset += vertex_db.size();
		fp_f_out << obj_db << std::endl;
	}
	vertex_db.clear();
	obj_db.clear();
}

void driver::write() {
	
	fp_f_out.close();
	std::ifstream fp_f_out_in(filename_out + ".obj_faces");
	fp_out << fp_f_out_in.rdbuf();
	fp_out.close();
	
	std::remove((filename_out + ".obj_faces").c_str());
	
	for ( size_t i=1; const auto &v : colour_db ) { 
		fp_mat << "newmtl color_" << i++ << std::endl;
		fp_mat << "Kd " << v[0] << " " << v[1] << " " << v[2] << " " << std::endl;
	}
	fp_mat.close();
}

}// NAMEPSACE dawn 
