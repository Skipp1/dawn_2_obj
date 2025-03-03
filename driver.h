#ifndef __DAWN_2_OBJ_MODEL__
#define __DAWN_2_OBJ_MODEL__
#include <string> 
#include <array> 
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <iomanip>

namespace dawn {
class driver {
	public:
	driver(const std::string &fn_in, const std::string &fn_out)
		: fp_in(fn_in, std::ios::in)
		, filename_out(fn_out)
	  , fp_out(fn_out + ".obj")
		, fp_mat(fn_out + ".mtl")
	  , fp_f_out(fn_out + ".obj_faces")
		{
			std::printf("constructing ...\n");
			if (!fp_in.is_open()) {
				std::fprintf(stderr, "Unable to open file %s.", fn_in.c_str());
				exit(1);
			}
			if (!fp_out.is_open()) {
				std::fprintf(stderr, "Unable to open file %s.obj", fn_out.c_str());
				exit(1);
			}
			fp_out << std::setprecision(8) << std::fixed;
			if (!fp_mat.is_open()) {
				std::fprintf(stderr, "Unable to open file %s..mat", fn_out.c_str());
				exit(1);
			}
			fp_mat << std::setprecision(8) << std::fixed;
		};
	
	driver(const std::string &fn_in)
		: driver(fn_in, fn_in.substr(0, fn_in.find_last_of('.')))
		{};
	
	void add_line(const std::vector<std::array<double, 3>> &line);
	void add_box(const std::array<double, 3> &oset, const std::array<double, 3> &s);
	void add_sphere( const std::array<double, 3> &o, double r, double res);
	void add_tubs(double rmin, double rmax, double dz, double sphi, double dphi);
	void add_cons(std::array<double, 2> rmin, std::array<double, 2> rmax, double dz, double sphi, double dphi);
	void add_polyhedron(const std::vector<std::array<double, 3>> &v, const std::vector<std::vector<long int>> &f);
	void add_trap(double dz, double theta, double phi, double h1, double bl1, double tl1, double alpha1, double h2, double bl2, double tl2, double alpha2 );
	void add_trd(double dx1, double dx2, double dy1, double dy2, double dz);
	void add_torus(double rmin, double rmax, double rtor, double phi_0, double phi);
	
	void set_basis(const std::array<double, 3> &e1, const std::array<double, 3> &e2);
	virtual bool filter_object() {return false;};
	void write_tmp();
	void write();
	
	std::string help() const;
	
	std::array<double, 3> origin;
	std::array<double, 3> rgb;
	std::string pv_name;
	size_t ndiv;
	
	const double markres  = 3.;
	const double marksize = 1/10.;
	
	std::ifstream fp_in;
	const std::string filename_out;
	
	private:
	size_t vertex_idx_offset = 0;
	
	std::ofstream fp_out;
	std::ofstream fp_mat;
	std::ofstream fp_f_out;
	
	std::pair<std::string, bool> o_name(const std::string &dflt, size_t colour_idx);
	std::string old_pv_name;
	
	std::stringstream add_object(const std::string &type);
	size_t add_unique_vertex(const std::array<double, 3> &p);
	void   eight_faces(std::stringstream &ss, const std::array<size_t, 8> &vs) const;
	size_t get_rgb(const std::array<double, 3> &v);
	
	std::vector<std::array<double, 3>> vertex_db, colour_db; 
	std::string obj_db;
	std::map<std::string, size_t> o_name_db;
	
	std::array<std::array<double, 3>, 3> xform = { std::array<double, 3>{1, 0, 0}
	                                             , std::array<double, 3>{0, 1, 0}
	                                             , std::array<double, 3>{0, 0, 1}};
	
};
} // NAMESPACE dawn

#endif
