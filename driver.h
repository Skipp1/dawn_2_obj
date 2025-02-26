#ifndef __DAWN_2_OBJ_MODEL__
#define __DAWN_2_OBJ_MODEL__
#include <string> 
#include <array> 
#include <vector>
#include <iostream>
#include <fstream>

namespace dawn {
class driver {
	public:
	driver(const std::string &fn_in, const std::string &fn_out)
		: fp_in(fn_in, std::ios::in)
		, filename_out(fn_out)
		{
			if (!fp_in.is_open()) {
				std::fprintf(stderr, "Unable to open file %s.", fn_in.c_str());
				exit(1);
			}
		};
	
	driver(const std::string &fn_in)
		: driver(fn_in, fn_in.substr(0, fn_in.find_last_of('.')))
		{};
	
	driver(int argc, char ** argv) 
		: driver( 
			  (argc==0) ? help() 
			: (argc==1) ? driver(argv[1])
			: (argc==2) ? driver(argv[1], argv[2])
			: help()
		) 
		{};
			
	void add_line(const std::vector<std::array<double, 3>> &line);
	void add_box(const std::array<double, 3> &oset, const std::array<double, 3> &s);
	void add_sphere( const std::array<double, 3> &o, double r, double res);
	void add_polyhedron(const std::vector<std::array<double, 3>> &v, const std::vector<std::vector<long int>> &f);
	void add_basis(const std::array<double, 3> &e1, const std::array<double, 3> &e2);
	void write();
	
	std::string help();
	
	std::array<double, 3> rgb;
	std::array<double, 3> origin;
	
	const double marker_res = 5;
	const double model_res = 5;
	
	std::ifstream fp_in;
	const std::string filename_out;
	
	private:
	size_t add_unique_vertex(std::vector<std::array<double, 3>> &v, const std::array<double, 3> &p);
	
	std::vector<std::array<double, 3>> vertex_db, colour_db; 
	std::vector<std::string> obj_db;
	
	std::array<std::array<double, 3>, 3> xform = { std::array<double, 3>{1, 0, 0}
	                                             , std::array<double, 3>{0, 1, 0}
	                                             , std::array<double, 3>{0, 0, 1}};
	
	size_t group_no = 0;
};
} // NAMESPACE dawn

#endif
