#ifndef __DAWN_2_OBJ_MODEL__
#define __DAWN_2_OBJ_MODEL__
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <string>
#include <vector>

namespace dawn {
class driver {
public:
	driver( const std::string& fn_in, const std::string& fn_out ) : fp_in( fn_in, std::ios::in ), filename_out( fn_out ) {
		std::printf( "constructing ...\n" );
		if ( !fp_in.is_open() ) {
			std::fprintf( stderr, "Unable to open file %s.", fn_in.c_str() );
			exit( 1 );
		}
	};

	driver( const std::string& fn_in ) : driver( fn_in, fn_in.substr( 0, fn_in.find_last_of( '.' ) ) ) {};

	void add_line( const std::vector<std::array<double, 3>>& line );
	void add_line( const std::vector<std::array<double, 4>>& line );
	void add_line_v1( const std::vector<std::array<double, 4>>& line );
	void add_box( const std::array<double, 3>& oset, const std::array<double, 3>& s );
	void add_sphere( const std::array<double, 3>& o, double r, double res );
	void add_tubs( double rmin, double rmax, double dz, double sphi, double dphi );
	void add_cons( std::array<double, 2> rmin, std::array<double, 2> rmax, double dz, double sphi, double dphi );
	void add_polyhedron( const std::vector<std::array<double, 3>>& v, const std::vector<std::vector<long int>>& f );
	void add_trap( double dz, double theta, double phi, double h1, double bl1, double tl1, double alpha1, double h2,
	               double bl2, double tl2, double alpha2 );
	void add_trd( double dx1, double dx2, double dy1, double dy2, double dz );
	void add_torus( double rmin, double rmax, double rtor, double phi_0, double phi );

	void set_basis( const std::array<double, 3>& e1, const std::array<double, 3>& e2 );
	void write();

	std::array<double, 3> origin;
	std::array<double, 3> rgb;
	std::string pv_name;
	size_t ndiv;

	const double markres = 3.;
	const double marksize = 1 / 10.;
	double time_max = std::numeric_limits<double>::max();
	double time_min = std::numeric_limits<double>::lowest();
	virtual double time_step( double t ) {
		// return  (t < 1.0) ? 0.0005
		//                   : 0.05;
		return 0.1 * ( 0.5 * std::tanh( t - 6 ) + 0.5 ) + 0.0001;
	};

	std::ifstream fp_in;
	const std::string filename_out;

	std::vector<std::regex> remove_pv;

private:
	struct object {
		std::vector<std::vector<size_t>> f;
		std::vector<std::vector<double>> v;
		size_t c;
	};

	object* add_object( const std::string& type );
	template <size_t N>
	size_t add_unique_vertex( object* obj, const std::array<double, N>& p );
	void six_faces( object* obj, const std::array<size_t, 8>& vs ) const;
	size_t get_rgb( const std::array<double, 3>& v );
	bool filter_object();

	std::vector<std::array<double, 3>> colour_db;
	std::map<std::string, object> obj_db;
	std::map<std::string, object> line_db;

	std::array<std::array<double, 3>, 3> xform = { std::array<double, 3>{ 1, 0, 0 }, //
	                                               std::array<double, 3>{ 0, 1, 0 }, //
	                                               std::array<double, 3>{ 0, 0, 1 } };
};
} // NAMESPACE dawn

#endif
