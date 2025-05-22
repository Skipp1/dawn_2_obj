#include "driver.h"

#include <memory>
#include <regex>
#include <stdarg.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace dawn {
template <size_t N>
size_t driver::add_unique_vertex( driver::object* obj, const std::array<double, N>& p ) {
	/* Adds a new vertex applying the origin offset and the basis transfomation.
	 * Returns the index of the vertex for use later in the obj file.
	 *
	 * We check to make sure the vertex is unique before adding it to prevent the obj file from getting too big.
	 * This is important when creating tracks as every child track will have shared verticies.
	 * It also helps blender as it makes all fances part of the same mesh.
	 */
	std::array<double, N> xp = p;
	for ( size_t i = 0; i < 3; i++ ) {
		xp[i] = std::inner_product( xform[i].begin(), xform[i].end(), p.begin(), static_cast<double>( 0. ) );
		xp[i] += origin[i];
	}

	auto it = std::find_if(
	    obj->v.crbegin(), ( obj->v.crend() - obj->v.crbegin() < 1000 ) ? obj->v.crend() : obj->v.crbegin() + 1000,
	    [&xp]( const auto& x ) -> bool {
		    return std::inner_product( xp.begin(), xp.end(), x.begin(), true, std::logical_and<>(), std::equal_to<>() );
	    } );

	if ( it == ( ( obj->v.crend() - obj->v.crbegin() < 1000 ) ? obj->v.crend() : obj->v.crbegin() + 1000 ) ) {
		obj->v.push_back( std::vector( xp.begin(), xp.end() ) );
		return obj->v.size();
	} else {
		return obj->v.size() + std::distance( it, obj->v.crbegin() ); // distance is negative?!?!
	}
}

size_t driver::get_rgb( const std::array<double, 3>& v ) {
	auto it = std::find_if( colour_db.begin(), colour_db.end(), [&v]( const std::array<double, 3>& x ) -> bool {
		return v[0] == x[0] && v[1] == x[1] && v[2] == x[2];
	} );
	if ( it == colour_db.end() ) {
		colour_db.push_back( v );
		return colour_db.size();
	} else {
		return std::distance( colour_db.begin(), it ) + 1;
	}
}

driver::object* driver::add_object( const std::string& dflt ) {
	size_t colour_idx = get_rgb( rgb );
	std::string name = ( pv_name.size() == 0 ) ? dflt + "_" + std::to_string( colour_idx )
	                                           : pv_name.substr( 0, pv_name.find_last_not_of( "0123456789.#" ) );
	pv_name.clear();

	obj_db[name].c = colour_idx;
	return &obj_db[name];
}

void driver::set_basis( const std::array<double, 3>& e1, const std::array<double, 3>& e2 ) {
	/* Defines the current rotation that the next object created will follow.
	 * \begin{verbatim}
	 * /BaseVector 1 0 0 0 1 0
	 * \end{verbatim}
	 *
	 * the obj file does not have an equivelent so all the transformations must be applied.
	 *
	 * DAWN defines it's 'default' basis vectors as $e_1=(1, 0, 0)$, $e_2=(0, 1, 0)$ and $e_3 = e_1 \wedge e_2$.
	 * any change to this basis defines a new transformed set of coordinates.
	 * As such we need to transform any object in the new coordinates ($e'_i$) back to the standard ones before writing to
	 * the obj file.
	 *
	 * To do this we construct a transformation matix $\mathbf{X}$ such that $X e_i = e'_i$
	 * As such a transformation is a special orthogonal transformation, we cn use $X X^{-1} = X X^\intercal = \mathbb{1}$
	 */
	const std::array<double, 3> e3 = { +e1[1] * e2[2] - e1[2] * e2[1],   //
	                                   -e1[0] * e2[2] + e1[2] * e2[0],   //
	                                   +e1[0] * e2[1] - e1[1] * e2[0] }; //
	for ( size_t i = 0; i < 3; i++ ) { xform[i] = { e1[i], e2[i], e3[i] }; }
}

void driver::six_faces( driver::object* obj, const std::array<size_t, 8>& vs ) const {
	obj->f.push_back( { vs[0], vs[1], vs[3], vs[2] } );
	obj->f.push_back( { vs[4], vs[5], vs[1], vs[0] } );
	obj->f.push_back( { vs[5], vs[7], vs[3], vs[1] } );
	obj->f.push_back( { vs[7], vs[6], vs[2], vs[3] } );
	obj->f.push_back( { vs[6], vs[4], vs[0], vs[2] } );
	obj->f.push_back( { vs[6], vs[7], vs[5], vs[4] } );
}

bool driver::filter_object() {
	return std::accumulate( remove_pv.begin(), remove_pv.end(), false,
	                        [&]( bool b, std::regex a ) { return b || std::regex_match( pv_name, a ); } );
}

void driver::write() {
	std::cerr << "writing..." << std::endl;
	std::ofstream fp_mat;
	FILE* fp_out = fopen( ( filename_out + ".obj" ).c_str(), "w" );
	fprintf( fp_out, "mtllib %s.mtl\n", filename_out.c_str() );

	size_t v_offset = 0;
	for ( const auto& [name, obj] : obj_db ) {
		fprintf( fp_out, "o %s\n", name.c_str() );
		fprintf( fp_out, "usemtl colour_%ld\n", obj.c );
		for ( const auto& v : obj.v ) { fprintf( fp_out, "v %.17g %.17g %.17g\n", v[0], v[1], v[2] ); }
		for ( const auto& fs : obj.f ) {
			fprintf( fp_out, "f" );
			for ( const auto& f : fs ) fprintf( fp_out, " %ld", f + v_offset );
			fprintf( fp_out, "\n" );
		}
		v_offset += obj.v.size();
	}

	for ( const auto& [name, line] : line_db ) {
		if ( line.v.size() == 0 ) continue;
		fprintf( fp_out, "o %s\n", name.c_str() );
		fprintf( fp_out, "usemtl colour_%ld\n", line.c );
		for ( const auto& v : line.v ) {
			fprintf( fp_out, "v" );
			for ( const double& vi : v ) fprintf( fp_out, " %.17g", vi );
			fprintf( fp_out, "\n" );
		}
		for ( const auto& fs : line.f ) {
			fprintf( fp_out, "l" );
			for ( size_t len = 0; const auto& f : fs ) {
				fprintf( fp_out, " %ld", f + v_offset );
				if ( len % 10000 == 0 && len > 1 ) fprintf( fp_out, "\nl" ); // don't make lines too long or blender chokes.
				len++;
			}
			fprintf( fp_out, "\n" );
		}
		v_offset += line.v.size();
	}
	fclose( fp_out );

	fp_mat.open( filename_out + ".mtl" );
	for ( size_t i = 1; const auto& v : colour_db ) {
		fp_mat << "newmtl colour_" << i++ << std::endl;
		fp_mat << "Ns 250" << std::endl;
		fp_mat << "Ka 1 1 1" << std::endl;
		fp_mat << "Kd " << v[0] << " " << v[1] << " " << v[2] << " " << std::endl;
		fp_mat << "Ks 0.5 0.5 0.5" << std::endl;
		fp_mat << "Ke 0 0 0" << std::endl;
		fp_mat << "Ni 1.45" << std::endl;
		fp_mat << "d 1" << std::endl;
		fp_mat << "illum 2" << std::endl;
		fp_mat << std::endl;
	}
	fp_mat.close();
}

void driver::add_line( const std::vector<std::array<double, 3>>& line ) {
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
	 *
	 * We write out all the lines at the end as there are normally a lot of them and its best if we group them.
	 */

	size_t colour_idx = get_rgb( rgb );
	pv_name = ( !pv_name.empty() ) ? "line_" + std::to_string( get_rgb( { rgb[0], rgb[1], rgb[2] + stod( pv_name ) } ) )
	                                     + "_" + ( ( stod( pv_name ) > 1000030060 ) ? std::string( "heavy" ) : pv_name )
	                               : "line";

	if ( filter_object() or line.size() < 2 ) {
		pv_name.clear();
		return;
	}

	auto l = &line_db[pv_name];
	l->c = colour_idx;

	std::vector<size_t> segments;
	for ( const auto& p : line ) { segments.push_back( add_unique_vertex( l, p ) ); }
	l->f.push_back( segments );

	pv_name.clear();
}

// Modified version of geant that writes out timing information of tracks too.
// This is a nonstandard dawn file, but it does allow us to make pretty animations.
/*
diff --git a/source/visualization/FukuiRenderer/include/G4FRSceneFunc.icc
b/source/visualization/FukuiRenderer/include/G4FRSceneFunc.icc index 0f836d1ee..050df4d39 100644
--- a/source/visualization/FukuiRenderer/include/G4FRSceneFunc.icc
+++ b/source/visualization/FukuiRenderer/include/G4FRSceneFunc.icc
@@ -78,11 +98,19 @@ void G4FRSCENEHANDLER::AddPrimitive (const G4Polyline& polyline)
  SendStr( FR_POLYLINE );

  //----- vertices on polyline
+	if ( nPoints==2
+		&& !std::isnan(polyline.GetVisAttributes()->GetEndTime())
+		&& !std::isnan(polyline.GetVisAttributes()->GetStartTime())
+	) {
+		// If we are drawing with time, each polyline will only have 2 verticies.
+		SendStrDouble4( FR_PL_VERTEX, polyline.GetVisAttributes()->GetStartTime(), polyline[0].x(), polyline[0].y(),
polyline[0].z() ); +		SendStrDouble4( FR_PL_VERTEX, polyline.GetVisAttributes()->GetEndTime(),   polyline[1].x(),
polyline[1].y(), polyline[1].z() ); +	} else { for ( i = 0; i < nPoints; i++ ) {
-	  SendStrDouble3( FR_PL_VERTEX   , \
-			  polyline[i].x(), \
-			  polyline[i].y(), \
-			  polyline[i].z()   );
+			SendStrDouble3( FR_PL_VERTEX,
+				polyline[i].x(), polyline[i].y(), polyline[i].z()
+			);
+		}
  }

    //----- send ending of polyline
*/
void driver::add_line( const std::vector<std::array<double, 4>>& line ) {
	size_t colour_idx = get_rgb( rgb );
	pv_name = ( !pv_name.empty() ) ? "line_" + ( ( stod( pv_name ) > 1000030060 ) ? std::string( "heavy" ) : pv_name )
	                               : "line";

	if ( filter_object() or line.size() < 2 ) {
		printf( "filtered line\n" );
		pv_name.clear();
		return;
	}

	auto l = &line_db[pv_name];
	l->c = colour_idx;

	if ( line.size() != 2 ) {
		std::cerr << "skipped polyline " << pv_name << std::endl
		          << "\tReason: polyline had timing information, but more than 2 vertices." << std::endl;
		return;
	}

	auto vv_op3 = []( auto op, const std::array<double, 3>& a, const std::array<double, 3>& b ) -> std::array<double, 3> {
		return { op( a[0], b[0] ), op( a[1], b[1] ), op( a[2], b[2] ) };
	};
	auto sv_op3 = []( auto op, double s, const std::array<double, 3>& a ) -> std::array<double, 3> {
		return { op( s, a[0] ), op( s, a[1] ), op( s, a[2] ) };
	};

	std::array<double, 3> start_pnt = { line[0][1], line[0][2], line[0][3] };
	std::array<double, 3> end_pnt = { line[1][1], line[1][2], line[1][3] };
	auto a = start_pnt;

	if ( l->f.size() == 0
	     || ( l->f.back().size() != 0
	          && l->f.back().back() != add_unique_vertex<6>( l, { a[0], a[1], a[2], line[0][0], rgb[0], rgb[1] } ) ) ) {
		l->f.push_back( { add_unique_vertex<6>( l, { line[0][1], line[0][2], line[0][3], line[0][0], rgb[0], rgb[1] } ) } );
	}

	for ( double t = line[0][0]; t < line[1][0]; t += time_step( t ) ) {
		// Geant only puts the start and end times at the point where the particle interacts.
		// If we want a smooth animation between the interactions then we need to interpolate.
		// This is just a simple linear interpolation between 2 points in 3d space.
		double f = ( t - line[0][0] ) / ( line[1][0] - line[0][0] );

		a = vv_op3( std::plus<double>{}, start_pnt,
		            sv_op3( std::multiplies<double>{}, f, vv_op3( std::minus<double>{}, end_pnt, start_pnt ) ) );
		size_t v = add_unique_vertex<6>( l, { a[0], a[1], a[2], t, rgb[0], rgb[1] } );
		if ( !( l->f.back().size() != 0 && l->f.back().back() == v ) ) { l->f.back().push_back( v ); }
	}
	l->f.back().push_back(
	    add_unique_vertex<6>( l, { end_pnt[0], end_pnt[1], end_pnt[2], line[1][0], rgb[0], rgb[1] } ) );
}

void driver::add_line_v1( const std::vector<std::array<double, 4>>& line ) {
	size_t colour_idx = get_rgb( rgb );
	pv_name = ( !pv_name.empty() ) ? "line_" + std::to_string( get_rgb( { rgb[0], rgb[1], rgb[2] + stod( pv_name ) } ) )
	                                     + "_" + ( ( stod( pv_name ) > 1000030060 ) ? std::string( "heavy" ) : pv_name )
	                               : "line";

	if ( filter_object() or line.size() < 2 ) {
		printf( "filtered line\n" );
		pv_name.clear();
		return;
	}

	auto l = &line_db[pv_name];
	l->c = colour_idx;

	if ( line.size() != 2 ) {
		std::cerr << "skipped polyline " << pv_name << std::endl
		          << "\tReason: polyline had timing information, but more than 2 vertices." << std::endl;
		return;
	}

	auto vv_op3 = []( auto op, const std::array<double, 3>& a, const std::array<double, 3>& b ) -> std::array<double, 3> {
		return { op( a[0], b[0] ), op( a[1], b[1] ), op( a[2], b[2] ) };
	};
	auto sv_op3 = []( auto op, double s, const std::array<double, 3>& a ) -> std::array<double, 3> {
		return { op( s, a[0] ), op( s, a[1] ), op( s, a[2] ) };
	};

	std::array<double, 3> a = { line[0][1], line[0][2], line[0][3] };
	std::array<double, 3> b = { line[1][1], line[1][2], line[1][3] };

	if ( time_min < line[0][0] && line[0][0] < time_max ) {
		if ( time_min < line[1][0] && line[1][0] < time_max ) {
			if ( l->f.size() != 0 && l->v.back()[0] == line[0][1] && l->v.back()[1] == line[0][2]
			     && l->v.back()[2] == line[0][3] ) {
				l->f.back().push_back( add_unique_vertex( l, b ) );
			} else {
				l->f.push_back( { add_unique_vertex( l, a ), add_unique_vertex( l, b ) } );
			}
		} else {
			// Geant only puts the start and end times at the point where the particle interacts.
			// If we want a smooth animation between the interactions then we need to interpolate.
			// This is just a simple linear interpolation between 2 points in 3d space.
			double f = ( line[1][0] < time_max ) ? 1. : ( time_max - line[0][0] ) / ( line[1][0] - line[0][0] );

			std::array<double, 3> new_pnt = vv_op3(
			    std::plus<double>{}, a, sv_op3( std::multiplies<double>{}, f, vv_op3( std::minus<double>{}, b, a ) ) );
			l->f.push_back( { add_unique_vertex( l, a ), add_unique_vertex( l, new_pnt ) } );
		}
		return;
	}
}

void driver::add_torus( double rmin, double rmax, double rtor, double phi_0, double phi ) {
	std::cerr << "skipped torus " << pv_name << std::endl << "\tReason: torus is not yet implemented" << std::endl;
}

void driver::add_cons( std::array<double, 2> rmin, std::array<double, 2> rmax, double dz, double phi_0, double phi ) {
	if ( filter_object() ) return;
	auto obj = add_object( "cons" );

	std::array<std::vector<size_t>, 4> rings;

	for ( size_t ring_idx = 0; ring_idx < rings.size(); ring_idx++ ) {
		for ( size_t i = 0; i < ndiv; i++ ) {
			double r = ( ring_idx == 0 ) ? rmax[0] : ( ring_idx == 1 ) ? rmin[0] : ( ring_idx == 2 ) ? rmin[1] : rmax[1];
			double x = r * sin( phi_0 + i * phi / ( (double)ndiv ) );
			double y = r * cos( phi_0 + i * phi / ( (double)ndiv ) );
			double z = ( ( ring_idx / 2 ) % 2 == 0 ) ? dz : -dz;
			rings[ring_idx].push_back( add_unique_vertex<3>( obj, { x, y, z } ) );
		}
	}

	for ( size_t j = 0; j < rings.size(); j++ )
		for ( size_t i = 0; i < ndiv; i++ ) {
			obj->f.push_back( { rings[j][i], rings[j][( i + 1 ) % ndiv], rings[( j + 1 ) % 4][( i + 1 ) % ndiv],
			                    rings[( j + 1 ) % 4][i] } );
		}
}

void driver::add_sphere( const std::array<double, 3>& o, double r, double res ) {
	/* Add a sphere of radius $rR at the point $\vec o$.
	 *
	 * We currently construct a uv-sphere as that is the easiest for me to code up.
	 * This is currently only used as a marker shape, but if I ever am bothered enough to make a geant simulation that has
	 * a sphere in it, I'll update the parsing code.
	 *
	 * The defintion for $(x, y, z)$ is taken from the conversion from spherical coordinates to cartesians.
	 */

	if ( filter_object() ) return;
	auto obj = add_object( "sphere" );
	std::vector<size_t> cur_ring( std::lrint( 2 * res ) + 1 );
	std::vector<size_t> old_ring( cur_ring.size() );

	for ( long i = 0; i <= std::lrint( res ); i++ ) {
		for ( size_t j = 0; j < cur_ring.size(); j++ ) {
			double z = r * std::cos( M_PI * (double)i / res );
			double y = r * std::sin( M_PI * (double)i / res ) * std::sin( 2 * M_PI * (double)j / ( 2 * res ) );
			double x = r * std::sin( M_PI * (double)i / res ) * std::cos( 2 * M_PI * (double)j / ( 2 * res ) );
			cur_ring[j] = add_unique_vertex<3>( obj, { x + o[0], y + o[1], z + o[2] } );
		}
		if ( i != 0 ) {
			for ( size_t j = 0; j < cur_ring.size(); j++ ) {
				obj->f.push_back( { cur_ring[j], cur_ring[( j + 1 ) % old_ring.size()], old_ring[( j + 1 ) % old_ring.size()],
				                    old_ring[j] } );
			}
		}
		old_ring = cur_ring;
	}
}

void driver::add_polyhedron( const std::vector<std::array<double, 3>>& vs,
                             const std::vector<std::vector<long int>>& f ) {
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
	 */

	if ( filter_object() ) return;
	auto obj = add_object( "mesh" );

	std::vector<size_t> m;
	for ( const auto& v : vs ) { m.push_back( add_unique_vertex( obj, v ) ); }

	for ( const auto& idx : f ) {
		std::vector<size_t> facet;
		for ( const auto& i : idx ) { facet.push_back( m[i - 1] ); }
		obj->f.push_back( facet );
		facet.clear();
	}
}

void driver::add_tubs( double rmin, double rmax, double dz, double phi_0, double phi ) {
	/*
	 * /Tubs describes a tube or its segment in azimuthal angle with constant minimum (inside) and maximum (outside)
	 * radii. Its height extends along the z axis. The top facet is on plane $z = +\dd z$, and the bottom facet on plane
	 * $z = -\dd z$. Centers of the top and bottom facets are $(0, 0, +\dd z)$ and (0, 0, -\dd z)$, respectively. /Tubs is
	 * equivalent to /Cons with rmin1 = rmin2 and rmax1 = rmax2. /Tubs corresponds to class G4Tubs of GEANT4.
	 * \begin{tabular}{ll}
	 * Format & /Tubs rmin rmax dz sphi dphi                     \\\hline
	 * rmin   & minimum (inside) radius                          \\
	 * rmax   & maximum (outside) radius                         \\
	 * dz     & half height along the $z$ axis                   \\
	 * sphi   & starting azimuthal angle, $[−2\pi, +2\pi]$       \\
	 * dphi   & extension of azimuthal angle, dphi = $[0, 2\pi]$ \\
	 * \end{tabular}
	 */

	if ( filter_object() ) return;
	auto obj = add_object( "tubs" );
	std::array<std::vector<size_t>, 4> rings;

	for ( size_t ring_idx = 0; auto& ring : rings ) {
		for ( size_t i = 0; i < ndiv; i++ ) {
			double r = ( ( ( ring_idx + 1 ) / 2 ) % 2 == 0 ) ? rmin : rmax;
			double x = r * sin( phi_0 + i * phi / ( (double)ndiv ) );
			double y = r * cos( phi_0 + i * phi / ( (double)ndiv ) );
			double z = ( ( ring_idx / 2 ) % 2 == 0 ) ? dz : -dz;
			ring.push_back( add_unique_vertex<3>( obj, { x, y, z } ) );
		}
		ring_idx++;
	}

	for ( size_t j = 0; j < rings.size(); j++ )
		for ( size_t i = 0; i < ndiv; i++ ) {
			obj->f.push_back( { rings[j][i], rings[j][( i + 1 ) % ndiv], rings[( j + 1 ) % 4][( i + 1 ) % ndiv],
			                    rings[( j + 1 ) % 4][i] } );
		}
}

void driver::add_trap( double dz, double theta, double phi, double h1, double bl1, double tl1, double alpha1, double h2,
                       double bl2, double tl2, double alpha2 ) {
	/*
	 * /Trap is a skewed version of /Trd, i.e., asymmetric pyramid with its upper part cut away.
	 * Its top and bottom facets are asymmetric trapezoids.
	 * Its skew is expressed with direction of a line joining the centers of the top and bottom trapezoids.
	 * (Here we deﬁne a center of a trapezoid by an intersection point of two lines passing through middle points of
	 * opponent edges.) This line should pass through the origin. The top trapezoid is on plane z = +dz, and the bottom
	 * trapezoid on plane z = -dz. Note that there are 11 parameters, but only 9 are really independent. Meanings of some
	 * parameters are similar to those of /Parallelepiped. /Trap corresponds to class G4Trap of GEANT4.
	 * \begin{tabular}{ll}
	 * dz     & half height of this shape along the z axis \\
	 * theta  & polar angle of the line expressing skew \\
	 * phi    & azimuthal angle of the line expressing skew \\
	 * h1     & half height of the bottom trapezoid along the y axis \\
	 * tl1    & half length along the x direction of the side at minumum y of the bottom trapezoid \\
	 * alpha1 & angle formed the y axis and a line joining middle points of the two x-directional edges of the bottom
	 * trapezoid \\
	 * h2     & half height of the top trapezoid along the y axis \\
	 * bl2    & half length along the x direction of the side at minimum y of the top trapezoid \\
	 * tl2    & half length along the x direction of the side at maximum y of the top trapezoid \\
	 * alpha2 & angle formed the y axis and a line joining middle points of the two x-directional edges of the bottom
	 * trapezoid \\
	 * \end{tabular}
	 */

	if ( filter_object() ) return;
	auto obj = add_object( "trap" );
	double cx = dz * tan( theta ) * cos( phi );
	double cy = dz * tan( theta ) * sin( phi );

	double dx1 = h1 * tan( alpha1 );
	double dx2 = h2 * tan( alpha2 );

	std::array<size_t, 8> vs;
	vs[0] = add_unique_vertex<3>( obj, { +tl2 - dx2 + cx, +h2 + cy, +dz } );
	vs[1] = add_unique_vertex<3>( obj, { -tl2 - dx2 + cx, +h2 + cy, +dz } );
	vs[2] = add_unique_vertex<3>( obj, { +bl2 + dx2 + cx, -h2 + cy, +dz } );
	vs[3] = add_unique_vertex<3>( obj, { -bl2 + dx2 + cx, -h2 + cy, +dz } );
	vs[4] = add_unique_vertex<3>( obj, { +tl1 - dx1 - cx, +h1 - cy, -dz } );
	vs[5] = add_unique_vertex<3>( obj, { -tl1 - dx1 - cx, +h1 - cy, -dz } );
	vs[6] = add_unique_vertex<3>( obj, { +bl1 + dx1 - cx, -h1 - cy, -dz } );
	vs[7] = add_unique_vertex<3>( obj, { -bl1 + dx1 - cx, -h1 - cy, -dz } );

	six_faces( obj, vs );
}

void driver::add_trd( double dx1, double dx2, double dy1, double dy2, double dz ) {
	if ( filter_object() ) return;
	auto obj = add_object( "trd" );
	std::array<size_t, 8> vs;
	vs[0] = add_unique_vertex<3>( obj, { +dx2, +dy2, +dz } );
	vs[1] = add_unique_vertex<3>( obj, { -dx2, +dy2, +dz } );
	vs[2] = add_unique_vertex<3>( obj, { +dx2, -dy2, +dz } );
	vs[3] = add_unique_vertex<3>( obj, { -dx2, -dy2, +dz } );
	vs[4] = add_unique_vertex<3>( obj, { +dx1, +dy1, -dz } );
	vs[5] = add_unique_vertex<3>( obj, { -dx1, +dy1, -dz } );
	vs[6] = add_unique_vertex<3>( obj, { +dx1, -dy1, -dz } );
	vs[7] = add_unique_vertex<3>( obj, { -dx1, -dy1, -dz } );

	six_faces( obj, vs );
}

void driver::add_box( const std::array<double, 3>& oset, const std::array<double, 3>& s ) {
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

	if ( filter_object() ) return;
	auto obj = add_object( "box" );
	std::array<size_t, 8> vs;
	vs[0] = add_unique_vertex<3>( obj, { +s[0] + oset[0], +s[1] + oset[1], +s[2] + oset[2] } );
	vs[1] = add_unique_vertex<3>( obj, { -s[0] + oset[0], +s[1] + oset[1], +s[2] + oset[2] } );
	vs[2] = add_unique_vertex<3>( obj, { +s[0] + oset[0], -s[1] + oset[1], +s[2] + oset[2] } );
	vs[3] = add_unique_vertex<3>( obj, { -s[0] + oset[0], -s[1] + oset[1], +s[2] + oset[2] } );
	vs[4] = add_unique_vertex<3>( obj, { +s[0] + oset[0], +s[1] + oset[1], -s[2] + oset[2] } );
	vs[5] = add_unique_vertex<3>( obj, { -s[0] + oset[0], +s[1] + oset[1], -s[2] + oset[2] } );
	vs[6] = add_unique_vertex<3>( obj, { +s[0] + oset[0], -s[1] + oset[1], -s[2] + oset[2] } );
	vs[7] = add_unique_vertex<3>( obj, { -s[0] + oset[0], -s[1] + oset[1], -s[2] + oset[2] } );

	six_faces( obj, vs );
}

} // namespace dawn
