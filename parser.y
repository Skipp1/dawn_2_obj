%require "3.2"
%language "c++"

%defines 
%define api.namespace {dawn}
%define api.parser.class {parser}
%define parse.trace
%verbose 


%code requires{
	#define YYDEBUG 1
	
	namespace dawn {
		class driver;
		class lexer;
	}
	
	// The following definitions is missing when %locations isn't used
	# ifndef YY_NULLPTR
	#  if defined __cplusplus && 201103L <= __cplusplus
	#   define YY_NULLPTR nullptr
	#  else
	#   define YY_NULLPTR 0
	#  endif
	# endif
	
	#include <array>
	#include <string>
	#include <vector>
	
	#include "driver.h"
}

%parse-param { lexer  &lex }
%parse-param { driver &drv }

%code {
	#include "lexer.h"
	#include "driver.h"
	
	#undef  yylex
	#define yylex lex.yylex
}

%define api.value.type variant

%token NEWLINE

%token DIV
%token BOUNDINGBOX
%token SETCAMERA
%token OPENDEVICE
%token BEGINMODELING
%token COLORRGB
%token FORCEWIREFRAME
%token ORIGIN
%token BASEVECTOR
%token BOX
%token MARKCIRCLE2D
%token MARKCIRCLE2DS
%token MARKSQUARE2D
%token MARKSQUARE2DS
%token NDIV 
%token TUBS
%token TRAP
%token TRD
%token CONS
%token TORUS
%token POLYHEDRON
%token VERTEX
%token FACET
%token ENDPOLYHEDRON
%token POLYLINE
%token PLVERTEX
%token ENDPOLYLINE
%token ENDMODELING
%token DRAWALL
%token CLOSEDEVICE


%token <std::string> PVNAME
%token <std::string> VERSION
%token <double> REAL_MOD_INT
%token <long int> INT

%start program

%%
program : program statement 
        | statement 
        ;

statement : ignored
          | bounding_box 
          | pvname
          | color_rgb
          | force_wireframe
          | origin
          | base_vector
          | ndiv
          | box
          | polyline
          | polyhedron
          | tubs
          | cons
          | trap
          | trd
          | torus
          | mark
          | drawall
          ;

ignored : VERSION NEWLINE
        | DIV NEWLINE 
        | SETCAMERA NEWLINE
        | OPENDEVICE NEWLINE
        | BEGINMODELING NEWLINE
        | ENDMODELING NEWLINE
        | CLOSEDEVICE NEWLINE
        ;


%nterm <double> real;
real : REAL_MOD_INT {$$=std::move($1);}
     | INT {$$=static_cast<double>($1);}
     ;

pvname : PVNAME NEWLINE 
       { drv.pv_name = std::move($1); }
       ;

bounding_box    : BOUNDINGBOX real real real real real real NEWLINE
                { std::printf("%f %f %f %f %f %f\n", $2, $3, $4, $5, $6, $7); } 
force_wireframe : FORCEWIREFRAME INT NEWLINE;
base_vector     : BASEVECTOR real real real real real real NEWLINE
                { drv.set_basis({std::move($2), std::move($3), std::move($4)}, {std::move($5), std::move($6), std::move($7)}); }; 

drawall : DRAWALL NEWLINE { drv.write(); };

color_rgb : COLORRGB real real real NEWLINE
          { drv.rgb = {std::move($2), std::move($3), std::move($4)}; }; 
origin    : ORIGIN real real real NEWLINE
          { drv.origin = {std::move($2), std::move($3), std::move($4)}; }; 
box       : BOX real real real NEWLINE
          { drv.add_box({0., 0., 0.}, {std::move($2), std::move($3), std::move($4)}); }; 

mark_sphere : MARKCIRCLE2D 
            | MARKCIRCLE2DS
            ;

mark_box : MARKSQUARE2D 
         | MARKSQUARE2DS
         ;


mark : mark_sphere real real real real NEWLINE 
     { drv.add_sphere({std::move($2), std::move($3), std::move($4)}, $5*drv.marksize, drv.markres ); }
     | mark_box real real real real NEWLINE
     { drv.add_box({std::move($2), std::move($3), std::move($4)}, {$5*drv.marksize, $5*drv.marksize, $5*drv.marksize} ); }
     ;

ndiv : NDIV INT NEWLINE 
     { drv.ndiv = std::move($2); }
     ;

tubs : TUBS real real real real real NEWLINE
     { drv.add_tubs(std::move($2), std::move($3), std::move($4), std::move($5), std::move($6)); }
     ;

cons : CONS real real real real real real real NEWLINE
     { drv.add_cons({std::move($2), std::move($3)}, {std::move($4), std::move($5)}
                   , std::move($6), std::move($7) ,  std::move($8) ); }
     ;

trap : TRAP real real real real real real real real real real real NEWLINE
     { drv.add_trap( std::move($2) , std::move($3) , std::move($4)
                   , std::move($5) , std::move($6) , std::move($7)
                   , std::move($8) , std::move($9) , std::move($10)
                   , std::move($11), std::move($12)); }
     ;

trd : TRD real real real real real NEWLINE
     { drv.add_trd( std::move($2) , std::move($3) , std::move($4), std::move($5) , std::move($6) ); }
     ;

torus : TORUS real real real real real NEWLINE
      { drv.add_torus(std::move($2), std::move($3), std::move($4), std::move($5), std::move($6)); }
      ;

%nterm <std::array<double, 3>> pl_3vertex;
pl_3vertex : PLVERTEX real real real NEWLINE
          { $$={std::move($2), std::move($3), std::move($4)}; }; 
          
%nterm <std::array<double, 4>> pl_4vertex;
pl_4vertex : PLVERTEX real real real real NEWLINE
        { $$={std::move($2), std::move($3), std::move($4), std::move($5)}; }
        ;

%nterm <std::vector<std::array<double, 4>>> multi_pl_4vertex;
multi_pl_4vertex : multi_pl_4vertex pl_4vertex { $1.push_back($2); $$ = std::move($1); }
                 | pl_4vertex { $$ = std::vector<std::array<double, 4>> {std::move($1)}; }
                 ;
                 
%nterm <std::vector<std::array<double, 3>>> multi_pl_3vertex;
multi_pl_3vertex : multi_pl_3vertex pl_3vertex { $1.push_back($2); $$ = std::move($1); }
                 | pl_3vertex { $$ = std::vector<std::array<double, 3>>{std::move($1)}; }
                 ;


polyline : POLYLINE NEWLINE multi_pl_3vertex ENDPOLYLINE NEWLINE 
         { drv.add_line(std::move($3)); }
         | POLYLINE NEWLINE multi_pl_4vertex ENDPOLYLINE NEWLINE
         { drv.add_line(std::move($3)); }

%nterm <std::array<double, 3>> polyhedron_vertex;
polyhedron_vertex : VERTEX real real real NEWLINE
                  { $$={std::move($2), std::move($3), std::move($4)}; }; 

%nterm <std::vector<long int>> n_idx;
n_idx : n_idx INT
      { $1.push_back($2); $$ = std::move($1); } 
      | INT
      { $$ = std::vector<long int>({$1}); }
      ;

%nterm <std::vector<long int>> polyhedron_facet;
polyhedron_facet : FACET n_idx NEWLINE
                 { $$ = std::move(std::move($2)); };

%nterm <std::vector<std::array<double, 3>>> multi_polyhedron_vertex;
multi_polyhedron_vertex : multi_polyhedron_vertex polyhedron_vertex
                        { $1.push_back(std::move($2)); $$ = std::move($1); } 
                        | polyhedron_vertex
                        { $$ = std::vector<std::array<double,3>>({std::move($1)}); } 
                        ;

%nterm <std::vector<std::vector<long int>>> multi_polyhedron_facet;
multi_polyhedron_facet : multi_polyhedron_facet polyhedron_facet
                       { $1.push_back(std::move($2)); $$ = std::move($1); } 
                       | polyhedron_facet
                       { $$ = std::vector<std::vector<long int>>({std::move($1)}); } 
                       ;

polyhedron : POLYHEDRON NEWLINE multi_polyhedron_vertex multi_polyhedron_facet ENDPOLYHEDRON NEWLINE
           { drv.add_polyhedron(std::move($3), std::move($4)); };

%%

void dawn::parser::error(const std::string &msg) {
	fprintf(stderr, "parser error: %s on line %ld \n", msg.c_str(), lex.lineno+1);
	
	std::string line;
	drv.fp_in.seekg(0);
	for (size_t i=0; i<lex.lineno+1; i++) {
		std::getline(drv.fp_in, line);
	}
	std::cerr << line << std::endl;
	
	exit(1);
	return;
}
