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
%token PVNAME
%token COLORRGB
%token FORCEWIREFRAME
%token ORIGIN
%token BASEVECTOR
%token BOX
%token MARKCIRCLE2D
%token MARKCIRCLE2DS
%token MARKSQUARE2D
%token MARKSQUARE2DS
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
          | color_rgb
          | force_wireframe
          | origin
          | base_vector
          | box
          | polyline
          | polyhedron
          | mark
          | drawall
          ;

ignored : VERSION NEWLINE
        | DIV NEWLINE 
        | SETCAMERA NEWLINE
        | OPENDEVICE NEWLINE
        | BEGINMODELING NEWLINE
        | PVNAME NEWLINE
        | ENDMODELING NEWLINE
        | CLOSEDEVICE NEWLINE
        ;


%nterm <double> real;
real : REAL_MOD_INT {$$=std::move($1);}
     | INT {$$=static_cast<double>($1);}
     ;

bounding_box    : BOUNDINGBOX real real real real real real NEWLINE
                { std::printf("%f %f %f %f %f %f\n", $2, $3, $4, $5, $6, $7); } 
force_wireframe : FORCEWIREFRAME INT NEWLINE;
base_vector     : BASEVECTOR real real real real real real NEWLINE
                { drv.add_basis({std::move($2), std::move($3), std::move($4)}, {std::move($5), std::move($6), std::move($7)}); }; 

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

%nterm <std::array<double, 3>> pl_vertex;
pl_vertex : PLVERTEX real real real NEWLINE
          { $$={std::move($2), std::move($3), std::move($4)}; }; 

%nterm <std::vector<std::array<double, 3>>> multi_pl_vertex;
multi_pl_vertex : multi_pl_vertex pl_vertex { $1.push_back($2); $$ = std::move($1); }
                | pl_vertex { $$ = std::vector<std::array<double, 3>>({std::move($1)}); }
                ;

polyline : POLYLINE NEWLINE multi_pl_vertex ENDPOLYLINE NEWLINE 
         { drv.add_line(std::move($3)); };

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

int main(int argc, char **argv) {
	
	
	setlocale(LC_ALL, "en_GB.UTF-8");
	
	if (argc != 3 || argc != 2) {
		fprintf(stderr, "usage: dawn_2_obj INPUT.prim OUTPUT.obj\n");
	}
	
	auto drv = dawn::driver(argc, argv);
	auto lex = dawn::lexer(&drv.fp_in);
	
	dawn::parser parse(lex, drv);
	//parse.set_debug_level(1);
	parse();
	
	return 0;
}

void dawn::parser::error(const std::string &msg) {
	fprintf(stderr, "parser error: %s %ld \n", msg.c_str(), lex.lineno);
	return;
}
