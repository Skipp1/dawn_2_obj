%{
#include <stdio.h>
#include <locale.h>
#include <stdlib.h>
#include <string.h>
#include "dawn_2_obj.h"

#ifndef YYSTYPE
# define YYSTYPE char*
#endif

//#define YYDEBUG 1

void yyerror(const char *c);
int yylex(void);
int yywrap(void);

static void yytext_gc_init(void);
static void yytext_gc_cleanup(void);
static void yytext_gc_run(void);

struct yyt_gcr {
	char **d;
	size_t i;
	size_t N;
} yyt_gcr;

void *self=NULL;
void *line=NULL;
void *shape=NULL;

%}

%token VERSION
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

%token NEWLINE
%token NUM

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
          ;

ignored : VERSION NEWLINE
        | DIV NEWLINE 
        | SETCAMERA NEWLINE
        | OPENDEVICE NEWLINE
        | BEGINMODELING NEWLINE
        | PVNAME NEWLINE
        | ENDMODELING NEWLINE
        | DRAWALL NEWLINE
        | CLOSEDEVICE NEWLINE
        ;

bounding_box    : BOUNDINGBOX NUM NUM NUM NUM NUM NUM NEWLINE
                { printf("%f %f %f %f %f %f\n", atof($2), atof($3), atof($4), atof($5), atof($6), atof($6) ); } 
force_wireframe : FORCEWIREFRAME NUM NEWLINE;
base_vector     : BASEVECTOR NUM NUM NUM NUM NUM NUM NEWLINE
                { set_basis(self, atof($2), atof($3), atof($4), atof($5), atof($6), atof($7)); }; 

color_rgb : COLORRGB NUM NUM NUM NEWLINE
          { set_rgb(self, atof($2), atof($3), atof($4) ); }; 
origin    : ORIGIN NUM NUM NUM NEWLINE
          { set_origin(self, atof($2), atof($3), atof($4) ); }; 
box       : BOX NUM NUM NUM NEWLINE
          { add_box(self, atof($2), atof($3), atof($4) ); yytext_gc_run(); }; 


beginpolyline : POLYLINE NEWLINE
              { line = polyline(); };

pl_vertex : PLVERTEX NUM NUM NUM NEWLINE
          { line_add_vertex(line, atof($2), atof($3), atof($4) ); }; 

endpolyline : ENDPOLYLINE NEWLINE
            { add_line(self, line); line=NULL, yytext_gc_run(); }

multi_pl_vertex : multi_pl_vertex pl_vertex
                | pl_vertex
                ;

polyline : beginpolyline multi_pl_vertex endpolyline;

beginpolyhedron : POLYHEDRON NEWLINE
                { shape = polyhedron(); };

polyhedron_vertex : VERTEX NUM NUM NUM NEWLINE 
                  { polyhedron_add_vertex(shape, atof($2), atof($3), atof($4) ); };

polyhedron_facet : FACET NUM NUM NUM NEWLINE 
                 { polyhedron_add_face(shape, atol($2), atol($3), atol($4) ); };

endpolyhedron : ENDPOLYHEDRON NEWLINE 
              { add_polyhedron(self, shape); shape=NULL, yytext_gc_run(); } 


multi_polyhedron_vertex : multi_polyhedron_vertex polyhedron_vertex
                        | polyhedron_vertex
                        ; 

multi_polyhedron_facet : multi_polyhedron_facet polyhedron_facet
                       | polyhedron_facet
                       ;

polyhedron : beginpolyhedron multi_polyhedron_vertex multi_polyhedron_facet endpolyhedron;
          
%%

int main(int argc, char **argv) {
	
	//extern int yydebug;
	//yydebug = 1;
	
	setlocale(LC_ALL, "en_GB.UTF-8");
	extern FILE *yyin;
	
	if (argc != 3 || argc != 2) {
		fprintf(stderr, "usage: dawn_2_obj INPUT.prim OUTPUT.obj\n");
	}
	
	yyin = fopen( argv[1], "r+" );
	
	if ( yyin == NULL ) {
		fprintf(stderr, "unable to open file %s\n", argv[1] ); 
		return 1;
	}
	
	self = construct_model();
	
	yytext_gc_init();
	yyparse();
	yytext_gc_cleanup();
	fclose(yyin);
	
	write_obj(self, argc, argv);
	return 0;
}

static void yytext_gc_init(void) {
	yyt_gcr.i = 0;
	yyt_gcr.N = 64;
	yyt_gcr.d = (char **)malloc( 64 * sizeof(char**));
	return;
}

static void yytext_gc_cleanup(void) {
	yytext_gc_run();
	free(	yyt_gcr.d );
}

void yytext_dup( const char *c ) {
	yylval = strdup(c);
	if ( yyt_gcr.i >= yyt_gcr.N ) {
		yyt_gcr.d = (char **)realloc( yyt_gcr.d, ( yyt_gcr.N + 64 ) * sizeof(char **) );
		yyt_gcr.N += 64;
	}
	yyt_gcr.d[yyt_gcr.i] = yylval;
	yyt_gcr.i++;
	return;
}

static void yytext_gc_run(void) {
	for(size_t i=0; i<yyt_gcr.i; i++) {
		free(yyt_gcr.d[i]);
	}
	yyt_gcr.i = 0;
	return;
}


void yyerror(const char *c) {
	extern char *yytext;
	fprintf(stderr, "yerror %s | %s |", c, yytext);
	for (size_t i=0; i<strlen(yytext); i++) {
		fprintf(stderr, " 0x%02x ", yytext[i]);
	}
	fprintf(stderr, "\n");
	return;
}
