%{
#include <stdio.h>
#include <locale.h>
#include <stdlib.h>
#include <string.h>
#include "dawn_2_obj.h"

#ifndef YYSTYPE
# define YYSTYPE char*
#endif

#define YYDEBUG 1

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

void *self;
void *line;

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

program : program statement NEWLINE
        | statement NEWLINE
        ;

statement : VERSION
          | DIV
          | SETCAMERA
          | OPENDEVICE
          | BEGINMODELING
          | bounding_box 
          | PVNAME
          | color_rgb
          | force_wireframe
          | origin
          | base_vector
          | box
          | polyline
          | pl_vertex
          | endpolyline
          | ENDMODELING
          | DRAWALL
          | CLOSEDEVICE
          ;

bounding_box    : BOUNDINGBOX NUM NUM NUM NUM NUM NUM
                { printf("%f %f %f %f %f %f\n", atof($2), atof($3), atof($4), atof($5), atof($6), atof($6) ); } 
force_wireframe : FORCEWIREFRAME NUM;
base_vector     : BASEVECTOR NUM NUM NUM NUM NUM NUM ; 

color_rgb : COLORRGB  NUM  NUM  NUM
          { set_rgb(self, atof($2), atof($3), atof($4) ); }; 
origin    : ORIGIN  NUM  NUM  NUM
          { set_origin(self, atof($2), atof($3), atof($4) ); }; 
pl_vertex : PLVERTEX  NUM  NUM  NUM
          { set_vertex(line, atof($2), atof($3), atof($4) ); }; 
box       : BOX  NUM  NUM  NUM
          { add_box(self, atof($2), atof($3), atof($4) ); yytext_gc_run(); }; 

polyline    : POLYLINE 
            { line = polyline(); };
endpolyline : ENDPOLYLINE
            { add_line(self, line); yytext_gc_run(); }
%%

int main(int argc, char **argv) {
	
	extern int yydebug;
	yydebug = 1;
	
	setlocale(LC_ALL, "en_GB.UTF-8");
	extern FILE *yyin;
	
	if (argc != 3) {
		fprintf(stderr, "usage: dawn_2_obj INPUT.prim OUTPUT.obj\n");
		return 1;
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
	
	write_obj(self, argv[2]);
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
