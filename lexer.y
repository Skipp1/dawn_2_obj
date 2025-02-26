%{

#include <cstdlib>

#include "lexer.h"

//#define ECHO printf("%s", yytext);
#define ECHO {};

#undef  YY_DECL
#define YY_DECL int dawn::lexer::yylex( dawn::parser::semantic_type * const lval )

using token = dawn::parser::token;

%}

%option yyclass="dawn::lexer"
%option debug 
%option noyywrap 
%option nodefault 
%option yylineno 
%option c++

%% 

%{
	yylval = lval;
%}


^##G4.PRIM-FORMAT-[0-9\.]+ { 
	ECHO;
	return token::VERSION;
}

"#--------------------" { ECHO; return token::DIV; }

"/BoundingBox"    { ECHO; return token::BOUNDINGBOX; }
"!SetCamera"      { ECHO; return token::SETCAMERA; }
"!OpenDevice"     { ECHO; return token::OPENDEVICE; }
"!BeginModeling"  { ECHO; return token::BEGINMODELING; }
"/ColorRGB"       { ECHO; return token::COLORRGB; }
"/ForceWireframe" { ECHO; return token::FORCEWIREFRAME; }
"/Origin"         { ECHO; return token::ORIGIN; }
"/BaseVector"     { ECHO; return token::BASEVECTOR; }
"/Box"            { ECHO; return token::BOX; }
"/MarkCircle2DS"  { ECHO; return token::MARKCIRCLE2DS; }
"/MarkCircle2D"   { ECHO; return token::MARKCIRCLE2D; }
"/MarkSquare2DS"  { ECHO; return token::MARKSQUARE2DS; }
"/MarkSquare2D"   { ECHO; return token::MARKSQUARE2D; }
"/Polyhedron"     { ECHO; return token::POLYHEDRON; }
"/Vertex"         { ECHO; return token::VERTEX; }
"/Facet"          { ECHO; return token::FACET; }
"/EndPolyhedron"  { ECHO; return token::ENDPOLYHEDRON; }
"/Polyline"       { ECHO; return token::POLYLINE; }
"/PLVertex"       { ECHO; return token::PLVERTEX; }
"/EndPolyline"    { ECHO; return token::ENDPOLYLINE; }
"!EndModeling"    { ECHO; return token::ENDMODELING; }
"!DrawAll"        { ECHO; return token::DRAWALL; }
"!CloseDevice"    { ECHO; return token::CLOSEDEVICE; }


[+-]?[0-9]+\.([0-9]*)?([eE][+-]?[0-9]+)? {
	ECHO;
	yylval->build<double>(std::atof(yytext));
	return token::REAL_MOD_INT;
}

-?[0-9]+ {
	ECHO;
	yylval->build<long int>(std::atol(yytext));
	return token::INT;
}

\/PVName.* { ECHO; return token::PVNAME;  }
"\n"       { ECHO; return token::NEWLINE; }

%%
