%option yylineno 
%option noyywrap
%{

#include "y.tab.h"
#include "dawn_2_obj.h"

//#define ECHO printf("%s", yytext);
#define ECHO {};
%}

%% 

^##G4.PRIM-FORMAT-[0-9\.]+ { 
	ECHO;
	yytext_dup(yytext); 
	return VERSION;
}

"#--------------------" { ECHO; return DIV; }

"/BoundingBox"    { ECHO; return BOUNDINGBOX; }
"!SetCamera"      { ECHO; return SETCAMERA; }
"!OpenDevice"     { ECHO; return OPENDEVICE; }
"!BeginModeling"  { ECHO; return BEGINMODELING; }
"/ColorRGB"       { ECHO; return COLORRGB; }
"/ForceWireframe" { ECHO; return FORCEWIREFRAME; }
"/Origin"         { ECHO; return ORIGIN; }
"/BaseVector"     { ECHO; return BASEVECTOR; }
"/Box"            { ECHO; return BOX; }
"/MarkCircle2DS"  { ECHO; return MARKCIRCLE2DS; }
"/MarkCircle2D"   { ECHO; return MARKCIRCLE2D; }
"/MarkSquare2DS"  { ECHO; return MARKSQUARE2DS; }
"/MarkSquare2D"   { ECHO; return MARKSQUARE2D; }
"/Polyhedron"     { ECHO; return POLYHEDRON; }
"/Vertex"         { ECHO; return VERTEX; }
"/Facet"          { ECHO; return FACET; }
"/EndPolyhedron"  { ECHO; return ENDPOLYHEDRON; }
"/Polyline"       { ECHO; return POLYLINE; }
"/PLVertex"       { ECHO; return PLVERTEX; }
"/EndPolyline"    { ECHO; return ENDPOLYLINE; }
"!EndModeling"    { ECHO; return ENDMODELING; }
"!DrawAll"        { ECHO; return DRAWALL; }
"!CloseDevice"    { ECHO; return CLOSEDEVICE; }


[+-]?[0-9]+(\.[0-9]*)?([eE][+-]?[0-9]+)? {
	ECHO; 
	yytext_dup(yytext);
	return NUM;
}

\/PVName.* { ECHO; return PVNAME; }
"\n" { ECHO; return NEWLINE; }

%%
