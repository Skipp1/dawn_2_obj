YACC=bison 
LEX=flex
CC=gcc 
CXX=g++
LD=ld

CFLAGS=-Ofast -march=native -mtune=native -flto -Wall -Wextra -c
LINK=-flto

dawn_2_obj: dawn_2_obj.cpp dawn_2_obj.h lex.yy.c y.tab.c 
	${CXX} ${CFLAGS} dawn_2_obj.cpp -o dawn_2_obj.o
	${CC}  ${CFLAGS} lex.yy.c -o lex.yy.o
	${CC}  ${CFLAGS} y.tab.c -o y.tab.o
	${CXX} ${LINK}  y.tab.o lex.yy.o dawn_2_obj.o -o dawn_2_obj

lex.yy.c: dawn_parser.l 
	${LEX} $< 

y.tab.c: dawn_parser.y 
	${YACC} -dyv $<
