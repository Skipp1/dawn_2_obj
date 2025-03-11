YACC=bison 
LEX=flex++
CXX=g++
LD=ld

CFLAGS=-Ofast -march=native -mtune=native -flto -Wall -Wextra -c
LINK=-flto -fwhole-program

FILES=$(shell ls)

dawn_2_obj: driver.cpp driver.h lex.yy.cc parser.tab.cc main.cpp
	${CXX} ${CFLAGS} driver.cpp -o driver.o
	${CXX} ${CFLAGS} main.cpp -o main.o
	${CXX} ${CFLAGS} lex.yy.cc -o lex.yy.o
	${CXX} ${CFLAGS} parser.tab.cc -o parser.tab.o
	${CXX} ${LINK}   parser.tab.o lex.yy.o driver.o main.o -o dawn_2_obj

lex.yy.cc: lexer.l lexer.h parser.tab.cc
	${LEX} lexer.l

parser.tab.cc: parser.y 
	${YACC} --no-lines -dv $<

clean: 
	-rm dawn_2_obj
	-rm lex.yy.cc lex.yy.o 
	-rm parser.tab.cc parser.tab.o parser.tab.hh
	-rm driver.o main.o 
	-rm parser.output 
