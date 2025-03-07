#include "driver.h"
#include "parser.tab.hh"
#include "lexer.h"

int main(int argc, char **argv) {
	if (argc != 1+1) {
		std::fprintf(stderr, "./dawn_2_obj <input.prim>");
		return 1;
	}
	
	dawn::driver   drv(argv[1]);
	//drv.remove_pv.push_back(std::regex("^lv.*"));
	
	dawn::lexer    lex(&drv.fp_in); 
	dawn::parser parse(lex, drv);
	parse();
	
	return 0;
}
