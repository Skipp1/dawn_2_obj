#include "driver.h"
#include "parser.tab.hh"
#include "lexer.h"


// Example of creating a volume filter:
/*
class filter_driver : public dawn::driver {
	bool filter_object() override {
		return pv_name.rfind("lv", 0) == 0;
	}
};
*/

int main(int argc, char **argv) {
	if (argc != 1+1) {
		std::fprintf(stderr, "./dawn_2_obj <input.prim>");
		return 1;
	}
	
	dawn::driver drv(argv[1]);
	dawn::lexer  lex(&drv.fp_in);
	dawn::parser parse(lex, drv);
	
	//parse.set_debug_level(1);
	parse();
	
	return 0;
}
