#ifndef __lexer_h__ 
#define __lexer_h__

#if ! defined(yyFlexLexerOnce)
	#include <FlexLexer.h>
#endif

#include "parser.tab.hh"

namespace dawn {
	class lexer : public yyFlexLexer {
	public:
		lexer(std::istream *in) : yyFlexLexer(in) {};
		
		using FlexLexer::yylex;
		
		virtual int yylex(
			parser::semantic_type * const lval
		);
	private:
		dawn::parser::semantic_type *yylval = nullptr;
	};
}

#endif 
