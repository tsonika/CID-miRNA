//**********************************************************************************
// File		: grammar_reader.h										
// Purpose	: Declaration of classes for reading and storing of MiRNA General SCFG
// Author	: Vipin Gupta (viping@gmail.com), Manish Kushwaha (manish.kushwaha@gmail.com)
// Revision	: 1.0	-	23-Feb-2005	Initial Draft (Vipin Gupta)
//		: 2.0	-	03-Jun-2005	Minor Unpdation for compatibility with the downstream code (Manish Kushwaha)
//		: 3.0	-	14-Jun-2005	Minor Unpdation for compatibility with the downstream code (Manish Kushwaha)
//**********************************************************************************

#ifndef _GRAMMAR__READER_H_
#define _GRAMMAR__READER_H_

#include "allincludes.h"
#include "cykmatrix.h"

enum TYPE_OF_RULE { _UNKNOWN_RULE = 0, _TERM_RULE = 1, _NONTERM_RULE = 2 };

class SCFGrule {
//Friend Functions required for direct member access to these functions
friend void FillFRCellUList(CykMatrixCell &,char,int);
friend bool FillFRCellDList(CykMatrixCell &,int);
friend void FillRestCellProcess1(const vector<CykMatrixCell> &,CykMatrixCell &,int,int);
friend bool FillRestCellProcess2(const vector<CykMatrixCell> &,CykMatrixCell &,int,int);
friend int main(int,char**);

private : string ruleLHS, ruleRHS;
		  double RuleScore;
		  bool UsedOnceFlag;
		  unsigned long RuleNumber;
		  unsigned int No_of_RHSsymbols;
		  TYPE_OF_RULE RuleType;

		  static unsigned long NTcount;		
		  //storing the total number of Non-terminals read
		  static unsigned long Tcount;
		  //storing the total number of Terminals read
		  static unsigned long RuleCount;
		  //storing the total number of Rules read
		  static string NT;
		  //string which stores all the non-terminal (NT) symbols in the form of *NT*
		  static string T;
		  //string which stores all the Terminal (T) symbols in the form of *T*
		  // It may be noted that Grammar assumes non-terminals to be in capital letters
		  // and terminals in small letters,
		  // and any two NT or T on RHS have a space between them
		  static string NTonLHS;
		  // This string stores the occurence of all unique NT which occur 
		  // on LHS of a rule, in form *NT* , this is to check the validity of Grammar
		  // as there should be no NT which does not give a production, i.e. only occurs on RHS

public  : SCFGrule() : RuleScore(MINUSINF), UsedOnceFlag(false), RuleNumber(0), No_of_RHSsymbols(0), RuleType(_UNKNOWN_RULE) { }
		  //Constructor Initializing basic stuff
		  
		  void InitializeRule( double & score, string & left, string & right);
		  //Function to initialize the basic rule parameters, also calls CountUpdate
			
		  void SetUsedOnceFlag(void);
		  // Simply sets the UsedOnce Flag to true to indicate the rule has been used
		  // in tree generation at least once

		  void CountUpdate(void);
		  // This function updates all the static members of this class like counts and lists 

		  static bool ValidateGrammar(void);
		  // returns true if Grammar is valid as described above
};

// This class basically instantiates the SCFGrule class
// and reads members from the Input Grammar File
vector<SCFGrule> ReadGrammarFile(ifstream * grammarfileptr = NULL);

#endif

