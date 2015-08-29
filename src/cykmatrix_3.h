//**********************************************************************************
// File		: cykmatrix.h
// Purpose	: Declaration of classes for creation of CYK Parser Matrix
// Author	: Vipin Gupta (viping@gmail.com), Manish Kushwaha (manish.kushwaha@gmail.com)
// Revision	: 1.0	-	20-Feb-2005	Initial Draft
//		: 2.0	-	03-Jun-2005	Minor Unpdation for compatibility with the downstream code (Manish Kushwaha)
//		: 3.0	-	14-Jun-2005	Minor Unpdation for compatibility with the downstream code (Manish Kushwaha)
//**********************************************************************************

#ifndef _CYK_MATRIX__H_
#define _CYK_MATRIX__H_


#include <fstream>
#include <iostream>
#include <list>
#include <vector>
#include <string>
#include <string.h>
#include <sstream>
#include <algorithm>
#include <map>

using namespace std;

#include "constants.h"



enum TYPE_OF_LIST { _UNKNOWN = 0, _UPLIST = 1, _DOWNLIST = 2, _PARENTREF =3, _TERMRULE =4 };


//**********************************************************************************
//This class simply defines the coordinates of an element in the matrix cell's list
//**********************************************************************************

struct ElementCoordinates {
	
	long row,col,index;
	TYPE_OF_LIST type;
	
	//Constructor
	ElementCoordinates(): row(-1), col(-1), index(-1), type(_UNKNOWN) { }
};

class CykMatrixCell;
//**********************************************************************************
// This class defines one element of a matrix cell's Uplist or Downlist
//**********************************************************************************

class CykListElement {
//Friend Function required for direct member access to these functions
friend void FillFRCellUList(CykMatrixCell &,char,int);
friend bool FillFRCellDList(CykMatrixCell &,int);
friend void FillRestCellProcess1(const vector<CykMatrixCell> &,CykMatrixCell &,int,int);
friend bool FillRestCellProcess2(const vector<CykMatrixCell> &,CykMatrixCell &,int,int);
friend int main(int,char**);
friend bool IsPossibleTree(const vector<CykMatrixCell> &);
friend string BuildTrace(const vector<CykMatrixCell> &);
friend void TraceFurther(ElementCoordinates,const vector<CykMatrixCell> &,ostringstream &);

private : string ruleLHS, ruleRHSpreDot, ruleRHSpostDot;
		  ElementCoordinates coordinates;
                  string Parents;
		  double score;

public  : CykListElement() : score(MINUSINF) {}

bool operator==(const CykListElement&) const;
bool operator!=(const CykListElement&) const;

};


//**********************************************************************************
// This class defines one cell of the upper triangular CYK Matrix
//**********************************************************************************

typedef vector<CykListElement> ElementList;
typedef vector<CykListElement>::iterator ElementListItr;


class CykMatrixCell {
//Friend Function required for direct member access to these functions
friend void FillFRCellUList(CykMatrixCell &,char,int);
friend bool FillFRCellDList(CykMatrixCell &,int);
friend void FillRestCellProcess1(const vector<CykMatrixCell> &,CykMatrixCell &,int,int);
friend bool FillRestCellProcess2(const vector<CykMatrixCell> &,CykMatrixCell &,int,int);
friend int main(int,char**);
friend bool IsPossibleTree(const vector<CykMatrixCell> &);
friend string BuildTrace(const vector<CykMatrixCell> &);
friend void TraceFurther(ElementCoordinates,const vector<CykMatrixCell> &,ostringstream &);

private : unsigned long UpListCount;
		  unsigned long DownListCount;
		  ElementList UpList;
		  ElementListItr UpListItr;
		  ElementList DownList;
		  ElementListItr DownListItr;

public	: CykMatrixCell() : UpListCount(0), DownListCount(0), UpListItr(UpList.begin()), DownListItr(DownList.begin()) {}
		  
};

#endif

