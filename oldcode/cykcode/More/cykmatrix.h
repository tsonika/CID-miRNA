#ifndef _CYK_MATRIX__H_
#define _CYK_MATRIX__H_


#include "allincludes.h"


enum TYPE_OF_LIST { _UNKNOWN = 0, _UPLIST = 1, _DOWNLIST = 2, _PARENTREF =3, _TERMRULE =4 };


//**********************************************************************************
//This class simply defines the coordinates of an element in the matrix cell's list
//**********************************************************************************

struct ElementCoordinates {
	
	int row,col,index;
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
friend int main();
friend bool IsPossibleTree(const vector<CykMatrixCell> &);
friend string BuildTrace(const vector<CykMatrixCell> &);
friend void TraceFurther(ElementCoordinates,const vector<CykMatrixCell> &,ostringstream &);
friend string make_parent(CykListElement temp);

private : string ruleLHS, ruleRHSpreDot, ruleRHSpostDot;
		  ElementCoordinates coordinates;
                  string Parents;
				  int s1,s2,s3;
				  int flag;
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
friend int main();
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

class cartesian1
{
	private :	string str;
    ElementCoordinates parent3;
	string parent1,parent2;
	double scoretemp;
friend	int chek(vector<cartesian1> ss,cartesian1 strtemp);
friend void FillRestCellProcess1(const vector<CykMatrixCell> &mat,CykMatrixCell &cykmc,int i,int j);


	public:
cartesian1()
{};

};

#endif

