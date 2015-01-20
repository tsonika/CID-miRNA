//**********************************************************************************
// File		: structstats.cpp										
// Purpose	: Definition of functions for reading and storing of MiRNA General SCFG	
// Authors	: Vipin Gupta (viping@gmail.com), Manish Kushwaha (manish.kushwaha@gmail.com)
// Revision	:1.0	-	26-Feb-2005	Initial Draft (Vipin Gupta)			
//		:2.0	-	10-Mar-2005	Grammar Reading and Error Handling perfected (Vipin Gupta)
//		:3.0	-	03-Jun-2005	Minor Upgrade in Grammer Reading and Validation functions,
//						Filling the first row of a CYK Upper Triangular Matrix (Manish Kushwaha)
//		:4.0	-	14-Jun-2005	Filling the rest of the rows of the CYK Upper Triangular Matrix and Tracing it to print a human-readable tree (Manish Kushwaha)
//		:5.0	-	15-Jun-2005	Included 'drawmirna.h' functionality to display graphical output
//		:5.5	-	25-Jun-2005	Removed the bias in asymmetric bulge display
//              :6.0    -       25-Jun-2005     Return statistics of a given structure
//**********************************************************************************

#include "grammar_reader.h"
#include "cykmatrix.h"
char matrix[5][81];//This matrix represents the screen of 5x80 (Height x Width) characters
#include "drawmirna.h"


//********************** The Code for Grammer Reading and Validation Begins ********
unsigned long SCFGrule::NTcount = 0;
unsigned long SCFGrule::Tcount = 0;
unsigned long SCFGrule::RuleCount = 0;
string SCFGrule::NT;
string SCFGrule::T;
string SCFGrule::NTonLHS;

//Function to initialize rule
void SCFGrule::InitializeRule(double & score, string & left, string & right)
{
  RuleScore		= score;
  ruleLHS		= left;
  ruleRHS		= right;
  RuleNumber		= RuleCount+1;
  RuleType		= _UNKNOWN_RULE;
  UsedOnceFlag		= false;

  No_of_RHSsymbols = 0;
}

//Function to set OnceUsedFlag if rule has been used at least once
void SCFGrule::SetUsedOnceFlag(void)
{
	UsedOnceFlag = true;
}

// Function for rule checking the grammar and update parameters has to be run for each object, optional
//optional
void SCFGrule::CountUpdate(void)
{
	string dummy, newRHS;
	string star = "*";
	dummy = star + ruleLHS + star;
	//Updating each terminal or non-terminal by enclosing it in *X*
	

	if(!isupper(ruleLHS[0]))
	{
		cout << endl << "Lowercase Symbol on LHS : " << ruleLHS 
			 << " at Rule No. " << RuleNumber << endl;
		myexit(1);
	}

    ruleLHS = dummy;

	if(NT.empty())
	{	
		NT = dummy; 
		++NTcount; 
		NTonLHS = NT;
	} 
	else if(NT.find(dummy) == string::npos) //dummy not found
	{ 
		NT += dummy; 
		NTonLHS += dummy; 
		++NTcount; 
	}
	else if(NTonLHS.find(dummy) == string::npos) 
	{
		NTonLHS += dummy; 
	}

	istringstream iss(ruleRHS);
	string temp;
    
	int count = 0;
	//This variable keeps track of how many symbols read from RHS
     

	while (iss >> temp)
	{
		++count;
		if(!isupper(temp[0]))
		{//If the rule does not begin with upper case
			if((count != 1) || (!isalpha(temp[0])))
			{
				cout << endl << "Lowercase Symbol on RHS : " << ruleRHS 
					 << " at Rule No. " << RuleNumber << endl;
				myexit(1);
			}
			else 
			{
				RuleType = _TERM_RULE;
				dummy = star + temp + star;
				if(T.empty())
				{	T = dummy; ++Tcount; } 
				else if((T.find(dummy) < 0) || (T.find(dummy) >= T.length()))
				{ T += dummy; ++Tcount; }

				newRHS += dummy;

			}
			break;

		}

		else
		{//If the rule begins with upper case= it is a non-terminal rule
			dummy = star + temp + star;
		
			//int x = NT.find(dummy);
			RuleType = _NONTERM_RULE;

			if(NT.find(dummy) == string::npos)
			{NT += dummy;
                        ++NTcount;//To confirm from Vipin- Was the updation of NTcount deliberately avoided here?
                        }
			newRHS += dummy;
		}
		temp.erase(temp.begin(),temp.end());
	}
	ruleRHS = newRHS;

	No_of_RHSsymbols = count;
	++RuleCount;
}

// validation of Grammar has to be run once (optional)

bool SCFGrule::ValidateGrammar(void)
{
	//validate Terminal strings
	for(int i = 0; i < T.length() ; ++i)
	{
		if((T[i] == '*') || (islower(T[i]))) continue;
		else
		{
			cout << endl <<"Some Error in Terminal string : " << T.substr(i-1,3) << endl << T;
			return false;
		}
	}	                              

	//validation of Non Terminal String
	//Checking if all NTs occur on LHS
/*
	if(NT.length() != NTonLHS.length())
	{
		cout << endl << "NT and NTonLHS are different in size" << endl
			 << "NT      : " << NT << endl
			 << "NTonLHS : " << NTonLHS << endl;
		return false;
	}
*/
	//char* temp = new char[NT.length() + 1];
	
	string sortedNT = NT;
	string sortedNTonLHS = NTonLHS;

	sort(sortedNT.begin(),sortedNT.end());
	sort(sortedNTonLHS.begin(),sortedNTonLHS.end());

	if(sortedNT != sortedNTonLHS)
	{
		cout << endl << "NT and NTonLHS are different in nature" << endl
			 << "NT      : " << NT << endl
			 << "NTonLHS : " << NTonLHS << endl
			 << "T       : " << T << endl;
		return false;
	}

	return true;
}


// This function reads the grammar and assumes the grammar according to following format
// Each Rule on a separate line : <probability> <Rule LHS> --> <Rule RHS>
// Each LHS can have a single NT name with no space as space is a delimiter for NTs
// RHS should have either a T or string of NTs separated by a space
// Also small letters are assumed to be Terminals (T) and Capital ones as Non-Terminals (NT)

vector<SCFGrule> ReadGrammarFile(ifstream * grammarfileptr)
{
vector<SCFGrule> TempRuleList;

	if(grammarfileptr == NULL)
	{
		grammarfileptr = new ifstream;
		if(grammarfileptr == NULL)
		{
			cout << endl << "Unable to Read or Allocate memory for reading Grammar File" << endl;
			myexit(1);
		}
		char filename[MAXPATH+1];
		cout << endl << "\n Enter Name of Grammar File : ";
		cin.getline(filename,MAXPATH);
	        grammarfileptr->open(filename, ios::in);
		if(!grammarfileptr->good())
		{
			cout << endl << "Unable to Read Grammar File : " << filename << endl;
			myexit(1);
		}
	}

	// Start Parsing the Grammar file

    	SCFGrule tempruleobj;

	while(!grammarfileptr->eof())
	{
		char temprule[MAXPATH];
		grammarfileptr->getline(temprule,MAXPATH-1);
		if(!strlen(temprule)) continue;
		
		double score;
		istringstream iss(temprule);
		string rule;

        if((iss >> score) && (score >= 0.0) && (score <= 1.0));
		else
		{
			cout << "\n Error Reading score in rule !!\n"
				 << temprule;
			myexit(1);
		}

		rule = iss.str();
		int tempitr = -1;
		while(!isspace(rule[++tempitr]));
		rule.erase(0,tempitr);
                tempitr = -1;
		while(!isalpha(rule[++tempitr]));
		rule.erase(0,tempitr);
		//this is done because to remove the initial prob value from input line
		int delimiter = rule.find("-->");

		if((delimiter < 0) || (delimiter >= rule.length()))
		{
			cout << "\n Error Reading --> in rule !!\n"
				 << temprule;
			myexit(1);
		}
		
		string tempLHS, tempRHS;
		tempLHS = rule.substr(0,delimiter);
		tempRHS = rule.substr(delimiter+3,rule.length()-delimiter-3);

		//PRUNING STARTING AND TRAILING WHITESPACES
		long spacecount = -1;
		while(!isalpha(tempLHS[++spacecount]));
		tempLHS.erase(0,spacecount);

		spacecount = tempLHS.length();
		while(!isalnum(tempLHS[--spacecount])); 
		tempLHS.erase(spacecount+1, tempLHS.length() - spacecount -1);

		spacecount = -1;
		while(!isalpha(tempRHS[++spacecount]));
		tempRHS.erase(0,spacecount);

		spacecount = tempRHS.length();
		while(!isalnum(tempRHS[--spacecount])); 
		tempRHS.erase(spacecount+1, tempRHS.length() - spacecount -1);
		//PRUNING ENDS

		tempruleobj.InitializeRule(score, tempLHS, tempRHS);
		tempruleobj.CountUpdate();
				
		TempRuleList.push_back(tempruleobj);
	}

	//IF OPTIONAL ERRORCHECK REQUIRED PUT HERE
	if(!SCFGrule::ValidateGrammar())
	{
		cout << "\n Error in Grammar.\n";
		myexit(1);
	}
return (TempRuleList);
}
//********************** The Code for Grammer Reading and Validation Ends **********


//************* Definition of operator overloading for CykListElement Begins *******
//The == and != operators overloaded for CykListElement keeping in view the frequent equality checks required

bool CykListElement::operator==(const CykListElement& InCle) const
{
if (this->ruleLHS != InCle.ruleLHS) return false;
if (this->ruleRHSpreDot != InCle.ruleRHSpreDot) return false;
if (this->ruleRHSpostDot != InCle.ruleRHSpostDot) return false;
if (this->score != InCle.score) return false;
if (this->coordinates.row != InCle.coordinates.row) return false;
if (this->coordinates.col != InCle.coordinates.col) return false;
if (this->coordinates.type != InCle.coordinates.type) return false;
if (this->Parents != InCle.Parents) return false;

return true;
}

bool CykListElement::operator!=(const CykListElement& InCle) const
{
if (*this == InCle) return false; else return true;
}
//************** Definition of operator overloading for CykListElement Ends ********

//********************** The Code for First Row Matrix Filling Begins **************
int SeqLen;//Global Containers
string Sequence;
vector<SCFGrule> RuleList;

int AbsPos(int i, int j)//This function returns the absolute position of the vector CykCell with coordinates i,j
{
int r,c,abspos=0;
for (r = 1;r <= SeqLen;r++)
    for (c = 1;c <= (SeqLen-r+1);c++)
        {
        ++abspos;
        if (r==i && c==j)
           return (abspos-1);
        }
}

void XYPos(int pos, int &i, int &j)//This function takes the absolute position of the vector CykCell and feeds the XY coordinates into i,j
{
int r,c,abspos=0;
for (r = 1;r <= SeqLen;r++)
    for (c = 1;c <= (SeqLen-r+1);c++)
        {
        ++abspos;
        if ((abspos-1) == pos)
           {
           i=r;
           j=c;
           return;
           }
        }
}

bool FillFRCellDList(CykMatrixCell &matcel, int pos)//This function fills the Lower List of a CykCell depending on the members of the Upper List
{
string tmpNT;
bool found;
bool NeedRepeat=false;
int itr1,itr2;

for (itr1=0; itr1< matcel.UpList.size(); itr1++)//Run through the Upper List of the cell
    {
    tmpNT=matcel.UpList[itr1].ruleLHS;//This is one such non terminal that gives the terminal represented by the current cell

    for (int x=0;x < RuleList[0].RuleCount;x++)//Run through each rule
        {
        if (RuleList[x].RuleType == _NONTERM_RULE && RuleList[x].ruleRHS == tmpNT)//If the non-terminal in the upperlist is directly given by a non-terminal in the RuleList
           {
           CykListElement tempcle;
           tempcle.ruleLHS=RuleList[x].ruleLHS;
           tempcle.ruleRHSpreDot=tmpNT;
           tempcle.ruleRHSpostDot="\0";//Since this is a terminal rule, RHS not required to be filled
           tempcle.score=RuleList[x].RuleScore;
           RuleList[x].SetUsedOnceFlag();
           tempcle.coordinates.row=1;
           tempcle.coordinates.col=pos;
           tempcle.coordinates.index=matcel.UpList.size();
           tempcle.coordinates.type=_UPLIST;
           ostringstream oss;
           oss<<"*"<<matcel.UpList[itr1].coordinates.row<<" "<<matcel.UpList[itr1].coordinates.col<<" "<<itr1<<"*";
           tempcle.Parents=oss.str();

           found=false;
           for (itr2=0; itr2< matcel.UpList.size(); itr2++)//This loop checks if the CykList Element is not already present
               {if (matcel.UpList[itr2] == tempcle) {found=true;break;} else found=false;}

           if (!found)
              {
              matcel.UpList.push_back(tempcle);
              ++(matcel.UpListCount);
              NeedRepeat=true;
              }
           }
           else {
           if(RuleList[x].RuleType == _NONTERM_RULE && RuleList[x].ruleRHS.find(tmpNT) == 0)//The non-terminal in the upperlist is NOT directly given by a non-terminal in the RuleList, but is the first non-terminal on RHS
           {
           CykListElement tempcle;
           ElementCoordinates tempec;
           tempcle.ruleLHS=RuleList[x].ruleLHS;
           tempcle.ruleRHSpreDot=tmpNT;
           tempcle.ruleRHSpostDot=RuleList[x].ruleRHS.substr(tmpNT.length());

           tempcle.score=RuleList[x].RuleScore;
           RuleList[x].SetUsedOnceFlag();
           tempcle.coordinates.row=1;
           tempcle.coordinates.col=pos;
           tempcle.coordinates.index=matcel.DownList.size();
           tempcle.coordinates.type=_DOWNLIST;
           ostringstream oss;
           oss<<"*"<<matcel.UpList[itr1].coordinates.row<<" "<<matcel.UpList[itr1].coordinates.col<<" "<<itr1<<" &*";
           tempcle.Parents=oss.str();

           found=false;
           for (itr2=0; itr2< matcel.DownList.size(); itr2++)//This loop checks if the CykList Element is not already present
               {if (matcel.DownList[itr2] == tempcle) {found=true;break;} else found=false;}

	   if (!found)//If the CykList Element is not already prsent
              {
              matcel.DownList.push_back(tempcle);
              ++(matcel.DownListCount);
              }
           }}

        }//End- Run thru each rule
    }//End- Run thru the Uplist of the cell
return (NeedRepeat);
}

void FillFRCellUList(CykMatrixCell &matcel,char term, int pos)//This function fills the UpperList of a CykCell
{
string dummy;
string strterm="x";
strterm[0]=term;
dummy = "*" + strterm + "*";

bool found;
int itr;

for (int x=0;x < RuleList[0].RuleCount;x++)//Run through each rule
    {
    if (RuleList[x].RuleType == _TERM_RULE && RuleList[x].ruleRHS == dummy)//If this is a terminal rule and the terminal is the same as our input terminal character (because terminal rules will always have 1 RHS symbol, a string search can be avoided.)
       {
       CykListElement tempcle;
       tempcle.ruleLHS=RuleList[x].ruleLHS;
       tempcle.ruleRHSpreDot=dummy;
       tempcle.ruleRHSpostDot="\0";//Since this is a terminal rule, RHSpostDot not required to be filled
       tempcle.score=RuleList[x].RuleScore;
       RuleList[x].SetUsedOnceFlag();
       tempcle.coordinates.row=1;
       tempcle.coordinates.col=pos;
       tempcle.coordinates.index=matcel.UpList.size();
       tempcle.coordinates.type=_UPLIST;
       ostringstream oss;
       oss<<"*1 "<<pos<<" "<<tempcle.coordinates.index<<" &*";
       tempcle.Parents=oss.str();

       found=false;
       for (itr=0; itr< matcel.UpList.size(); itr++)//This loop checks if the CykList Element is not already present
           {if (matcel.UpList[itr] == tempcle) {found=true;break;} else found=false;}

       if (!found)
          {
          matcel.UpList.push_back(tempcle);
          ++(matcel.UpListCount);
          }
       }
    }//End - Run through each rule
}

CykMatrixCell FillFRCell(char term, int pos)//Fill a Cell of the first row
{//This function accepts one terminal character, and the position of the cell and returns a filled in CykMatrixCell
CykMatrixCell cykmatcel;

FillFRCellUList(cykmatcel, term, pos);//To fill the upper list of the CykCell
while (FillFRCellDList(cykmatcel, pos));//To fill the lower list of the CykCell
return (cykmatcel);
}

void FillFirstRow(vector<CykMatrixCell> &mat)//This function accepts a reference to a vector of type CykMatrixCell
{
//The filling of the first row is dependent on the input string for validation
CykMatrixCell tempcykcel;

for (int x=1;x <= SeqLen;x++)
    {
    tempcykcel=FillFRCell(Sequence[x-1],x);
    mat.push_back(tempcykcel);
    }
}
//********************** The Code for First Row Matrix Filling Ends *****************


//*************** The Code for Filling Rest of the Matrix Begins ********************
bool FillRestCellProcess2(const vector<CykMatrixCell> &mat,CykMatrixCell &cykmc,int i,int j)//This function runs the defined Process-II for filling  the non-first-row matrix
{
string dummy;
bool found;
int itr1,itr2;
bool NeedRepeat=false;

for (itr1=0; itr1< cykmc.UpList.size(); itr1++)
    {//Run through the UpperList of cell i,j
    dummy=cykmc.UpList[itr1].ruleLHS;

    for (int x=0;x < RuleList[0].RuleCount;x++)//Run through each rule
        {
        if (RuleList[x].RuleType == _NONTERM_RULE && RuleList[x].ruleRHS == dummy)//If the non-terminal we are looking for is directly given by a non-terminal in the RuleList
           {
           CykListElement tempcle;
           tempcle.ruleLHS=RuleList[x].ruleLHS;
           tempcle.ruleRHSpreDot=dummy;
           tempcle.ruleRHSpostDot="\0";//Since this non-terminal is given directly

           tempcle.score=RuleList[x].RuleScore * cykmc.UpList[itr1].score;
           RuleList[x].SetUsedOnceFlag();
           tempcle.coordinates.row=i;
           tempcle.coordinates.col=j;
           tempcle.coordinates.index=cykmc.UpList.size();
           tempcle.coordinates.type=_UPLIST;
           ostringstream oss;
           oss<<"*"<<cykmc.UpList[itr1].coordinates.row<<" "<<cykmc.UpList[itr1].coordinates.col<<" "<<itr1<<"*";
           tempcle.Parents=oss.str();

           found=false;
           for (itr2=0; itr2< cykmc.UpList.size(); itr2++)//This loop checks if the CykList Element is not already present
               {if (cykmc.UpList[itr2] == tempcle) {found=true;break;} else found=false;}

      	   if (!found)//If the CykList Element is not already prsent
              {
              cykmc.UpList.push_back(tempcle);
              ++(cykmc.UpListCount);

              NeedRepeat = true;//Since the upper list has been updated the process needs to be repeated
              }
           }
        else if (RuleList[x].RuleType == _NONTERM_RULE && RuleList[x].ruleRHS.find(dummy) == 0)//The non-terminal we are looking for is NOT directly given by a non-terminal in the RuleList, but is the first non-terminal on RHS
           {
           CykListElement tempcle;
           tempcle.ruleLHS=RuleList[x].ruleLHS;
           tempcle.ruleRHSpreDot=dummy;
           tempcle.ruleRHSpostDot=RuleList[x].ruleRHS.substr(dummy.length());//Since this non-terminal is NOT given directly, we extract the post-dot portion

           tempcle.score=RuleList[x].RuleScore * cykmc.UpList[itr1].score;
           RuleList[x].SetUsedOnceFlag();
           tempcle.coordinates.row=i;
           tempcle.coordinates.col=j;
           tempcle.coordinates.index=cykmc.DownList.size();
           tempcle.coordinates.type=_DOWNLIST;
           ostringstream oss;
           oss<<"*"<<cykmc.UpList[itr1].coordinates.row<<" "<<cykmc.UpList[itr1].coordinates.col<<" "<<itr1<<"*";
           tempcle.Parents=oss.str();

           found=false;
           for (itr2=0; itr2< cykmc.UpList.size(); itr2++)//This loop checks if the CykList Element is not already present
               {if (cykmc.UpList[itr2] == tempcle) {found=true;break;} else found=false;}

      	   if (!found)//If the CykList Element is not already prsent
              {
              cykmc.DownList.push_back(tempcle);
              ++(cykmc.DownListCount);
              }

           }
        }//End- Run through each rule
    }//End- Run through the UpperList of cell i,j
return (NeedRepeat);
}//FillRestCellProcess2- function ends

void FillRestCellProcess1(const vector<CykMatrixCell> &mat,CykMatrixCell &cykmc,int i,int j)//This function runs the defined Process-I for filling  the non-first-row matrix
{
int maxk=i-1;
int DwnListPos,UprListPos;
string dummy;
bool found;
int dwnitr,upitr,itr2;

for (int k=1;k<=maxk;k++)
    {	
    DwnListPos=AbsPos(k,j);
    UprListPos=AbsPos(i-k,j+k);

    for (dwnitr=0;dwnitr < mat[DwnListPos].DownList.size();dwnitr++)//Run through the entire DownList of the cell at DwnListPos
        {
        for (upitr=0;upitr < mat[UprListPos].UpList.size();upitr++)//Run through the entire UpList of the cell at UprListPos
            {
            dummy=mat[DwnListPos].DownList[dwnitr].ruleRHSpreDot + mat[UprListPos].UpList[upitr].ruleLHS;
            //Using the specific Upper List element and Lower List element we now fill the current cell

            for (int x=0;x < RuleList[0].RuleCount;x++)//Run through each rule
                {
                if (RuleList[x].RuleType == _NONTERM_RULE && RuleList[x].ruleRHS == dummy)//If the non-terminal set we are looking for is directly given by a non-terminal in the RuleList
                   {
                   CykListElement tempcle;
                   tempcle.ruleLHS=RuleList[x].ruleLHS;
                   tempcle.ruleRHSpreDot=dummy;
                   tempcle.ruleRHSpostDot="\0";//Since this non-terminal set is given directly

                   tempcle.score=RuleList[x].RuleScore * mat[UprListPos].UpList[upitr].score;
                   RuleList[x].SetUsedOnceFlag();
                   tempcle.coordinates.row=i;
                   tempcle.coordinates.col=j;
                   tempcle.coordinates.index=cykmc.UpList.size();
                   tempcle.coordinates.type=_UPLIST;
                   tempcle.Parents=mat[DwnListPos].DownList[dwnitr].Parents;
                   ostringstream oss;
                   oss<<"*"<<mat[UprListPos].UpList[upitr].coordinates.row<<" "<<mat[UprListPos].UpList[upitr].coordinates.col<<" "<<upitr;
                   if (mat[UprListPos].UpList[upitr].coordinates.row == 1) oss<<"&";
                   oss<<"*";
                   tempcle.Parents+=oss.str();

                   found=false;
                   for (itr2=0;itr2<cykmc.UpList.size();itr2++)//This loop checks if the CykList Element is not already present
                       {if (cykmc.UpList[itr2] == tempcle) {found=true;break;} else found=false;}

        	   if (!found)//If the CykList Element is not already prsent
                      {
                      cykmc.UpList.push_back(tempcle);
                      ++(cykmc.UpListCount);
                      }
                   }
                else if (RuleList[x].RuleType == _NONTERM_RULE && RuleList[x].ruleRHS.find(dummy) == 0)//The non-terminal set we are looking for is NOT directly given by a non-terminal in the RuleList, but is the first non-terminal set on RHS
                   {
                   CykListElement tempcle;
                   tempcle.ruleLHS=RuleList[x].ruleLHS;
                   tempcle.ruleRHSpreDot=dummy;
                   tempcle.ruleRHSpostDot=RuleList[x].ruleRHS.substr(dummy.length());//Since this non-terminal set is NOT given directly

                   tempcle.score=RuleList[x].RuleScore * mat[UprListPos].UpList[upitr].score;
                   RuleList[x].SetUsedOnceFlag();
                   tempcle.coordinates.row=i;
                   tempcle.coordinates.col=j;
                   tempcle.coordinates.index=cykmc.DownList.size();
                   tempcle.coordinates.type=_DOWNLIST;
                   tempcle.Parents=mat[DwnListPos].DownList[dwnitr].Parents;
                   ostringstream oss;
                   oss<<"*"<<mat[UprListPos].UpList[upitr].coordinates.row<<" "<<mat[UprListPos].UpList[upitr].coordinates.col<<" "<<upitr;
                   if (mat[UprListPos].UpList[upitr].coordinates.row == 1) oss<<"&";
                   oss<<"*";
                   tempcle.Parents+=oss.str();

                   found=false;
                   for (itr2=0;itr2<cykmc.DownList.size();itr2++)//This loop checks if the CykList Element is not already present
                       {if (cykmc.DownList[itr2] == tempcle) {found=true;break;} else found=false;}

        	   if (!found)//If the CykList Element is not already prsent
                      {
                      cykmc.DownList.push_back(tempcle);
                      ++(cykmc.DownListCount);
                      }
                   }//End- if else
                }//End- Run through each rule
            }//End- Run through the entire UpList of the UprListPos
        }//End- Run through the entire DownList of the DwnListPos

    }//end of loop for 'k'
}//FillRestCellProcess1- function ends

void FillRestCell(const vector<CykMatrixCell> &mat,CykMatrixCell &cykmc, int &row, int &col)//This function accepts one CykMatrixCell and its position as row, column (all by reference) and calls Processes-I and II for filling  the non-first-row matrix
{
FillRestCellProcess1(mat,cykmc,row,col);
while (FillRestCellProcess2(mat,cykmc,row,col));
}//FillRestCell()- function ends

void FillRestMatrix(vector<CykMatrixCell> &mat)//This function accepts a reference to a vector of type CykMatrixCell
{
//The filling of the other rows is dependent on the cells of the rows already filled above

int row=2;//Since row filling now begins at row 2

for (int maxcol=SeqLen-1; maxcol>=1; maxcol--)//The Maximum number of columns for each row will keep falling at every step
    {
    for (int col=1; col<=maxcol; col++)//Keep incrementing the column number
        {
        CykMatrixCell tempcykcel;
        FillRestCell(mat,tempcykcel,row,col);
        mat.push_back(tempcykcel);
        }//col for loop ends
    ++row;
    }//maxcol for loop ends
}//FillRestMatrix()- function ends

//**************** The Code for Filling Rest of the Matrix Ends *********************

//************** Code for checking the Parse Tree Feasibility Begins ****************
bool IsPossibleTree(const vector<CykMatrixCell> &UTM)
{
int CheckPos=UTM.size()-1,itr;

for (itr=0; itr < UTM[CheckPos].UpList.size(); itr++)
    {
    if (UTM[CheckPos].UpList[itr].ruleLHS == "*ROOT*") return true;
    }
return false;
}
//*************** Code for checking the Parse Tree Feasibility Ends *****************

//******************** Code for Printing the Parse Tree Begins **********************
ElementCoordinates GetCoordinates(string inp)//This function splits up the sting-encoded coordinates into the structure 'ElementCoordinates' form
{
ElementCoordinates out;
istringstream iss(inp);
iss>>out.row>>out.col>>out.index>>inp;
if (inp == "&")
   {
   out.type=_TERMRULE;
   }
return out;
}//End- GetCoordinates

void TraceFurther(ElementCoordinates Pos,const vector<CykMatrixCell> &UTM,ostringstream &oss)//This function traces one step of a matrix, prints it to oss, and recalls itself for following the rest of the tree
{
CykListElement temp=UTM[AbsPos(Pos.row,Pos.col)].UpList[Pos.index];

if (temp.ruleLHS.length()>0)
   {
   oss<<temp.ruleLHS.substr(1,temp.ruleLHS.length()-2);//This is to remove the leading and trailing '*' (stars)
   oss<<" --> ";

   string tempParents=temp.Parents;
   while (tempParents.length() > 0)
         {
         string PosPart,tempstr;
         PosPart=tempParents.substr(1,tempParents.find("*",1)-1);

         Pos=GetCoordinates(PosPart);

         if (Pos.type == _TERMRULE)
            tempstr=UTM[AbsPos(Pos.row,Pos.col)].UpList[Pos.index].ruleRHSpreDot;
            else
            tempstr=UTM[AbsPos(Pos.row,Pos.col)].UpList[Pos.index].ruleLHS;
         oss<<tempstr.substr(1,tempstr.length()-2)<<" ";

         if (tempParents.length() > PosPart.length()+2)
            tempParents=tempParents.substr(PosPart.length()+2);
            else
            break;
         }//End- while
   oss<<endl;

   tempParents=temp.Parents;
   while (tempParents.length() > 0)
         {
         string PosPart,tempstr;
         PosPart=tempParents.substr(1,tempParents.find("*",1)-1);
         Pos=GetCoordinates(PosPart);

         if (Pos.type != _TERMRULE)
            TraceFurther(Pos,UTM,oss);

         if (tempParents.length() > PosPart.length()+2)
            tempParents=tempParents.substr(PosPart.length()+2);
            else
            break;
         }//End- while
   }//End- if
}//End- TraceFurther()


string BuildTrace(const vector<CykMatrixCell> &UTM)//This function determines the ROOT with highest score, and uses that UpList to trace back the tree by calling TraceFurther()
{
int CheckPos=UTM.size()-1;//This is the last element in the matrix, the ROOT must lie here
int BestRootPos=-1,itr;

for (itr=0; itr < UTM[CheckPos].UpList.size(); itr++)
    {
    if (UTM[CheckPos].UpList[itr].ruleLHS == "*ROOT*" && BestRootPos == -1) BestRootPos=itr;
    if (UTM[CheckPos].UpList[itr].ruleLHS == "*ROOT*" && UTM[CheckPos].UpList[itr].score > UTM[CheckPos].UpList[BestRootPos].score)
       BestRootPos = itr;//Thus 'BestRootPos' gets the position of the greatest score amongst all the ROOTs
    }

if (BestRootPos == -1)
   {
   cout<<"ERROR tracing the tree!! No ROOT found!!\n";
   exit(1);
   }
   else
   {
   ostringstream oss;
   TraceFurther(UTM[CheckPos].UpList[BestRootPos].coordinates, UTM, oss);
   oss<<endl;
   return oss.str();
   }
}//End- BuildTrace()
//********************* Code for Printing the Parse Tree Ends ***********************

//************************ Code for Building nodes Begins ***************************
node* BuildNodes(string &trace)
{
node *finaltree= new node;
finaltree->initialize(ROOT);
node *currnode=finaltree;
node *a;
istringstream iss(trace);
string currline,onlynode;
int findpos;
bool LoopOccured=false;

while (!iss.eof() && !LoopOccured)
      {
      getline(iss,currline);
      if (currline.length() == 0) continue;
      a=new node;

      onlynode=currline.substr(currline.find("-->")+4);

      findpos=onlynode.find("STEM");
      if (findpos != string::npos)
         {//Make node for STEM
         a->initialize(STEM);
         stem *qstem=new stem;
         qstem->l=onlynode[0];

         findpos=onlynode.find(" ",findpos);
         qstem->r=onlynode[findpos+1];
         a->storedata(qstem);
         onlynode="";
         }//End- STEM

      findpos=onlynode.find("ASYM");
      if (findpos != string::npos)
         {//Make node for ASYM
         a->initialize(ASYM);
         asym *qasym=new asym;
         qasym->l=onlynode.substr(0,findpos);

         findpos=onlynode.find(" ",findpos);
         qasym->r=onlynode.substr(findpos+1);
         
         findpos=(qasym->l).find(" ");
         while (findpos != string::npos)
               {
               (qasym->l).erase(findpos,1);
               findpos=(qasym->l).find(" ");
               }

         findpos=(qasym->r).find(" ");
         while (findpos != string::npos)
               {
               (qasym->r).erase(findpos,1);
               findpos=(qasym->r).find(" ");
               }

         a->storedata(qasym);
         onlynode="";
         }//End- ASYM

      findpos=onlynode.find("SYM");
      if (findpos != string::npos)
         {//Make node for SYM
         a->initialize(SYM);
         sym *qsym=new sym;
         qsym->l=onlynode.substr(0,findpos);

         findpos=onlynode.find(" ",findpos);
         qsym->r=onlynode.substr(findpos+1);
         
         findpos=(qsym->l).find(" ");
         while (findpos != string::npos)
               {
               (qsym->l).erase(findpos,1);
               findpos=(qsym->l).find(" ");
               }

         findpos=(qsym->r).find(" ");
         while (findpos != string::npos)
               {
               (qsym->r).erase(findpos,1);
               findpos=(qsym->r).find(" ");
               }

         a->storedata(qsym);
         onlynode="";
         }//End- SYM

      while (isspace(onlynode[0]))  onlynode.erase(0,1);//Remove all leading spaces

      findpos=onlynode.find("LOOP");
      if (findpos != string::npos && findpos>0)
         {//Make node for LOOP, if this point is like the STEM
         a->initialize(STEM);
         stem *qstem=new stem;
         qstem->l=onlynode[0];
                                                                                                                                                                 
         findpos=onlynode.find(" ",findpos);
         qstem->r=onlynode[findpos+1];
         a->storedata(qstem);

         //This portion fills in this special point into the node tree
         if (currnode == finaltree)
            finaltree->storep(a);
         else
            currnode->storep(a);
         currnode=a;
         a=new node;
         }//End- LOOP, like the STEM

      findpos=onlynode.find("LOOP");
      if (findpos != string::npos)
         {//Make node for LOOP, if this point has the LOOP string
         a->initialize(LOOP);
         loop *qloop=new loop;
         string qloops;
         getline(iss,qloops);
         qloops=qloops.substr(qloops.find("-->")+4);
         qloop->s=qloops;

         findpos=(qloop->s).find(" ");
         while (findpos != string::npos)
               {
               (qloop->s).erase(findpos,1);
               findpos=(qloop->s).find(" ");
               }
         a->storedata(qloop);
         LoopOccured=true;
         }//End- LOOP

      //Now 'a' has the current node
      if (currnode == finaltree)
         finaltree->storep(a);
      else
         currnode->storep(a);

      currnode=a;
      }//End- while()

return finaltree;
}

//************************* Code for Building nodes Ends ****************************
int main(int argc, char *argv[])
{
//Recieve input from command line
if (argc < 3)
   {
   cout<<"ERROR!! Insufficient arguments.\n";
   cout<<"Usage Format:\n";
   cout<<"<Program Name> <Input File Name with structures> <Input File Name with details>\n";
   exit(1);
   }

char *InFileName=argv[1];
char *OutFileName=argv[2];

//Open input and output files for reading
ifstream * infile;
infile= new ifstream;
infile->open(InFileName, ios::in);
if (!infile->good())
   {
   cout<<"ERROR!! Could not open INPUT file.\n\n";
   exit(1);
   }
ofstream * outfile;
outfile= new ofstream;
outfile->open(OutFileName, ios::out);
if (!outfile->good())
   {
   cout<<"ERROR!! Could not open OUTPUT file.\n\n";
   exit(1);
   }

string inpline,inpline2,miRNAstruct;
bool tracking=false;
miRNAstats structstats;

while(!infile->eof())
     {//Run through the entire input file
     getline(*infile,inpline);//Read one line of the file in 'inpline'
     inpline2=inpline;

     while (isspace(inpline2[0])) inpline2.erase(0,1);//Removes all leading spaces

     if (inpline2.find(">")==0)
        {//miRNA name line
        if (tracking)
           {//End and write an miRNA structure
           structstats=getstats(miRNAstruct);
           *outfile<<miRNAstruct;

           *outfile<<"Number of Continuous Stems :"<<structstats.n_stems;
           *outfile<<"\nNumber of Asymmetric Bulges :"<<structstats.n_asyms;
           *outfile<<"\nNumber of Symmetric Bulges :"<<structstats.n_syms;
           *outfile<<"\nNumber of Bases :"<<structstats.n_bases;
           *outfile<<"\nNumber of Base Pairs :"<<structstats.n_base_pairs;
           *outfile<<"\nNumber of Bases in the Terminal Loop :"<<structstats.n_loop_bases<<"\n\n"; 
           }//End- End and write an miRNA structure //Initiate an miRNA structure

        //Initiate an miRNA structure
        miRNAstruct="";
        tracking=true;
        *outfile<<inpline2<<endl;
        }//End- miRNA name line
        else
        {
        miRNAstruct+=inpline+"\n";
        }
     }//End- Run through the entire input file - while

//Write the last miRNA structure
if (miRNAstruct.length() != 0)
   {
   structstats=getstats(miRNAstruct);
   *outfile<<miRNAstruct;
                                                                                                                  
   *outfile<<"Number of Continuous Stems :"<<structstats.n_stems;
   *outfile<<"\nNumber of Asymmetric Bulges :"<<structstats.n_asyms;
   *outfile<<"\nNumber of Symmetric Bulges :"<<structstats.n_syms;
   *outfile<<"\nNumber of Bases :"<<structstats.n_bases;
   *outfile<<"\nNumber of Base Pairs :"<<structstats.n_base_pairs;
   *outfile<<"\nNumber of Bases in the Terminal Loop :"<<structstats.n_loop_bases<<"\n\n"; 
   }       
//End- Write the last miRNA structure


outfile->close();
infile->close();
delete outfile;
delete infile;

return 0;
}
