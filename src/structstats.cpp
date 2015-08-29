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
   cout<<"<Program Name> <Grammar File Name> <Sequence to Parse>\n";
   exit(1);
   }

char *GFileName=argv[1];
Sequence = argv[2];

//Read Grammar from file into the vector RuleList
ifstream * gfile;
gfile= new ifstream;
gfile->open(GFileName, ios::in);
RuleList=ReadGrammarFile(gfile);
gfile->close();
delete gfile;

//The following portion of the program displays the grammer as read from the file
//A sample output can be seen in RuleReference.txt (Grammar file used was 'result.txt')
/*
cout<<"Printing the values of various variables in RuleList...\n\n";

cout<<"STATIC MEMBERS:\n";
cout<<"NTcount: "<<RuleList[0].NTcount<<endl;
cout<<"Tcount: "<<RuleList[0].Tcount<<endl;
cout<<"RuleCount: "<<RuleList[0].RuleCount<<endl;
cout<<"NT: "<<RuleList[0].NT<<endl;
cout<<"T: "<<RuleList[0].T<<endl;
cout<<"NTonLHS: "<<RuleList[0].NTonLHS<<endl<<endl;
void FillFirstRow(vector<CykMatrixCell> &mat)//This function accepts a reference to a vector of type CykMatrixCell
cout<<"VARYING MEMBERS:\n";
for (int x=0;x<RuleList[0].RuleCount;x++)
    {
    cout<<"RuleNumber: "<<RuleList[x].RuleNumber<<endl;
    cout<<"RuleScore: "<<RuleList[x].RuleScore<<endl;
    cout<<"Rule:\truleLHS --> ruleRHS\n"<<RuleList[x].ruleLHS<<" --> "<<RuleList[x].ruleRHS<<endl;
    cout<<"UsedOnceFlag: "<<RuleList[x].UsedOnceFlag<<endl;
    cout<<"No_of_RHSsymbols: "<<RuleList[x].No_of_RHSsymbols<<endl;
    cout<<"RuleType: "<<RuleList[x].RuleType<<"\n\n\n";
    }//*/


//Make the matrix as a vector
SeqLen=Sequence.length();
vector<CykMatrixCell> UTMatrix;
FillFirstRow(UTMatrix);
FillRestMatrix(UTMatrix);

//The following portion of the program displays the results of the filled Matrix
//A sample output can be seen in MatrixRef.txt (Grammar file used: smade02.txt, sequence 'aauggau')
/*
int x,y,abspos,itr;
for (abspos=0;abspos < UTMatrix.size(); abspos++)
    {
    XYPos(abspos,x,y);

    cout<<"\n\nCell Position: "<<x<<","<<y<<" = Absolute Position "<<abspos<<endl;
    cout<<"Sub String Represented: "<<Sequence.substr(y-1,x)<<"\n\n";

    cout<<"\tUpListCount :"<<UTMatrix[abspos].UpListCount<<endl;
    cout<<"\tDownListCount :"<<UTMatrix[abspos].DownListCount<<endl;

    cout<<"\tUpList:\n";
    for (itr=0;itr<UTMatrix[abspos].UpList.size();itr++)
        {
        cout<<"\t\truleLHS: "<<UTMatrix[abspos].UpList[itr].ruleLHS<<endl;
        cout<<"\t\truleRHSpreDot: "<<UTMatrix[abspos].UpList[itr].ruleRHSpreDot<<endl;
        cout<<"\t\truleRHSpostDot: "<<UTMatrix[abspos].UpList[itr].ruleRHSpostDot<<endl;
        cout<<"\t\tscore: "<<UTMatrix[abspos].UpList[itr].score<<endl;
        cout<<"\t\tcoordinates:\n";
        cout<<"\t\t\trow: "<<UTMatrix[abspos].UpList[itr].coordinates.row<<endl;
        cout<<"\t\t\tcol: "<<UTMatrix[abspos].UpList[itr].coordinates.col<<endl;
        cout<<"\t\t\tindex: "<<UTMatrix[abspos].UpList[itr].coordinates.index<<endl;
        cout<<"\t\t\ttype: "<<UTMatrix[abspos].UpList[itr].coordinates.type<<endl;
        cout<<"\t\tParents: "<<UTMatrix[abspos].UpList[itr].Parents<<"\n\n";
        }

    cout<<"\tDownList:\n";
    for (itr=0;itr<UTMatrix[abspos].DownList.size();itr++)
        {
        cout<<"\t\truleLHS: "<<UTMatrix[abspos].DownList[itr].ruleLHS<<endl;
        cout<<"\t\truleRHSpreDot: "<<UTMatrix[abspos].DownList[itr].ruleRHSpreDot<<endl;
        cout<<"\t\truleRHSpostDot: "<<UTMatrix[abspos].DownList[itr].ruleRHSpostDot<<endl;
        cout<<"\t\tscore: "<<UTMatrix[abspos].DownList[itr].score<<endl;
        cout<<"\t\tcoordinates:\n";
        cout<<"\t\t\trow: "<<UTMatrix[abspos].DownList[itr].coordinates.row<<endl;
        cout<<"\t\t\tcol: "<<UTMatrix[abspos].DownList[itr].coordinates.col<<endl;
        cout<<"\t\t\tindex: "<<UTMatrix[abspos].DownList[itr].coordinates.index<<endl;
        cout<<"\t\t\ttype: "<<UTMatrix[abspos].DownList[itr].coordinates.type<<endl;
        cout<<"\t\tParents: "<<UTMatrix[abspos].DownList[itr].Parents<<"\n\n";

        }
    }//End- Printing of Filled Matrix*/

//Checking if the input sequence can be parsed using the grammar, and printing the tree if it can be parsed through the grammar
if (IsPossibleTree(UTMatrix))
   {
   cout<<"\nThe sequence CAN be parsed using the grammar!!\n\n";
   string trace;
   trace=BuildTrace(UTMatrix);//'trace' contains the output stream for the BuildTrace(); This is the same as the output seen on screen in the previous version
   cout<<trace;

   node *RootNode= new node;//Create a new node
   RootNode = BuildNodes(trace);//Calling local function BuildNodes() to build the nodes that can be recognised by the fillmatrix() code of 'drawmirna.h'

   initmatrix();
   fillmatrix(RootNode,0);

   string miRNAstruct;
   miRNAstruct=returnmatrix();
   cout<<miRNAstruct;
                                                                                                                  
   miRNAstats structstats;
   structstats=getstats(miRNAstruct);
                                                                                                                  
   cout<<"Stats for the above structure:\n";
   cout<<"Total number of Continuous Stems :"<<structstats.n_stems;
   cout<<"\nTotal number of Asymmetric Bulges :"<<structstats.n_asyms;
   cout<<"\nTotal number of Symmetric Bulges :"<<structstats.n_syms;
   cout<<"\nTotal number of Bases :"<<structstats.n_bases;
   cout<<"\nTotal number of Base Pairs :"<<structstats.n_base_pairs;
   cout<<"\nTotal number of Bases in Loop :"<<structstats.n_loop_bases<<"\n\n";



//   displaymatrix();
   }
   else
   cout<<"\n\nThe sequence CANNOT be parsed using the grammar!!\n\n";

return 0;
}
