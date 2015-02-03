//**********************************************************************************
// File		: SplitFile.cpp										
// Purpose	: To split any fasta type file into two different fastas depending on an index list file
// Author	: Manish Kushwaha (manish.kushwaha@gmail.com)
// Revision	: 1.0	-	27-Jun-2005	Working  (Manish Kushwaha)			
//**********************************************************************************

#include "allincludes.h"

int main(int argc, char *argv[])
{
//Recieve input from command line
if (argc < 5)
   {
   cout<<"ERROR!! Insufficient arguments.\n";
   cout<<"Usage Format:\n";
   cout<<"<Program Name> <Input Fasta-type File Name> <List File for Segregation> <Output File Name with the items listed in List File> <Output File Name with the items NOT listed in List File>\n";
   exit(1);
   }

char *InFileName=argv[1];
char *ListFileName=argv[2];
char *ListOutFileName=argv[3];
char *ListNotOutFileName=argv[4];

//Open input and output files for reading
ifstream * InFile;
InFile= new ifstream;
InFile->open(InFileName, ios::in);
if (!InFile->good())
   {
   cout<<"ERROR!! Could not open 'Input Fasta-type File'.\n\n";
   exit(1);
   }

ifstream * ListFile;
ListFile= new ifstream;
ListFile->open(ListFileName, ios::in);
if (!ListFile->good())
   {
   cout<<"ERROR!! Could not open 'List File for Segregation'.\n\n";
   exit(1);
   }

ofstream * ListOutFile;
ListOutFile= new ofstream;
ListOutFile->open(ListOutFileName, ios::out);
if (!ListOutFile->good())
   {
   cout<<"ERROR!! Could not open 'Output File Name with the items listed in List File'.\n\n";
   exit(1);
   }

ofstream * ListNotOutFile;
ListNotOutFile= new ofstream;
ListNotOutFile->open(ListNotOutFileName, ios::out);
if (!ListNotOutFile->good())
   {
   cout<<"ERROR!! Could not open 'Output File Name with the items NOT listed in List File'.\n\n";
   exit(1);
   }

string ListNames[1000],inpline;
int ListCount=0;

//Read in the titles in the List File
while (!ListFile->eof())
      {
      getline(*ListFile,inpline);
      while (isspace(inpline[0]) || inpline[0]=='>') inpline.erase(0,1);//Remove all leading spaces and '>'
      inpline=inpline.substr(0,inpline.find(" "));//Keep only the part before a space on line
      if (inpline.length()>0) ListNames[ListCount++]=inpline;
      }
ListFile->close();
delete ListFile;

//Run through the Input Fasta-type file
string inpline2,middlecontent="";
bool tracking=false;
bool islist;
ofstream *ToWrite=new ofstream;

while (!InFile->eof())
      {//Run through the Input Fasta-type file
      getline(*InFile,inpline);
      inpline2=inpline;

      while (isspace(inpline2[0])) inpline2.erase(0,1);//Removes all leading spaces
      if (inpline2.find(">")==0)
         {//Title Line
         if (tracking)
            {//End and write middle content
            *ToWrite<<middlecontent;//This line prints the original middle content to the output file
            middlecontent="";
            }//End- End and write middle content
         //Initiate a middle content
         if (middlecontent.length() != 0) *ToWrite<<middlecontent;
         middlecontent="";
         tracking=true;

         //This portion checks where the output must be written
         islist=false;
         for (int temp=0;temp<ListCount;temp++)
             {
             if (inpline.find(ListNames[temp]+" ") != string::npos)
                {
                islist=true;
                break;
                }
             else if (inpline.find(ListNames[temp]) == (inpline.length()-ListNames[temp].length()))
                {
                islist=true;
                break;
                }
             }//End- for
         if (islist) ToWrite=ListOutFile; else ToWrite=ListNotOutFile;
         //End- This portion checks and assigns where the output must be written

         *ToWrite<<inpline<<endl;
         }//End- Title Line
        else
        {
        middlecontent+=inpline+"\n";
        }
      }//End- while

if (middlecontent.length() != 0) *ToWrite<<middlecontent;

InFile->close();
delete InFile;
ListOutFile->close();
delete ListOutFile;
ListNotOutFile->close();
delete ListNotOutFile;

return 0;
}
