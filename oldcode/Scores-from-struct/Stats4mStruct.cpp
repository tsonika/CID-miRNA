//**********************************************************************************
// File		: Stats4mStruct.cpp										
// Purpose	: Reading an miRNA structure file and outputting the characteristics of the structure in terms of some handy structural details	
// Author	: Manish Kushwaha (manish.kushwaha@gmail.com)
// Revision	: 1.0	-	25-Jun-2005	Working  (Manish Kushwaha)			
//**********************************************************************************

char matrix[5][81];//This matrix represents the screen of 5x80 (Height x Width) characters
#include "drawmirna.h"

int main(int argc, char *argv[])
{
//Recieve input from command line
if (argc < 3)
   {
   cout<<"ERROR!! Insufficient arguments.\n";
   cout<<"Usage Format:\n";
   cout<<"<Program Name> <Input File Name with structures> <Output File Name with details>\n";
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
           structstats=getstats(miRNAstruct,true);
           *outfile<<miRNAstruct;//This line prints the original structure to the output file

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
   *outfile<<miRNAstruct;//This line prints the original structure to the output file
                                                                                                                  
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
