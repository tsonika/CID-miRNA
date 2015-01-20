#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <cstdlib>
#include <math.h>
 
using namespace std;
 
#define MAXLEN 200
#define MAXPATH 100
 
char infile[MAXPATH],outfile[MAXPATH];
double cscore;
 
void passfile();
 
int main(int argc,char *argv[])
{
if (argc < 3)
        {
        cout<<"\nInsufficient number of arguments!!\n\n";
        cout<<"SYNTAX:\n<Program name> <Input file name> <cutoff score> <Output file name>\n";
        return -1;
        }
 
strcpy(infile,argv[1]);
cscore=atof(argv[2]);
strcpy(outfile,argv[3]);
 
//cout<<cyintercept<<endl<<cslope<<endl;
passfile();
 
return 0;
}
 
void passfile()
{
fstream filein, fileout;
 
//Open the input file for reading
filein.open(infile,ios::in);
if(!filein.good())
        {
        cout << "\n Error opening Input File!!\n";
        exit(1);
        }
 
//Open the output file for writing
fileout.open(outfile,ios::out);
if(!fileout.good())
        {
        cout << "\n Error opening Output File!!\n";
        exit(1);
        }
 
cout<<"Going to segregate files.";
//Run through the input file
char inputline1[MAXLEN],inputline2[MAXLEN],inputline3[MAXLEN];
inputline1[0]=inputline2[0]=inputline3[0]='\0';
char tmp[20];
double seqlen,seqnscore;
while (!filein.eof())
        {
        filein.getline(inputline1,MAXLEN-1);
        if (strstr(inputline1,"Sequence :"))
                {//Thus the first line with the pattern "Sequence :" has been found
                filein.getline(inputline2,MAXLEN-1);
                if (strlen(inputline2))
                        {//Thus the second line with some non empty entry found
                        filein.getline(inputline3,MAXLEN-1);
                        if (strlen(inputline3))
                                {//Thus the third line with some non empty entry found
                                istringstream iss(inputline3);
                                iss >> tmp >> tmp >> seqlen >> tmp >> tmp >> tmp >> seqnscore >> tmp >> tmp >> tmp;
 
                                //Now checking if the score is valid and then writing to outputfile
                                if (seqnscore >= cscore && seqlen >= 60)
                                        {
//                                      cout<<inputline1<<endl<<inputline2<<endl<<inputline3<<endl<<((cslope*seqlen)+cyintercept)<<endl;
                                        fileout<<inputline1<<endl<<inputline2<<endl<<inputline3<<endl;
                                        }
                                }
                        }
                }
        }
//End- Run through the input file
 
cout<<"\nDone segregating the passed files!!\n";
}
 


