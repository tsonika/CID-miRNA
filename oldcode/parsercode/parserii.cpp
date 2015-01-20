#include <iostream>
#include <fstream>
#include <string.h>
#include <stdlib.h>

using namespace std;

const int MAXPATH = 100;
const int MAXLEN = 200;
char verboseflag;

class selector {

private : 
char infile[MAXPATH];
char outfile[MAXPATH];
char window[MAXLEN];
char searchwindow[MAXLEN];

char sequence[MAXLEN];
char header[MAXLEN];
char temp[MAXLEN];
int searchlen;
int seqlen;
int headerlen;
int windowlen;

int searchingat;
int foundat;

long incount;
long outcount;
void complement(char *str);
//void memallocator(void);
void parser(void);
int check(void);
char altcomplement(char ichar);
int mystrstr(char* targetstr, char* searchstr);
void filler(int start);

public:
void initializer(void);

};

void selector::initializer(void)
{
 cout << "\n SECOND STAGE PARSER \n";
 cout << "\n Enter the name of input file in Fasta Format : ";
 cin.getline(infile,MAXPATH-1);
 cout << "\n Enter The Name of Output file : ";
 cin.getline(outfile, MAXPATH-1);
 
 cout << "\n Enter the Window length at beginning : ";
 cin >> windowlen;
 cout << "\n Enter the length to search inside : ";
 cin >> searchlen;
 cout <<"\n Do You want verbose file output (y/n) : ";
 cin >> verboseflag;
 verboseflag = tolower(verboseflag);
 parser();
}


/* void selector::memallocator()
{

	window = new char[seqlen-2*windowlen+1];
	searchwindow = new char[searchlen + 1];

	if(!window || !searchwindow)
	{
		cout << "\n Memory Allocation Error \n";
		exit(1);
	}
	
} */
 
void selector::parser(void)
{
	
	int ctr1, ctr2;
	incount = 0;
	outcount = 0;
	fstream fin;
	fin.open(infile,ios::in);


	fstream fout, foutfa;
	fout.open(outfile,ios::out|ios::trunc);

	strcat(outfile,".fa");
	foutfa.open(outfile,ios::out|ios::trunc);

	if(!fin.good() || !fout.good() || !foutfa.good())
	{
		cout << "\n Error Opening files "<<fin.good() <<" "<< fout.good() <<" " << foutfa.good() <<"\n";
		exit(1);
	}

	while(!fin.eof())
	{
		fin.getline(temp,MAXLEN-1);
		if(!strlen(temp)) continue;

		if(temp[0] == '>')
		{
			++incount;
			strcpy(header,temp);
			fin.getline(temp, MAXLEN-1);
			
			strcpy(sequence,temp);
			seqlen = strlen(sequence);

			//memallocator();

			for(ctr1 = windowlen, ctr2 = 0; ctr1 < (seqlen-windowlen); ++ctr1, ++ctr2)
			  window[ctr2] = sequence[ctr1];
			
			window[ctr2] = '\0';
			
			//cout<<"\n WINDOW : "<<window << " \\ "<<strlen(window);
			//cin.get();
			
			if(check())
			{
				++outcount;
				if(verboseflag == 'y')
				{
				 fout<<sequence<<"\n";
				 foutfa<<header<<" Search :"<<searchingat<<" Found :"<<foundat
					   <<"\n"<<sequence<<"\n";
				}
				 
			}
		}
		else continue;
	}
	fin.close();
	
	
	
	cout << "\n Second stage Parsing Complete !!"
		 << "\n RESULT : Parsed "<<outcount<<" out of " <<incount<< " sequences "
		 << "\n See " << outfile << " For output \n";
	fout.close();
	foutfa.close();

}


char selector::altcomplement(char ichar)
{
 if (ichar == 'A') return 'G';
 else if(ichar == 'C') return 'T';
 return ichar;
}


int selector::mystrstr(char* targetstr, char* searchstr)
{
 int searchlen = strlen(searchstr);
 int targetlen = strlen(targetstr);
 
 if(targetlen < searchlen) return 0;
 
 char* temp = 0;
 short matchflag = 1;
 int ctr, cur = 0;
 
 while(cur <= (targetlen-searchlen))
 {
 matchflag = 1;
 for(ctr = 0; (ctr < searchlen) && matchflag; ++ctr)
 {
 if((targetstr[cur+ctr] != searchstr[ctr]) && (targetstr[cur+ctr] != altcomplement(searchstr[ctr])))
 {
  matchflag = 0;
  break;
 }
 }
 if(!matchflag)
 cur = cur+ctr+1;
 else
 {
  temp = targetstr + cur;
  //cout << "\nRESULT\n " << targetstr << "\n" << searchstr <<"\t" << cur;
  foundat = searchingat + searchlen + cur;
  return 1;
 }
 }
 return 0;
 }
 

void selector::filler(int start)
{
	int i;
	//cout<< "\n Hello Filler";
	//cin.get();
	for( i = 0; (i < searchlen) && (start <= (strlen(window) - searchlen)); ++i)
		searchwindow[i] = window[start+i];

	searchwindow[i] = '\0';
	//cout <<"\n inside filler : "<<start<<searchwindow;
	//cin.get();
	
}

void selector::complement(char *str)
{
	int i;
	for(i = 0; i < strlen(str); ++i)
	{
		if (str[i] == 'A') str[i] = 'U';
		else if (str[i] == 'U') str[i] = 'A';
		else if (str[i] == 'G') str[i] = 'C';
		else if (str[i] == 'C') str[i] = 'G';
		else str[i] = 'X';
	}
}

int selector::check(void)
{
	int wlen = strlen(window);
	//cout << "\n slen : "<<searchlen;
	//cin.get();
	int i;
	char* temp;

	for(i = 0; i < (wlen-2*searchlen); ++i)
	{
		filler(i);
		
	        //cout << "\n SEARCHWINDOW : " << searchwindow << "\n"<<wlen;
		complement(searchwindow);
		//cout<<"\n Comp "<< searchwindow;
		temp = window+i+searchlen;
		//cout<<"\n Temp "<<temp;
		//cin.get();
		searchingat = i+windowlen;

		if(mystrstr(temp, searchwindow))
			return 1;
	}

	return 0;
}

int main()
{
selector s;
s.initializer();
return 0;
}








