//**********************************************************************************
// File		: drawmirna.h
// Purpose	: Declaration and definition of classes and functions for the generation of a graphical output of an miRNA
// Author	: Manish Kushwaha (manish.kushwaha@gmail.com)
// Revision	: 1.0	-	28-May-2005	Initial Draft that used a Stub to randomly generate trees and printed them
//		: 2.0	-	15-Jun-2005	Modification to accept input from the miRNA program as a string stream
//		: 3.0	-	25-Jun-2005	Modification to return the output as a string, rather than a character type matrix, and calculate statistical details about a structure
//**********************************************************************************

#ifndef _DRAW__MIRNA_H_
#define _DRAW__MIRNA_H_

#include "allincludes.h"

enum TYPE_OF_NODE { ROOT = 0, STEM = 1, SYM = 2, ASYM =3, LOOP =4 };

struct stem
       {
       char l,r;
       };//This structure implements a STEM. STEM can have one base on the left and another on the right (characters 'l' and 'r' respectively)

struct sym
       {
       string l,r;
       };//This structure implements a SYM (=symmetrical bulge). SYM can have any number of characters on the right and the left (strings 'l' and 'r' respectively)

struct asym
       {
       string l,r;
       };//This structure implements an ASYM (=asymmetrical bulge). ASYM can have any number of characters on the right and the left (strings 'l' and 'r' respectively)

struct loop
       {
       string s;
       };//This structure implements a LOOP. SYM can have one stretch of string that folds like a cap at the terminus of the loop

struct miRNAstats
       {
       double n_stems,n_asyms,n_syms,n_bases,n_base_pairs,n_loop_bases;
       vector <double> stem_len,asym_len,asym_nbases,sym_len;
       };//This structure includes information about the structure of a given miRNA sequence

class node
      {
      public:
      TYPE_OF_NODE type;   //Type of the current node
      node *p;    //Pointer to the next node
      void *data; //Pointer to data in the current node

      public:
      node(TYPE_OF_NODE);
      void initialize(TYPE_OF_NODE);
      void storep(node*);
      void storedata(void*);
      TYPE_OF_NODE gettype();
      node* getp();
      void* getdata();
      };

node::node(TYPE_OF_NODE t=STEM)//By default, assume that the node is STEM type and initialise it accordingly
{
initialize(t);//Calling the initialiser
}

void node::initialize(TYPE_OF_NODE t)//This function initialises a node
{
p=NULL; //Needs to be defined after the initialisation, unless type is LOOP

if (t==ROOT || t==ASYM || t==SYM || t==LOOP) {type =t;} else {type=STEM;}

switch (t)//Memory is allocated according to the type of node
       {
       case ASYM:
	    if (!(data = new asym))
	       {
	       perror ("Insufficient memory for ASYM\n");
	       exit (1);
	       }
            break;

       case SYM:
	    if (!(data = new sym))
	       {
	       perror ("Insufficient memory for ASYM\n");
	       exit (1);
	       }
            break;
            
       case LOOP:
	    if (!(data = new loop))
	       {
	       perror ("Insufficient memory for ASYM\n");
	       exit (1);
	       }
            break;

       case ROOT:
            data=NULL;
            break;


       case STEM:
	    if (!(data = new stem))
	       {
	       perror ("Insufficient memory for ASYM\n");
	       exit (1);
	       }
       }
}

void node::storep(node *a)//Stores the pointer to the next node
{p=a;}

void node::storedata(void *inp)//Stores the data associated with current node
{data=inp;}

TYPE_OF_NODE node::gettype()//Returns the current node-type
{
return (type);
}

node* node::getp()//Returns the pointer to the next node
{
return (p);
}

void* node::getdata()//Returns the data stored at the current node
{
return (data);
}

void fillmatrix(node*,int);//For Matrix
void initmatrix();//For Matrix
void displaymatrix();//For Matrix
string returnmatrix();//For Matrix
miRNAstats getstats(string,bool);//For Matrix

void fillmatrix(node *inpnode,int column)//This function fills in the appropriate characters at the right place in the screen-matrix, the global 2-d array 'matrix'. 'column' tells the function which column of the matrix to print this node to
{
TYPE_OF_NODE inpnodetype;
inpnodetype=inpnode->gettype();//Get the current node type

//Just declaring some variables required later
int p,q,r;
stem *qstem=new stem;
asym *qasym=new asym;
sym *qsym=new sym;
loop *qloop=new loop;
string newshrtstr,shrtstr,lngstr,qsymr,qsyml;

//string newshrtstr;

switch(inpnodetype)//Filling done depending on the node type
      {
      case ROOT:
           fillmatrix(inpnode->getp(),column);//Nothing needs to be done at the ROOT. Recurse using the next node
           break;

      case STEM:
          qstem=(stem*)(inpnode->getdata());
          matrix[1][column]=qstem->l;
          matrix[2][column]='|';
          matrix[3][column]=qstem->r;
          fillmatrix(inpnode->getp(),column+1);//Recurse using the next node
          break;

     case ASYM:
          qasym=(asym*)(inpnode->getdata());
          p= (qasym->r).length();
          q= (qasym->l).length();
          if (q>p) {r=q;} else {r=p;}//r gets the length of the bigger

          if (q==r) {shrtstr=qasym->r;lngstr=qasym->l;} else {shrtstr=qasym->l;lngstr=qasym->r;p=q;}//So, 'shrtstr' gets the shorter string and 'p' gets the length of the shorter string
          //newshrtstr=new char(r+1);
          newshrtstr.assign(r-p,'-');
          newshrtstr.append(shrtstr);

          if (shrtstr == qasym->r) {qasym->r = newshrtstr;} else {qasym->l = newshrtstr;}
       
          for (q=0;q<r;q++)
              {
              matrix[0][column+q]=qasym->l[q];
              matrix[4][column+q]=qasym->r[q];
              }
          fillmatrix(inpnode->getp(),column+r);//Recurse using the next node
          break;

     case SYM:
          qsym=(sym*)(inpnode->getdata());
          p= (qsym->r).length();
          qsymr=qsym->r;
          qsyml=qsym->l;

          for (q=0;q<p;q++)
              {
              matrix[0][column+q]=qsyml[q];
              matrix[4][column+q]=qsymr[q];
              }
          fillmatrix(inpnode->getp(),column+p);//Recurse using the next node
          break;

     case LOOP:
          qloop=(loop*)(inpnode->getdata());
          p=(qloop->s).length();
          string seq=qloop->s;

          if (p>3)
             {
             q=int(p/2)-1;
             for (r=0;r<q;r++)
                 {
                 matrix[0][column+r]=seq[r];
                 matrix[4][column+r]=seq[p-1-r];
                 }
             column=column+q;
             seq=seq.substr(q,p-(q*2));
             }
   
          p=seq.length();
          switch (p)
                 {
                 case 1:
                      matrix[2][column]=seq[0];
                      break;
                 case 2:
                      matrix[1][column]=seq[0];
                      matrix[3][column]=seq[1];
                      break;
                 case 3:
                      matrix[1][column]=seq[0];
                      matrix[2][column]=seq[1];
                      matrix[3][column]=seq[2];
                 }
          break; 
      }
}

void initmatrix()//This function initialises the 'matrix' with spaces
{
int r,c;
for (r=0;r<5;r++)
    {
    for (c=0;c<55;c++)
        matrix[r][c]=' ';
    matrix[r][c]='\0';
    }
}

void displaymatrix()//This function prints the matrix on screen
{
printf ("%s",matrix[0]);
printf ("\n%s",matrix[1]);
printf ("\n%s",matrix[2]);
printf ("\n%s",matrix[3]);
printf ("\n%s\n\n",matrix[4]);
}

string returnmatrix()//This function returns the matrix into a string
{
ostringstream oss;
oss<<matrix[0]<<"\n"<<matrix[1]<<"\n"<<matrix[2]<<"\n"<<matrix[3]<<"\n"<<matrix[4]<<"\n\n";
return oss.str();
}

miRNAstats getstats(string inp, bool considertruncated=false)
{
//       int n_stems,n_asyms,n_syms,n_bases,n_base_pairs,n_loop_bases
//cout<<inp<<"\n\n";
miRNAstats tempstats;
//Initialise 'tempstats'
tempstats.n_stems=0;
tempstats.n_asyms=0;
tempstats.n_syms=0;
tempstats.n_bases=0;
tempstats.n_base_pairs=0;
tempstats.n_loop_bases=0;

istringstream iss(inp);
string mat[5];
string temp="";

//To initialise all strings in the matrix 'mat' with "" 
int r,c;
for (r=0;r<5;r++) mat[r]="";
r=-1;
int maxlen=0;

//To distribute each non-zero length line into the matrix 'mat'
do
  {
  getline(iss,temp);
  if (temp.length() > 0 && r<4)
     {
     mat[++r]=temp;
     if (mat[r].length() > maxlen) maxlen=mat[r].length();
     }
  } while (!iss.eof());

if (maxlen==0 || mat[4].length()==0)
   {
   cout<<"ERROR!! Invalid or incomplete input string.\n";
   exit(1);
   }

/*cout<<"The string as matrix\n\n";
cout<<"0:"<<mat[0]<<endl;
cout<<"1:"<<mat[1]<<endl;
cout<<"2:"<<mat[2]<<endl;
cout<<"3:"<<mat[3]<<endl;
cout<<"4:"<<mat[4]<<endl<<endl;//*/

//Now analyse each column
enum KEEPING_TRACK {NONE=0, ROOT = 1, STEM = 2, SYM = 3, ASYM =4, LOOP =5 };
KEEPING_TRACK track=NONE;
int temploopbases=0;
bool foundbp=false;
double temp_stem_len, temp_asym_len, temp_asym_nbases, temp_sym_len;
//       int n_asyms,n_syms,n_bases,n_loop_bases
for (c=0;c<maxlen;c++)
    {
    if (track == NONE)
       {//If a track is OFF

       if (considertruncated && !foundbp)
          {//If considertruncated is switched on, the function will ignore anything before a basepair
          if (mat[2][c] == '|') foundbp=true;
          if (!foundbp) continue;
          }

       if (mat[2][c] == '|')
          {//This is part of STEM
          track=STEM;
          temp_stem_len = 1;
          ++(tempstats.n_stems);//Since a new STEM has been encountered, increase STEM count
          tempstats.n_bases+=2;//Since each STEM point has two bases, increase base count by 2
          ++(tempstats.n_base_pairs);//Since a new STEM has been encountered, increase base pair count by 1
          continue;
          }

       if (mat[2][c] != '|' && isalpha(mat[0][c]) && isalpha(mat[4][c]))
          {//This is part of SYM
          track=SYM;
          temp_sym_len = 1;
          tempstats.n_bases+=2;//Since each SYM point has two bases, increase base count by 2
          temploopbases=2;//Since each SYM point can later become a loop point and has two bases
          continue;
          }
       else if (mat[2][c] != '|' && (mat[0][c]=='-' || mat[4][c]=='-'))
          {//This is part of ASYM
          track=ASYM;
          temp_asym_len = 1;
          temp_asym_nbases = 1;
          ++(tempstats.n_bases);//Since each ASYM point has one base, increase base count by 1
          continue;
          }

       if ((isalpha(mat[1][c]) && isalpha(mat[3][c])) || isalpha(mat[2][c]))
          {//This is part of LOOP
          track=LOOP;
          --c;
          continue;
          }
       }//End- If a track is OFF
    else
       {//If a track is ON
//       int n_bases,n_loop_bases

       if (track == STEM)
          {//That is, if a STEM is being monitored
          if (mat[2][c] == '|')
             {//This point is also part of STEM
             tempstats.n_bases+=2;//Since each STEM point has two bases, increase base count by 2
             temp_stem_len++;
             ++(tempstats.n_base_pairs);//Since a new STEM has been encountered, increase base pair count by1
             }
          else
             {//This point no more represents a STEM; so switch OFF the track
             track=NONE;
             tempstats.stem_len.push_back(temp_stem_len);
             --c;
             continue;
             }
          }//End- if a STEM is being monitored

       if (track == SYM)
          {//That is, if a SYM is being monitored
          if (mat[2][c] != '|' && isalpha(mat[0][c]) && isalpha(mat[4][c]))
             {//This point is also part of SYM
             tempstats.n_bases+=2;//Since each SYM point has two bases, increase base count by 2
             temp_sym_len++;
             temploopbases+=2;//Since each SYM point can later become a loop point and has two bases
             }//End- This point is also part of SYM
          else
             {//This point no more represents a SYM; so switch OFF or modify the track
             if (mat[2][c] == '|')
                {//A STEM has begun, switch OFF track and increment SYM count by 1
                track=NONE;
                tempstats.sym_len.push_back(temp_sym_len);
                --c;
                ++(tempstats.n_syms);
                temploopbases=0;
                continue;
                }//End- A STEM has begun

             if (mat[2][c] != '|' && (mat[0][c]=='-' || mat[4][c]=='-'))
                {//This is ASYM, modify track and repeat
                track=ASYM;
                temp_asym_len = temp_sym_len;
                temp_asym_nbases = temp_sym_len*2;
                --c;
                temploopbases=0;
                continue;
                }//End- This is ASYM

             if ((isalpha(mat[1][c]) && isalpha(mat[3][c])) || isalpha(mat[2][c]))
                {//This is LOOP, modify track and go ahead
                track=LOOP;
                --c;
                continue;
                }
             }//End- This point no more represents a SYM
          }//End- if a SYM is being monitored

       if (track == ASYM)
          {//That is, if an ASYM is being monitored
          if (mat[2][c] != '|')
             {//This point is also part of ASYM
             ++(tempstats.n_bases);//Since each ASYM point has one base, increase base count by 1
             temp_asym_len++;
             temp_asym_nbases++;
             }//End- This point is also part of ASYM
          else
             {//This point no more represents an ASYM; so switch OFF the track and increment ASYM count by 1
             track=NONE;
             tempstats.asym_len.push_back(temp_asym_len);
             tempstats.asym_nbases.push_back(temp_asym_nbases);
             --c;
             ++(tempstats.n_asyms);
             continue;
             }//End- This point no more represents an ASYM
          }//End- if an ASYM is being monitored

       if (track == LOOP)
          {//That is, if a LOOP is being monitored
//(isalpha(mat[1][c]) && isalpha(mat[3][c])) || isalpha(mat[2][c])
          if (isalpha(mat[1][c]) && isalpha(mat[2][c]) && isalpha(mat[3][c]))
             {//Last position with three bases
             tempstats.n_bases+=3;//Since this Last LOOP point has three bases, increase base count by 3
             temploopbases+=3;//Since this Last LOOP point has three bases
             }
          else if (isalpha(mat[1][c]) && isalpha(mat[3][c]))
             {//Last position with two bases
             tempstats.n_bases+=2;//Since this Last LOOP point has two bases, increase base count by 2
             temploopbases+=2;//Since this Last LOOP point has two bases
             }
          else if (isalpha(mat[2][c]))
             {//Last position with one base
             ++(tempstats.n_bases);//Since this Last LOOP point has one base, increase base count by 1
             ++temploopbases;//Since this Last LOOP point has one base
             }
          }//End- if a LOOP is being monitored

       }//End- If a track is ON
    }//End- for

if (temploopbases !=0) tempstats.n_loop_bases=temploopbases;

return tempstats;
}

#endif
