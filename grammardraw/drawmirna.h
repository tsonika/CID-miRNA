//**********************************************************************************
// File		: drawmirna.h
// Purpose	: Declaration and definition of classes and functions for the generation of a graphical output of an miRNA
// Author	: Manish Kushwaha (manish.kushwaha@gmail.com)
// Revision	: 1.0	-	28-May-2005	Initial Draft that used a Stub to randomly generate trees and printed them
//		: 2.0	-	15-Jun-2005	Modification to accept input from the miRNA program as a string stream
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
int rnd(int);

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

#endif
