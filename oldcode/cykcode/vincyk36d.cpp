#include <fstream>
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <string.h>
#include <time.h>

using namespace std;

#define NSYM 58
// No. of Non terminal symbols
#define TSYM 4
// No. of Terminal symbols

#define MAXLEN 120 
// Max length of a sequence

#define MINSCORE -0.62

#define MINUSINF -1e90

#define MAXPATH 100
/*
 double ***a; 		//NSYM X NSYM X NSYM
 double **b;  		//NSYM X TSYM
 double ***gammaseq;	//MAXLEN X MAXLEN X NSYM
 bool ***valida;	//NSYM X NSYM X NSYM
 bool **validb;		//NSYM X TSYM
 bool ***validgamma;
 bool skipsequenceflag;
 */
 
 double a[NSYM][NSYM][NSYM]; 		//NSYM X NSYM X NSYM
 double b[NSYM][TSYM];  			//NSYM X TSYM
 double gammaseq[MAXLEN][MAXLEN][NSYM];	//MAXLEN X MAXLEN X NSYM
 bool valida[NSYM][NSYM][NSYM];	//NSYM X NSYM X NSYM
 bool validb[NSYM][TSYM];		//NSYM X TSYM
 bool validgamma[MAXLEN][MAXLEN][NSYM];
 bool skipsequenceflag;
  
 time_t t1, t2;
 fstream fout;
 
 void givetime (long int t_in_sec)
{
 cout << "\t || "<<(t_in_sec / 3600)<<" HOURS - "<< ((t_in_sec % 3600)/60) <<" MINUTES - "<<(t_in_sec % 60)<<" SECONDS ||\n";
 fout << "\t || "<<(t_in_sec / 3600)<<" HOURS - "<< ((t_in_sec % 3600)/60) <<" MINUTES - "<<(t_in_sec % 60)<<" SECONDS ||\n";
}

 
class cyk{

 
 char seq[MAXLEN];
 
 short seqctr;
 short seqlength;
 short code_generator(char x);
 
 void filevaluereader(void);
 void calcgamma(void);
 void initializer(void);
 //void memoryallocator(void);
 //void deallocator(void);
 
 public :
 void invoke(void);
};

void cyk::invoke(void)
{
//memoryallocator();
initializer();
filevaluereader();
calcgamma();
//deallocator();
}


void cyk::calcgamma(void)
{
fstream gm;
char infile[MAXPATH+1];
char outfile[MAXPATH+1];
char sequence[MAXLEN+1];
short base;
short seqlen;
long seqctr = 0, parsable = 0, notparsable = 0;
double maxscore,score,pul;
short i,j,k,v,y,z;

cout <<"\n Enter Sequence file : ";
cin.getline(infile,MAXPATH);

cout << "\n Enter Output file : ";
cin.getline(outfile,MAXPATH);

fout.open(outfile,ios::out);
gm.open(infile,ios::in);
if(!gm.good() || !fout.good())
{
 cout << "\n Error Opening Sequence // Output file \n";
 exit(1);
 }
 
(void)time(&t1); 
// start clocking

gm >> sequence;

do
{
skipsequenceflag = false;
//gm.getline(sequence,MAXLEN);

seqlen = strlen(sequence);
cout << "\n sequence : " << seqlen << "\t"<<sequence;
if(seqlen <= 0)
continue;

++seqctr;

// clear gammaseq and validgamma
for(i = 0; i < MAXLEN; ++i)
 for(j = 0; j < MAXLEN; ++j)
  for(v = 0; v < NSYM; ++v)
  {
  gammaseq[i][j][v] = MINUSINF;
  validgamma[i][j][v] = false;
  }

// actual work begins

for(i = 0; ((i < seqlen) && (!skipsequenceflag)); ++i)
 for(v = 0; ((v < NSYM) && (!skipsequenceflag)); ++v)
 {
  //cout << "\nBASE"<<(base = code_generator(sequence[i]))<<" ## ";
  base = code_generator(sequence[i]);
  
  if(skipsequenceflag) continue;
  if(validb[v][base])
  {
   //cout << "\n HELLO ";
   gammaseq[i][i][v] = b[v][base];
   //cout<<"\nGAMMA "<<i<<"\t"<<i<<"\t"<<v <<"\t"<<exp(gammaseq[i][i][v]);
   validgamma[i][i][v] = true;
   }
   
   }
   
   if(skipsequenceflag)
   {
    cout << "\n Unknow character found in sequence # " << seqctr;
    continue;
   }
 
 
  for(i= 1;i < (seqlen-1); ++i)
  {
  for(j = 0; j < (seqlen-i); ++j)
  {
  for(v = NSYM-1; v >= 0; --v)
  { 
  // calcing gammaseq j,j+i,v
  maxscore = MINUSINF;
  
   for(y = NSYM-1; y >= 0; --y)
    for(z = NSYM-1; z >= 0; --z)
  {
  if(valida[v][y][z])
  {
  
   for(k=0; k < i; ++k)
   {
    if(validgamma[j][j+k][y] && validgamma[j+k+1][j+i][z])
    {
    
    score = gammaseq[j][j+k][y] + gammaseq[j+k+1][j+i][z] + a[v][y][z];
    if(score > maxscore)
     {
      maxscore = score;
      validgamma[j][j+i][v] = true;
      //cout<<"\nGAMMA "<<j<<"\t"<<j+i<<"\t"<<v <<"\t"<<exp(maxscore);
     }
    }
   } // k loop
    
   } // if valid a v y z
   
   } // y-z loop
   if(validgamma[j][j+i][v])
   {
   gammaseq[j][j+i][v] = maxscore;
   //cout<<"\n FINAL SCORE "<<exp(gammaseq[j][j+i][v]);
   
   }
   } // v loop
  } // j loop
   } // i loop
   
 i = seqlen - 1;
 j = 0;
 v = 0;
  
  // calcing gammaseq 0,len-1,0
  maxscore = MINUSINF;
  
   for(y = NSYM-1; y >= 0; --y)
    for(z = NSYM-1; z >= 0; --z)
  {
  if(valida[v][y][z])
  {
  
   for(k=0; k < i; ++k)
   {
    if(validgamma[j][j+k][y] && validgamma[j+k+1][j+i][z])
    {
    
    score = gammaseq[j][j+k][y] + gammaseq[j+k+1][j+i][z] + a[v][y][z];
    if(score > maxscore)
     {
      maxscore = score;
      validgamma[j][j+i][v] = true;
      //cout<<"\nGAMMA "<<j<<"\t"<<j+i<<"\t"<<v <<"\t"<<exp(maxscore);
     }
    }
   } // k loop
    
   } // if valid a v y z
   
   } // y-z loop
   if(validgamma[j][j+i][v])
   {
   gammaseq[j][j+i][v] = maxscore;
   //cout<<"\n FINAL SCORE "<<exp(gammaseq[j][j+i][v]);
   
   }
    // v loop ends

   pul = gammaseq[0][seqlen-1][0]/seqlen;  // score per unit length
  
 
 if(!(validgamma[0][seqlen-1][0]) || pul < MINSCORE)
  {
   	//fout << "\n Sequence # "<<seqctr << " CANNOT BE PARSED "<<(gammaseq[0][seqlen-1][0]);
   cout << "\n Sequence # "<<seqctr << " CANNOT BE PARSED "<<(gammaseq[0][seqlen-1][0]);
   notparsable++;
   }
 else
  {
   cout << "\n Sequence # "<<seqctr << "\t SCORE = " << (gammaseq[0][seqlen-1][0]) <<"\t p.u.l = "<< pul;
   
   fout << "\n Sequence # "<<seqctr << "\t LENGTH = " << seqlen<< "\t SCORE = " << (gammaseq[0][seqlen-1][0]) <<"\t p.u.l = "<< pul
         << "\n" <<sequence;
   parsable++;
   
   }
   gm >> sequence;
   }while(!gm.eof());
   
   
   (void)time(&t2);
   cout << "\n Finally "<<parsable<<" sequences parsed and "<< notparsable <<" sequences rejected";
   cout << "\n End of Parsing of total "<< seqctr << " sequences. \n";
   fout << "\n Finally "<<parsable<<" sequences parsed and "<< notparsable <<" sequences rejected";
   fout << "\n End of Parsing of total "<< seqctr << " sequences. \n";

   cout << "\n TOTAL TIME TAKEN BY PROGRAM : "; givetime((long)(t2-t1)); 

   fout.close();
   gm.close();
}
   
 
 

short cyk::code_generator(char x)
{
 if((x == 'a')||(x == 'A')) return 0;
 if((x == 'u')||(x == 'U')) return 1;
 if((x == 'g')||(x == 'G')) return 2;
 if((x == 'c')||(x == 'C')) return 3;
 return -1;
 cout <<"\n This sequence contains characters which cannot be parsed";
 skipsequenceflag = true;
}


void cyk::initializer(void)
{

short i,j,k;

for( i = 0; i < NSYM; ++i)
 for( j = 0; j < NSYM; ++j)
  for( k = 0; k < NSYM; ++k)
   {
    a[i][j][k] = MINUSINF;
    valida[i][j][k] = false;
   } 
for( i = 0; i < NSYM; ++i)
 for( j = 0; j < TSYM; ++j)
 {
 b[i][j] = MINUSINF;
 validb[i][j] = false;
 }
} 


void cyk::filevaluereader(void)
{
 
// This function reads the a and b values in the file format :
// a 1 2 3 0.67
// b 4 5 0.55

// Basically the first line says that a[1][2][3] = 0.67
// and second line says b[4][5] = 0.55

 char filename[MAXPATH];
 fstream fvr;
 char ch = 0;
 short pm[3] = {0,0,0};
 double val = 0.0;

 cout << "\n Enter the name of the Input Probability File : ";
 cin.getline(filename, MAXPATH - 1);

 fvr.open(filename,ios::in);
 if(!fvr.good())
 {
  cout << "\n Error Opening File !! \n";
  exit(1);
 }
 while(!fvr.eof())
 {
  fvr>>ch;
  if((ch == 'a') || (ch == 'A'))
  {
   fvr>>pm[0];
   fvr>>pm[1];
   fvr>>pm[2];
   fvr>>val;
   if((pm[0] < NSYM) &&(pm[1] < NSYM) &&(pm[2] < NSYM))
    {
  a[pm[0]][pm[1]][pm[2]] = log10(val);
  valida[pm[0]][pm[1]][pm[2]] = true;
  cout << "\n a : "<< pm[0] <<"\t" << pm[1] << "\t"<< pm[2] << "\t"<<exp(a[pm[0]][pm[1]][pm[2]]);

  }
  }
  else
  { if((ch == 'b') || (ch == 'B'))
  {
   fvr>>pm[0];
   fvr>>pm[1];
   fvr>>val;

   if((pm[0] < NSYM) &&(pm[1] < TSYM))
    {
    b[pm[0]][pm[1]] = log10(val);
    cout << "\n b : "<< pm[0] <<"\t" << pm[1]  << "\t"<< exp(b[pm[0]][pm[1]]);
    validb[pm[0]][pm[1]] = true;
  }
  }
  
  else
{
cout << "\n Unable to parse Input Probability File"
     << "\n Check the Format";
fvr.close();
exit(1);
}
}
 }
 fvr.close();
}

/*
void cyk::memoryallocator(void)
{
short i,j,k;

gammaseq = new double**[MAXLEN];
if(!gammaseq)
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }
 
for(i = 0; i < MAXLEN; ++i)
{
gammaseq[i] = new double*[MAXLEN];
if(!gammaseq[i])
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }
for(j = 0; j < MAXLEN; ++j)
{
gammaseq[i][j] = new double[NSYM];
if(!gammaseq[i][j])
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }
}
}


a = new double**[NSYM];
if(!a)
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }
 
for(i = 0; i < NSYM; ++i)
{
a[i] = new double*[NSYM];
if(!a[i])
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }
for(j = 0; j < NSYM; ++j)
{
a[i][j] = new double[NSYM];
if(!a[i][j])
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }
}
}

b = new double*[NSYM];
if(!b)
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }

for(i = 0; i < NSYM; ++i)
{
b[i] = new double[TSYM];
if(!b[i])
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }
 
}

// bools

validgamma = new bool**[MAXLEN];
if(!validgamma)
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }
 
for(i = 0; i < MAXLEN; ++i)
{
validgamma[i] = new bool*[MAXLEN];
if(!validgamma[i])
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }
for(j = 0; j < MAXLEN; ++j)
{
validgamma[i][j] = new bool[NSYM];
if(!valida[i][j])
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }
}
}


valida = new bool**[NSYM];
if(!valida)
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }
 
for(i = 0; i < NSYM; ++i)
{
valida[i] = new bool*[NSYM];
if(!valida[i])
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }
for(j = 0; j < NSYM; ++j)
{
valida[i][j] = new bool[NSYM];
if(!valida[i][j])
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }
}
}

validb = new bool*[NSYM];
if(!validb)
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }

for(i = 0; i < NSYM; ++i)
{
validb[i] = new bool[TSYM];
if(!validb[i])
{
 cout << "\n Unable to allocate Memory\n";
 exit(1);
 }
 
}

}


void cyk::deallocator(void)
{
short i,j,k;

for(i = 0; i < MAXLEN; ++i)
{
for(j = 0; j < MAXLEN; ++j)
{
delete [] gammaseq[i][j];
}

delete [] gammaseq[i];
}
delete [] gammaseq;


for(i = 0; i < NSYM; ++i)
{
for(j = 0; j < NSYM; ++j)
{
delete [] a[i][j];
}

delete [] a[i];
}
delete [] a;

for(i = 0; i < NSYM; ++i)
{
delete [] b[i];
}
delete [] b;

// bools



for(i = 0; i < MAXLEN; ++i)
{
for(j = 0; j < MAXLEN; ++j)
{
delete [] validgamma[i][j];
}

delete [] validgamma[i];
}
delete [] validgamma;

for(i = 0; i < NSYM; ++i)
{
for(j = 0; j < NSYM; ++j)
{
delete [] valida[i][j];
}

delete [] valida[i];
}
delete [] valida;

for(i = 0; i < NSYM; ++i)
{
delete [] validb[i];
}
delete [] validb;

}
*/

int main()
{

cyk c;
c.invoke();
return 0;
}
