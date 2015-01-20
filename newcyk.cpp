#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <cmath>
#include <ctime>
using namespace std;

#define MAXLEN 130
#define NSYM 60
#define TSYM 4
#define MAXSEQ 300
#define MINUSINF -1e90
#define START_NT 0
#define MAXPATH 200
#define MAXNUM 240002 // maximum no. of sequence that can be picked up from a file

struct entry {
       long nterm;
       double score;
};

struct cell {

       long members;           //basically stores the number of members in the array
       entry node[NSYM];

       void operator= (cell &xyz);
};

long eqi;

void cell::operator=(cell &xyz)
{
       members = xyz.members;
       for(eqi = 0 ; eqi < NSYM; ++eqi)
               node[eqi] = xyz.node[eqi];

}

ifstream inprobfile, inseqfile;
ofstream resultfile;

static char inprobfilename[MAXPATH];
static char inseqfilename[MAXPATH];
static char resultfilename[MAXPATH];
static char allseq[MAXNUM][MAXLEN+1];

//bool existrule[NSYM][NSYM];
cell hashrules[NSYM][NSYM];
cell hashterminals[TSYM];

cell matrix[MAXLEN][MAXLEN];

long modsequence[MAXLEN];
long seqctr;

long i;
long j;
long p;
long n1;
long n2;
long n3;
long n4;
long nrule;

bool flag = false;
//cell &x, &y, &z;
cell *x, *y, *z;
double score;

void initial_value_reader();
// This function has to first reset everything and then
// Fill up Hashrule 2-D matrix and
// Fill up the hash matrix for Terminals

bool modsequence_creator(const char* seq);
// This will convert the character sequence longo an equivalent
// numeric sequence by converting each Terminal longo equivalent numeric code
// returns true if successful conversion else returns false

double cyk_calculator(const char* sequence);
// calculates the score

void read_sequences(void);
// read all sequences

void givetime (long int t_in_sec);
// give time in human readable format

time_t t1,t2, t3, t4;

int main(int argc, char* argv[])
{
 //const char sequence[MAXLEN] = "acauugcuacuuacaauuaguuuugcagguuugcauuucagcguauauauguauauguggcugugcaaauccaugcaaaacugauugugauaaugu";
 //strcpy(inprobfilename, "finalprobs.txt");

 //const char sequence[MAXLEN] = "acucaacuc";
 //strcpy(inprobfilename, "testprobs.txt");

       if(argc < 4)
       { cout  << endl << "INSUFFICIENT ARGMENTS" << endl
                       << "Format is : Program <input sequence file> <input prob. file> <output result file> "
                       << endl;
       return 1;
       }

       strcpy(inseqfilename,argv[1]);
       strcpy(inprobfilename,argv[2]);
       strcpy(resultfilename,argv[3]);

 initial_value_reader();
 read_sequences();

 resultfile.open(resultfilename);

 if(!resultfile.good())
 {
         cout << endl << " ERROR OPENING RESULT FILE !! " << endl;
         exit(1);
 }

 (void)time(&t1);
       //Start Clocking;

 double myscore;
 long seqlen, mi;
 for( mi = 0; mi < seqctr; ++mi)
 {
       myscore = cyk_calculator(allseq[mi]);
       seqlen  = strlen(allseq[mi]);

       cout << "\n Sequence : " << mi+1 << endl << allseq[mi] ;
       cout << "\n Length : " << seqlen;
       cout << "\t Normal SCORE = " << myscore/seqlen;
       cout << "\t SCORE = " << myscore;

       resultfile << "\n Sequence : " << mi+1 << endl << allseq[mi] ;
       resultfile << "\n Length : " << seqlen;
       resultfile << "\t Normal SCORE = " << myscore/seqlen;
       resultfile << "\t SCORE = " << myscore;

 }

 (void)time(&t2);
       //End Clocking;
       cout << "\n TOTAL TIME TAKEN BY PROGRAM : \n";
       givetime((long int)(t2-t1));

 resultfile.close();
 //cin.get();
 return 0;

}

double cyk_calculator(const char* sequence)
{

       long seqlen = strlen(sequence);
       if(seqlen <= 0) return MINUSINF;

       if(!modsequence_creator(sequence)) return MINUSINF;

       for(j = 0; j < seqlen; ++j)
               matrix[0][j] = hashterminals[modsequence[j]];

       // I have initialized the first row according to each terminal

       for(i = 1; i < seqlen; ++i)
       {
               for(j = 0; j < (seqlen - i); ++j)
               {
                       matrix[i][j].members = 0;
                       for(p = 0; p < i; ++p)
                       {
                               // Taking cartesian product of matrix[p][j] and matrix[i-p-1][j+p+1]

                               x = &matrix[p][j];
                               y = &matrix[i-p-1][j+p+1];

                               //cell &x = matrix[p][j];
                               //cell &y = matrix[i-p-1][j+p+1];

                               //for(n1 = 0; n1 < matrix[p][j].members; ++n1)
                                 for(n1 = 0; n1 < x->members; ++n1)
                               //  for(n1 = 0; n1 < x.members; ++n1)

                               {
                                       //for(n2 = 0; n2 < matrix[i-p-1][j+p+1].members; ++n2)
                                       for(n2 = 0; n2 < y->members; ++n2)
                                       //  for(n2 = 0; n2 < y.members; ++n2)
                                       {
                                               z = &hashrules[matrix[p][j].node[n1].nterm][matrix[i-p-1][j+p+1].node[n2].nterm];
                                               //cell &z = hashrules[matrix[p][j].node[n1].nterm][matrix[i-p-1][j+p+1].node[n2].nterm];

                                               //for(n3 = 0; n3 < hashrules[matrix[p][j].node[n1].nterm][matrix[i-p-1][j+p+1].node[n2].nterm].members; ++n3)
                                                 for(n3 = 0; n3 < z->members; ++n3)
                                               // for(n3 = 0; n3 < z.members; ++n3)
                                               {
                                                       flag = false;
                                                       //score = matrix[p][j].node[n1].score + matrix[i-p-1][j+p+1].node[n2].score + hashrules[matrix[p][j].node[n1].nterm][matrix[i-p-1][j+p+1].node[n2].nterm].node[n3].score;
                                                         score = x->node[n1].score + y->node[n2].score + z->node[n3].score;
                                                       // score = x.node[n1].score + y.node[n2].score + z.node[n3].score;
                                                       //nrule = hashrules[matrix[p][j].node[n1].nterm][matrix[i-p-1][j+p+1].node[n2].nterm].node[n3].nterm;
                                                         nrule = z->node[n3].nterm;
                                                       //   nrule = z.node[n3].nterm;

                                                               for(n4 = 0; n4 < matrix[i][j].members; ++n4)
                                                               {
                                                                       if(nrule == matrix[i][j].node[n4].nterm)
                                                                       {
                                                                               flag = true;
                                                                               if (score > matrix[i][j].node[n4].score)
                                                                               {
                                                                                       matrix[i][j].node[n4].score = score;
                                                                                       break;
                                                                               }

                                                                       }
                                                               }

                                                               if(!flag)
                                                               {
                                                                       // adding new entry
                                                                       matrix[i][j].node[matrix[i][j].members].nterm = nrule;
                                                                       matrix[i][j].node[matrix[i][j].members++].score = score;
                                                               }

                                               } // n3 loop
                                       }
                               }
                       }
               }
       }

       flag = false;
       for(n1 = 0; n1 < matrix[seqlen-1][0].members; ++n1)
       {
               if(matrix[seqlen-1][0].node[n1].nterm == START_NT)
               {
                       score = matrix[seqlen-1][0].node[n1].score;
                       flag = true;
                       break;
               }
       }

       if(flag) return score;

       return MINUSINF;
}

// This function has to first reset everything and then
// Fill up Hashrules 2-D matrix and
// Fill up the hash matrix for Terminals
//cell hashrules[NSYM][NSYM];
//cell hashterminals[TSYM];

//cell matrix[MAXLEN][MAXLEN];

void initial_value_reader()
{

       long vi , j;

       for(vi = 0; vi < MAXLEN; ++vi)
               for(j = 0 ; j < MAXLEN; ++j)
                       matrix[vi][j].members = 0;

       for(vi = 0; vi < NSYM; ++vi)
               for(j = 0; j < NSYM; ++j)
                       hashrules[vi][j].members = 0;

       for(vi = 0; vi < TSYM; ++vi)
               hashterminals[vi].members = 0;

       inprobfile.open(inprobfilename, ios::in);

       char ch = 0;
       long temp;
       long pm[3] = {0,0,0};
       long double val = 0.0;
       char tempstream[MAXPATH+1];

       if(!inprobfile.good())
       {
               cout << "\n Error Opening Input Probability File !! \n";
               exit(1);
       }

       while(!inprobfile.eof())
       {
               inprobfile.getline(tempstream,MAXPATH);
               if(5 > strlen(tempstream)) continue;
               //just checking for basic length of one input line (can't be less than 5)

               istringstream iss(tempstream);

               if(!(iss>>ch)) continue;

               if((ch == 'a') || (ch == 'A'))
               {
                       iss>>pm[0];
                       iss>>pm[1];
                       iss>>pm[2];
                       iss>>val;
                       if((pm[0] < NSYM) &&(pm[1] < NSYM) &&(pm[2] < NSYM))
                       {
                               temp = hashrules[pm[1]][pm[2]].members++;
                               hashrules[pm[1]][pm[2]].node[temp].nterm = pm[0];
                               hashrules[pm[1]][pm[2]].node[temp].score = log10(val);

                               cout << "\n a : "<< pm[0] <<"\t" << pm[1] << "\t"<< pm[2] << "\t"<< val; // (a[pm[0]][pm[1]][pm[2]]);
                       }
               }
               else
               {
                       if((ch == 'b') || (ch == 'B'))
                       {
                               iss>>pm[0];
                               iss>>pm[1];
                               iss>>val;

                               if((pm[0] < NSYM) &&(pm[1] < TSYM))
                               {
                                       //b[pm[0]][pm[1]] = val;
                                       //validb[pm[0]][pm[1]] = true;

                                       temp = hashterminals[pm[1]].members++;
                                       hashterminals[pm[1]].node[temp].nterm = pm[0];
                                       hashterminals[pm[1]].node[temp].score = log10(val);

                                       cout << "\n b : "<< pm[0] <<"\t" << pm[1] << "\t"<< val; //(b[pm[0]][pm[1]]);
                               }
                       }

                       else if(ch == '#') continue;
                       //Comment Symbol

                       else if(!inprobfile.eof())
                       {
                               cout << "\n Unable to parse Input Probability File" << "\n Check the Format";
                               inprobfile.close();
                               exit(1);
                       }
               }
       }
       inprobfile.close();
}

bool modsequence_creator(const char* sequence)
{
       long len = strlen(sequence);
       if (len <= 0) return false;

       for(int j = 0; j < len; ++j)
       {
                        if((sequence[j] == 'a') || (sequence[j] == 'A')) modsequence[j] = 0;
                       else if((sequence[j] == 'u') || (sequence[j] == 'U')) modsequence[j] = 1;
                       else if((sequence[j] == 'g') || (sequence[j] == 'G')) modsequence[j] = 2;
                       else if((sequence[j] == 'c') || (sequence[j] == 'C')) modsequence[j] = 3;
                       else
                       {
                               cout << endl << "Seq # " << sequence << endl << " Illegal char found. Rejected";
                               return false;

                       }
       }
       return true;
}

void givetime (long int t_in_sec)
{
cout << "\t || "<<(t_in_sec / 3600)<<" HOURS - "
         << ((t_in_sec % 3600)/60) <<" MINUTES - "<<(t_in_sec % 60)<<" SECONDS ||\n";
}

void read_sequences(void)
{
       inseqfile.open(inseqfilename);
       if(!inseqfile.good())
       {
               cout << "\n Error Opening Input Sequence File !! \n";
               exit(1);
       }

       seqctr = 0;
       static char tempstream[MAXLEN+1];

       while(!inseqfile.eof() && (seqctr < MAXNUM))
       {
               //inseqfile.getline(tempstream,MAXLEN+1,'\n');
               inseqfile >> tempstream;
               if((strlen(tempstream) <= 1) || strstr(tempstream, ">")) continue;

               strcpy(allseq[seqctr],tempstream);
               ++seqctr;
       }

       inseqfile.close();
}



