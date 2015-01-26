//______________________________________________________________
//
// INSIDE OUTSIDE ALGORITHM v4.0 (2004)
// Author       : Vipin Gupta
//______________________________________________________________

 
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <ctime>
    using namespace std;

 
#define MAXPATH 100
static char inprobfile[MAXPATH];

static char outprobfile[MAXPATH];

static char inseqfile[MAXPATH];

 
#define MAXLEN 120
#define NSYM 58
#define TSYM 4
#define MAXNUM 100
    
#define APPROX_ZERO 1e-7
#define THRESH_LIMIT 1e-5
//threshhold percentage limit to break reestimation

static long double a[NSYM][NSYM][NSYM];

static long double b[NSYM][TSYM];

static bool valida[NSYM][NSYM][NSYM];

static bool validb[NSYM][TSYM];

 
static bool validalpha[MAXLEN][MAXLEN][NSYM];

static bool validbeta[MAXLEN][MAXLEN][NSYM];

 
static long double Na[NSYM][NSYM][NSYM];

static long double Nb[NSYM][TSYM];

static long double Da[NSYM][NSYM][NSYM];

static long double Db[NSYM][TSYM];

 
static long double alpha[MAXLEN][MAXLEN][NSYM];

static long double beta[MAXLEN][MAXLEN][NSYM];

 
static long double newsum = 0;

static long double oldsum = 0;

 
static char seqholder[MAXNUM][MAXLEN + 1];

static int num_of_seq = 0;

static char sequence[MAXLEN + 1];

static int modsequence[MAXLEN + 1];

static int seqlen;

 
time_t t1, t2, t3, t4;

 
void givetime (long int);

//function to measure time taken by program
    
class ALPHA {
  
 
public:
void initialize_alpha (void);
    
void calc_alpha (void);

} _alpha;


 
 
 
void
ALPHA::initialize_alpha (void) 
{
    
int i, j, v;
    
 
for (i = 0; i < MAXLEN; ++i)
        
 {
        
for (j = 0; j < MAXLEN; ++j)
            
 {
            
for (v = 0; v < NSYM; ++v)
                
 {
                
alpha[i][j][v] = 0;
                
validalpha[i][j][v] = false;
                
}
            
}
        
}
    
 
for (i = 0; i < seqlen; ++i)
        
 {
        
for (v = 0; v < NSYM; ++v)
            
 {
            
if (validb[v][modsequence[i]])
                
 {
                
alpha[i][i][v] = b[v][modsequence[i]];
                
validalpha[i][i][v] = true;
                
}
            
}
        
}

}


 
void
ALPHA::calc_alpha (void) 
{
    
int i, j, v, k, y, z, diff;
    
 
for (diff = 1; diff < seqlen; ++diff)
        
 {
        
for (i = 0; i < (seqlen - 1); ++i)
            
 {
            
j = i + diff;
            
if (j < seqlen)
                
 {
                
for (v = 0; v < NSYM; ++v)
                    
 {
                    
for (k = i; k < j; ++k)
                        
 {
                        
for (z = 0; z < NSYM; ++z)
                            
 {
                            
for (y = 0; y < NSYM; ++y)
                                
 {
                                
if (valida[v][y][z]
                                     && validalpha[i][k][y]
                                     && validalpha[k + 1][j][z])
                                    
 {
                                    
alpha[i][j][v] +=
                                        (alpha[i][k][y] *
                                         alpha[k + 1][j][z] * a[v][y][z]);
                                    
validalpha[i][j][v] = true;
                                    
}
                                
}
                            
}
                        
}
                    
                        //if(validalpha[i][j][v])
                        //cout << "\n alpha "<< i <<" "<<j<<" "<<v<<" " << alpha[i][j][v];              
                    }
                
}
            
}
        
}

 
}


 
 
 
class BETA {
  
public:
void initialize_beta (void);
    
void calc_beta (void);

} _beta;


 
void
BETA::initialize_beta (void) 
{
    
int i, j, v;
    
 
for (i = 0; i < MAXLEN; ++i)
        
 {
        
for (j = 0; j < MAXLEN; ++j)
            
 {
            
for (v = 0; v < NSYM; ++v)
                
 {
                
beta[i][j][v] = 0;
                
validbeta[i][j][v] = false;
                
}
            
}
        
}
    
 
beta[0][seqlen - 1][0] = 1;
    
validbeta[0][seqlen - 1][0] = true;

}


 
void
BETA::calc_beta (void) 
{
    
int i, j, v, y, z, k, diff;
    
long double sum1, sum2;
    
 
for (diff = seqlen - 2; diff >= 0; --diff)
        
for (i = 0; i < (seqlen - diff); ++i)
            
 {
            
j = i + diff;
            
for (v = 0; v < NSYM; ++v)
                
 {
                
sum1 = 0;
                
sum2 = 0;
                
 
if (i > 0)
                    
 {
                    
for (k = 0; k < i; ++k)
                        
 {
                        
for (z = 0; z < NSYM; ++z)
                            
 {
                            
for (y = 0; y < NSYM; ++y)
                                
 {
                                
if (valida[y][z][v]
                                     && validalpha[k][i - 1][z]
                                     && validbeta[k][j][y])
                                    
 {
                                    
sum1 +=
                                        (alpha[k][i - 1][z] *
                                         beta[k][j][y] * a[y][z][v]);
                                    
}
                                
}
                            
}
                        
}
                    
 
}
                
 
if (j < (seqlen - 1))
                    
 {
                    
for (k = j + 1; k < seqlen; ++k)
                        
 {
                        
for (z = 0; z < NSYM; ++z)
                            
 {
                            
for (y = 0; y < NSYM; ++y)
                                
 {
                                
if (valida[y][v][z]
                                     && validalpha[j + 1][k][z]
                                     && validbeta[i][k][y])
                                    
 {
                                    
sum2 +=
                                        (alpha[j + 1][k][z] *
                                         beta[i][k][y] * a[y][v][z]);
                                    
}
                                
}
                            
}
                        
}
                    
}
                
 
beta[i][j][v] = sum1 + sum2;
                
if (beta[i][j][v] > 0.0)
                    
 {
                    
validbeta[i][j][v] = true;
                    
                        //cout << "\n beta "<< i <<" "<<j<<" "<<v<<" " << beta[i][j][v];
                    }
                
}
            
}

}


 
 
class reestimation {
  
public:
void updateA (void);
    
void updateB (void);

} _reestimate;


 
void
reestimation::updateB (void) 
{
    
long double sum1, sum2;
    
int i, j, v, t;
    
 
for (v = 0; v < NSYM; ++v)
        
for (t = 0; t < TSYM; ++t)
            
 {
            
 
if (validb[v][t])
                
 {
                
sum1 = 0;
                
sum2 = 0;
                
 
for (i = 0; i < seqlen; ++i)
                    
 {
                    
if ((modsequence[i] == t) && validbeta[i][i][v])
                        
sum1 += (beta[i][i][v] * b[v][t]);
                    
}
                
 
for (i = 0; i < seqlen; ++i)
                    
 {
                    
for (j = i; j < seqlen; ++j)
                        
 {
                        
if (validbeta[i][j][v]
                             && validalpha[i][j][v])
                            
sum2 += (beta[i][j][v] * alpha[i][j][v]);
                        
}
                    
}
                
 
Nb[v][t] += (sum1 / alpha[0][seqlen - 1][0]);
                
Db[v][t] += (sum2 / alpha[0][seqlen - 1][0]);
                
}
            
 
}

}


 
void
reestimation::updateA (void) 
{
    
int i, j, k, v, y, z;
    
long double sum1, sum2;
    
 
for (v = 0; v < NSYM; ++v)
        
for (y = 0; y < NSYM; ++y)
            
for (z = 0; z < NSYM; ++z)
                
 {
                
 
if (valida[v][y][z])
                    
 {
                    
sum1 = 0;
                    
sum2 = 0;
                    
 
for (i = 0; i < (seqlen - 1); ++i)
                        
for (j = i + 1; j < seqlen; ++j)
                            
if (validbeta[i][j][v])
                                
 {
                                
for (k = i; k < j; ++k)
                                    
if (validalpha[i][k][y]
                                         && validalpha[k + 1][j][z])
                                        
sum1 +=
                                            (beta[i][j][v] *
                                             a[v][y][z] *
                                             alpha[i][k][y] * alpha[k + 1]
                                             [j][z]);
                                
}
                    
 
for (i = 0; i < seqlen; ++i)
                        
for (j = i; j < seqlen; ++j)
                            
if (validbeta[i][j][v]
                                 && validalpha[i][j][v])
                                
sum2 += (beta[i][j][v] * alpha[i][j][v]);
                    
 
Na[v][y][z] += (sum1 / alpha[0][seqlen - 1][0]);
                    
Da[v][y][z] += (sum2 / alpha[0][seqlen - 1][0]);
                    
}
                
 
}

}


 
 
class IOcontrol {
  
private:
void initializer (void);
    // set Na Nb Da Db to zero
    void probreader (void);     // and put up a and b values
    void seqreader (void);
    
void estimator (void);
    
        // for each sequence read,make modsequence, run alpha and beta calcs, and reestimate
  public:
void finalcalculator (void);
    
        // calculate final reestimations and output to file
};


 
void
IOcontrol::initializer (void) 
{
    
 
int i, j, k;
    
 
for (i = 0; i < NSYM; ++i)
        
 {
        
for (j = 0; j < NSYM; ++j)
            
for (k = 0; k < NSYM; ++k)
                
 {
                
Na[i][j][k] = 0;
                
Da[i][j][k] = 0;
                
}
        
 
for (j = 0; j < TSYM; ++j)
            
 {
            
Nb[i][j] = 0;
            
Db[i][j] = 0;
            
}
        
}

 
}


 
void
IOcontrol::probreader () 
{
    
 
int i, j, k;
    
 
for (i = 0; i < NSYM; ++i)
        
 {
        
for (j = 0; j < NSYM; ++j)
            
for (k = 0; k < NSYM; ++k)
                
 {
                
a[i][j][k] = 0;
                
valida[i][j][k] = false;
                
}
        
 
for (j = 0; j < TSYM; ++j)
            
 {
            
b[i][j] = 0;
            
validb[i][j] = false;
            
}
        
}
    
 
 
ifstream fvr (inprobfile);
    
 
 
        // This function reads the a and b values in the file format :
        // a 1 2 3 0.67
        // b 4 5 0.55
        
        // Basically the first line says that a[1][2][3] = 0.67
        // and second line says b[4][5] = 0.55
    
char ch = 0;
    
long pm[3] = { 0, 0, 0 };
    
long double val = 0.0;
    
char tempstream[MAXPATH + 1];
    
 
if (!fvr.good ())
        
 {
        
cout << "\n Error Opening Input Probability File !! \n";
        
exit (1);
        
}
    
 
while (!fvr.eof ())
        
 {
        
fvr.getline (tempstream, MAXPATH);
        
if (5 > strlen (tempstream))
            continue;
        
            //just checking for basic length of one input line (can't be less than 5)
            
istringstream iss (tempstream);
        
 
if (!(iss >> ch))
            continue;
        
 
if ((ch == 'a') || (ch == 'A'))
            
 {
            
iss >> pm[0];
            
iss >> pm[1];
            
iss >> pm[2];
            
iss >> val;
            
if ((pm[0] < NSYM) && (pm[1] < NSYM) && (pm[2] < NSYM))
                
 {
                
a[pm[0]][pm[1]][pm[2]] = val;
                
valida[pm[0]][pm[1]][pm[2]] = true;
                
cout << "\n a : " << pm[0] << "\t" << pm[1] << "\t" <<
                    pm[2] << "\t" << (a[pm[0]][pm[1]][pm[2]]);
                
}
            
}
        
        else
            
 {
            
if ((ch == 'b') || (ch == 'B'))
                
 {
                
iss >> pm[0];
                
iss >> pm[1];
                
iss >> val;
                
 
if ((pm[0] < NSYM) && (pm[1] < TSYM))
                    
 {
                    
b[pm[0]][pm[1]] = val;
                    
validb[pm[0]][pm[1]] = true;
                    
cout << "\n b : " << pm[0] << "\t" << pm[1] <<
                        "\t" << (b[pm[0]][pm[1]]);
                    
}
                
}
            
 
            else if (ch == '#')
                continue;
            
                //Comment Symbol
                
            else if (!fvr.eof ())
                
 {
                
cout << "\n Unable to parse Input Probability File" <<
                    "\n Check the Format";
                
fvr.close ();
                
exit (1);
                
}
            
}
        
}
    
fvr.close ();

}


 
 
void
IOcontrol::seqreader (void) 
{
    
ifstream finseq (inseqfile);
    
 
 
if (!finseq.good ())
        
 {
        
cout << "\n Unable to parse Input Sequence File \n";
        
exit (1);
        
}
    
 
while (!finseq.eof ())
        
 {
        
finseq.getline (seqholder[num_of_seq], MAXLEN);
        
seqlen = strlen (seqholder[num_of_seq]);
        
if (!seqlen)
            continue;
        
 
++num_of_seq;
        
 
}
    
 
finseq.close ();

}


 
 
void
IOcontrol::estimator (void) 
{
    
// for each sequence read,make modsequence, run alpha and beta calcs, and reestimate
    int i = 0, j = 0, reject = 0;
    
bool skipsequenceflag = false;
    
 
for (i = 0; i < num_of_seq; ++i)
        
 {
        
strcpy (sequence, seqholder[i]);
        
seqlen = strlen (sequence);
        
cout << endl << " Analyzing Seq. " << "# " << i +
            1 << "\t Length : " << seqlen;
        
 
skipsequenceflag = false;
        
 
for (j = 0; j < seqlen; ++j)
            
 {
            
if ((sequence[j] == 'a') || (sequence[j] == 'A'))
                modsequence[j] = 0;
            
            else if ((sequence[j] == 'u') || (sequence[j] == 'U'))
                modsequence[j] = 1;
            
            else if ((sequence[j] == 'g') || (sequence[j] == 'G'))
                modsequence[j] = 2;
            
            else if ((sequence[j] == 'c') || (sequence[j] == 'C'))
                modsequence[j] = 3;
            
            else
                
 {
                
cout << endl << "Seq # " << i +
                    1 << " Illegal char found. Rejected";
                
skipsequenceflag = true;
                
++reject;
                
break;
                
}
            
}
        
 
if (skipsequenceflag)
            continue;
        
 
_alpha.initialize_alpha ();
        
_alpha.calc_alpha ();
        
 
if (alpha[0][seqlen - 1][0] <= 0)
            
 {
            
cout << endl << "Seq # " << i +
                1 << " sequence prob zero. Rejected";
            
++reject;
            
skipsequenceflag = true;
            
}
        
 
if (skipsequenceflag)
            continue;
        
 
newsum += -log (alpha[0][seqlen - 1][0]);    //checking for threshhold
        
_beta.initialize_beta ();
        
_beta.calc_beta ();
        
_reestimate.updateB ();
        
_reestimate.updateA ();
        
}
    
 
cout << endl << reject << " Sequences Rejected !! ";

}


 
void
IOcontrol::finalcalculator (void) 
{
    
probreader ();
    
seqreader ();
    
 
ofstream fout (outprobfile);
    
if (!fout.good ())
        
 {
        
cout << "\n Unable to open Output Probability File \n";
        
exit (1);
        
}
    
 
 
long double delta = 100, temp;
    
 
int iter = 0;
    
int v, y, z, t;
    
bool terminatorflag = false;
    
 
    do
        
 {
        
            //initialize, estimate and calc newsum and oldsum, check delta percentage
            
(void) time (&t3);
        
            // start clocking current iteration
            
++iter;
        
initializer ();
        
 
oldsum = newsum;
        
newsum = 0;
        
 
estimator ();        //calculates newsum = SIGMA -log(seq. probs)
        
cout << endl <<
            "______________________________________________________________"
            << endl;
        
cout << "\n OLD SCORE = " << oldsum << "\t NEW SCORE = " << newsum;
        
cout << endl <<
            "______________________________________________________________"
            << endl;
        
 
for (v = 0; v < NSYM; ++v)
            
for (y = 0; y < NSYM; ++y)
                
for (z = 0; z < NSYM; ++z)
                    
if (valida[v][y][z])
                        
 {
                        
temp = Na[v][y][z] / Da[v][y][z];
                        
 
if (temp < APPROX_ZERO)
                            
temp = 0;
                        
                            /*
                               if(temp > a[v][y][z])
                               newsum += (temp - a[v][y][z]);
                               else
                               newsum += (a[v][y][z] - temp);
                             */ 
                            
a[v][y][z] = temp;
                        
cout << "\n a " << v << " " << y << " " << z <<
                            " " << a[v][y][z];
                        
}
        
 
for (v = 0; v < NSYM; ++v)
            
for (t = 0; t < TSYM; ++t)
                
if (validb[v][t])
                    
 {
                    
temp = Nb[v][t] / Db[v][t];
                    
 
if (temp < APPROX_ZERO)
                        
temp = 0;
                    
 
                        /*
                           if(temp > b[v][t])
                           newsum += (temp - b[v][t]);
                           else
                           newsum += (b[v][t] - temp);
                         */ 
                        
b[v][t] = temp;
                    
cout << "\n b " << v << " " << t << " " << b[v][t];
                    
}
        
 
cout << endl <<
            "______________________________________________________________";
        
cout << endl << "\n ITERATION # " << iter << "\t Score = " <<
            (newsum);
        
 
(void) time (&t4);
        
            //End clocking current Iteration
            
if ((iter == 1) || (iter == 2))
            
 {
            
cout << "\n TIME TAKEN : ";
            
givetime ((long int) (t4 - t3));
            
cout <<
                "______________________________________________________________"
                << endl;
            
 
continue;
            
}
        
 
 
if (newsum == oldsum)
            
delta = 0;
        
        else
            
 {
            
delta = ((newsum - oldsum) / oldsum) * 100.0;
            
if (delta < 0)
                delta = -delta;
            
}
        
 
cout << "\t Delta = " << delta << " % ";
        
 
cout << "\n TIME TAKEN : ";
        
givetime ((long int) (t4 - t3));
        
cout <<
            "______________________________________________________________"
            << endl;
        
 
if (delta <= THRESH_LIMIT)
            
terminatorflag = true;
        
 
}                    //while(iter < 10);
    while (!terminatorflag);
    
 
        // Writing Output to file
        
for (v = 0; v < NSYM; ++v)
        
for (y = 0; y < NSYM; ++y)
            
for (z = 0; z < NSYM; ++z)
                
if (valida[v][y][z])
                    
 {
                    
fout << "a " << v << " " << y << " " << z << " " <<
                        a[v][y][z] << endl;
                    
}
    
 
for (v = 0; v < NSYM; ++v)
        
for (t = 0; t < TSYM; ++t)
            
if (validb[v][t])
                
 {
                
fout << "b " << v << " " << t << " " << b[v][t] << endl;
                
}
    
 
 
fout.close ();
    
cout << endl << " || END OF PROGRAM || " << endl;

}


 
 
//int main()
    int
main (int argc, char *argv[]) 
{
    
 
(void) time (&t1);
    
        //Start Clocking;
        
 
if (argc < 4)
        
 {
        cout << endl << "INSUFFICIENT ARGMENTS" << endl 
            <<"Format is : Program <input sequence file> <input prob. file> <output prob. file> "
            
 <<endl;
        
return 1;
        
}
    
 
strcpy (inseqfile, argv[1]);
    
strcpy (inprobfile, argv[2]);
    
strcpy (outprobfile, argv[3]);
    
 
/*
	strcpy(inseqfile,"seq_test.txt");
	strcpy(inprobfile,"inprob_test.txt");
	strcpy(outprobfile,"outprob_test.txt");
*/ 
        
IOcontrol c;
    
c.finalcalculator ();
    
 
(void) time (&t2);
    
        //End Clocking;
        cout << "\n TOTAL TIME TAKEN BY PROGRAM : \n";
    
givetime ((long int) (t2 - t1));
    
 
return 0;

}


 
 
void
givetime (long int t_in_sec) 
{
    
cout << "\t || " << (t_in_sec / 3600) << " HOURS - " 
        <<((t_in_sec % 3600) / 60) << " MINUTES - " << (t_in_sec %
                                                        60) <<
        " SECONDS ||\n";

} 
 
