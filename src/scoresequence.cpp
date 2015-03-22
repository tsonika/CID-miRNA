/*
Scores RNA sequences based on a stochastic context-free grammar model
Scoring is done by a CYK algorithm (Cocke-Younger-Kasami)
*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <cmath>
#include <ctime>

using namespace std;

#define MAX_SEQUENCE_LENGTH 130
#define MAX_INPUT_LINE_LENGTH 50000
#define NONTERMINAL_SYMBOLS 60
#define TERMINAL_SYMBOLS 4
#define MINUSINF -1e90
#define START_NONTTERMINAL 0
#define MAX_PATH 5000

struct entry {
    long non_terminal;
    double score;
};

struct cell {

    long members;           //basically stores the number of members in the array
    entry node[NONTERMINAL_SYMBOLS];

    void operator= (cell &xyz);
};

void cell::operator=(cell &xyz)
{
    long eqi;
    members = xyz.members;
    for(eqi = 0 ; eqi < NONTERMINAL_SYMBOLS; ++eqi)
        node[eqi] = xyz.node[eqi];
}

cell rule_transitions[NONTERMINAL_SYMBOLS][NONTERMINAL_SYMBOLS];
cell terminal_emissions[TERMINAL_SYMBOLS];

cell matrix[MAX_SEQUENCE_LENGTH][MAX_SEQUENCE_LENGTH];



// This will convert the character sequence longo an equivalent
// numeric sequence by converting each Terminal longo equivalent numeric code
// returns true if successful conversion else returns false
bool convert_sequence(const char* sequence, long* converted_sequence)
{
    long len = strlen(sequence);
    if (len <= 0) return false;

    for(int j = 0; j < len; ++j)
    {
            if((sequence[j] == 'a') || (sequence[j] == 'A')) converted_sequence[j] = 0;
            else if((sequence[j] == 'u') || (sequence[j] == 'U')) converted_sequence[j] = 1;
            else if((sequence[j] == 'g') || (sequence[j] == 'G')) converted_sequence[j] = 2;
            else if((sequence[j] == 'c') || (sequence[j] == 'C')) converted_sequence[j] = 3;
            else
            {
                cerr << endl << "Seq # " << sequence << endl << " Illegal char found. Rejected";
                return false;

            }
    }
    return true;
}

/* 
Implement the CYK (Cocke-Younger-Kasami) for parsing stochastic context-free grammars as described in "Stochastic Context-Free Grammars and RNA Secondary Structure Prediction" by Martin Knudsen - http://cs.au.dk/~cstorm/students/Knudsen_Jun2005.pdf, page 28.
*/

double cyk_calculator(const char* sequence)
{

    long i, j, p;
    long n1, n2, n3, n4;
    long nrule;
    bool full_sequence_scored = false;
    bool found = false;
    double score = 0.0;
    cell *x, *y, *z;    

    long sequence_length = strlen(sequence);
    if(sequence_length <= 0) return MINUSINF;

    long converted_sequence[MAX_SEQUENCE_LENGTH];
    if(!convert_sequence(sequence, converted_sequence)) return MINUSINF;

    for(j = 0; j < sequence_length; ++j)
        matrix[0][j] = terminal_emissions[converted_sequence[j]];

    // I have initialized the first row according to each terminal


    // From the comment in the algorithm desciption:
    // i is the length of the span, j the start and p where to split into two subspans

    for(i = 1; i < sequence_length; ++i)
    {
        for(j = 0; j < (sequence_length - i); ++j)
        {
            matrix[i][j].members = 0;
            for(p = 0; p < i; ++p)
            {
                // Taking cartesian product of matrix[p][j] and matrix[i-p-1][j+p+1]

                x = &matrix[p][j];
                y = &matrix[i-p-1][j+p+1];

                for(n1 = 0; n1 < x->members; ++n1)
                {
                    for(n2 = 0; n2 < y->members; ++n2)
                    {
                        z = &rule_transitions[x->node[n1].non_terminal][y->node[n2].non_terminal];

                        for(n3 = 0; n3 < z->members; ++n3)
                        {
                            found = false;
                            score = x->node[n1].score + y->node[n2].score + z->node[n3].score;
                            nrule = z->node[n3].non_terminal;

                            for(n4 = 0; n4 < matrix[i][j].members; ++n4)
                            {
                                if(nrule == matrix[i][j].node[n4].non_terminal)
                                {
                                    found = true;
                                    if (score > matrix[i][j].node[n4].score)
                                    {
                                        matrix[i][j].node[n4].score = score;
                                    }
                                    break;
                                }
                            }

                            if(!found)
                            {
                                // adding new entry
                                matrix[i][j].node[matrix[i][j].members].non_terminal = nrule;
                                matrix[i][j].node[matrix[i][j].members++].score = score;
                            }

                        } // n3 loop
                    }
                }
            }
        }
    }

    full_sequence_scored = false;
    for(n1 = 0; n1 < matrix[sequence_length-1][0].members; ++n1)
    {
        if(matrix[sequence_length-1][0].node[n1].non_terminal == START_NONTTERMINAL)
        {
            score = matrix[sequence_length-1][0].node[n1].score;
            full_sequence_scored = true;
            break;
        }
    }

    if (full_sequence_scored) return score;

    return MINUSINF;
}





// This function has to first reset everything and then
// Fill up Hashrule 2-D matrix and
// Fill up the hash matrix for Terminals
void initial_value_reader(const char* input_probability_filename)
{

    long vi , j;

    for(vi = 0; vi < MAX_SEQUENCE_LENGTH; ++vi)
        for(j = 0 ; j < MAX_SEQUENCE_LENGTH; ++j)
            matrix[vi][j].members = 0;

    for(vi = 0; vi < NONTERMINAL_SYMBOLS; ++vi)
        for(j = 0; j < NONTERMINAL_SYMBOLS; ++j)
            rule_transitions[vi][j].members = 0;

    for(vi = 0; vi < TERMINAL_SYMBOLS; ++vi)
        terminal_emissions[vi].members = 0;


    ifstream inprobfile;
    inprobfile.open(input_probability_filename, ios::in);

    char ch = 0;
    long temp;
    long pm[3] = {0,0,0};
    long double val = 0.0;
    char tempstream[MAX_PATH+1];

    if(!inprobfile.good())
    {
        cerr << "\n Error Opening Input Probability File !! \n";
        exit(1);
    }

    while(!inprobfile.eof())
    {
        inprobfile.getline(tempstream,MAX_PATH);
        if(5 > strlen(tempstream)) continue;
        //just checking for basic length of one input line (can't be less than 5)

        istringstream iss(tempstream);

        if(!(iss>>ch)) continue;

        if((ch == 'a') || (ch == 'A'))
        {
            // 'a' lines describe probability of transitions from non-terminals to non-terminals
            // eg a 57 51 37 0.00120985 means
            // the probability of going from non-terminal 37 to non-terminals 57 and 51 is 0.001209...
            iss>>pm[0];
            iss>>pm[1];
            iss>>pm[2];
            iss>>val;
            if((pm[0] < NONTERMINAL_SYMBOLS) &&(pm[1] < NONTERMINAL_SYMBOLS) &&(pm[2] < NONTERMINAL_SYMBOLS))
            {
                temp = rule_transitions[pm[1]][pm[2]].members++;
                rule_transitions[pm[1]][pm[2]].node[temp].non_terminal = pm[0];
                rule_transitions[pm[1]][pm[2]].node[temp].score = log10(val);
            }
        }
        else if((ch == 'b') || (ch == 'B'))
        {
            // 'b' lines describe probabilities of emitting terminals
            // eg b 51 0 0.243964
            // is the probability of emitting terminal '0' from non-terminal 51 (0.243...)

            iss>>pm[0];
            iss>>pm[1];
            iss>>val;

            if((pm[0] < NONTERMINAL_SYMBOLS) &&(pm[1] < TERMINAL_SYMBOLS))
            {
                temp = terminal_emissions[pm[1]].members++;
                terminal_emissions[pm[1]].node[temp].non_terminal = pm[0];
                terminal_emissions[pm[1]].node[temp].score = log10(val);
            }
        }

        else if(ch == '#') continue;
            //Comment Symbol

        else if(!inprobfile.eof())
        {
            cerr << "\n Unable to parse Input Probability File" << "\n Check the Format";
            inprobfile.close();
            exit(1);
        }

    }
    inprobfile.close();
}


// give time in human readable format
void print_time (long int t_in_sec)
{
    cout << "\t || "<<(t_in_sec / 3600)<<" HOURS - "
     << ((t_in_sec % 3600)/60) <<" MINUTES - "<<(t_in_sec % 60)<<" SECONDS ||\n";
}


int main(int argc, char* argv[])
{

    time_t t1,t2;
    char input_probability_filename[MAX_PATH];
    char input_sequence_filename[MAX_PATH];
    char result_filename[MAX_PATH];


    if(argc < 4) { 
        cerr  << endl << "INSUFFICIENT ARGMENTS" << endl
            << "Format is : Program <input sequence file> <input prob. file> <output result file> "
            << endl;
        return 1;
    }

    if (strlen(argv[1]) >= MAX_PATH - 1) {
        cerr << "The sequence file name is too long" << endl;
        return 1;
    }

    if (strlen(argv[2]) >= MAX_PATH - 1) {
        cerr << "The input probability file name is too long" << endl;
        return 1;
    }

    if (strlen(argv[3]) >= MAX_PATH - 1) {
        cerr << "The output result file name is too long" << endl;
        return 1;
    }

    strcpy(input_sequence_filename,argv[1]);
    strcpy(input_probability_filename,argv[2]);
    strcpy(result_filename,argv[3]);

    initial_value_reader(input_probability_filename);

    ofstream resultfile;
    resultfile.open(result_filename);

    if(!resultfile.good())
    {
        cerr << endl << " ERROR OPENING RESULT FILE !! " << endl;
        exit(1);
    }

    ifstream inseqfile;
    inseqfile.open(input_sequence_filename);
    if(!inseqfile.good())
    {
        cerr << "\n Error Opening Input Sequence File !! \n";
        exit(1);
    }


    (void)time(&t1);
    //Start Clocking;

    long counter = 0;
    long line_length;
    char input_line[MAX_INPUT_LINE_LENGTH+1];
    double myscore;

    while(!inseqfile.eof())
    {
        inseqfile.getline(input_line, MAX_INPUT_LINE_LENGTH+1);
        counter++;

        if (inseqfile.fail()) {
            if (!inseqfile.eof()) {
                // Something went wrong. Most likely that line was too long
                cerr << "Something went wrong reading line " << counter << ". Probably too long. Bailing" << endl;
            }
            break;
        }

        line_length  = strlen(input_line);
        if((line_length <= 1) || strstr(input_line, ">")) continue;

        if (line_length > MAX_SEQUENCE_LENGTH) {
            cerr << "Line " << counter << " too long. Skipping" << endl;            
            continue;
        }
        
        myscore = cyk_calculator(input_line);
        
        resultfile << "\n Sequence : " << counter << endl << input_line
        << "\n Length : " << line_length
        << "\t Normal SCORE = " << myscore/line_length
        << "\t SCORE = " << myscore;

    }

    inseqfile.close();

    (void)time(&t2);
    //End Clocking;
    cout << "\n TOTAL TIME TAKEN BY PROGRAM : \n";
    print_time((long int)(t2-t1));

    resultfile.close();
    return 0;

}
