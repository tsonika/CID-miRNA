//**********************************************************************************
// File         : Stats4mStruct.cpp                                                                             
// Purpose      : Reading an miRNA structure file and outputting the characteristics of the structure in terms of some handy structural details 
// Author       : Manish Kushwaha (manish.kushwaha@gmail.com)
// Revision     : 1.0   -       25-Jun-2005     Working  (Manish Kushwaha)                      
//**********************************************************************************
#include <fstream>
#include <iostream>
#include <list>
#include <vector>

#include <sstream>
#include <stdlib.h>


#include "constants.h"
#include "mirnastats.h"

using namespace std;

char matrix[5][81];             //This matrix represents the screen of 5x80 (Height x Width) characters


double score_mirna(string& miRNAstruct) {
    miRNAstats structstats = getstats (miRNAstruct);
    double structscore = 0;

    structscore = structstats.n_base_pairs * 2;     //(Number of paired bases)*(+1)

    for (size_t x = 0; x < structstats.sym_len.size (); x++)
        if (structstats.sym_len[x] > 2)
            structscore += (structstats.sym_len[x] * 2) * (-2);     //(Number of unpaired bases in SYM of size>2)*(-2)
        else
            structscore += structstats.sym_len[x];  //(Number of unpaired bases in SYM of size<=2)*(+0.5)

    for (size_t x = 0; x < structstats.asym_nbases.size (); x++)
        structscore += structstats.asym_nbases[x] * (-2);   //(Number of bases in ASYM)*(-2)

    if (structstats.n_loop_bases > 5)
        structscore += (5 - structstats.n_loop_bases);      //(Number of bases in terminal loop of size >5)*(-1)

    return structscore;
}

int main (int argc, char *argv[])
{
//Recieve input from command line
    if (argc < 3) {
        cout << "ERROR!! Insufficient arguments.\n";
        cout << "Usage Format:\n";
        cout <<
            "<Program Name> <Input File Name with structures> <Output File Name with details>\n";
        exit (1);
    }

    char *InFileName = argv[1];
    char *OutFileName = argv[2];

    //Open input and output files for reading
    ifstream *infile;
    infile = new ifstream;
    infile->open (InFileName, ios::in);
    if (!infile->good ()) {
        cout << "ERROR!! Could not open INPUT file.\n\n";
        exit (1);
    }
    ofstream *outfile;
    outfile = new ofstream;
    outfile->open (OutFileName, ios::out);
    if (!outfile->good ()) {
        cout << "ERROR!! Could not open OUTPUT file.\n\n";
        exit (1);
    }

    string inpline, inpline2, inpline3, miRNAstruct;
    bool tracking = false;

    while (!infile->eof ()) {
        //Run through the entire input file
        getline (*infile, inpline);     //Read one line of the file in 'inpline'
        inpline2 = inpline;

        while (isspace (inpline2[0]))
            inpline2.erase (0, 1);      //Removes all leading spaces

        if (inpline2.empty () && !tracking)
            continue;           // Skip empty lines, but only until we are into the structures section

        if (inpline2.find (">") == 0) { //miRNA name line
            if (tracking) {     //End and write an miRNA structure
                *outfile << miRNAstruct;        //This line prints the original structure to the output file
                *outfile << "SCORE: " << score_mirna(miRNAstruct) << "\n";    
            }                   

            //Initiate an miRNA structure
            miRNAstruct = "";
            tracking = true;
            *outfile << inpline2 << endl;
            getline (*infile, inpline3);
            *outfile << inpline3 << endl;
        }                       //End- miRNA name line
        else {
            miRNAstruct += inpline + "\n";
        }
    }                           //End- Run through the entire input file - while

    //Write the last miRNA structure
    if (miRNAstruct.length () != 0) {
        *outfile << miRNAstruct;        //This line prints the original structure to the output file
        *outfile << "SCORE: " << score_mirna(miRNAstruct) << "\n";    
    }


    outfile->close ();
    infile->close ();
    delete outfile;
    delete infile;

    return 0;
}
