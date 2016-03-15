//**********************************************************************************
// File         : structstats.cpp                                                                               
// Purpose      : Definition of functions for reading and storing of MiRNA General SCFG 
// Authors      : Vipin Gupta (viping@gmail.com), Manish Kushwaha (manish.kushwaha@gmail.com)
// Revision     :1.0    -       26-Feb-2005     Initial Draft (Vipin Gupta)                     
//              :2.0    -       10-Mar-2005     Grammar Reading and Error Handling perfected (Vipin Gupta)
//              :3.0    -       03-Jun-2005     Minor Upgrade in Grammer Reading and Validation functions,
//                                              Filling the first row of a CYK Upper Triangular Matrix (Manish Kushwaha)
//              :4.0    -       14-Jun-2005     Filling the rest of the rows of the CYK Upper Triangular Matrix and Tracing it to print a human-readable tree (Manish Kushwaha)
//              :5.0    -       15-Jun-2005     Included 'drawmirna.h' functionality to display graphical output
//              :5.5    -       25-Jun-2005     Removed the bias in asymmetric bulge display
//              :6.0    -       25-Jun-2005     Return statistics of a given structure
//**********************************************************************************

#include "grammar_reader.h"
#include "cykmatrix.h"
#include "drawmirna.h"



int main (int argc, char *argv[])
{
//Recieve input from command line
    if (argc < 3) {
        cout << "ERROR!! Insufficient arguments.\n";
        cout << "Usage Format:\n";
        cout <<
            "<Program Name> <Input File Name with structures> <Input File Name with details>\n";
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

    string inpline, inpline2, miRNAstruct;
    bool tracking = false;
    miRNAstats structstats;

    while (!infile->eof ()) {   //Run through the entire input file
        getline (*infile, inpline);     //Read one line of the file in 'inpline'
        inpline2 = inpline;

        while (isspace (inpline2[0]))
            inpline2.erase (0, 1);      //Removes all leading spaces

        if (inpline2.find (">") == 0) { //miRNA name line
            if (tracking) {     //End and write an miRNA structure
                structstats = getstats (miRNAstruct);
                *outfile << miRNAstruct;

                *outfile << "Number of Continuous Stems :" << structstats.
                    n_stems;
                *outfile << "\nNumber of Asymmetric Bulges :" << structstats.
                    n_asyms;
                *outfile << "\nNumber of Symmetric Bulges :" << structstats.
                    n_syms;
                *outfile << "\nNumber of Bases :" << structstats.n_bases;
                *outfile << "\nNumber of Base Pairs :" << structstats.
                    n_base_pairs;
                *outfile << "\nNumber of Bases in the Terminal Loop :" <<
                    structstats.n_loop_bases << "\n\n";
            }                   //End- End and write an miRNA structure //Initiate an miRNA structure

            //Initiate an miRNA structure
            miRNAstruct = "";
            tracking = true;
            *outfile << inpline2 << endl;
        }                       //End- miRNA name line
        else {
            miRNAstruct += inpline + "\n";
        }
    }                           //End- Run through the entire input file - while

//Write the last miRNA structure
    if (miRNAstruct.length () != 0) {
        structstats = getstats (miRNAstruct);
        *outfile << miRNAstruct;

        *outfile << "Number of Continuous Stems :" << structstats.n_stems;
        *outfile << "\nNumber of Asymmetric Bulges :" << structstats.n_asyms;
        *outfile << "\nNumber of Symmetric Bulges :" << structstats.n_syms;
        *outfile << "\nNumber of Bases :" << structstats.n_bases;
        *outfile << "\nNumber of Base Pairs :" << structstats.n_base_pairs;
        *outfile << "\nNumber of Bases in the Terminal Loop :" << structstats.
            n_loop_bases << "\n\n";
    }
//End- Write the last miRNA structure


    outfile->close ();
    infile->close ();
    delete outfile;
    delete infile;

    return 0;
}
