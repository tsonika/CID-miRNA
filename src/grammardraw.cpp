//**********************************************************************************
// File         : grammardraw.cpp                                                                               
// Purpose      : Definition of functions for reading and storing of MiRNA General SCFG 
// Authors:                   Vipin Gupta(viping@gmail.com),
//(Alphabetically,Last Name)  Vishal Kapoor(vkapoor@cse.iitd.ernet.in),
//                                                        Manish Kushwaha(manishkushwaha@gmail.com),
//                                                        Vasudha Verma(vasudha.verma@gmail.com)
//      :1.0    -       26-Feb-2005     Initial Draft (Vipin Gupta)                     
//              :2.0    -       10-Mar-2005     Grammar Reading and Error Handling perfected (Vipin Gupta)
//              :3.0    -       03-Jun-2005     Minor Upgrade in Grammer Reading and Validation functions,
//                                              Filling the first row of a CYK Upper Triangular Matrix (Manish Kushwaha)
//              :4.0    -       14-Jun-2005     Filling the rest of the rows of the CYK Upper Triangular Matrix and Tracing it to print a human-readable tree (Manish Kushwaha)
//              :5.0    -       15-Jun-2005     Included 'drawmirna.h' functionality to display graphical output
//              :5.5    -       25-Jun-2005     Removed the bias in asymmetric bulge display
//      :6.0    -   29-Jun-2005 Optimization with introducing associative containers,removal of dead code
//      :8.0    -       02-Jul-2005 Final draft with HUGE bug in osstreamdtring removed,Code now runs faster thatn before
//**********************************************************************************
#pragma warning (push)          // store the warning level
#pragma warning (disable:4786)  // identifier was truncated to '255' characters in the debug information
#include "grammar_reader.h"
#include "cykmatrix.h"
#include "drawmirna.h"


int main ()
{

    char *GFileName = "vinfinalgrammar.txt";
    Sequence =
        "ucauugguccagaggggagauagguuccugugauuuuuccuucuucucuauagaauaaauga";
//Read Grammar from file into the vector RuleList
    ifstream *gfile;
    gfile = new ifstream;
    gfile->open (GFileName, ios::in);
    RuleList = ReadGrammarFile (gfile);
    gfile->close ();
    delete gfile;

    string s1, s2;
    int len;
    for (int x = 0; x < RuleList[0].RuleCount; x++) {
        s1 = RuleList[x].ruleRHS;
        if (RuleList[x].RuleType == 1) {
            Tmap.insert (pair < string,
                         int >(RuleList[x].ruleRHS, RuleList[x].RuleNumber));

        } else if (RuleList[x].RuleType == 2) {
            NTmap.insert (pair < string,
                          int >(RuleList[x].ruleRHS, RuleList[x].RuleNumber));
            for (int ss = 1; ss < RuleList[x].No_of_RHSsymbols; ss++) {

                len = s1.length ();
                s1.erase (len - 1, 1);
                int delimiter = s1.rfind ("*");
                s2 = s1.substr (0, delimiter);
                s1 = s2;
                NTmapfirst.insert (pair < string,
                                   int >(s2, RuleList[x].RuleNumber));
            }
        }
    }

//Make the matrix as a vector
    SeqLen = Sequence.length ();
    vector < CykMatrixCell > UTMatrix;
    FillFirstRow (UTMatrix);
    FillRestMatrix (UTMatrix);


//Checking if the input sequence can be parsed using the grammar, and printing the tree if it can be parsed through the grammar
    if (IsPossibleTree (UTMatrix)) {
        cout << "\nThe sequence CAN be parsed using the grammar!!\n\n";
        string trace;
        trace = BuildTrace (UTMatrix);  //'trace' contains the output stream for the BuildTrace(); This is the same as the output seen on screen in the previous version
        cout << trace;

        node *RootNode = new node (STEM);       //Create a new node
        RootNode = BuildNodes (trace);  //Calling local function BuildNodes() to build the nodes that can be recognised by the fillmatrix() code of 'drawmirna.h'

        initmatrix ();
        fillmatrix (RootNode, 0);
        displaymatrix ();
    } else
        cout << "\n\nThe sequence CANNOT be parsed using the grammar!!\n\n";

    return 0;
}
