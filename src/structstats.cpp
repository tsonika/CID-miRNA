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




//********************* Code for Printing the Parse Tree Ends ***********************

//************************ Code for Building nodes Begins ***************************
node *BuildNodes (string & trace)
{
    node *finaltree = new node;
    finaltree->initialize (ROOT);
    node *currnode = finaltree;
    node *a;
    istringstream iss (trace);
    string currline, onlynode;
    int findpos;
    bool LoopOccured = false;

    while (!iss.eof () && !LoopOccured) {
        getline (iss, currline);
        if (currline.length () == 0)
            continue;
        a = new node;

        onlynode = currline.substr (currline.find ("-->") + 4);

        findpos = onlynode.find ("STEM");
        if (findpos != string::npos) {  //Make node for STEM
            a->initialize (STEM);
            stem *qstem = new stem;
            qstem->l = onlynode[0];

            findpos = onlynode.find (" ", findpos);
            qstem->r = onlynode[findpos + 1];
            a->storedata (qstem);
            onlynode = "";
        }                       //End- STEM

        findpos = onlynode.find ("ASYM");
        if (findpos != string::npos) {  //Make node for ASYM
            a->initialize (ASYM);
            asym *qasym = new asym;
            qasym->l = onlynode.substr (0, findpos);

            findpos = onlynode.find (" ", findpos);
            qasym->r = onlynode.substr (findpos + 1);

            findpos = (qasym->l).find (" ");
            while (findpos != string::npos) {
                (qasym->l).erase (findpos, 1);
                findpos = (qasym->l).find (" ");
            }

            findpos = (qasym->r).find (" ");
            while (findpos != string::npos) {
                (qasym->r).erase (findpos, 1);
                findpos = (qasym->r).find (" ");
            }

            a->storedata (qasym);
            onlynode = "";
        }                       //End- ASYM

        findpos = onlynode.find ("SYM");
        if (findpos != string::npos) {  //Make node for SYM
            a->initialize (SYM);
            sym *qsym = new sym;
            qsym->l = onlynode.substr (0, findpos);

            findpos = onlynode.find (" ", findpos);
            qsym->r = onlynode.substr (findpos + 1);

            findpos = (qsym->l).find (" ");
            while (findpos != string::npos) {
                (qsym->l).erase (findpos, 1);
                findpos = (qsym->l).find (" ");
            }

            findpos = (qsym->r).find (" ");
            while (findpos != string::npos) {
                (qsym->r).erase (findpos, 1);
                findpos = (qsym->r).find (" ");
            }

            a->storedata (qsym);
            onlynode = "";
        }                       //End- SYM

        while (isspace (onlynode[0]))
            onlynode.erase (0, 1);      //Remove all leading spaces

        findpos = onlynode.find ("LOOP");
        if (findpos != string::npos && findpos > 0) {   //Make node for LOOP, if this point is like the STEM
            a->initialize (STEM);
            stem *qstem = new stem;
            qstem->l = onlynode[0];

            findpos = onlynode.find (" ", findpos);
            qstem->r = onlynode[findpos + 1];
            a->storedata (qstem);

            //This portion fills in this special point into the node tree
            if (currnode == finaltree)
                finaltree->storep (a);
            else
                currnode->storep (a);
            currnode = a;
            a = new node;
        }                       //End- LOOP, like the STEM

        findpos = onlynode.find ("LOOP");
        if (findpos != string::npos) {  //Make node for LOOP, if this point has the LOOP string
            a->initialize (LOOP);
            loop *qloop = new loop;
            string qloops;
            getline (iss, qloops);
            qloops = qloops.substr (qloops.find ("-->") + 4);
            qloop->s = qloops;

            findpos = (qloop->s).find (" ");
            while (findpos != string::npos) {
                (qloop->s).erase (findpos, 1);
                findpos = (qloop->s).find (" ");
            }
            a->storedata (qloop);
            LoopOccured = true;
        }                       //End- LOOP

        //Now 'a' has the current node
        if (currnode == finaltree)
            finaltree->storep (a);
        else
            currnode->storep (a);

        currnode = a;
    }                           //End- while()

    return finaltree;
}

//************************* Code for Building nodes Ends ****************************
int main (int argc, char *argv[])
{
//Recieve input from command line
    if (argc < 3) {
        cout << "ERROR!! Insufficient arguments.\n";
        cout << "Usage Format:\n";
        cout << "<Program Name> <Grammar File Name> <Sequence to Parse>\n";
        exit (1);
    }

    char *GFileName = argv[1];
    Sequence = argv[2];

//Read Grammar from file into the vector RuleList
    ifstream *gfile;
    gfile = new ifstream;
    gfile->open (GFileName, ios::in);
    RuleList = ReadGrammarFile (gfile);
    gfile->close ();
    delete gfile;

//The following portion of the program displays the grammer as read from the file
//A sample output can be seen in RuleReference.txt (Grammar file used was 'result.txt')
/*
cout<<"Printing the values of various variables in RuleList...\n\n";

cout<<"STATIC MEMBERS:\n";
cout<<"NTcount: "<<RuleList[0].NTcount<<endl;
cout<<"Tcount: "<<RuleList[0].Tcount<<endl;
cout<<"RuleCount: "<<RuleList[0].RuleCount<<endl;
cout<<"NT: "<<RuleList[0].NT<<endl;
cout<<"T: "<<RuleList[0].T<<endl;
cout<<"NTonLHS: "<<RuleList[0].NTonLHS<<endl<<endl;
void FillFirstRow(vector<CykMatrixCell> &mat)//This function accepts a reference to a vector of type CykMatrixCell
cout<<"VARYING MEMBERS:\n";
for (int x=0;x<RuleList[0].RuleCount;x++)
    {
    cout<<"RuleNumber: "<<RuleList[x].RuleNumber<<endl;
    cout<<"RuleScore: "<<RuleList[x].RuleScore<<endl;
    cout<<"Rule:\truleLHS --> ruleRHS\n"<<RuleList[x].ruleLHS<<" --> "<<RuleList[x].ruleRHS<<endl;
    cout<<"UsedOnceFlag: "<<RuleList[x].UsedOnceFlag<<endl;
    cout<<"No_of_RHSsymbols: "<<RuleList[x].No_of_RHSsymbols<<endl;
    cout<<"RuleType: "<<RuleList[x].RuleType<<"\n\n\n";
    }//*/


//Make the matrix as a vector
    SeqLen = Sequence.length ();
    vector < CykMatrixCell > UTMatrix;
    FillFirstRow (UTMatrix);
    FillRestMatrix (UTMatrix);

//The following portion of the program displays the results of the filled Matrix
//A sample output can be seen in MatrixRef.txt (Grammar file used: smade02.txt, sequence 'aauggau')
/*
int x,y,abspos,itr;
for (abspos=0;abspos < UTMatrix.size(); abspos++)
    {
    XYPos(abspos,x,y);

    cout<<"\n\nCell Position: "<<x<<","<<y<<" = Absolute Position "<<abspos<<endl;
    cout<<"Sub String Represented: "<<Sequence.substr(y-1,x)<<"\n\n";

    cout<<"\tUpListCount :"<<UTMatrix[abspos].UpListCount<<endl;
    cout<<"\tDownListCount :"<<UTMatrix[abspos].DownListCount<<endl;

    cout<<"\tUpList:\n";
    for (itr=0;itr<UTMatrix[abspos].UpList.size();itr++)
        {
        cout<<"\t\truleLHS: "<<UTMatrix[abspos].UpList[itr].ruleLHS<<endl;
        cout<<"\t\truleRHSpreDot: "<<UTMatrix[abspos].UpList[itr].ruleRHSpreDot<<endl;
        cout<<"\t\truleRHSpostDot: "<<UTMatrix[abspos].UpList[itr].ruleRHSpostDot<<endl;
        cout<<"\t\tscore: "<<UTMatrix[abspos].UpList[itr].score<<endl;
        cout<<"\t\tcoordinates:\n";
        cout<<"\t\t\trow: "<<UTMatrix[abspos].UpList[itr].coordinates.row<<endl;
        cout<<"\t\t\tcol: "<<UTMatrix[abspos].UpList[itr].coordinates.col<<endl;
        cout<<"\t\t\tindex: "<<UTMatrix[abspos].UpList[itr].coordinates.index<<endl;
        cout<<"\t\t\ttype: "<<UTMatrix[abspos].UpList[itr].coordinates.type<<endl;
        cout<<"\t\tParents: "<<UTMatrix[abspos].UpList[itr].Parents<<"\n\n";
        }

    cout<<"\tDownList:\n";
    for (itr=0;itr<UTMatrix[abspos].DownList.size();itr++)
        {
        cout<<"\t\truleLHS: "<<UTMatrix[abspos].DownList[itr].ruleLHS<<endl;
        cout<<"\t\truleRHSpreDot: "<<UTMatrix[abspos].DownList[itr].ruleRHSpreDot<<endl;
        cout<<"\t\truleRHSpostDot: "<<UTMatrix[abspos].DownList[itr].ruleRHSpostDot<<endl;
        cout<<"\t\tscore: "<<UTMatrix[abspos].DownList[itr].score<<endl;
        cout<<"\t\tcoordinates:\n";
        cout<<"\t\t\trow: "<<UTMatrix[abspos].DownList[itr].coordinates.row<<endl;
        cout<<"\t\t\tcol: "<<UTMatrix[abspos].DownList[itr].coordinates.col<<endl;
        cout<<"\t\t\tindex: "<<UTMatrix[abspos].DownList[itr].coordinates.index<<endl;
        cout<<"\t\t\ttype: "<<UTMatrix[abspos].DownList[itr].coordinates.type<<endl;
        cout<<"\t\tParents: "<<UTMatrix[abspos].DownList[itr].Parents<<"\n\n";

        }
    }//End- Printing of Filled Matrix*/

//Checking if the input sequence can be parsed using the grammar, and printing the tree if it can be parsed through the grammar
    if (IsPossibleTree (UTMatrix)) {
        cout << "\nThe sequence CAN be parsed using the grammar!!\n\n";
        string trace;
        trace = BuildTrace (UTMatrix);  //'trace' contains the output stream for the BuildTrace(); This is the same as the output seen on screen in the previous version
        cout << trace;

        node *RootNode = new node;      //Create a new node
        RootNode = BuildNodes (trace);  //Calling local function BuildNodes() to build the nodes that can be recognised by the fillmatrix() code of 'drawmirna.h'

        initmatrix ();
        fillmatrix (RootNode, 0);

        string miRNAstruct;
        miRNAstruct = returnmatrix ();
        cout << miRNAstruct;

        miRNAstats structstats;
        structstats = getstats (miRNAstruct);

        cout << "Stats for the above structure:\n";
        cout << "Total number of Continuous Stems :" << structstats.n_stems;
        cout << "\nTotal number of Asymmetric Bulges :" << structstats.
            n_asyms;
        cout << "\nTotal number of Symmetric Bulges :" << structstats.n_syms;
        cout << "\nTotal number of Bases :" << structstats.n_bases;
        cout << "\nTotal number of Base Pairs :" << structstats.n_base_pairs;
        cout << "\nTotal number of Bases in Loop :" << structstats.
            n_loop_bases << "\n\n";



//   displaymatrix();
    } else
        cout << "\n\nThe sequence CANNOT be parsed using the grammar!!\n\n";

    return 0;
}
