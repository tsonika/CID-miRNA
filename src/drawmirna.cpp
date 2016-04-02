#include <fstream>
#include <iostream>
#include <list>
#include <vector>

#include <sstream>
#include <algorithm>
#include <map>

#include "constants.h"
#include "drawmirna.h"



using namespace std;

char matrix[5][81];             //This matrix represents the screen of 5x80 (Height x Width) characters



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


void fillmatrix (node * inpnode, int column)    //This function fills in the appropriate characters at the right place in the screen-matrix, the global 2-d array 'matrix'. 'column' tells the function which column of the matrix to print this node to
{
    TYPE_OF_NODE inpnodetype;
    inpnodetype = inpnode->gettype ();  //Get the current node type

//Just declaring some variables required later
    int p, q, r;
    stem *qstem = new stem;
    asym *qasym = new asym;
    sym *qsym = new sym;
    loop *qloop = new loop;
    string newshrtstr, shrtstr, lngstr, qsymr, qsyml;

//string newshrtstr;

    switch (inpnodetype)        //Filling done depending on the node type
    {
    case ROOT:
        fillmatrix (inpnode->getp (), column);  //Nothing needs to be done at the ROOT. Recurse using the next node
        break;

    case STEM:
        qstem = (stem *) (inpnode->getdata ());
        matrix[1][column] = qstem->l;
        matrix[2][column] = '|';
        matrix[3][column] = qstem->r;
        fillmatrix (inpnode->getp (), column + 1);      //Recurse using the next node
        break;

    case ASYM:
        qasym = (asym *) (inpnode->getdata ());
        p = (qasym->r).length ();
        q = (qasym->l).length ();
        if (q > p) {
            r = q;
        } else {
            r = p;
        }                       //r gets the length of the bigger

        if (q == r) {
            shrtstr = qasym->r;
            lngstr = qasym->l;
        } else {
            shrtstr = qasym->l;
            lngstr = qasym->r;
            p = q;
        }                       //So, 'shrtstr' gets the shorter string and 'p' gets the length of the shorter string
        //newshrtstr=new char(r+1);
        newshrtstr.assign (r - p, '-');
        newshrtstr.append (shrtstr);

        if (shrtstr == qasym->r) {
            qasym->r = newshrtstr;
        } else {
            qasym->l = newshrtstr;
        }

        for (q = 0; q < r; q++) {
            matrix[0][column + q] = qasym->l[q];
            matrix[4][column + q] = qasym->r[q];
        }
        fillmatrix (inpnode->getp (), column + r);      //Recurse using the next node
        break;

    case SYM:
        qsym = (sym *) (inpnode->getdata ());
        p = (qsym->r).length ();
        qsymr = qsym->r;
        qsyml = qsym->l;

        for (q = 0; q < p; q++) {
            matrix[0][column + q] = qsyml[q];
            matrix[4][column + q] = qsymr[q];
        }
        fillmatrix (inpnode->getp (), column + p);      //Recurse using the next node
        break;

    case LOOP:
        qloop = (loop *) (inpnode->getdata ());
        p = (qloop->s).length ();
        string seq = qloop->s;

        if (p > 3) {
            q = int (p / 2) - 1;
            for (r = 0; r < q; r++) {
                matrix[0][column + r] = seq[r];
                matrix[4][column + r] = seq[p - 1 - r];
            }
            column = column + q;
            seq = seq.substr (q, p - (q * 2));
        }

        p = seq.length ();
        switch (p) {
        case 1:
            matrix[2][column] = seq[0];
            break;
        case 2:
            matrix[1][column] = seq[0];
            matrix[3][column] = seq[1];
            break;
        case 3:
            matrix[1][column] = seq[0];
            matrix[2][column] = seq[1];
            matrix[3][column] = seq[2];
        }
        break;
    }
}

void initmatrix ()              //This function initialises the 'matrix' with spaces
{
    int r, c;
    for (r = 0; r < 5; r++) {
        for (c = 0; c < 55; c++)
            matrix[r][c] = ' ';
        matrix[r][c] = '\0';
    }
}

void displaymatrix ()           //This function prints the matrix on screen
{
    printf ("%s", matrix[0]);
    printf ("\n%s", matrix[1]);
    printf ("\n%s", matrix[2]);
    printf ("\n%s", matrix[3]);
    printf ("\n%s\n\n", matrix[4]);
}

string returnmatrix ()          //This function returns the matrix into a string
{
    ostringstream oss;
    oss << matrix[0] << "\n" << matrix[1] << "\n" << matrix[2] << "\n" <<
        matrix[3] << "\n" << matrix[4] << "\n\n";
    return oss.str ();
}

