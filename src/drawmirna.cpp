#include "drawmirna.h"

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
