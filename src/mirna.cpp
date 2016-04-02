#include <stdio.h>

#include "constants.h"
#include "mirna.h"



using namespace std;


node::node (TYPE_OF_NODE t = STEM)      //By default, assume that the node is STEM type and initialise it accordingly
{
    initialize (t);             //Calling the initialiser
}

void node::initialize (TYPE_OF_NODE t)  //This function initialises a node
{
    p = NULL;                   //Needs to be defined after the initialisation, unless type is LOOP

    if (t == ROOT || t == ASYM || t == SYM || t == LOOP) {
        type = t;
    } else {
        type = STEM;
    }

    switch (t)                  //Memory is allocated according to the type of node
    {
    case ASYM:
        if (!(data = new asym)) {
            perror ("Insufficient memory for ASYM\n");
            exit (1);
        }
        break;

    case SYM:
        if (!(data = new sym)) {
            perror ("Insufficient memory for ASYM\n");
            exit (1);
        }
        break;

    case LOOP:
        if (!(data = new loop)) {
            perror ("Insufficient memory for ASYM\n");
            exit (1);
        }
        break;

    case ROOT:
        data = NULL;
        break;


    case STEM:
        if (!(data = new stem)) {
            perror ("Insufficient memory for ASYM\n");
            exit (1);
        }
    }
}

void node::storep (node * a)    //Stores the pointer to the next node
{
    p = a;
}

void node::storedata (void *inp)        //Stores the data associated with current node
{
    data = inp;
}

TYPE_OF_NODE node::gettype ()   //Returns the current node-type
{
    return (type);
}

node *node::getp ()             //Returns the pointer to the next node
{
    return (p);
}

void *node::getdata ()          //Returns the data stored at the current node
{
    return (data);
}
