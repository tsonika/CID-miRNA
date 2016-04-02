#ifndef _MIRNA_H_
#define _MIRNA_H_

// miRNA structure definitions

#include <string>
#include <vector>

enum TYPE_OF_NODE
{ ROOT = 0, STEM = 1, SYM = 2, ASYM = 3, LOOP = 4 };

struct stem
{
    char l, r;
};                              //This structure implements a STEM. STEM can have one base on the left and another on the right (characters 'l' and 'r' respectively)

struct sym
{
    std::string l, r;
};                              //This structure implements a SYM (=symmetrical bulge). SYM can have any number of characters on the right and the left (strings 'l' and 'r' respectively)

struct asym
{
    std::string l, r;
};                              //This structure implements an ASYM (=asymmetrical bulge). ASYM can have any number of characters on the right and the left (strings 'l' and 'r' respectively)

struct loop
{
    std::string s;
};                              //This structure implements a LOOP. SYM can have one stretch of string that folds like a cap at the terminus of the loop


class node
{
  public:
    TYPE_OF_NODE type;          //Type of the current node
    node *p;                    //Pointer to the next node
    void *data;                 //Pointer to data in the current node

  public:
    node (TYPE_OF_NODE);
    void initialize (TYPE_OF_NODE);
    void storep (node *);
    void storedata (void *);
    TYPE_OF_NODE gettype ();
    node *getp ();
    void *getdata ();
};


#endif
