//**********************************************************************************
// File         : drawmirna.h
// Purpose      : Declaration and definition of classes and functions for the generation of a graphical output of an miRNA
// Author       : Manish Kushwaha (manish.kushwaha@gmail.com)
// Revision     : 1.0   -       28-May-2005     Initial Draft that used a Stub to randomly generate trees and printed them
//              : 2.0   -       15-Jun-2005     Modification to accept input from the miRNA program as a string stream
//              : 3.0   -       25-Jun-2005     Modification to return the output as a string, rather than a character type matrix, and calculate statistical details about a structure
//**********************************************************************************

#ifndef _DRAW__MIRNA_H_
#define _DRAW__MIRNA_H_

#include <string>

#include "mirna.h"


void fillmatrix (node *, int);  //For Matrix
void initmatrix ();             //For Matrix
void displaymatrix ();          //For Matrix
std::string returnmatrix ();         //For Matrix
node *BuildNodes (std::string & trace);

#endif
