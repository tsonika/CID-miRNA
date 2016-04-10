#include <fstream>
#include <iostream>

#include <sstream>
#include <stdlib.h>

#include "mirnastats.h"


using namespace std;

miRNAstats getstats (string inp, bool considertruncated)
{
    miRNAstats tempstats;
    //Initialise 'tempstats'
    tempstats.n_stems = 0;
    tempstats.n_asyms = 0;
    tempstats.n_syms = 0;
    tempstats.n_bases = 0;
    tempstats.n_base_pairs = 0;
    tempstats.n_loop_bases = 0;

    istringstream iss (inp);
    string mat[5];
    string temp = "";

    //To initialise all strings in the matrix 'mat' with "" 
    int r;
    size_t c;
    for (r = 0; r < 5; r++)
        mat[r] = "";
    r = -1;
    size_t maxlen = 0;

    //To distribute each non-zero length line into the matrix 'mat'
    do {
        getline (iss, temp);
        if (temp.length () > 0 && r < 4) {
            mat[++r] = temp;
            if (mat[r].length () > maxlen)
                maxlen = mat[r].length ();
        }
    }
    while (!iss.eof ());

    if (maxlen == 0 || mat[4].length () == 0) {
        cout << "ERROR!! Invalid or incomplete input string.\n";
        exit (1);
    }

    
    //Now analyse each column
    enum KEEPING_TRACK
    { NONE = 0, ROOT = 1, STEM = 2, SYM = 3, ASYM = 4, LOOP = 5 };
    KEEPING_TRACK track = NONE;
    int temploopbases = 0;
    bool foundbp = false;
    double temp_stem_len, temp_asym_len, temp_asym_nbases, temp_sym_len;

    for (c = 0; c < maxlen; c++) {
        if (track == NONE) {    //If a track is OFF

            if (considertruncated && !foundbp) {        //If considertruncated is switched on, the function will ignore anything before a basepair
                if (mat[2][c] == '|') {
                    foundbp = true;
                } else continue;
            }

            if (mat[2][c] == '|') {     //This is part of STEM
                track = STEM;
                temp_stem_len = 1;
                ++(tempstats.n_stems);  //Since a new STEM has been encountered, increase STEM count
                tempstats.n_bases += 2; //Since each STEM point has two bases, increase base count by 2
                ++(tempstats.n_base_pairs);     //Since a new STEM has been encountered, increase base pair count by 1
            } else if (isalpha (mat[0][c]) && isalpha (mat[4][c])) {       //This is part of SYM
                track = SYM;
                temp_sym_len = 1;
                tempstats.n_bases += 2; //Since each SYM point has two bases, increase base count by 2
                temploopbases = 2;      //Since each SYM point can later become a loop point and has two bases
            } else if (mat[0][c] == '-' || mat[4][c] == '-') {    //This is part of ASYM
                track = ASYM;
                temp_asym_len = 1;
                temp_asym_nbases = 1;
                ++(tempstats.n_bases);  //Since each ASYM point has one base, increase base count by 1
            } else if ((isalpha (mat[1][c]) && isalpha (mat[3][c])) || isalpha (mat[2][c])) {       
                //This is part of LOOP
                track = LOOP;
                --c;
            }
        }                       //End- If a track is OFF
        else {                  //If a track is ON

            if (track == STEM) {        //That is, if a STEM is being monitored
                if (mat[2][c] == '|') { //This point is also part of STEM
                    tempstats.n_bases += 2;     //Since each STEM point has two bases, increase base count by 2
                    temp_stem_len++;
                    ++(tempstats.n_base_pairs); //Since a new STEM has been encountered, increase base pair count by1
                } else {        //This point no more represents a STEM; so switch OFF the track
                    track = NONE;
                    tempstats.stem_len.push_back (temp_stem_len);
                    --c;
                }
            }                   //End- if a STEM is being monitored
            else if (track == SYM) { //That is, if a SYM is being monitored
                if (mat[2][c] != '|' && isalpha (mat[0][c])
                    && isalpha (mat[4][c])) {   //This point is also part of SYM
                    tempstats.n_bases += 2;     //Since each SYM point has two bases, increase base count by 2
                    temp_sym_len++;
                    temploopbases += 2; //Since each SYM point can later become a loop point and has two bases
                }               //End- This point is also part of SYM
                else {          //This point no more represents a SYM; so switch OFF or modify the track
                    if (mat[2][c] == '|') {     //A STEM has begun, switch OFF track and increment SYM count by 1
                        track = NONE;
                        tempstats.sym_len.push_back (temp_sym_len);
                        --c;
                        ++(tempstats.n_syms);
                        temploopbases = 0;
                    }           //End- A STEM has begun
                    else if (mat[0][c] == '-' || mat[4][c] == '-') {   //This is ASYM, modify track and repeat
                        track = ASYM;
                        temp_asym_len = temp_sym_len;
                        temp_asym_nbases = temp_sym_len * 2;
                        --c;
                        temploopbases = 0;
                    }           //End- This is ASYM

                    else if ((isalpha (mat[1][c]) && isalpha (mat[3][c])) || isalpha (mat[2][c])) {       
                    //This is LOOP, modify track and go ahead
                        track = LOOP;
                        --c;
                    }
                }               //End- This point no more represents a SYM
            }                   //End- if a SYM is being monitored
            else if (track == ASYM) {        //That is, if an ASYM is being monitored
                if (mat[2][c] != '|') { //This point is also part of ASYM
                    ++(tempstats.n_bases);      //Since each ASYM point has one base, increase base count by 1
                    temp_asym_len++;
                    temp_asym_nbases++;
                }               //End- This point is also part of ASYM
                else {          //This point no more represents an ASYM; so switch OFF the track and increment ASYM count by 1
                    track = NONE;
                    tempstats.asym_len.push_back (temp_asym_len);
                    tempstats.asym_nbases.push_back (temp_asym_nbases);
                    --c;
                    ++(tempstats.n_asyms);
                }               //End- This point no more represents an ASYM
            }                   //End- if an ASYM is being monitored
            else if (track == LOOP) {        //That is, if a LOOP is being monitored
//(isalpha(mat[1][c]) && isalpha(mat[3][c])) || isalpha(mat[2][c])
                if (isalpha (mat[1][c]) && isalpha (mat[2][c])
                    && isalpha (mat[3][c])) {   //Last position with three bases
                    tempstats.n_bases += 3;     //Since this Last LOOP point has three bases, increase base count by 3
                    temploopbases += 3; //Since this Last LOOP point has three bases
                } else if (isalpha (mat[1][c]) && isalpha (mat[3][c])) {        //Last position with two bases
                    tempstats.n_bases += 2;     //Since this Last LOOP point has two bases, increase base count by 2
                    temploopbases += 2; //Since this Last LOOP point has two bases
                } else if (isalpha (mat[2][c])) {       //Last position with one base
                    ++(tempstats.n_bases);      //Since this Last LOOP point has one base, increase base count by 1
                    ++temploopbases;    //Since this Last LOOP point has one base
                }
            }                   //End- if a LOOP is being monitored

        }                       //End- If a track is ON
    }                           //End- for

    if (temploopbases != 0)
        tempstats.n_loop_bases = temploopbases;

    return tempstats;
}


