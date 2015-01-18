#include <iostream>
#include <string.h>
#include <fstream>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

using namespace std;

// Creating a Parser for the Data File

const long MAXPATH = 100;
//const unsigned long long MAXMEM = 100000000;
const unsigned long int MAXFILESEQ = 2000000;
const int MAXFILENUM = 625;
//Max sequences per file

char cykoutputflag = 'Y';

class parser
{
  private:unsigned int windowlen;
    unsigned int minrange;
    unsigned int maxrange;
    unsigned long int cachelen;
    unsigned long int totalseqlen;
    unsigned long int totalacclen;
    unsigned long cachecur;
    unsigned long int filelen;
    unsigned long int seqcount;
    unsigned long int *seqlentracker;
    char filecount[3];
    unsigned long int fileseqcounter;
    char flagverbose;

    char *cache;
    char *window;
    char *compwindow;
    char *searchwindow;
    char *sequence;
    char infile[MAXPATH];
    char outfile[MAXPATH];
    char currentfile[MAXPATH + 5];

    void mem_allocator (void);
    // allocate memory to cache and window
    int complementor (void);
    // Complements the current window and returns the validity of window
    void filler (char *str, unsigned long start, unsigned long end);
    void cacheparser (void);
    // The actual parser
    char *mystrstr (char *targetstr, char *searchstr);
    char altcomplement (char ichar);
    char DNAtoRNAcomplement (char ichar);
    bool checkalienchar (char *checkstr);

  public:void initializer (void);
    // get all the parameters manually
    void initializer2 (char filein[MAXPATH], char fileout[MAXPATH],
                       unsigned long winlen, unsigned long minr,
                       unsigned long maxr, unsigned long int cachelength = 0);
    void freememory (void);


};

void
parser::cacheparser (void)
{

    fstream fout[MAXFILENUM];
    int filectr;


    filecount[0] = filecount[1] = filecount[2] = '\0';
    filectr = 0;

    filecount[0] = 'A';
    filecount[1] = 'A';
    filecount[2] = '\0';
    strcpy (currentfile, outfile);
    strcat (currentfile, filecount);

    fout[filectr].open (currentfile, ios::out | ios::trunc);

    if (!fout[filectr].good ()) {
        cout << "\n Error Opening Output File";
        exit (1);
    }

    cachecur = 0;
    seqcount = 0;
    fileseqcounter = 0;
    char *ptr = 0;
    short foundmatch = 0;
    unsigned long seqlen = 0;
    totalseqlen = 0;
    totalacclen = 0;
    unsigned long int curwindowstart = 0, curwindowend = 0;

    cout << "\n PARSING CACHE SIZE " << cachelen << "\n Please Wait....";

    cout << "\n Do you want output in CYK Format (y/n) : y";

    //cin >> cykoutputflag;
    //cykoutputflag = toupper(cykoutputflag);

    cykoutputflag = 'Y';

    if (cykoutputflag == 'N')
        fout[filectr] << "\n !!  SCANNING STEM LENGTH [" << windowlen <<
            "] IN CACHE SIZE [" << cachelen << "]  !!\n\n";

    while (cachecur < (cachelen - 2 * windowlen - minrange)) {
        foundmatch = 0;
        filler (window, cachecur, windowlen);
        //cout <<"\n" << window <<"\t"<<cachecur<<" of " <<cachelen;

        if (!complementor ())
            continue;

        filler (searchwindow, cachecur + windowlen + minrange,
                maxrange - minrange + windowlen);

        ptr = mystrstr (searchwindow, compwindow);
        while (ptr) {
            foundmatch = 1;
            // *ptr = '\0';
            seqlen =
                windowlen + minrange + strlen (searchwindow) -
                strlen (ptr) + windowlen;
            filler (sequence, cachecur, seqlen);
            //cout << "\n Sequence Found !! COUNT = " <<
            ++seqcount;
            ++fileseqcounter;
            if (fileseqcounter > MAXFILESEQ) {
                fileseqcounter = 0;

                // incrementing file extension

                if (filecount[1] == 'Z') {
                    filecount[1] = 'A';
                    filecount[0]++;
                } else
                    filecount[1]++;

                fout[filectr].close ();

                ++filectr;      // increasing file number

                strcpy (currentfile, outfile);
                strcat (currentfile, filecount);

                fout[filectr].open (currentfile, ios::out | ios::trunc);
            }
            if (!fout[filectr].good ()) {
                cout << "\n Error Writing to Output File !!\n";
                exit (1);
            }

            if ((flagverbose == 'y') && checkalienchar (sequence))      //if TRUE
            {
                if (cykoutputflag == 'N')
                    fout[filectr] << "\n> Position : " << cachecur <<
                        "\t Length : " << seqlen << "\n";
                fout[filectr] << sequence << "\n";
            }

            totalseqlen += seqlen;
            // Just keeping and incrementing stats of all lengths found
            ++seqlentracker[seqlen - 2 * windowlen - minrange];

            ptr++;
            ptr = mystrstr (ptr, compwindow);
        }
        if (foundmatch) {
            // This particular portion is calculating the acceptance rate
            if ((cachecur > curwindowstart) && (cachecur <= curwindowend)
                && ((cachecur + seqlen) > curwindowend))
                curwindowend = cachecur + seqlen;
            else if (cachecur > curwindowend) {
                totalacclen += (curwindowend - curwindowstart);
                curwindowstart = cachecur;
                curwindowend = cachecur + seqlen;
            }
            // Calculated the accepted portion
        }
        ++cachecur;
    }

    // Compensating Acceptance for last window accepted
    totalacclen += (curwindowend - curwindowstart);

    // Calculate Acceptance Rate
    double rejrate;
    rejrate = ((double) totalacclen / (double) cachelen) * 100.0;

    // Calculate Avg Sequence length
    double avgseqlen;
    avgseqlen = (double) totalseqlen / (double) seqcount;

    // Calculating the Mode Sequence Length
    unsigned long int maxmode = 0, posmaxmode = 0, modectr = 0;
    for (modectr = 0; modectr <= (maxrange - minrange); ++modectr) {
        if (maxmode < seqlentracker[modectr]) {
            maxmode = seqlentracker[modectr];
            posmaxmode = modectr;
        }
    }

    // Calculating the Median Sequence length
    unsigned long int medianpos = 0;
    long int medianctr = (seqcount + 1) / 2;

    while ((medianctr > 0) && (medianpos <= (maxrange - minrange))) {
        medianctr -= seqlentracker[medianpos];
        ++medianpos;
    }
    --medianpos;                // compensating one extra increment in while loop

    // Calculating Standard Deviation from Mode
    double stddevsum = 0.0, stddev = 0.0, stddevpc = 0.0;
    for (modectr = 0; modectr <= (maxrange - minrange); ++modectr)
        stddevsum +=
            seqlentracker[modectr] * (posmaxmode - modectr) * (posmaxmode -
                                                               modectr);

    stddev = sqrt (stddevsum / seqcount);
    stddevpc = stddev / (posmaxmode + 2 * windowlen + minrange) * 100;


    cout << "\n\n End of Parsing !!"
        << "\n\n " << seqcount << " SEQUENCES FOUND !!\n"
        << "\n ACCEPTANCE RATE : " << rejrate << " % \n"
        // Just a tester
        //<<(double)totalseqlen/(double)cachelen*100.0 <<" % \n"
        << "\n MINIMUM SEQUENCE LENGTH \t: " << 2 * windowlen + minrange
        << "\n MAXIMUM SEQUENCE LENGTH \t: " << 2 * windowlen + maxrange
        << "\n AVERAGE SEQUENCE LENGTH \t: " << avgseqlen
        << "\n MEDIAN  SEQUENCE LENGTH \t: " << medianpos + 2 * windowlen +
        minrange << "\n MODE    SEQUENCE LENGTH \t: " << posmaxmode +
        2 * windowlen +
        minrange << "\n STD. DEVIATION FROM MODE\t: " << stddev << " ( " <<
        stddevpc << " % ) \n" << "\n Please see file [" << outfile <<
        " AA to " << filecount << "] for output !!\n";

    // Printing Distribution
    cout << "\n Press any key to view distribution.....";
    //cin.get();
    cout << "\n | LENGTH\tFREQ\t| LENGTH\tFREQ\t| LENGTH\tFREQ\t|";
    short int colctr;
    for (modectr = 0; modectr <= (maxrange - minrange - 3); modectr += 3) {
        cout << "\n ";
        for (colctr = 0; colctr < 3; ++colctr)
            cout << "| " << modectr + colctr + 2 * windowlen +
                minrange << "\t\t" << seqlentracker[modectr + colctr] << "\t";
        cout << "|";
    }
    cout << "\n ";
    for (colctr = 0; colctr <= (maxrange - minrange - modectr); ++colctr)
        cout << "| " << modectr + colctr + 2 * windowlen +
            minrange << "\t\t" << seqlentracker[modectr + colctr] << "\t";
    cout << "| ";
    // Distribution Printed
    /* if(cykoutputflag == 'N')
       {
       fout[filectr] << "\n\n PARSING RESULT :"
       << "\n\n " << seqcount <<" SEQUENCES FOUND !!\n"
       << "\n ACCEPTANCE RATE : "<< rejrate << " % \n"
       << "\n MINIMUM SEQUENCE LENGTH \t: "<<2*windowlen+minrange
       << "\n MAXIMUM SEQUENCE LENGTH \t: "<<2*windowlen+maxrange
       << "\n AVERAGE SEQUENCE LENGTH \t: "<< avgseqlen
       << "\n MEDIAN  SEQUENCE LENGTH \t: "<< medianpos + 2*windowlen + minrange
       << "\n MODE    SEQUENCE LENGTH \t: "<< posmaxmode + 2*windowlen + minrange
       << "\n STD. DEVIATION FROM MODE\t: "<< stddev << " ( "<<stddevpc<<" % ) \n"; */

    // Same Distribution in file
    /* fout[filectr] <<"\n SEQUENCE LENGTH DISTRIBUTION \n";
       fout[filectr] <<"\n | LENGTH\tFREQ\t| LENGTH\tFREQ\t| LENGTH\tFREQ";

       for(modectr = 0; modectr <= (maxrange-minrange-3); modectr += 3)
       {
       fout[filectr]<<"\n ";
       for(colctr = 0; colctr < 3; ++colctr)
       fout[filectr] <<"| "<<modectr + colctr + 2*windowlen + minrange <<"\t\t"<< seqlentracker[modectr+colctr]<<"\t";
       fout[filectr] <<"|";
       }
       fout[filectr]<<"\n ";
       for(colctr = 0; colctr <= (maxrange - minrange - modectr); ++colctr)
       fout[filectr] <<"| "<<modectr + colctr + 2*windowlen + minrange <<"\t\t"<<seqlentracker[modectr+colctr]<<"\t";
       fout[filectr] <<"| \n\n";
       // Distribution in file printed
       } */
    fout[filectr] << "\n";
    fout[filectr].close ();

}

void
parser::filler (char *str, unsigned long start, unsigned long len)
{

    unsigned long ctr;

    for (ctr = 0; ((ctr < len) && (cache[start + ctr] != '\0')); ++ctr)
        str[ctr] = cache[start + ctr];

    str[ctr] = '\0';
}

int
parser::complementor (void)
{
    unsigned long int ctr = 0;
    strcpy (compwindow, window);
    while ((ctr < windowlen) && (window[ctr] != '\0')) {
        if (window[ctr] == 'A')
            compwindow[windowlen - 1 - ctr] = 'U';
        else if (window[ctr] == 'U')
            compwindow[windowlen - 1 - ctr] = 'A';
        else if (window[ctr] == 'G')
            compwindow[windowlen - 1 - ctr] = 'C';
        else if (window[ctr] == 'C')
            compwindow[windowlen - 1 - ctr] = 'G';
        else {
            //cout<<"\n Unknown Character found near " << cachecur+ctr;
            cachecur = cachecur + ctr + 1;
            return 0;
        }
        ++ctr;
    }
    return 1;
}


void
parser::initializer (void)
{
    //cout << "\n Enter the name of Input File : ";
    //cin.getline(infile,MAXPATH-1);
    strcpy (infile, "ch22mask");

    //cout << "\n Enter the Name of output file : ";
    //cin.getline(outfile,MAXPATH-1);
    strcpy (outfile, "incyk");
    //cout << "\n Enter Window Length (Stem size to match) : ";
    //cin >> windowlen;
    windowlen = 3;

    //cout << "\n Enter Min Distance to scan from window : ";
    //cin >> minrange;
    minrange = 54;

    //cout << "\n Enter Max Distance to scan from window : ";
    //cin >> maxrange;
    maxrange = 104;

    do {
        //cout << "\n Do you want verbose file output ? (y/n) : ";
        //cin >> flagverbose;
        //flagverbose = tolower(flagverbose);
        flagverbose = 'y';

    }
    while ((flagverbose != 'y') && (flagverbose != 'n'));
    //cout << "\n Enter cachelength (Min : "<< 2*windowlen+maxrange <<", Enter 0 for Max) : ";
    //cin >> cachelen;
    //setting cachelen permanently at zero for MAX CACHE SIZE
    cachelen = 0;
    cout << "\n Allocating Memory. Please wait......";
    mem_allocator ();
}


void
parser::initializer2 (char filein[MAXPATH], char fileout[MAXPATH],
                      unsigned long int winlen, unsigned long int minr,
                      unsigned long int maxr, unsigned long int cachelength)
{

    strcpy (infile, filein);
    strcpy (outfile, fileout);
    windowlen = winlen;
    cachelen = cachelength;
    minrange = minr;
    maxrange = maxr;
    mem_allocator ();
}



void
parser::mem_allocator (void)
{
    filelen = 0;
    char ch;

    fstream f, g;
    f.open (infile, ios::in);
    if (f.good ()) {
        while (!f.eof ()) {
            f.get ();
            ++filelen;
        }
    } else {
        cout << "\n Error Opening Input File !!\n";
        exit (1);
    }

    f.close ();

    if (cachelen == 0) {
        cachelen = filelen + 2;
        //cout << "\n CacheLen : " << cachelen;
    }
    // Allocating Memory

    cache = 0;
    cache = new char[cachelen];
    //cout << "\n CacheLen : " << strlen(cache);

    unsigned long int fillctr;
    for (fillctr = 0; fillctr < strlen (cache); ++fillctr)
        cache[fillctr] = 0;

    //cin.get();
    window = 0;
    compwindow = 0;
    window = new char[windowlen];
    compwindow = new char[windowlen];
    searchwindow = new char[maxrange - minrange + windowlen];
    sequence = new char[2 * windowlen + maxrange];
    seqlentracker = new unsigned long int[maxrange - minrange + 2];

    if (!cache || !window || !compwindow || !searchwindow || !sequence
        || !seqlentracker) {
        cout << "\n Error Allocating Memory !!\n";
        exit (1);
    } else {

        // Setting count of all lengths found as zero
        for (fillctr = 0; fillctr <= (maxrange - minrange + 1); ++fillctr)
            seqlentracker[fillctr] = 0;

        cachecur = 0;

        g.open (infile, ios::in);
        if (!g.good ()) {
            cout << "\n Error Opening Input File !!\n";
            exit (1);
        }

        while ((!g.eof ()) && (cachecur < (cachelen - 1))) {
            ch = g.get ();
            if (isalpha (ch)) {
                //cout<<ch;
                cache[cachecur] = DNAtoRNAcomplement (toupper (ch));
                ++cachecur;
            }
        }
        cache[cachecur] = '\0';

        cout << "\n Cache Fill Size : " << cachecur
            << "\n File Fill Size  : " << filelen;
        cachelen = strlen (cache);
        cout << "\n CACHE LEN : " << cachelen;
        //cin.get();
    }

    g.close ();
    cachecur = 0;
    cacheparser ();
}

void
parser::freememory (void)
{
    delete[]cache;
    delete[]window;
    delete[]compwindow;
    delete[]searchwindow;
    delete[]sequence;
}


// my own strstr

char
parser::altcomplement (char ichar)
{
    if (ichar == 'A')
        return 'G';
    else if (ichar == 'C')
        return 'U';
    return ichar;
}

char *
parser::mystrstr (char *targetstr, char *searchstr)
{
    int searchlen = strlen (searchstr);
    int targetlen = strlen (targetstr);

    if (targetlen < searchlen)
        return 0;

    char *temp = 0;
    short matchflag = 1;
    int ctr, cur = 0;

    while (cur <= (targetlen - searchlen)) {
        matchflag = 1;
        for (ctr = 0; (ctr < searchlen) && matchflag; ++ctr) {
            if ((targetstr[cur + ctr] != searchstr[ctr])
                && (targetstr[cur + ctr] != altcomplement (searchstr[ctr]))) {
                matchflag = 0;
                break;
            }
        }
        if (!matchflag)
            cur = cur + ctr + 1;
        else {
            temp = targetstr + cur;
            return temp;
        }
    }
    return 0;
}

char
parser::DNAtoRNAcomplement (char ichar)
{
    if (ichar == 'A')
        return 'U';
    else if (ichar == 'T')
        return 'A';
    else if (ichar == 'G')
        return 'C';
    else if (ichar == 'C')
        return 'G';
    else
        return ichar;
}

bool parser::checkalienchar (char *checkstr)
{
    int
        i,
        len;
    len = strlen (checkstr);


    for (i = 0; i < len; ++i)
        if ((checkstr[i] != 'A') && (checkstr[i] != 'U')
            && (checkstr[i] != 'G') && (checkstr[i] != 'C')) {
            seqcount--;
            fileseqcounter--;
            return false;
        }

    return true;
}



int
main ()
{
    parser
        p;
    p.initializer ();
    p.freememory ();
    //cin.get();
    cout << "\n";
    return 0;
}
