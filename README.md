
#CID-miRNA: A web server for prediction of novel miRNA precursors in human genome 
http://melb.agrf.org.au:8888/

CIDmiRNA is a tool for computer-assisted identification of micro-RNA using an
SCFG model and has been designed to analyze either a single sequence or complete
genome.

If you use the results of this software in a published paper, please cite the
following article:

    Tyagi S. et. al., Biochemical and Biophysical Research Communications                       
    Volume 372, Issue 4, 8 August 2008, Pages 831-834                                           

For help or comments, please use the following group:
https://groups.google.com/forum/#!forum/cid-mirna

It is written in Python. It runs on Python 2.7.x.

It also depends on RNAfold, which can be downloaded from
http://www.tbi.univie.ac.at/RNA/RNAfold.html  It has been tested with version
2.1.9, but it is very likely to work with earlier versions too.


###Input sequence

One can either upload or paste FASTA sequence(s). The sequence should not
contain any ambiguous characters (only A, C, G, T and U are allowed).  If DNA
sequence is submitted all the Ts are converted to Us and processed.


###Running

Run 'make' before running the program the first time to compile cutoffpassscore
and scorestructure


You can get a complete list of options and defaults by running:

```
./bin/cidmirna.py --help
```

By default, all output is directed to the current directory. You can change that
by specifying the preferred directory with the -o option.

If you are running on an SGE cluster, you can pass --sge to make use of it. You
can specify what queue to run on by passing it as --sge-queue <queuename>, and
direct the SGE logs with --sge-logs <directory>


Try Yourself:

```
./bin/cidmirna.py testin.fa
```


### Organism

Although current version supports human miRNA detection, miRNA of closely
related species can also be identified using grammar derived from human miRNAs.

###Minimum stem length

Minimum number of complementary base pairs at 5` and 3` end of potential miRNA
sequence is 3.


###Window

Length of known human miRNA, used for training, varies from 60 to 125. Window
length can be selected from this range and the difference  between the two
should be at least equal to the stem length chosen.


###Score cut-off

The optimal cutoff for human miRNAs is -0.609999.
