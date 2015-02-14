
#CID-miRNA:CID-miRNA: A web server for prediction of novel miRNA precursors in human genome 
Web site appearing soon

Please cite this article as:

    Tyagi S. et. al., Biochemical and Biophysical Research Communications                       
    Volume 372, Issue 4, 8 August 2008, Pages 831-834                                           


CIDmiRNA is the tool for computer-assisted identification of micro-RNA using an SCFG model and has been designed to analyze either a single sequence or complete genome.

It is written in Python. It runs on Python 2.7.x.

It also depends on RNAfold, which can be downloaded from http://www.tbi.univie.ac.at/RNA/RNAfold.html It has been tested 
with version 2.1.8, but it is very likely to work with earlier versions too.


###Input sequence

One can either upload or paste FASTA sequence(s). The sequence should not contain any ambiguous characters (only A, C, G, T and U are allowed). 
If DNA sequence is submitted all the Ts are converted to Us and processed. 


###Running

Run 'make' before running the program the first time to compile cutoffpassscore and Scores4mStruct


You can get a complete list of options and defaults by running:

```
./cidmirna.py --help
```

By default, all output is directed to the current directory. You can change that by specifying
the preferred directory with the -o option.

If you are running on an SGE cluster, you can pass --sge to make use of it. You can specify
what queue to run on by passing it as --sge-queue <queuename>, and direct the SGE logs with --sge-logs <directory>


Try Yourself:

```
./cidmirna.py testin.fa
```


### Organism

Although current version supports human miRNA detection, miRNA of closely related species can also be identified using grammar derived from human miRNAs.

###Minimum stem length

Minimum number of complementary base pairs at 5` and 3` end of potential miRNA sequence is 3.


###Window

Length of known human miRNA, used for training, varies from 60 to 125. Window length can be selected from this range and the difference 
between the two should be at least equal to the stem length chosen.


###Score cut-off

The optimal cutoff for human miRNAs is -0.609999.


Please contact Dr Sonika Tyagi (sonika.tyagi@gmail.com) for any queries
regarding the usage of the tool.
