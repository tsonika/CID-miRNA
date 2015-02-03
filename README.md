
#CID-miRNA:CID-miRNA: A web server for prediction of novel miRNA precursors in human genome 
http://mirna.jnu.ac.in/cidmirna

Please cite this article as:                                                                 
Tyagi S. et. al., Biochemical and Biophysical Research Communications                       
Volume 372, Issue 4, 8 August 2008, Pages 831-834                                           


CIDmiRNA is the tool for computer-assisted identification of micro-RNA using an SCFG model and has been designed to analyze either a single sequence or complete genome.

It depends on RNAfold, which can be downloaded from http://www.tbi.univie.ac.at/RNA/RNAfold.html It has been tested 
with version 2.1.8, but it is very likely that earlier versions will work too.


###Input sequence

One can either upload or paste FASTA sequence(s). The sequence should not contain any ambiguous characters (only A, C, G, T and U are allowed). 
If DNA sequence is submitted all the Ts are converted to Us and processed. 

### Organism

Although current version supports human miRNA detection, miRNA of closely related species can also be identified using grammar derived from human miRNAs.

###Minimum stem length

Minimum number of complementary base pairs at 5` and 3` end of potential miRNA sequence is 3.


###Window

Length of known human miRNA, used for training, varies  from 60 to 125. Window length can be selected from this range and the difference 
between the two should be at least equal to the stem length chosen.


###Score cut-off

The optimal cutoff for human miRNAs is -0.609999.


You can get a complete list of options and defaults by running:

```
./cidmirna.py --help
```

###Example usage of the CID-miRNA

```
./cidmirna.py <filename>
```


Try Yourself:

```
./cidmirna.py testin.fa
```

Please contact Dr Sonika Tyagi (sonika.tyagi@gmail.com) for any queries
regarding the usage of the tool.
