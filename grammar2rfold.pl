#!/usr/bin/perl


$PGM = $0;
$PGM =~ s#.*/##;              
$usage = <<USAGE;               # usage message
  USAGE:
        $PGM -i  <input file name>

       -i <input file name>  It reads in the  grammar output from a file.
       
       	Example of input format:
	>Predicted      miRNA Position: 23209   Length: 100     Score: -60.1264 Normalised Score: -0.601264
	CCTTCTAGTGGCAAGAGTGACGTAAGTGATATGCGGAAATTTCTTTCCAAGCCTGCTTGGAGAAGCTTCCTCTGCCTGCTTCTCTTTGGCCACCTCCAGG
	>Predicted      miRNA Position: 23509   Length: 75      Score: -45.5868 Normalised Score: -0.607824
	GTTTTCCCTCTTATGTCCAGCAAATGCTGCATGGAGCCCTGGAATTCTATGTGGAAAGCTAGGAAGAGGGAGAGC
	>Predicted      miRNA Position: 38122   Length: 62      Score: -37.5368 Normalised Score: -0.605432
	GGGTCTTTGTGTCAATCTGAGCTCTGATGTCCACCTAGAGATTGGGTATCCACCTAAGGCCC

	This program converts the grammar out format to the one
	recognized by RNAfold with constraints options.

	Writes to STDOUT.

USAGE

$nargs = 2;                     # number of required args
                                                                                
#print $#ARGV+1";
if ($#ARGV+1 <$nargs) { print "$usage"; exit(1); }

# get input arguments
while ($#ARGV >= 0) {
  $_ = shift;
  if($_ eq "-h") {
    print  "$usage"; exit(1);
  } elsif ($_ eq "-i") {
    $inputfile = shift;
  } else {
     print "$usage"; exit(1);
  }
}

open(GM, $inputfile) || die ("Could not open the grammar output file\n");

while (<GM>) {
  chomp;
  if ($_ !~ /^#/ or $_ !~ /^$/) {
    if(/^>/) { # read the header line
      @temparray = split(/\s+/);
      $pos = $temparray[3];
      $score = $temparray[7];
      $normscore = $temparray[10];
      undef @temparray;
      print ">$pos($score)($normscore)\n";
    } else {
      $seq = $_;
      print "$seq\n";
      for($i=0; $i<int(length($seq)/2); $i++) {
	print "<";
      } 
      for($i=0; $i<(length($seq)-int(length($seq)/2)); $i++) {
	print ".";
      } 
    print "\n";
    } #end else 
  } # if not a blank or comment line

} # end while 
close GM;
