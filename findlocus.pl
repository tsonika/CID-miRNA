#!/usr/bin/perl


#This program looks for the given sequences returned from the CutOff selector program in the source Contig file, and
#writes another file that can serve as input for the Overlap Remover

#Check the number of arguments
if (scalar(@ARGV) < 3)
	{
	die "\nInsufficient number of arguments!!\n\nProgramName <input file from Cutoff Selector> <Source Contig File> <output file for Overlap Remover>\n";
	}
	else
	{
	$InpFile=$ARGV[0];
	$ContigFile=$ARGV[1];
	$OutFile=$ARGV[2];
	}


#Open the files for reading
print "\nOpening input files for reading...\n";
open INPFILE, "<$InpFile" || die "Could not open file: $InpFile\n";
open CONTIGFILE, "<$ContigFile" || die "Could not open file: $ContigFile\n";


#Open the file for writing
print "\nOpening output file...\n";
open OUTFILE, ">$OutFile" || die "Could not open file: $OutFile\n";

#Pickup the Entire Contig
$Contig="";
while ($ContigLine=<CONTIGFILE>)
	{
	chomp($ContigLine);
	if ($ContigLine !~ /^\>/)
		{# If the file line does not have the line beginning with '>'
		$Contig=$Contig.$ContigLine;
		}
	}
close CONTIGFILE;
$Contig =~ s/\s+//g; #Remove all whitespaces
$Contig =~ tr/tT/uU/;
$Contig = uc($Contig);

#open  CONTIGFILE, '>Modified_ContigFile.fa';
#print CONTIGFILE $Contig;
#close CONTIGFILE;

#Run through the Input File and look for its position in the ContigFile
print "\nRunning through the Input file and determining the positions...\n\nPlease wait...\n";

while ($InputLine1=<INPFILE>)
	{
	chomp($InputLine1);
	if ($InputLine1 =~ /Sequence \:/)
		{#Thus the first line with the pattern "Sequence :" has been found
		$InputLine2=<INPFILE>;
		chomp($InputLine2);
		if (length($InputLine2))
			{#Thus the second line with some non empty entry found
			$InputLine3=<INPFILE>;
			if (length($InputLine3))
				{#Thus the third line with some non empty entry found

				($tmp,$tmp,$tmp,$SeqLen,$tmp,$tmp,$tmp,$NSeqScore,$tmp,$tmp,$SeqScore)=split (/\s+/,$InputLine3);

				$InputLine2Flipped=$InputLine2;
                                $InputLine2Flipped=uc($InputLine2Flipped);
				#$InputLine2Flipped =~ tr/uU/tT/;
				#$InputLine2Flipped =~ s/U/T/g;
				#$InputLine2Flipped =~ s/u/t/g;
                                $SeqPos=0;$LastSeqPos=-1;

                                while ($SeqPos >= 0)
                                      {
     				      $SeqPos=index($Contig,$InputLine2Flipped,$SeqPos);
				      if ($SeqPos == -1)
				     	     {#The sequence is not found in the Contig
					     print STDERR "Error!! A sequence in the input file was not found in the Source fasta.\n\($InputLine1\)\nContinuing with other sequences...\n\n";
					     }
					     elsif ($SeqPos >= 0)
					     {
     					     #Now writing the file in the output format required for the Overlap Remover
                                             $LastSeqPos = $SeqPos;
					     $SeqPos++;
					     print OUTFILE " Predicted miRNA\tPosition: $SeqPos   Length = $SeqLen\n";
					     print OUTFILE " Score = $SeqScore\tNormalized Score= $NSeqScore\n";
                                             $InputLine2 =~ tr/tT/uU/;
                                             $InputLine2 = uc($InputLine2);
					     print OUTFILE "$InputLine2\n";
					     }
					     else
					     {
					     if ($LastSeqPos == -1) {print "Error!! Negative Position!!\n$InputLine2Flipped\n$SeqPos\n";}
#					     print "$SeqPos"
					     }
                                      }#end while


				}
			}
		}
	}

print "DONE!!\n";
close OUTFILE;
close INPFILE;

