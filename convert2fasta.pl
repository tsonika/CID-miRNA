#!/usr/bin/perl

if (scalar(@ARGV) <= 1)
   {
   die ("Insufficient arguments!!\nSYNTAX:\nconvert2fasta.pl <File in the final results format> <Output FASTA file name>\n\n");
   }

$UnorderedLoci = $ARGV[0];
$FastaF = $ARGV[1];


#Read the gene Loci
open (INWITHLOCI, $UnorderedLoci) ||  die ("ERROR!! Could not open Input File in the final results format.\n");
open (OUTFASTA, ">$FastaF") || die ("ERROR!! Could not open output file $FastaF\n");

$x=0;
print "Reading and converting lines...\n";
#print "Doing sequence... \n";
while ($GeneLociLine = <INWITHLOCI>)
      {
      if ($GeneLociLine =~ /Predicted miRNA/)
         {#Start Extraction
         ($tmp,$tmp,$tmp,$tmp,$Locus,$tmp,$tmp,$Len) = split /\s+/, $GeneLociLine;
#         print "$OnlyLoci[$x]\t$OnlyLen[$x]\n"; 
         chomp($GeneLociLine);

         $FastaLine= ">Predicted\tmiRNA Position: ".$Locus."\tLength: ".$Len;
         $GeneLociLine = <INWITHLOCI>;
#         chomp($GeneLociLine);
         ($tmp,$tmp,$tmp,$Score,$tmp,$tmp,$NScore,$tmp)= split /\s+/,$GeneLociLine;
         $FastaLine .= "\tScore: ".$Score."\tNormalised Score: ".$NScore."\n";
         $GeneLociLine = <INWITHLOCI>;
         $GeneLociLine =~ s/U/T/g;
         $FastaLine .= $GeneLociLine;
         print OUTFASTA $FastaLine;
         $x++;
         }#End Extraction

      print "$x\t";
      }
close INWITHLOCI;
close OUTFASTA;
print "\n\n$x sequences done!!\n\n";
