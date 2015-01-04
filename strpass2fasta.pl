#!/usr/bin/perl

if (scalar(@ARGV) <= 1)
   {
   die ("Insufficient arguments!!\nSYNTAX:\nstrpass2fasta.pl <File in the final results format> <Output FASTA file name>\n\n");
   }

$UnorderedLoci = $ARGV[0];
$FastaF = $ARGV[1];


#Read the gene Loci
open (INWITHSTR, $UnorderedLoci) ||  die ("ERROR!! Could not open Input File in the final results format.\n");
open (OUTFASTA, ">$FastaF") || die ("ERROR!! Could not open output file $FastaF\n");

$x=0;
print "Reading and converting lines...\n";
#print "Doing sequence... \n";
while ($GreaterThanLine = <INWITHSTR>)
      {
      if ($GreaterThanLine =~ /^>/)
         {#Start Extraction
         ($tmp,$Locus,$Score,$tmp,$NScore) = split /[>()\s+]/, $GreaterThanLine;
#print "$GreaterThanLine\t$Locus\t$Score\n";
#exit;

#         print "$OnlyLoci[$x]\t$OnlyLen[$x]\n"; 
#         chomp($GeneLociLine);

         $SeqLine = <INWITHSTR>;
         chomp($SeqLine);
         $Len = length($SeqLine);
#         $NScore = $Score / $Len;
         $FastaLine= ">Predicted-miRNA Position:".$Locus." Length:".$Len." GrammarScore:".$Score." NormalisedGrammarScore:".$NScore;

         $StructLine1 = <INWITHSTR>;
         $StructLine2 = <INWITHSTR>;
         $StructLine3 = <INWITHSTR>;
         $StructLine4 = <INWITHSTR>;
         $StructLine5 = <INWITHSTR>;

         $StructScoreLine = <INWITHSTR>;
         chomp($StructScoreLine);
         $StructScore = 0;
         if ($StructScoreLine =~ /SCORE: /)
            {
            $StructScore = $';
            }

         $FastaLine .= " StructuralScore:".$StructScore."\n";
         $SeqLine =~ s/U/T/g;

         print OUTFASTA $FastaLine;
         print OUTFASTA "$SeqLine\n";
         print OUTFASTA $StructLine1;
         print OUTFASTA $StructLine2;
         print OUTFASTA $StructLine3;
         print OUTFASTA $StructLine4;
         print OUTFASTA $StructLine5;

         $x++;
         }#End Extraction

      print "$x\t";
      }
close INWITHSTR;
close OUTFASTA;
print "\n\n$x sequences done!!\n\n";
