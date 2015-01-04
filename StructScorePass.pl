#!/usr/bin/perl

if (scalar(@ARGV) <= 2)
   {
   die ("Insufficient arguments!!\nSYNTAX:\nStructScorePass.pl <Input file in the final results format> <Structural Cutoff Score> <Output file in final results format>\n\n");
   }
$CSS = 0;
$BeforePass = $ARGV[0];
$CSS = $ARGV[1];
$AfterPass = $ARGV[2];


#Read the gene Loci
open (BEFOREPASS, $BeforePass) ||  die ("ERROR!! Could not open Input File in the final results format.\n");
open (AFTERPASS, ">$AfterPass") || die ("ERROR!! Could not open output file $AfterPass\n");

$x=0;
print "Reading and checking sequence scores...\n";
#print "Doing sequence... \n";
while ($GreaterThanLine = <BEFOREPASS>)
      {
      if ($GreaterThanLine =~ /^>/)
         {#Start Extraction
         $SeqLine = <BEFOREPASS>;
         $StructLine1 = <BEFOREPASS>;
         $StructLine2 = <BEFOREPASS>;
         $StructLine3 = <BEFOREPASS>;
         $StructLine4 = <BEFOREPASS>;
         $StructLine5 = <BEFOREPASS>;
         $BlankLine = <BEFOREPASS>;

         $BlankLine_Stripped = $BlankLine;
	 $BlankLine_Stripped =~ s/\s+//g;
         if (length($BlankLine_Stripped) != 0) {$StructScoreLine = $BlankLine;} else {$StructScoreLine = <BEFOREPASS>;}

#         $StructScoreLine = <BEFOREPASS>;
         chomp($StructScoreLine);
         $StructScore = 0;
         if ($StructScoreLine =~ /^SCORE: /)
            {
            $StructScore = $';
            }

         print "\nSCORE:$StructScore\n";
         if ($StructScore >= $CSS)
            {
            print AFTERPASS $GreaterThanLine;
            print AFTERPASS $SeqLine;
            print AFTERPASS $StructLine1;
            print AFTERPASS $StructLine2;
            print AFTERPASS $StructLine3;
            print AFTERPASS $StructLine4;
            print AFTERPASS $StructLine5;
            print AFTERPASS "$StructScoreLine\n";
            }

         $x++;
         }#End Extraction

      print "$x\t";
      }
close BEFOREPASS;
close AFTERPASS;
print "\n\n$x sequences done!!\n\n";
