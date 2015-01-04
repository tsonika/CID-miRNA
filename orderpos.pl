#!/usr/bin/perl

if (scalar(@ARGV) <= 0)
   {
   die ("No arguments specified!!\nSYNTAX:\norderpos <File with Unordered positions>\n\n");
   }

$UnorderedLoci = $ARGV[0];


#Read the gene Loci
open (INWITHLOCI, $UnorderedLoci) ||  die ("ERROR!! Could not open Input File with Unordered positions.\n");
$x=0;
print "Reading sequences...\n";
print "Doing sequence... \n";
while ($GeneLociLine = <INWITHLOCI>)
      {
      if ($GeneLociLine =~ /Predicted miRNA/)
         {#Start Extraction
         ($tmp,$tmp,$tmp,$tmp,$OnlyLoci[$x],$tmp) = split /\s+/, $GeneLociLine;
#         print "$tmp1\t$OnlyLoci[$x]\n"; 

         $AtLoci[$x]=$GeneLociLine;
         $GeneLociLine = <INWITHLOCI>;
         $AtLoci[$x] .= $GeneLociLine;
         $GeneLociLine = <INWITHLOCI>;
         $AtLoci[$x] .= $GeneLociLine;
         }#End Extraction

      $x++;
      print "$x\t";
      }
close INWITHLOCI;
print "\n\n$x sequences read!!\n\n";


#Ordering Process
print "Now ordering them...\n";
$on=0;
$full_len=scalar(@OnlyLoci);
for ($a=0; $a < $full_len ; $a++)
    {
    print "$a of $full_len\n";
    $LowPos = $a;

    for ($b=$a; $b < $full_len ; $b++)
        {
        if ($OnlyLoci[$b] < $OnlyLoci[$LowPos])
           {#Finding the position with lowest value
           $LowPos=$b;
           }
         }#end $b

           #Swapping
           $tmp=$OnlyLoci[$a];
           $OnlyLoci[$a]=$OnlyLoci[$LowPos];
           $OnlyLoci[$LowPos]=$tmp;

           $tmp=$AtLoci[$a];
           $AtLoci[$a]=$AtLoci[$LowPos];
           $AtLoci[$LowPos]=$tmp;
    }#end $a

print "\n DONE ORDERING!!\n";

print "\n\n Now Writing them to file....\n\n";
open (OUTLOCI, ">$UnorderedLoci.ordered") || die ("ERROR!! Could not open output file $UnorderedLoci.ordered\n");

for ($a=0; $a < scalar(@OnlyLoci) ; $a++)
    {
    print OUTLOCI $AtLoci[$a];
    }

close OUTLOCI;
print "Done writing them to file.\n\n";

