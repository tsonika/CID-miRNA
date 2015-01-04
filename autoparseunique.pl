#!/usr/bin/perl

#This program automates the process parsing sequences from the parser

$MaxListSize=850;

#Set the Start Benchmark
use Benchmark;
$prgstart = new Benchmark;

# Read the argument that is the FASTA file
$InputFasta=$ARGV[0]; $NForcePairs=$ARGV[1]; $MinBases=$ARGV[2]; $MaxBases=$ARGV[3];
if (scalar(@ARGV)<3 || $InputFasta eq "--help" || $InputFasta eq "-h"){
die "SYNTAX:\nautoparseunique \<Input File name with only sequences to be parsed\> \<The number of Bases to Force Pair at the ends\> \<Lowest Length allowed\> \<Highest Length allowed\>\n";
}

print "\n\nInput File: $InputFasta\n";

########### Each Input FASTA File will be of the type ###########
# >hsa-mir-abc							#
# augcaugcaugcaugc.........					#
# >hsa-mir-pqr							#
# gcaugcaugcaugcau.........		 			#
#################################################################


# Look for a temporary work folder name, that does not already exist
do {
$TempFolder = rand(999999);
$TempFolder = "autoparse.$TempFolder";
} while (opendir DIRHANDLE, $TempFolder);
close DIRHANDLE;

# Make the Temporary Work Folder
mkdir $TempFolder;

# Open the Fasta file for input, take contents and close it
if (! open INPUTFASTA, "<$InputFasta") {
   rmdir $TempFolder;
   die "\nERROR: Could not open Input File.";
   }
@FastaContents=<INPUTFASTA>;
$FastaLineCount=@FastaContents;
$SeqCount=$FastaLineCount;
close INPUTFASTA;

#Make the name for Output file and open it
$rpos = rindex $InputFasta,"\.";
if ($rpos != -1) {$OutputFile=substr($InputFasta,0,$rpos);}
else {
$OutputFile=$InputFasta;
}

$OutputFile2="$OutputFile.pars.performance";
$OutputFile= "$OutputFile.pars.pass";

$OutputFileCount=1;
print "Opening file $OutputFile.$OutputFileCount\n";
if (! open OUTPUTFILE, ">$OutputFile.$OutputFileCount") {
   rmdir $TempFolder;
   die "\nERROR: Could not open Output file $OutputFile.$OutputFileCount\n";
   }


#Change the folder for Temporary work
chdir $TempFolder;
$LParserLineCount=0;
$ParserLineCount=0;
@ParsedSequencesWaiting=('INIT');

#Run through the FASTA file contents
for ($i=0; $i < $FastaLineCount; $i++){#Fasta runthrough
$fileline1="";
$fileline1=$FastaContents[$i];

if ($fileline1 ne "") {#Non-empty line availability

   open TEMPFILE,">TempFile";#Opening Temporary File for MFold Run
   print TEMPFILE $fileline1;
   close TEMPFILE;#Temporary File made for MFold Run

   $CurrentSeq= $i+1;
   print "\n\n\n\nProcessing Sequence \($CurrentSeq of $FastaLineCount\)\: $fileline1";
   print "Calling parser4auto  ...\n";

   # Call MFold below, wait till the output line about dG is received, extract dG value and, then close the execution
   my $pid = open (AUTOPARSE,"\.\.\/parser4auto TempFile $NForcePairs $MinBases $MaxBases|") || die "\!\! COULD NOT EXECUTE parser4auto\!\! (Temporary Folder not removed.)\n";
   while (<AUTOPARSE>){#Begin trapping parser4auto results
         $parseroutput = $_;
  
         $aftermatch= "";
         if ($parseroutput =~ /PARSED SEQUENCE\:/) {
         $aftermatch=$';
         chomp($aftermatch);
         IfUniqueWrite();
         }
   }#End trapping MFold results
   close AUTOPARSE;# MFold Done
print "parser4auto successful!!\n";
   }#End of non-empty line availability
}#End of Fasta runthrough
FlushWaitList();
close OUTPUTFILE;

#Delete all files in Temporary Folder and Remove Folder

#Reset to the old folder, and Remove the temporary folder
chdir '..';
use File::Path;
rmtree $TempFolder;

#Print the total Time Details
print "\n\n\n\nNumber of Sequences scanned through the Parser: $FastaLineCount\n";
print "\n\n\n\nNumber of Sequences Output from the Parser: $LParserLineCount\n";
$prgfinish = new Benchmark;
$difference=timediff($prgstart,$prgfinish);
$perf=timestr($difference);
print "Performance:\n$perf\n\n";

#Write Performance
open OUTPUTFILE, ">$OutputFile2";
print OUTPUTFILE "Performance:\n$perf\n\n";
close OUTPUTFILE;
exit;

#Make IfUniqueWrite() subroutine
sub IfUniqueWrite {
if ($ParsedSequencesWaiting[0] eq 'INIT') {
$ParsedSequencesWaiting[0]=$aftermatch;
} else { #That is, the first sequence is not being written

#Check if this sequence is already there
$PSWCount=@ParsedSequencesWaiting;
$SeqFoundInPSW=0;
for ($tmp=0;$tmp<$PSWCount;$tmp++) {
    if ($ParsedSequencesWaiting[$tmp] eq $aftermatch) {$SeqFoundInPSW=1;last;}
    }#Check for duplicate entry over

if ($SeqFoundInPSW == 0) {#If the sequence is UNIQUE
   if ($PSWCount < $MaxListSize) {#Store this into the Waiting List
      $ParsedSequencesWaiting[$PSWCount]=$aftermatch;
      } else { #If $MaxListSize members are already in the waitlist, Put the first sequence to file
      print OUTPUTFILE "$ParsedSequencesWaiting[0]\n";
      $LParserLineCount++;
      $ParserLineCount++;
      if ($ParserLineCount >= 240002) {NextFile();}
      
      for ($tmp=0;$tmp<($MaxListSize-1);$tmp++) {#Shift all values up by one cell
          $ParsedSequencesWaiting[$tmp]=$ParsedSequencesWaiting[$tmp+1];
          }
      $ParsedSequencesWaiting[$MaxListSize-1]=$aftermatch; #Add this sequence to the end of the Waiting List
      }#Saving the sequence over
   }#If the sequence is unique- over
}#End of- The first sequence is not being written
}#End Sub


#Make FlushWaitList() subroutine
sub FlushWaitList(){
$PSWCount=@ParsedSequencesWaiting;
for ($tmp=0;$tmp<$PSWCount;$tmp++) {
    print OUTPUTFILE "$ParsedSequencesWaiting[$tmp]\n";
    $LParserLineCount++;
    $ParserLineCount++;
    if ($ParserLineCount >= 240002) {NextFile();}
                      }
}

#Open next file
sub NextFile(){
close OUTPUTFILE;

$OutputFileCount++;
#print "Next file\n";
#print "Opening file $OutputFile.$OutputFileCount\n";
if (! open OUTPUTFILE, ">../$OutputFile.$OutputFileCount"){
   rmdir $TempFolder;
   die "\nERROR: Could not open Output file $OutputFile.$OutputFileCount\n";
   }
$ParserLineCount=0;
}

