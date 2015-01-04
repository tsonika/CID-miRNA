#!/usr/bin/perl

#use strict;

#Check the number of arguments
if (scalar(@ARGV) < 2)
	{
	die "\nInsufficient number of arguments!!\n\nProgramName <Source Contig File> <Destination Contig File>\nThis program converts a broken FASTA sequence to continuous FASTA sequece\n\n";
	}
	else
	{
	$InpFile=$ARGV[0];
	$OutFile=$ARGV[1];
	}


#Open the files for reading
print "\nOpening input files ($InpFile) for reading...\n";
open INPCONTIG, "<$InpFile" || die "Could not open file: $InpFile\n";
open OUTCONTIG, ">$OutFile" || die "Could not open file: $OutFile\n";


#Pickup the Entire Contig
$Contig="";
while ($ContigLine=<INPCONTIG>)
	{
	chomp $ContigLine;
#print "$ContigLine";
	if ($ContigLine !~ /^\>/)
		{# If the file line does not have the line beginning with '>'
		$Contig .= $ContigLine;
		}
        else
                {
                if ($Contig ne "") {print OUTCONTIG "$Contig\n";}
                $Contig = "";

                print OUTCONTIG "$ContigLine\n";
                }
	}
close INPCONTIG;
#open  OUTCONTIGFILE, '>Modified_ContigFile.fa';
if ($Contig ne "") {print OUTCONTIG "$Contig\n";}

#print OUTCONTIG ">$InpFile\n",$Contig,"\n";
close OUTCONTIG;

