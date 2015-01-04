#!/usr/bin/perl


$PGM = $0;
$PGM =~ s#.*/##;
$usage = <<USAGE;               # usage message
  USAGE:
        $PGM -i  <input file name>


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
     	 		print "$_\n";
    		} elsif ($_ !~ /^[ACGTUacgtu]/) {
      			$string = $_;
      			$count = 0;       # Count of the pattern
      			$pos = 0;
      			$beforematch=$aftermatch=$pattern = "";
      			while ($string =~ /\).*.\(/g) {
				$count++;
      				#  print STDERR "before match = $`, After match = $'\nPattern matched:$&\n";
				$pattern = $&;
				$beforematch = $`;
				$aftermatch = $';
      			} # while looks for the multiple loops 
    			#  print STDERR "count:$count\n";	
      			$subpatbegin = $subpatend =  0;
      			if ($count >0) {  # if multiple loops found
#        			if ($pattern =~ /^\)+/) {
					$subpatbegin = $pattern =~ s/\)/\)/g;
#                                        print  "$pattern has $subpatbegin )\n";
 #         				$subpatbegin = $&;
  #        				$subpatbegin =~ s/\)/\(/g;
   #     			}
	
#				if($pattern =~ /\(+$/) {
                                        $subpatend = $pattern =~ s/\(/\(/g;
 #         				$subpatend = $&;
  #        				$subpatend =~ s/\(/\)/g;
   #     			}
        			if($subpatbegin>0 and $subpatend>0) {
	  				$beforematch = reverse $beforematch;
#print "$beforematch\t$subpatbegin\n";
	  				for($x=0; $x<$subpatbegin; $x++) {
 	     					$beforematch =~ s/\(/\./;
	  				} #end for
#print "$beforematch\n";
	  				$beforematch = reverse $beforematch;
#print "$aftermatch\t$subpatend\n";
	  				for($x=0; $x<$subpatend; $x++) {
 	     					$aftermatch =~ s/\)/\./;
	  				} #end for
#print "$aftermatch\n";
          				$pattern =~ s/./\./g;
          				$finalResult = $beforematch.$pattern.$aftermatch;
	  				print STDOUT "$finalResult\n";
 				} #if length
      			} else {    # if single loop, print as it is
        			print STDOUT "$_\n";
      			}
    		} else { 
			print STDOUT "$_\n";   #if not start with [acgt]
		}
  	} # if not a blank or comment line

} # end while 
close GM;
