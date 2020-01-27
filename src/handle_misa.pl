# Author: Junyang Yue
# Program name: handle_misa.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl handle_misa.pl [MISA file name]
## 
##    [MISA file name]     The file contains the SSR results produced by the MISA program.
## 
## USAGE: perl handle_misa.pl Hongyang.fasta.misa
## ________________________________________________________________________________
##


# Check for arguments. If none display syntax #

if (@ARGV == 0) {
	open (IN,"<$0");
	while (<IN>) {if (/^\#\# (.*)/) {$message .= "$1\n"}};
	close (IN);
	die $message;
};

# delete the annotation line #

open (IN, "./$ARGV[0]") || die ("\nError in $0: MISA result file doesn't exist: $! !\n\n");
open (OUT, ">./data/$ARGV[0]");

while (<IN>) {
	chomp;
	if (m/^\#/) {
		#nothing to do
	} else {
		print OUT $_."\n";
	}
}

close (IN);
close (OUT);
