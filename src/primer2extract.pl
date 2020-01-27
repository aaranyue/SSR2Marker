# Author: Junyang Yue
# Program name: primer2extract.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:	 perl primer2extract.pl
## 
## 
## USAGE: perl primer2extract.pl
## ________________________________________________________________________________
##


# extract primer sequences #

open (IN, "./data/primer_bank.txt") || die ("\nError in $0: The Primer3 result file with doesn't exist: $! !\n\n");
open (OUT, ">./data/primer_extract.txt");

$/ = "=\n";

while (<IN>) {
	my ($id) = (/SEQUENCE_ID=(.+?\n)/);
	$id =~ s/\r//;
	$id =~ s/\n//;
	
	
	# ------------------------- Primer 1 ------------------------- #
	/PRIMER_LEFT_0_SEQUENCE=(.*)/ || do {next};
	my $info = $id."1\t$1\t";
	/PRIMER_LEFT_0_TM=(.*)/; $info .= "$1\t";

	/PRIMER_RIGHT_0_SEQUENCE=(.*)/;	$info .= "$1\t";
	/PRIMER_RIGHT_0_TM=(.*)/; $info .= "$1\n";
	
	
	# ------------------------- Primer 2 ------------------------- #
	/PRIMER_LEFT_1_SEQUENCE=(.*)/; $info .= $id."2\t$1\t";
	/PRIMER_LEFT_1_TM=(.*)/; $info .= "$1\t";
		
	/PRIMER_RIGHT_1_SEQUENCE=(.*)/;	$info .= "$1\t";
	/PRIMER_RIGHT_1_TM=(.*)/; $info .= "$1\n";
	
	
	# ------------------------- Primer 3 ------------------------- #
	/PRIMER_LEFT_2_SEQUENCE=(.*)/; $info .= $id."3\t$1\t";
	/PRIMER_LEFT_2_TM=(.*)/; $info .= "$1\t";
		
	/PRIMER_RIGHT_2_SEQUENCE=(.*)/;	$info .= "$1\t";
	/PRIMER_RIGHT_2_TM=(.*)/; $info .= "$1\n";
	
	
	# ------------------------- Primer 4 ------------------------- #
	/PRIMER_LEFT_3_SEQUENCE=(.*)/; $info .= $id."4\t$1\t";
	/PRIMER_LEFT_3_TM=(.*)/; $info .= "$1\t";
		
	/PRIMER_RIGHT_3_SEQUENCE=(.*)/;	$info .= "$1\t";
	/PRIMER_RIGHT_3_TM=(.*)/; $info .= "$1\n";
	
	
	# ------------------------- Primer 5 ------------------------- #
	/PRIMER_LEFT_4_SEQUENCE=(.*)/; $info .= $id."5\t$1\t";
	/PRIMER_LEFT_4_TM=(.*)/; $info .= "$1\t";
		
	/PRIMER_RIGHT_4_SEQUENCE=(.*)/;	$info .= "$1\t";
	/PRIMER_RIGHT_4_TM=(.*)/; $info .= "$1\n";
	
	
	print OUT $info;
};

close (IN);
close (OUT);
