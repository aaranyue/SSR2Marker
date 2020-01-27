# Author: Junyang Yue
# Program name: match2select.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl match2select.pl
## 
## 
## USAGE: perl match2select.pl
## ________________________________________________________________________________
##


# select single copy #

open (IN, "./data/output/merge_blast_match.tab") || die ("\nError in $0: The matched file with count information doesn't exist: $! !\n\n");
open (OUT, ">./data/output/merge_blast_match_select.tab");

while (<IN>) {
	chomp;
	@line = split m/\t/, $_;
	chomp ($line[0]);
	chomp ($line[1]);
	chomp ($line[2]);
	chomp ($line[3]);
	chomp ($line[4]);
	if (($line[3] == 1) && ($line[4] == 1)) {
		print OUT $line[0]."\t".$line[1]."\t".$line[2]."\n";
	}
}

close (IN);
close (OUT);
