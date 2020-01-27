# Author: Junyang Yue
# Program name: blast2split.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl blast2split.pl
## 
## 
## USAGE: perl blast2split.pl
## ________________________________________________________________________________
##


# handle BLAST result #

open (IN, "./data/output/merge_blast.tab") || die ("\nError in $0: BLAST result file doesn't exist: $! !\n\n");
open (OUT, ">./data/output/blast_row1_row2.tab");
open (OUT1, ">./data/output/blast_row1.txt");
open (OUT2, ">./data/output/blast_row2.txt");

while (<IN>) {
	chomp;
	@line = split m/\t/, $_;
	chomp ($line[0]);
	chomp ($line[1]);
	chomp ($line[3]);
	chomp ($line[8]);
	chomp ($line[9]);
	if ($line[3] > 100) { #half of the select cut length, see $get_length in misa2grab.pl
		$id1 = substr ($line[0],0,-7);
		$id2 = substr ($line[1],0,-7);
		if ($line[8] < $line[9]) {
			print OUT $line[0]."\t".$line[1]."\t"."F"."\n"; #Forward
			print OUT1 $line[0]."\n";
			print OUT2 $line[1]."\n";
		} elsif ($line[8] > $line[9]) {
			print OUT $line[0]."\t".$line[1]."\t"."R"."\n"; #Reverse
			print OUT1 $line[0]."\n";
			print OUT2 $line[1]."\n";
		}
	}
}

close (IN);
close (OUT);
close (OUT1);
close (OUT2);
