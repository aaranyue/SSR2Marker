# Author: Junyang Yue
# Program name: class2compare.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl class2compare.pl
## 
## 
## USAGE: perl class2compare.pl
## ________________________________________________________________________________
##


# obtain SSR class #

open (IN, "./data/output/merge_blast_match_select_class.tab") || die ("\nError in $0: The BLAST result file with class information doesn't exist: $! !\n\n");
open (OUT, ">./data/output/merge_blast_match_select_class_compare.tab");

while (<IN>) {
	chomp;
	@line = split m/\t/, $_;
	chomp ($line[0]);
	chomp ($line[1]);
	chomp ($line[2]);
	chomp ($line[3]);
	chomp ($line[4]);
	if ($line[3] eq $line[4]) { #only reserve the same class
		$id1 = substr ($line[0],0,-1);
		$stream1 = substr ($line[0],-1,);
		$id2 = substr ($line[1],0,-1);
		$stream2 = substr ($line[1],-1,);
		print OUT $id1."-".$id2."_".$line[2]."-".$line[3]."\t".$stream1."-".$stream2."\n";
	}
}

close (IN);
close (OUT);