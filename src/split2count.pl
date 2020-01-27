# Author: Junyang Yue
# Program name: split2count.pl


## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl split2count.pl filename
## 
## 
## USAGE: perl split2count.pl blast_row1
## ________________________________________________________________________________
##


# count BLAST result #

$fileinput = $ARGV[0].".txt";
$fileoutput = $ARGV[0].".count";
$species = $ARGV[0];
$species =~ s/\.fasta//;

open (IN, "./data/output/$fileinput") || die ("\nError in $0: The BLAST result file of $species doesn't exist: $! !\n\n");
open (OUT, ">./data/output/$fileoutput");

my %result;

while (<IN>) {
	chomp;
	$line = $_;
	if (not exists $result{$line}){
		$result{$line} = 1;
	}
	else {
		$result{$line}++;
	}
}

foreach (keys %result) {
	if ($result{$_} == 1){
		print OUT $_."\t"."1"."\n";
	}
	else {
		print OUT $_."\t".$result{$_}."\n";
	}
}

close (IN);
close (OUT);
