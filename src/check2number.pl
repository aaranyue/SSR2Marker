# Author: Junyang Yue
# Program name: check2number.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl check2number.pl
## 
## 
## USAGE: perl check2number.pl
## ________________________________________________________________________________
##


# count the number of available primer #

open (IN, "./data/primer_check.txt") || die ("\nError in $0: The checked primer file doesn't exist: $! !\n\n");
open (OUT, ">./data/primer_number.txt");

my %result;

while (<IN>) {
	chomp;
	@line = split m/\t/, $_;
	chomp ($line[0]);
	if ($line[0] =~ m/(.+\d{6}_.+\d{6})_[A-Z]\d/) {
		$id = $1;
		if (not exists $result{$id}){
			$result{$id} = 1;
		} else {
			$result{$id}++;
		}	
	}
}

foreach (sort keys %result) {
	print OUT $_."\t".$result{$_}."\n";
}

close (IN);
close (OUT);

open (NO, "./data/primer_number.txt");
open (ORD, ">./data/primer_order.txt");

while (<NO>) {
	chomp;
	my $order = sprintf "%05d", $i+1;
	$ssr2marker_id = "P4SSR".$order;
	print ORD $_."\t".$ssr2marker_id."\n";
	$i++;
}

close (NO);
close (ORD);
