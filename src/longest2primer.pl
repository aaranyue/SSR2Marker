# Author: Junyang Yue
# Program name: longest2primer.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl longest2primer.pl
## 
## 
## USAGE: perl longest2primer.pl
## ________________________________________________________________________________
##


# run Primer3 program #

open (SEQ, "./data/longest_identical_region.tab") || die ("\nError in $0: The result file with information of the longest identical regions doesn't exist: $! !\n\n");

@seq = <SEQ>;
my %hash;
foreach (@seq) {
	chomp;
	if ($_) {
		@line = split m/\t/, $_;
		chomp ($line[0]);
		chomp ($line[1]);
		$hash{$line[0]} = $line[1];
	}
}
close (SEQ);

open (IN, "./data/primer_prepare.id") || die ("\nError in $0: The BLAST result file with information of two species ID doesn't exist: $! !\n\n");
open (OUT, ">./data/primer_prepare.txt");

while (<IN>) {
	chomp;
	$id = $_;
	$up_stream = $_."U";
	$target = "N" x 100; #mimic the target sequences
	$down_stream = $_."D";
	$seq = $hash{$up_stream}.$target.$hash{$down_stream};
	$start = length($hash{$up_stream});
	print OUT "SEQUENCE_ID=".$id."\n"; #needed
	print OUT "SEQUENCE_TEMPLATE=".$seq."\n"; #needed
	print OUT "PRIMER_TASK=generic"."\n";
	print OUT "PRIMER_PICK_LEFT_PRIMER=1"."\n";
	print OUT "PRIMER_PICK_INTERNAL_OLIGO=0"."\n";
	print OUT "PRIMER_PICK_RIGHT_PRIMER=1"."\n";
	print OUT "PRIMER_OPT_SIZE=20"."\n";
	print OUT "PRIMER_MIN_SIZE=16"."\n";
	print OUT "PRIMER_MAX_SIZE=28"."\n";
	print OUT "PRIMER_PRODUCT_SIZE_RANGE=100-500"."\n";
	print OUT "PRIMER_EXPLAIN_NLAG=1"."\n";
	print OUT "SEQUENCE_TARGET=".$start.","."100"."\n"; #needed
	print OUT "="."\n"; #needed
}

close (IN);
close (OUT);
