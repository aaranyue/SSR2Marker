# Author: Junyang Yue
# Program name: count2match.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl count2match.pl
## 
## 
## USAGE: perl count2match.pl
## ________________________________________________________________________________
##


# match count result #

open (CT1, "./data/output/blast_row1.count") || die ("\nError in $0: The count file of species IDs doesn't exist: $! !\n\n");
open (CT2, "./data/output/blast_row2.count") || die ("\nError in $0: The count file of species IDs doesn't exist: $! !\n\n");

@ct1 = <CT1>;
my %hash1;
foreach (@ct1) {
	chomp ($_);
	if ($_) {
		my @line1 = split m/\t/, $_;
		chomp ($line1[0]);
		chomp ($line1[1]);
		$hash1{$line1[0]} = $line1[1];
	}
}
close (CT1);

@ct2 = <CT2>;
my %hash2;
foreach (@ct2) {
	chomp ($_);
	if ($_) {
		my @line2 = split m/\t/, $_;
		chomp ($line2[0]);
		chomp ($line2[1]);
		$hash2{$line2[0]} = $line2[1];
	}
}
close (CT2);

open (IN, "./data/output/blast_row1_row2.tab") || die ("\nError in $0: The BLAST result file of species IDs doesn't exist: $! !\n\n");
open (OUT, ">./data/output/merge_blast_match.tab");

my @in = <IN>;

foreach (@in) {
	chomp ($_);
	if ($_) {
		my $this_line = $_;
		my @list = split m/\t/, $this_line;
		chomp ($list[0]);
		chomp ($list[1]);
		print OUT $this_line."\t".$hash1{$list[0]}."\t".$hash2{$list[1]}."\n";
	}
}

close (IN);
close (OUT);
