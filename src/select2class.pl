# Author: Junyang Yue
# Program name: select2class.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl select2class.pl
## 
## 
## USAGE: perl select2class.pl
## ________________________________________________________________________________
##


# add SSR class #

open (CT1, "./src/ssr_class.ini") || die ("\nError in $0: The SSR class file doesn't exist: $! !\n\n");
open (CT2, "./data/merge.tab") || die ("\nError in $0: The merged file doesn't exist: $! !\n\n");

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
		chomp ($line2[4]);
		$line2[4] =~ s/\(([A-Z]+)\)\d+/$1/;
		$hash2{$line2[0]} = $line2[4];
	}
}
close (CT2);


open (IN, "./data/output/merge_blast_match_select.tab") || die ("\nError in $0: The selected file doesn't exist: $! !\n\n");
open (OUT, ">./data/output/merge_blast_match_select_class.tab");

my @in = <IN>;

foreach(@in) {
	chomp ($_);
	if ($_) {
		my $this_line = $_;
		my @list = split m/\t/, $this_line;
		chomp($list[0]);
		chomp($list[1]);
		print OUT $this_line."\t".$hash1{$hash2{@list[0]}}."\t".$hash1{$hash2{@list[1]}}."\n";
	}
}

close (IN);
close (OUT);
