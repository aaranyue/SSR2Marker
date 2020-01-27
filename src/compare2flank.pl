# Author: Junyang Yue
# Program name: compare2flank.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl compare2flank.pl
## 
## 
## USAGE: perl compare2flank.pl
## ________________________________________________________________________________
##


# obtain SSR flanking sequences #

open (IN, "./data/output/merge_blast_match_select_class_compare.tab") || die ("\nError in $0: The BLAST result file with the same class information doesn't exist: $! !\n\n");
open (OUT, ">./data/output/merge_blast_match_select_class_compare_flank.tab");

my @in = <IN>;
my %hash;
my @name;
my @value;

foreach (@in) {
	chomp;
	if ($_) {
		my @line = split m/\t/, $_;
		chomp ($line[0]);
		chomp ($line[1]);
		push @name, $line[0];
		push @value, $line[1];
	}
}

foreach (0..$#value) {
    $hash{$name[$_]} = $hash{$name[$_]}."_".$value[$_];
}

my @name_new;
my @value_new;

foreach (sort keys %hash) {
    push @name_new, $_;
    push @value_new, $hash{$_};
}

for ($i=0; $i<=$#value_new; $i++) {
	$value_new[$i] =~ s/^_//;
	if ($value_new[$i] =~ m/_/) {
		print OUT $name_new[$i]."\t".$value_new[$i]."\n";
	}
}

close (IN);
close (OUT);
