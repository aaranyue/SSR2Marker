# Author: Junyang Yue
# Program name: eliminate2check.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl eliminate2check.pl
## 
## 
## USAGE: perl eliminate2check.pl
## ________________________________________________________________________________
##


# check primer sequences #

open (IN, "./data/primer_eliminate.txt") || die ("\nError in $0: The Primer3 result file with single-copy primers doesn't exist: $! !\n\n");
open (OUT, ">./data/primer_check.txt"); #reserved
open (OUT1, ">./data/primer_multimatch.txt"); #eliminated, just for statistics

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
		#chomp ($line[2]);
		chomp ($line[3]);
		chomp ($line[4]);
		chomp ($line[5]);
		chomp ($line[6]);
		chomp ($line[7]);
		chomp ($line[8]);
		chomp ($line[9]);
		$long = $line[9];
		$long =~ s/\/50\-1000//;
		$list = $line[1].":".$line[3].":".$line[4].":".$line[5].":".$line[6].":".$line[7].":".$line[8].":".$long;
		push @name, $line[0];
		push @value, $list;
	}
}

foreach (0..$#value) {
    $hash{$name[$_]} = $hash{$name[$_]}."##".$value[$_];
}

my @name_new;
my @value_new;

foreach (sort keys %hash) {
    push @name_new, $_;
    push @value_new, $hash{$_};
}

for ($i=0; $i<=$#value_new; $i++) {
	$value_new[$i] =~ s/^\#\#//;
	$count = () = ($value_new[$i] =~ m/\#\#/g);
	if ($count == 1) {
		print OUT $name_new[$i]."\t".$value_new[$i]."\n";
	} else {
		print OUT1 $name_new[$i]."\t".$value_new[$i]."\n";
	}
}

close (IN);
close (OUT);
close (OUT1);
