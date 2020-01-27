# Author: Junyang Yue
# Program name: overall2statistics.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSRs among close species.
## 
## SYNTAX:   perl overall2statistics.pl [FASTA file 1] [FASTA file 2]
## 
##    [FASTA file 1]     The file contains the sequences from species A in FASTA format.
##    [FASTA file 2]     The file contains the sequences from species B in FASTA format.
## 
## USAGE: perl overall2statistics.pl Hongyang.fasta White.fasta
## ________________________________________________________________________________
##


# Check for arguments. If none display syntax #

if (@ARGV == 0) {
	open (IN,"<$0");
	while (<IN>) {if (/^\#\# (.*)/) {$message .= "$1\n"}};
	close (IN);
	die $message;
};


# prepare file name #

$species1 = $ARGV[0];
$species2 = $ARGV[1];
$species1 =~ s/\.fasta//;
$species2 =~ s/\.fasta//;


# carries on the statistics #

$lastline = `tail -n 1 ./data/primer_overall.txt`;

@lastline = split m/\t/, $lastline;
chomp ($lastline[0]);
chomp ($lastline[1]);

$primer_count = $lastline[0];
$ssr2marker_count = $lastline[1];

$primer_count =~ s/Primer0+//; #line 1 #title a
$primer_count =~ s/Primer//; #line 1 #title a
$ssr2marker_count =~ s/P4SSR0+//; #line 2 #title b
$ssr2marker_count =~ s/P4SSR//; #line 2 #title b

open (IN, "./data/primer_overall.txt") || die ("\nError in $0: The overview file doesn't exist: $! !\n\n");
open (IN1, "./data/primer_multimatch.txt") || die ("\nError in $0: The multi-match primer file doesn't exist: $! !\n\n");
open (IN2, "./data/primer_shortlong.txt") || die ("\nError in $0: The short long primer file doesn't exist: $! !\n\n");
open (OUT, ">./data/SSR2Marker.stat");
open (OUT1, ">./data/SSR_Class.stat");
open (OUT2, ">./data/Target_Sequence_Difference.stat");
open (OUT3, ">./data/SSR_Length_$species1.stat");
open (OUT4, ">./data/SSR_Chromosome_$species1.stat");
open (OUT5, ">./data/Target_Sequence_Length_$species1.stat");
open (OUT6, ">./data/SSR_Length_$species2.stat");
open (OUT7, ">./data/SSR_Chromosome_$species2.stat");
open (OUT8, ">./data/Target_Sequence_Length_$species2.stat");

my @in1 = <IN1>;
my $size_primer_multimatch = @in1; #title x

my @in2 = <IN2>;
my $size_primer_shortlong = @in2; #title x

print OUT "\#Here is the statistical results\.\n";
print OUT "\#The employed program is ssr2marker.pl\.\n\n";

print OUT "The total number of designed primer pairs is $primer_count\.\n"; #title a

print OUT "The eliminated number of primer pairs is $size_primer_multimatch\, due to multiple matches in at least one of these two genomes\.\n"; #title x
print OUT "The eliminated number of primer pairs is $size_primer_shortlong\, due to short matches in at least one of these two genomes\.\n"; #title x

print OUT "The total number of identified SSR regions is $ssr2marker_count\.\n"; #title b

print OUT1 "\#"."SSR classes shared by $species1 and $species2"."\t"."Count"."\n";
print OUT2 "\#"."Length difference of target sequences between $species1 and $species2"."\t"."Count"."\n";
print OUT3 "\#"."SSR length in $species1"."\t"."Count"."\n";
print OUT4 "\#"."Chromosome distribution of SSRs in $species1"."\t"."Count"."\n";
print OUT5 "\#"."Target sequence length in $species1"."\t"."Count"."\n";
print OUT6 "\#"."SSR length in $species2"."\t"."Count"."\n";
print OUT7 "\#"."Chromosome distribution of SSRs in $species2"."\t"."Count"."\n";
print OUT8 "\#"."Target sequence length in $species2"."\t"."Count"."\n";

my %ssr_class; #line 3
my %target_sequence_difference; #line 12
my %ssr_length_species1; #line 16
my %ssr_chromosome_species1; #line 17
my %target_sequence_length_species1; #line 24
my %ssr_length_species2; #line 28
my %ssr_chromosome_species2; #line 29
my %target_sequence_length_species2; #line 36

while (<IN>) {
	chomp;
	if (m/^\#/) {
		#nothing to do
	} elsif ($_) {
		@line = split m/\t/, $_;
		
		chomp ($line[2]); #SSR_Class
		chomp ($line[11]); #Target_Sequence_Difference
		chomp ($line[15]); #SSR_Length_species1
		chomp ($line[16]); #SSR_Chromosome_species1
		chomp ($line[23]); #Target_Sequence_Length_species1
		chomp ($line[27]); #SSR_Length_species2
		chomp ($line[28]); #SSR_Chromosome_species2
		chomp ($line[35]); #Target_Sequence_Length_species2
		
		if (not exists $ssr_class{$line[2]}) {
			$ssr_class{$line[2]} = 1;
		} else {
			$ssr_class{$line[2]}++;
		}
		
		if (not exists $target_sequence_difference{$line[11]}) {
			$target_sequence_difference{$line[11]} = 1;
		} else {
			$target_sequence_difference{$line[11]}++;
		}
		
		if (not exists $ssr_length_species1{$line[15]}) {
			$ssr_length_species1{$line[15]} = 1;
		} else {
			$ssr_length_species1{$line[15]}++;
		}
		
		if (not exists $ssr_chromosome_species1{$line[16]}) {
			$ssr_chromosome_species1{$line[16]} = 1;
		} else {
			$ssr_chromosome_species1{$line[16]}++;
		}
		
		if (not exists $target_sequence_length_species1{$line[23]}) {
			$target_sequence_length_species1{$line[23]} = 1;
		} else {
			$target_sequence_length_species1{$line[23]}++;
		}
		
		if (not exists $ssr_length_species2{$line[27]}) {
			$ssr_length_species2{$line[27]} = 1;
		} else {
			$ssr_length_species2{$line[27]}++;
		}
		
		if (not exists $ssr_chromosome_species2{$line[28]}) {
			$ssr_chromosome_species2{$line[28]} = 1;
		} else {
			$ssr_chromosome_species2{$line[28]}++;
		}
		
		if (not exists $target_sequence_length_species2{$line[35]}) {
			$target_sequence_length_species2{$line[35]} = 1;
		} else {
			$target_sequence_length_species2{$line[35]}++;
		}
	}
}

# only count #
foreach (sort keys %ssr_class) {
	if ($_) {
		print OUT1 $_."\t".$ssr_class{$_}."\n";
	}
}

my @keys_ssr_class = keys %ssr_class;
my $size_ssr_class = @keys_ssr_class; #title c

print OUT "The total number of SSR classes is $size_ssr_class\.\n"; #title c

# count difference #
foreach (sort {$a <=> $b} keys %target_sequence_difference) {
	print OUT2 $_."\t".$target_sequence_difference{$_}."\n";
}

if (values %target_sequence_difference > 0) {
	my @keys_target_sequence_difference = keys %target_sequence_difference;
	my $size_target_sequence_difference = @keys_target_sequence_difference; #title d
	
	my $count_target_sequence_difference = 0;
	
	foreach ($i=0;$i<=$#keys_target_sequence_difference;$i++) {
		$count_target_sequence_difference += $keys_target_sequence_difference[$i];
	}
	
	print OUT "The total number of target sequences with different length between $species1 and $species2 is $count_target_sequence_difference\.\n"; #title d
}

# only min and max #
foreach (sort {$a <=> $b} keys %ssr_length_species1) {
	print OUT3 $_."\t".$ssr_length_species1{$_}."\n";
}

$min_ssr_length_species1 = (sort {$a <=> $b} keys %ssr_length_species1)[0]; #title e
$max_ssr_length_species1 = (sort {$b <=> $a} keys %ssr_length_species1)[0]; #title e

my @keys_ssr_length_species1 = keys %ssr_length_species1;
my $size_ssr_length_species1 = @keys_ssr_length_species1;

# only count #
foreach (sort keys %ssr_chromosome_species1) {
	print OUT4 $_."\t".$ssr_chromosome_species1{$_}."\n";
}

my @keys_ssr_chromosome_species1 = keys %ssr_chromosome_species1;
my $size_ssr_chromosome_species1 = @keys_ssr_chromosome_species1; #title f

# only min and max #
foreach (sort {$a <=> $b} keys %target_sequence_length_species1) {
	print OUT5 $_."\t".$target_sequence_length_species1{$_}."\n";
}

$min_target_sequence_length_species1 = (sort {$a <=> $b} keys %target_sequence_length_species1)[0]; #title g
$max_target_sequence_length_species1 = (sort {$b <=> $a} keys %target_sequence_length_species1)[0]; #title g

my @keys_target_sequence_length_species1 = keys %target_sequence_length_species1;
my $size_target_sequence_length_species1 = @keys_target_sequence_length_species1;

# only min and max #
foreach (sort {$a <=> $b} keys %ssr_length_species2) {
	print OUT6 $_."\t".$ssr_length_species2{$_}."\n";
}

$min_ssr_length_species2 = (sort {$a <=> $b} keys %ssr_length_species2)[0]; #title h
$max_ssr_length_species2 = (sort {$b <=> $a} keys %ssr_length_species2)[0]; #title h

my @keys_ssr_length_species2 = keys %ssr_length_species2;
my $size_ssr_length_species2 = @keys_ssr_length_species2;

# only count #
foreach (sort keys %ssr_chromosome_species2) {
	print OUT7 $_."\t".$ssr_chromosome_species2{$_}."\n";
}

my @keys_ssr_chromosome_species2 = keys %ssr_chromosome_species2;
my $size_ssr_chromosome_species2 = @keys_ssr_chromosome_species2; #title i

# only min and max #
foreach (sort {$a <=> $b} keys %target_sequence_length_species2) {
	print OUT8 $_."\t".$target_sequence_length_species2{$_}."\n";
}

$min_target_sequence_length_species2 = (sort {$a <=> $b} keys %target_sequence_length_species2)[0]; #title j
$max_target_sequence_length_species2 = (sort {$b <=> $a} keys %target_sequence_length_species2)[0]; #title g

my @keys_target_sequence_length_species2 = keys %target_sequence_length_species2;
my $size_target_sequence_length_species2 = @keys_target_sequence_length_species2;

print OUT "The length of SSRs in $species1 ranging from $min_ssr_length_species1 to $max_ssr_length_species1\.\n"; #title e
print OUT "The total number of chromosomes or sequences with identified SSRs in $species1 is $size_ssr_chromosome_species1\.\n"; #title f
print OUT "The length of target sequences in $species1 ranging from $min_target_sequence_length_species1 to $max_target_sequence_length_species1\.\n"; #title g
print OUT "The length of SSRs in $species2 ranging from $min_ssr_length_species2 to $max_ssr_length_species2\.\n"; #title h
print OUT "The total number of chromosomes with identified SSRs in $species2 is $size_ssr_chromosome_species2\.\n"; #title i
print OUT "The length of target sequences in $species2 ranging from $min_target_sequence_length_species2 to $max_target_sequence_length_species2\.\n"; #title j

close (IN);
close (OUT);
close (OUT1);
close (OUT2);
close (OUT3);
close (OUT4);
close (OUT5);
close (OUT6);
close (OUT7);
close (OUT8);
