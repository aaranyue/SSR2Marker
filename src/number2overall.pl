# Author: Junyang Yue
# Program name: number2overall.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl number2overall.pl [FASTA file 1] [FASTA file 2]
## 
##    [FASTA file 1]     The file contains the sequences from species A in FASTA format.
##    [FASTA file 2]     The file contains the sequences from species B in FASTA format.
## 
## USAGE: perl number2overall.pl Hongyang.fasta White.fasta
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


# build hash #

open (HA0, "./data/primer_order.txt") || die ("\nError in $0: The primer number file doesn't exist: $! !\n\n");

@ha0 = <HA0>;

my %hash0; #SSR primer count
my %hash0id; #SSR primer id

foreach (@ha0) {
	chomp ($_);
	if ($_) {
		my @line0 = split m/\t/, $_;
		chomp ($line0[0]);
		chomp ($line0[1]);
		chomp ($line0[2]);
		$hash0{$line0[0]} = $line0[1];
		$hash0id{$line0[0]} = $line0[2];
	}
}

close (HA0);

open (HA1, "./data/merge.unit") || die ("\nError in $0: The merged file doesn't exist: $! !\n\n");

@ha1 = <HA1>;

my %hash1; #SSR start
my %hash2; #SSR end
my %hash3; #SSR unit
my %hash4; #SSR repeat
my %hash3ori; #original SSR unit

foreach (@ha1) {
	chomp ($_);
	if ($_) {
		my @line1 = split m/\t/, $_;
		chomp ($line1[0]);
		chomp ($line1[1]);
		chomp ($line1[2]);
		chomp ($line1[3]);
		chomp ($line1[4]);
		
		my $unit = $line1[3];
		my $repeat = $line1[3];
		$unit =~ s/\(([A-Z]+)\)\d+/$1/;
		$repeat =~ s/\([A-Z]+\)(\d+)/$1/;
		
		$hash1{$line1[0]} = $line1[1];
		$hash2{$line1[0]} = $line1[2];
		$hash3{$line1[0]} = $unit;
		$hash4{$line1[0]} = $repeat;
		$hash3ori{$line1[0]} = $line1[4];
	}
}

close (HA1);

open (HA5, "./src/ssr_class.ini") || die ("\nError in $0: The SSR class file doesn't exist: $! !\n\n");

@ha5 = <HA5>;

my %hash5; #SSR class

foreach (@ha5) {
	chomp ($_);
	if ($_) {
		my @line5 = split m/\t/, $_;
		chomp ($line5[0]);
		chomp ($line5[1]);
		$hash5{$line5[0]} = $line5[1];
	}
}

close (HA5);

open (HA6, "./data/primer_extract.txt") || die ("\nError in $0: The extracted primer file doesn't exist: $! !\n\n");

@ha6 = <HA6>;

my %hash6; #SSR forward primer
my %hash7; #SSR forward Tm
my %hash8; #SSR reverse primer
my %hash9; #SSR reverse Tm

foreach (@ha6) {
	chomp ($_);
	if ($_) {
		my @line6 = split m/\t/, $_;
		chomp ($line6[0]);
		chomp ($line6[1]);
		chomp ($line6[2]);
		chomp ($line6[3]);
		chomp ($line6[4]);
		$hash6{$line6[0]} = $line6[1];
		$hash7{$line6[0]} = $line6[2];
		$hash8{$line6[0]} = $line6[3];
		$hash9{$line6[0]} = $line6[4];
	}
}

close (HA6);

open (HAA, "./data/$ARGV[0]") || die ("\nError in $0: FASTA file doesn't exist: $! !\n\n");

$/ = ">";
my @haa = <HAA>;
my %hasha; #species 1 or A

foreach (@haa) {
	chomp;
	@linea = split m/\n/, $_;
	chomp ($linea[0]);
	$linea[0] =~ s/\r//;
	chomp ($linea[1]);
	$linea[1] =~ s/\r//;
	$hasha{$linea[0]} = $linea[1];
}

close (HAA);

open (HAB, "./data/$ARGV[1]") || die ("\nError in $0: FASTA file doesn't exist: $! !\n\n");

$/ = ">";
my @hab = <HAB>;
my %hashb; #species 2 or B

foreach (@hab) {
	chomp;
	@lineb = split m/\n/, $_;
	chomp ($lineb[0]);
	$lineb[0] =~ s/\r//;
	chomp ($lineb[1]);
	$lineb[1] =~ s/\r//;
	$hashb{$lineb[0]} = $lineb[1];
}

close (HAB);


# organize into a file including all necessary results #

$/ = "\n";

open (IN, "./data/primer_check.txt") || die ("\nError in $0: The checked primer file doesn't exist: $! !\n\n");
open (OUT, ">./data/primer_overall.txt"); # totally 36 rows

print OUT "\#"."Primer_ID"."\t"."SSR2Marker_ID"."\t"."SSR_Class"."\t"."Primer_Pair_Count"."\t"."Primer_Pair_Number"."\t"."Primer_Forward"."\t"."Primer_Forward_Length"."\t"."Primer_Forward_Tm"."\t"."Primer_Reverse"."\t"."Primer_Reverse_Length"."\t"."Primer_Reverse_Tm"."\t"."Target_Sequence_Difference"."\t"."SSR_ID_".$species1."\t"."SSR_Unit_".$species1."\t"."SSR_Repeat_".$species1."\t"."SSR_Length_".$species1."\t"."SSR_Chromosome_".$species1."\t"."SSR_Chromosome_Direction_".$species1."\t"."SSR_Chromosome_Start_".$species1."\t"."SSR_Chromosome_End_".$species1."\t"."Target_Chromosome_Start_".$species1."\t"."Target_Chromosome_End_".$species1."\t"."Target_Seqence_".$species1."\t"."Target_Sequence_Length_".$species1."\t"."SSR_ID_".$species2."\t"."SSR_Unit_".$species2."\t"."SSR_Repeat_".$species2."\t"."SSR_Length_".$species2."\t"."SSR_Chromosome_".$species2."\t"."SSR_Chromosome_Direction_".$species2."\t"."SSR_Chromosome_Start_".$species2."\t"."SSR_Chromosome_End_".$species2."\t"."Target_Chromosome_Start_".$species2."\t"."Target_Chromosome_End_".$species2."\t"."Target_Seqence_".$species2."\t"."Target_Sequence_Length_".$species2."\n"; #title

while (<IN>) {
	chomp;
	my $order = sprintf "%06d", $i+1;
	
	$primer_id = "Primer".$order; # row 1
	
	my @line = split m/\t/, $_;
	chomp ($line[0]);
	chomp ($line[1]);
	
	$primer_forward = $hash6{$line[0]}; # row 6
	$primer_forward_tm = $hash7{$line[0]}; # row 8
	$primer_reverse = $hash8{$line[0]}; # row 9
	$primer_reverse_tm = $hash9{$line[0]}; # row 11
	
	$primer_forward_length = length($primer_forward); # row 7
	$primer_reverse_length = length($primer_reverse); # row 10
	
	
	if ($line[0] =~ m/(.+\d{6})_(.+\d{6})_([A-Z])(\d)/) {
		
		my $id = $1."_".$2;
		
		$ssr2marker_id = $hash0id{$id}; # row 2
		$primer_pair_count = $hash0{$id}; # row 4
		
		$ssr_id_species1 = $2; # row 13
		$ssr_id_species2 = $1; # row 25
		
		$primer_pair_number = $4; # row 5
		
		@list = split m/\#\#/, $line[1];
		chomp ($list[0]);
		chomp ($list[1]);
		
		@a = split m/:/, $list[0];
		@b = split m/:/, $list[1];
		
		chomp ($a[0]);
		chomp ($a[1]);
		chomp ($a[2]);
		chomp ($a[3]);
		chomp ($a[4]);
		chomp ($a[7]);
		
		chomp ($b[0]);
		chomp ($b[1]);
		chomp ($b[2]);
		chomp ($b[3]);
		chomp ($b[4]);
		chomp ($b[7]);
		
		if ($species1 eq $a[0]) {
			$ssr_chromosome_species1 = $a[1]; # row 17
			$ssr_chromosome_direction_species1 = $a[2]; # row 18
			$target_chromosome_start_species1 = $a[3]; # row 21
			$target_chromosome_end_species1 = $a[4]; # row 22
			$target_sequence_length_species1 = $a[7]; # row 24
		
			$ssr_chromosome_species2 = $b[1]; # row 29
			$ssr_chromosome_direction_species2 = $b[2]; # row 30
			$target_chromosome_start_species2 = $b[3]; # row 33
			$target_chromosome_end_species2 = $b[4]; # row 34
			$target_sequence_length_species2 = $b[7]; # row 36
			
			$target_seqence_species1 = substr($hasha{$ssr_chromosome_species1},$target_chromosome_start_species1-1,$target_chromosome_end_species1-$target_chromosome_start_species1+1); # row 23
			$target_seqence_species2 = substr($hashb{$ssr_chromosome_species2},$target_chromosome_start_species2-1,$target_chromosome_end_species2-$target_chromosome_start_species2+1); # row 35
		} else {
			$ssr_chromosome_species1 = $b[1]; # row 17
			$ssr_chromosome_direction_species1 = $b[2]; # row 18
			$target_chromosome_start_species1 = $b[3]; # row 21
			$target_chromosome_end_species1 = $b[4]; # row 22
			$target_sequence_length_species1 = $b[7]; # row 24
		
			$ssr_chromosome_species2 = $a[1]; # row 29
			$ssr_chromosome_direction_species2 = $a[2]; # row 30
			$target_chromosome_start_species2 = $a[3]; # row 33
			$target_chromosome_end_species2 = $a[4]; # row 34
			$target_sequence_length_species2 = $a[7]; # row 36
			
			$target_seqence_species1 = substr($hasha{$ssr_chromosome_species1},$target_chromosome_start_species1-1,$target_chromosome_end_species1-$target_chromosome_start_species1+1); # row 23
			$target_seqence_species2 = substr($hashb{$ssr_chromosome_species2},$target_chromosome_start_species2-1,$target_chromosome_end_species2-$target_chromosome_start_species2+1); # row 35
		}
		
		$target_sequence_difference = abs ($target_sequence_length_species1 - $target_sequence_length_species2); # row 12
		
		$ssr_chromosome_start_species1 = $hash1{$ssr_id_species1}; # row 19
		$ssr_chromosome_end_species1 = $hash2{$ssr_id_species1}; # row 20
		
		$ssr_chromosome_start_species2 = $hash1{$ssr_id_species2}; # row 31
		$ssr_chromosome_end_species2 = $hash2{$ssr_id_species2}; # row 32
		
		$ssr_unit_species1 = $hash3{$ssr_id_species1}; # row 14
		
		if ($ssr_unit_species1 eq 'N') {
			$ssr_unit_species1 = $hash3ori{$ssr_id_species1}; # row 14
			$ssr_repeat_species1 = 1; # row 15
			$ssr_length_species1 = abs ($ssr_chromosome_end_species1 - $ssr_chromosome_start_species1); # row 16
			$ssr_class = "N"; # row 3, only the same selected
		} else {
			$ssr_repeat_species1 = $hash4{$ssr_id_species1}; # row 15
			$ssr_length_species1 = length($ssr_unit_species1) * $ssr_repeat_species1; # row 16. Difference can also be used, but here use count.
			
			$ssr_class = $hash5{$ssr_unit_species1}; # row 3, only the same selected, the same of species1 and species2, only assign once
		}
		
		$ssr_unit_species2 = $hash3{$ssr_id_species2}; # row 26
		
		if ($ssr_unit_species2 eq 'N') {
			$ssr_unit_species2 = $hash3ori{$ssr_id_species2}; # row 26
			$ssr_repeat_species2 = 1; # row 27
			$ssr_length_species2 = abs ($ssr_chromosome_end_species2 - $ssr_chromosome_start_species2) # row 28
		} else {
			$ssr_repeat_species2 = $hash4{$ssr_id_species2}; # row 27
			$ssr_length_species2 = length($ssr_unit_species2) * $ssr_repeat_species2; # row 28. Difference can also be used, but here use count.
		}
		
		#$ssr_unit_species2 = $hash3{$ssr_id_species2}; # row 26
		#$ssr_repeat_species2 = $hash4{$ssr_id_species2}; # row 27
		#$ssr_length_species2 = length($ssr_unit_species2) * $ssr_repeat_species2; # row 28
		
		#$ssr_class = $hash5{$ssr_unit_species1}; # row 3, only the same selected
	}
	
	print OUT $primer_id."\t".$ssr2marker_id."\t".$ssr_class."\t".$primer_pair_count."\t".$primer_pair_number."\t".$primer_forward."\t".$primer_forward_length."\t".$primer_forward_tm."\t".$primer_reverse."\t".$primer_reverse_length."\t".$primer_reverse_tm."\t".$target_sequence_difference."\t".$ssr_id_species1."\t".$ssr_unit_species1."\t".$ssr_repeat_species1."\t".$ssr_length_species1."\t".$ssr_chromosome_species1."\t".$ssr_chromosome_direction_species1."\t".$ssr_chromosome_start_species1."\t".$ssr_chromosome_end_species1."\t".$target_chromosome_start_species1."\t".$target_chromosome_end_species1."\t".$target_seqence_species1."\t".$target_sequence_length_species1."\t".$ssr_id_species2."\t".$ssr_unit_species2."\t".$ssr_repeat_species2."\t".$ssr_length_species2."\t".$ssr_chromosome_species2."\t".$ssr_chromosome_direction_species2."\t".$ssr_chromosome_start_species2."\t".$ssr_chromosome_end_species2."\t".$target_chromosome_start_species2."\t".$target_chromosome_end_species2."\t".$target_seqence_species2."\t".$target_sequence_length_species2."\n"; #result
	
	$i++;
	if ($i > 999999) {
		die ("\nError in $0: The system broke down under excessive loads due to the huge SSR sequences !\nYou could split the input file and run this program again !\n\n");
	}
}

close (IN);
close (OUT);
