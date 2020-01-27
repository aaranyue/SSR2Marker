# Author: Junyang Yue
# Program name: misa2grab.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl misa2grab.pl [FASTA file name]
## 
##    [FASTA file name]     The file contains the sequences from a species in FASTA format.
## 
## USAGE: perl misa2grab.pl Hongyang.fasta
## ________________________________________________________________________________
##


# Check for arguments. If none display syntax #

if (@ARGV == 0) {
	open (IN,"<$0");
	while (<IN>) {if (/^\#\# (.*)/) {$message .= "$1\n"}};
	close (IN);
	die $message;
};

# build hash #
open (FA, "./data/$ARGV[0]") || die ("\nError in $0: FASTA file doesn't exist: $! !\n\n");

$/ = ">";
my @fa = <FA>;
my %hash;

foreach (@fa) {
	chomp;
	@line = split m/\n/, $_;
	chomp ($line[0]);
	$line[0] =~ s/\r//;
	chomp ($line[1]);
	$line[1] =~ s/\r//;
	$hash{$line[0]} = $line[1];
}

close (FA);

# grab sequences #

$filename = $ARGV[0];
$misa_filename = $ARGV[0].".misa";
$upstream_filename = $ARGV[0].".up";
$fulllength_filename = $ARGV[0].".full";
$downstream_filename = $ARGV[0].".down";
$unit_filename = $ARGV[0].".unit";
$filename =~ s/\.fasta//;

open (MI, "./data/$misa_filename") || die ("\nError in $0: MISA result file doesn't exist: $! !\n\n");
open (UP, ">./data/$upstream_filename");
open (FULL, ">./data/$fulllength_filename");
open (DOWN, ">./data/$downstream_filename");
open (UNIT, ">./data/$unit_filename");

$/ = "\n";

my @mi = <MI>;

my $get_length = 200;

foreach (@mi) {
	chomp;
	my $id = sprintf "%06d", $i+1;
	@list = split m/\t/, $_;
	chomp ($list[0]);
	chomp ($list[3]);
	chomp ($list[5]);
	chomp ($list[6]);
	chomp ($list[7]);
	if ($list[5] > $get_length) {
		print UP $filename.$id."U"."\t".$list[0]."\t".($list[5]-$get_length)."\t".($list[5]-1)."\t".$list[3]."\t".substr($hash{$list[0]},$list[5]-$get_length-1,$get_length)."\n";
		print FULL $filename.$id."\t".$list[0]."\t".$list[5]."\t".$list[6]."\t".$list[3]."\t"."$get_length+1"."\t".($list[6]-$list[5]+$get_length+1)."\t".substr($hash{$list[0]},$list[5]-$get_length-1,$list[6]-$list[5]+$get_length*2+1)."\n";
		print DOWN $filename.$id."D"."\t".$list[0]."\t".($list[6]+1)."\t".($list[6]+$get_length)."\t".$list[3]."\t".substr($hash{$list[0]},$list[6],$get_length)."\n";
	} elsif ($list[5] == 1) {
		print UP $filename.$id."U"."\t".$list[0]."\t"."0"."\t"."0"."\t".$list[3]."\t"."N"."\n";
		print FULL $filename.$id."\t".$list[0]."\t".$list[5]."\t".$list[6]."\t".$list[3]."\t".$list[5]."\t".$list[6]."\t".substr($hash{$list[0]},0,$list[6]+$get_length)."\n";
		print DOWN $filename.$id."D"."\t".$list[0]."\t".($list[6]+1)."\t".($list[6]+$get_length+1)."\t".$list[3]."\t".substr($hash{$list[0]},$list[6],$get_length)."\n";
	} else {
		print UP $filename.$id."U"."\t".$list[0]."\t"."1"."\t".($list[5]-1)."\t".$list[3]."\t".substr($hash{$list[0]},0,$list[5]-1)."\n";
		print FULL $filename.$id."\t".$list[0]."\t".$list[5]."\t".$list[6]."\t".$list[3]."\t".$list[5]."\t".$list[6]."\t".substr($hash{$list[0]},0,$list[6]+$get_length)."\n";
		print DOWN $filename.$id."D"."\t".$list[0]."\t".($list[6]+1)."\t".($list[6]+$get_length)."\t".$list[3]."\t".substr($hash{$list[0]},$list[6],$get_length)."\n";
	}
	print UNIT $filename.$id."\t".$list[5]."\t".$list[6]."\t".$list[3]."\t".$list[7]."\n";
	$i++;
	if ($i > 999999) {
		die ("\nError in $0: The system broke down under excessive loads due to the huge sequences !\nYou could split the input file and run this program again !\n\n");
	}
}

close (MI);
close (UP);
close (FULL);
close (DOWN);
close (UNIT);
