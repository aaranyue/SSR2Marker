# Author: Junyang Yue
# Program name: extract2epcr.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl extract2epcr.pl [FASTA file name]
## 
##    [FASTA file name]     The name of file contains the sequences from a species in FASTA format.
## 
## USAGE: perl extract2epcr.pl Hongyang
## ________________________________________________________________________________
##


# Check for arguments. If none display syntax #

if (@ARGV == 0) {
	open (IN,"<$0");
	while (<IN>) {if (/^\#\# (.*)/) {$message .= "$1\n"}};
	close (IN);
	die $message;
};


# run ePCR program #

open (IN, "./data/primer_extract.txt") || die ("\nError in $0: The Primer3 result file with extracted information doesn't exist: $! !\n\n");

while (<IN>) {
	chomp;
	@line = split m/\t/, $_;
	chomp ($line[0]);
	chomp ($line[1]);
	chomp ($line[2]);
	chomp ($line[3]);
	chomp ($line[4]);
	$filename = $line[0].".".$ARGV[0];
	system "re-PCR -s ./data/epcr/$ARGV[0].hash -n 1 -g 1 $line[1] $line[3] 50-1000 > ./data/epcr/output/$filename"
}

close (IN);
