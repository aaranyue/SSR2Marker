# Author: Junyang Yue
# Program name: fasta_fasta.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl fasta_fasta.pl [FASTA file name]
## 
##    [FASTA file name]     The file contains the sequences from a species in FASTA format.
## 
## USAGE: perl fasta_fasta.pl Hongyang.fasta
## ________________________________________________________________________________
##


# Check for arguments. If none display syntax #

if (@ARGV == 0) {
	open (IN,"<$0");
	while (<IN>) {if (/^\#\# (.*)/) {$message .= "$1\n"}};
	close (IN);
	die $message;
};

# change the format of FASTA file #
open (IN, "$ARGV[0]") || die ("\nError in $0: The FASTA file doesn't exist: $! !\n\n");
open (OUT, ">./data/$ARGV[0]");

$/ = ">";

while (<IN>) {
	chomp;
	if ($_ =~ m/(.+?)\n(.+)/sx) {
		$id = $1;
		$seq = $2;
		$seq =~ s/\r//g;
		$seq =~ s/\n//g;
		$sequence = uc ($seq);
		print OUT ">".$id."\n".$sequence."\n";
	}
}

close (IN);
close (OUT);
