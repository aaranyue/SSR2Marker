# Author: Junyang Yue
# Program name: handle2misa.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl handle2misa.pl [MISA file name]
## 
##    [MISA file name]     The file contains the SSR results produced by the MISA program.
## 
## USAGE: perl handle2misa.pl Hongyang.fasta.misa
## ________________________________________________________________________________
##


# Check for arguments. If none display syntax #

if (@ARGV == 0) {
	open (IN,"<$0");
	while (<IN>) {if (/^\#\# (.*)/) {$message .= "$1\n"}};
	close (IN);
	die $message;
};

# change the unit format in MISA result file #

system "mv ./data/$ARGV[0] ./data/$ARGV[0].bak";

open (IN, "./data/$ARGV[0].bak") || die ("\nError in $0: MISA result file doesn't exist: $! !\n\n");
open (OUT, ">./data/$ARGV[0]");

while (<IN>) {
	chomp;
	$line = $_;
	$line =~ s/\*//g;
	@list = split m/\t/, $line;
	chomp ($list[0]); #chromosome id
	chomp ($list[1]); #serial number
	chomp ($list[2]); #class
	chomp ($list[3]); #unit
	chomp ($list[4]); #length
	chomp ($list[5]); #start
	chomp ($list[6]); #end
	if ($list[2] eq "c") {
		print OUT $list[0]."\t".$list[1]."\t".$list[2]."\t"."\(N\)$list[4]"."\t".$list[4]."\t".$list[5]."\t".$list[6]."\t".$list[3]."\n";
	} else {
		print OUT $line."\t".$list[3]."\n";
	}
}

close (IN);
close (OUT);
