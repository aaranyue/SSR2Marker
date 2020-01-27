# Author: Junyang Yue
# Program name: grab2blast.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl grab2blast.pl [FASTA file 1] [FASTA file 2]
## 
##    [FASTA file 1]     The file contains the sequences from species A in FASTA format.
##    [FASTA file 2]     The file contains the sequences from species B in FASTA format.
## 
## USAGE: perl grab2blast.pl Hongyang.fasta White.fasta
## ________________________________________________________________________________
##


# Check for arguments. If none display syntax #

if (@ARGV == 0) {
	open (IN,"<$0");
	while (<IN>) {if (/^\#\# (.*)/) {$message .= "$1\n"}};
	close (IN);
	die $message;
};


# merge sequences #

system "cat ./data/*.up ./data/*.down > ./data/merge.tab";
system "cat ./data/*.unit > ./data/merge.unit";

# prepare FASTA file #

$filename1 = $ARGV[0];
$filename2 = $ARGV[1];
$filename1 =~ s/\.fasta//;
$filename2 =~ s/\.fasta//;
$filename_new1 = $filename1."Flank";
$filename_new2 = $filename2."Flank";

open (IN, "./data/merge.tab") || die ("\nError in $0: The merged file doesn't exist: $! !\n\n");
open (OUT1, ">./data/$filename_new1.fasta");
open (OUT2, ">./data/$filename_new2.fasta");

while (<IN>) {
	chomp;
	@line = split m/\t/, $_;
	chomp ($line[0]);
	chomp ($line[5]);
	if ($line[5] =~ m/^N+$/i) {
		#nothing to do
	} else {
		$line[5] =~ s/N//gi;
		$line[5] =~ s/(B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|U|V|W|X|Y|Z)//gi;
		if ($line[0] =~ m/$filename1/) {
			print OUT1 ">".$line[0]."\n".$line[5]."\n";
		} elsif ($line[0] =~ m/$filename2/) {
			print OUT2 ">".$line[0]."\n".$line[5]."\n";
		}
	}
}

close (IN);
close (OUT);

# run FASTA program #

system "rm -rf ./data/db/";
system "rm -rf ./data/input/";
system "rm -rf ./data/output/";

mkdir "./data/db/";
mkdir "./data/input/";
mkdir "./data/output/";

use File::Copy;
copy("./data/$filename_new1.fasta","./data/db/");
copy("./data/$filename_new2.fasta","./data/input/");

$cpu_count = `cat /proc/cpuinfo| grep "processor"| wc -l`;
chomp ($cpu_count);

system "makeblastdb -in ./data/db/$filename_new1.fasta -parse_seqids -hash_index -dbtype nucl";
system "blastn -db ./data/db/$filename_new1.fasta -query ./data/input/$filename_new2.fasta -out ./data/output/merge_blast.tab -evalue 1e-5 -outfmt 6 -perc_identity 80 -max_target_seqs 5 -num_threads $cpu_count";
