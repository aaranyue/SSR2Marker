#!/usr/bin/perl -w
# Program name: SSR2Marker.pl
# Author: Junyang Yue
# Email: aaran.yue@gmail.com
# Thank you for using SSR2Marker. If you encount any problems or find some bugs, 
# please contact us through the Email address.


### _______________________________________________________________________________
### 
### Program name: SSR2Marker.pl
### Author:       Junyang Yue
### Release date: 18/11/2019
### Version:      Version 1.0
### _______________________________________________________________________________
### 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among 
##              close species.
## 
## SYNTAX:      perl SSR2Marker.pl [FASTA file 1] [FASTA file 2] [MISA file 1] 
##              [MISA file 2]
## 
##    [FASTA file 1]    Single file in FASTA format containing the sequence(s) from 
##                      species A.
##    [FASTA file 2]    Single file in FASTA format containing the sequence(s) from 
##                      species B.
##    [MISA file 1]     The MISA result of species A produced by the MISA program. 
##    [MISA file 2]     The MISA result of species B produced by the MISA program.
##    [-help|-h]        Further information.
##    [-version|-v]     Version descriptions.
## 
## NOTICE:      The FASTA files are necessary. The MISA files are optional. If 
##              users have existing MISA results, they are suggested to be 
##              provided in the current step. This will save about 17% operation 
##              times to obtain the SSR primers. While the MISA results were not 
##              provided, SSR2Marker will begin with performing the MISA analysis 
##              until obtaining the SSR primers. What the users need to know is 
##              that the MISA results may present some difference due to the 
##              different parameter settings. Anyway, both the FASTA and MISA 
##              (if provided) files from two species should be provided together.
## 
## USAGE:       perl SSR2Marker.pl Hongyang.fasta White.fasta
##              or
##              perl SSR2Marker.pl Hongyang.fasta White.fasta Hongyang.fasta.misa 
##              White.fasta.misa
## ________________________________________________________________________________
##


use File::Copy;


# Check for arguments. If none display syntax #

if (@ARGV == 0) {
	$message = "\nError: More information is needed for running this program !\nYou may type '-help' or '-h' for further information\n\n";
	die $message;
};


# Check if help is required #

if (($ARGV[0] eq '-help')|($ARGV[0] eq '-h')) {
	open (IN,"<$0");
	while (<IN>) {if (/^\#\# (.*)/) {$message .= "$1\n\n"}};
	close (IN);
	die $message;
};


# Check if version is required #

if (($ARGV[0] eq '-version')|($ARGV[0] eq '-v')) {
	open (IN,"<$0");
	while (<IN>) {if (/^\#\#\# (.*)/) {$message .= "$1\n\n"}};
	close (IN);
	die $message;
};

print "Beginning to check the dependent softwares needed in this program \.\.\. \n"; #label

system "perl ./src/soft_soft.pl";

my $file = "./SSR2Marker.log";
if (-e $file) {
	die ("\nThe dependent software is missing. Please install it firstly $! !\n\n");
} else {
	print "All dependent softwares could be correctly called here\. \n"; #label
}


# Open input files, including FASTA or MISA files #

$arg = @ARGV;

if (($arg == 2) && ($ARGV[0] =~ m/\.fasta/) && ($ARGV[1] =~ m/\.fasta/)) {
	if ($ARGV[0] eq $ARGV[1]) {
		$message = "\nError: Two different FASTA files are needed for running this program !\n\n";
		die $message;
	} else {
		system "rm -rf ./data";
		print "Preparing the working directory \.\.\. \n"; #label
	
		mkdir "data";
		print "The working directory is successfully built\. \n"; #label
	
		$species1 = $ARGV[0];
		$species1 =~ s/\.fasta//;
		$species2 = $ARGV[1];
		$species2 =~ s/\.fasta//;
		
		print "Copy the FASTA files of species $species1 and $species2 to the \'data\' working directory \.\.\. \n"; #label
		system "perl ./src/fasta_fasta.pl $ARGV[0]";
		system "perl ./src/fasta_fasta.pl $ARGV[1]";
		print "Finishing the FASTA files copy\. \n"; #label
		
# ----------------------------- This is different line -----------------------------
		print "Beginning to identify the SSRs in $species1 species \.\.\. \n"; #label
		system "perl ./src/misa.pl $ARGV[0]";
		print "Finishing the identification of SSRs in $species1 species\. \n"; #label
	
		print "Beginning to identify the SSRs in $species2 species \.\.\. \n"; #label
		system "perl ./src/misa.pl $ARGV[1]";
		print "Finishing the identification of SSRs in $species2 species\. \n"; #label
# ----------------------------- This is different line -----------------------------
		
		system "perl ./src/handle2misa.pl $ARGV[0].misa";
		system "perl ./src/handle2misa.pl $ARGV[1].misa";
		
		print "Grabbing the SSR sequences as well as their upstream and downstream sequences of species $species1 and $species2 \.\.\. \n"; #label
		system "perl ./src/misa2grab.pl $ARGV[0]";
		system "perl ./src/misa2grab.pl $ARGV[1]";
		print "SSR, upstream and downstream sequences are obtained\. \n"; #label
	
		print "Beginning to run the BLASTN program \.\.\. \n"; #label
		system "perl ./src/grab2blast.pl $ARGV[0] $ARGV[1]";
		print "Finishing the BLASTN program\. \n"; #label
	
		print "Beginning to identify single-copy SSRs \.\.\. \n"; #label
		system "perl ./src/blast2split.pl";
		system "perl ./src/split2count.pl blast_row1";
		system "perl ./src/split2count.pl blast_row2";
		system "perl ./src/count2match.pl";
		system "perl ./src/match2select.pl";
		print "Finishing the identification of single-copy SSRs\. \n"; #label
	
		print "Beginning to identify the identical regions between $species1 and $species2 species within upstream and downstream sequences, respectively \.\.\. \n"; #label
		system "perl ./src/select2class.pl";
		system "perl ./src/class2compare.pl";
		system "perl ./src/compare2flank.pl";
		
		print "Beginning to perform BLAST program \.\.\. \n"; #label
		system "mkdir ./data/fasta";
		system "perl ./src/flank2fasta.pl";
		
		print "Beginning to run MAFFT program for the identical regions between $species1 and $species2 species \.\.\. \n"; #label
		system "mkdir ./data/mafft";
		system "sh ./src/fasta2mafft.sh" || die ("\nError in $0: MAFFT running is broken: $! !\n\n");
		system "perl ./src/mafft2identical.pl";
		print "Finishing the identification of identical regions\. \n"; #label
		
		print "Keeping the longest identical regions for the design of common primers \.\.\. \n"; #label
		system "perl ./src/identical2longest.pl";
		system "perl ./src/longest2primer.pl";
		
		print "Beginning to run Primer3 program for the design of common primers \.\.\. \n"; #label
		system "primer3_core ./data/primer_prepare.txt > ./data/primer_bank.txt" || die ("\nError in $0: Primer3 running is broken: $! !\n\n");
		print "Finishing the design of common primers\. \n"; #label
		
		print "Using NCBI-ePCR to verify sequence tagged sites (STSs) within genome sequences of $species1 and $species2 species \.\.\. \n"; #label
		system "perl ./src/primer2extract.pl";
		system "mkdir ./data/epcr";
		system "mkdir ./data/epcr/output";
		system "mv ./data/$ARGV[0] ./data/epcr/";
		system "mv ./data/$ARGV[1] ./data/epcr/";
		system "famap -t n -b ./data/epcr/$species1.famap ./data/epcr/$species1.fasta" || die ("\nError in $0: ePCR running is broken: $! !\n\n");
		system "famap -t n -b ./data/epcr/$species2.famap ./data/epcr/$species2.fasta" || die ("\nError in $0: ePCR running is broken: $! !\n\n");
		system "fahash -b ./data/epcr/$species1.hash -w 12 -f 3 ./data/epcr/$species1.famap" || die ("\nError in $0: ePCR running is broken: $! !\n\n");
		system "fahash -b ./data/epcr/$species2.hash -w 12 -f 3 ./data/epcr/$species2.famap" || die ("\nError in $0: ePCR running is broken: $! !\n\n");
		system "perl ./src/extract2epcr.pl $species1";
		system "perl ./src/extract2epcr.pl $species2";
		system "perl ./src/epcr2eliminate.pl";
		system "perl ./src/eliminate2check.pl";
		print "Five pairs of high-quality common primers are current obtained\. \n"; #label
		
		system "mv ./data/epcr/$ARGV[0] ./data/";
		system "mv ./data/epcr/$ARGV[1] ./data/";
		
		print "Beginning to organize and count all results \.\.\. \n"; #label
		system "perl ./src/check2number.pl";
		system "perl ./src/number2overall.pl $ARGV[0] $ARGV[1]";
		system "perl ./src/overall2statistics.pl $ARGV[0] $ARGV[1]";
		
		my $date = &getTime();
		my $worktime = $date -> {date}; #obtain localtime yyyymmdd
		my $resultfile = "result".$worktime;
		system "rm -rf ./$resultfile";
		system "mkdir ./$resultfile";
		
		system "cp ./data/primer_overall.txt ./$resultfile/SSR2Marker.txt";
		system "cp ./data/*.stat ./$resultfile/";
		
		open (RD, "./$resultfile/SSR2Marker.stat") || die ("\nError in $0: $! !\n\n");
		
		while (<RD>) {
			if (m/^\#/) {
				#nothing to do
			} else {
				print $_;
			}
		}
		
		close (RD);
		
		print "\nAll SSR Primers results as well as their corresponding statistics have already been generated by SSR2Marker\. \n"; #label
		print "You can view them in the folder of $resultfile now\. \n"; #label
	}
# ----------------------------- This is demarcation line -----------------------------
# ----------------------------- This is demarcation line -----------------------------
# ----------------------------- This is demarcation line -----------------------------
} elsif (($arg == 4) && ($ARGV[0] =~ m/\.fasta/) && ($ARGV[1] =~ m/\.fasta/) && ($ARGV[2] =~ m/\.misa/) && ($ARGV[3] =~ m/\.misa/)) {
	if ($ARGV[0] eq $ARGV[1]) {
		$message = "\nError: Two different FASTA files are needed for running this program !\n\n";
		die $message;
	} elsif ($ARGV[2] eq $ARGV[3]) {
		$message = "\nError: Two different MISA files are needed for running this program !\n\n";
		die $message;
	} else {
		system "rm -rf ./data";
		print "Preparing the working directory \.\.\. \n"; #label
	
		mkdir "data";
		print "The working directory is successfully built\. \n"; #label
	
		$species1 = $ARGV[0];
		$species1 =~ s/\.fasta//;
		$species2 = $ARGV[1];
		$species2 =~ s/\.fasta//;
		
		print "Copy the FASTA files of species $species1 and $species2 to the \'data\' working directory \.\.\. \n"; #label
		system "perl ./src/fasta_fasta.pl $ARGV[0]";
		system "perl ./src/fasta_fasta.pl $ARGV[1]";
		print "Finishing the FASTA files copy\. \n"; #label
		
# ----------------------------- This is different line -----------------------------
		print "Beginning to handle the MISA result files in $species1 and $species2 species \.\.\. \n"; #label
		system "perl ./src/handle_misa.pl $ARGV[2]";
		system "perl ./src/handle_misa.pl $ARGV[3]";
		print "Finishing the processing of MISA results in $species1 and $species2 species\. \n"; #label
# ----------------------------- This is different line -----------------------------
		
		system "perl ./src/handle2misa.pl $ARGV[2]";
		system "perl ./src/handle2misa.pl $ARGV[3]";
	
		print "Grabbing the SSR sequences as well as their upstream and downstream sequences of species $species1 and $species2 \.\.\. \n"; #label
		system "perl ./src/misa2grab.pl $ARGV[0]";
		system "perl ./src/misa2grab.pl $ARGV[1]";
		print "SSR, upstream and downstream sequences are obtained\. \n"; #label
	
		print "Beginning to run the BLASTN program \.\.\. \n"; #label
		system "perl ./src/grab2blast.pl $ARGV[0] $ARGV[1]";
		print "Finishing the BLASTN program\. \n"; #label
	
		print "Beginning to identify single-copy SSRs \.\.\. \n"; #label
		system "perl ./src/blast2split.pl";
		system "perl ./src/split2count.pl blast_row1";
		system "perl ./src/split2count.pl blast_row2";
		system "perl ./src/count2match.pl";
		system "perl ./src/match2select.pl";
		print "Finishing the identification of single-copy SSRs\. \n"; #label
	
		print "Beginning to identify the identical regions between $species1 and $species2 species within upstream and downstream sequences, respectively \.\.\. \n"; #label
		system "perl ./src/select2class.pl";
		system "perl ./src/class2compare.pl";
		system "perl ./src/compare2flank.pl";
		
		print "Beginning to perform BLAST program \.\.\. \n"; #label
		system "mkdir ./data/fasta";
		system "perl ./src/flank2fasta.pl";
		
		print "Beginning to run MAFFT program for the identical regions between $species1 and $species2 species \.\.\. \n"; #label
		system "mkdir ./data/mafft";
		system "sh ./src/fasta2mafft.sh" || die ("\nError in $0: MAFFT running is broken: $! !\n\n");
		system "perl ./src/mafft2identical.pl";
		print "Finishing the identification of identical regions\. \n"; #label
		
		print "Keeping the longest identical regions for the design of common primers \.\.\. \n"; #label
		system "perl ./src/identical2longest.pl";
		system "perl ./src/longest2primer.pl";
		
		print "Beginning to run Primer3 program for the design of common primers \.\.\. \n"; #label
		system "primer3_core ./data/primer_prepare.txt > ./data/primer_bank.txt" || die ("\nError in $0: Primer3 running is broken: $! !\n\n");
		print "Finishing the design of common primers\. \n"; #label
		
		print "Using NCBI-ePCR to verify sequence tagged sites (STSs) within genome sequences of $species1 and $species2 species \.\.\. \n"; #label
		system "perl ./src/primer2extract.pl";
		system "mkdir ./data/epcr";
		system "mkdir ./data/epcr/output";
		system "mv ./data/$ARGV[0] ./data/epcr/";
		system "mv ./data/$ARGV[1] ./data/epcr/";
		system "famap -t n -b ./data/epcr/$species1.famap ./data/epcr/$species1.fasta" || die ("\nError in $0: ePCR running is broken: $! !\n\n");
		system "famap -t n -b ./data/epcr/$species2.famap ./data/epcr/$species2.fasta" || die ("\nError in $0: ePCR running is broken: $! !\n\n");
		system "fahash -b ./data/epcr/$species1.hash -w 12 -f 3 ./data/epcr/$species1.famap" || die ("\nError in $0: ePCR running is broken: $! !\n\n");
		system "fahash -b ./data/epcr/$species2.hash -w 12 -f 3 ./data/epcr/$species2.famap" || die ("\nError in $0: ePCR running is broken: $! !\n\n");
		system "perl ./src/extract2epcr.pl $species1";
		system "perl ./src/extract2epcr.pl $species2";
		system "perl ./src/epcr2eliminate.pl";
		system "perl ./src/eliminate2check.pl";
		print "Five pairs of high-quality common primers are current obtained\. \n"; #label
		
		system "mv ./data/epcr/$ARGV[0] ./data/";
		system "mv ./data/epcr/$ARGV[1] ./data/";
		
		print "Beginning to organize and count all results \.\.\. \n"; #label
		system "perl ./src/check2number.pl";
		system "perl ./src/number2overall.pl $ARGV[0] $ARGV[1]";
		system "perl ./src/overall2statistics.pl $ARGV[0] $ARGV[1]";
		
		my $date = &getTime();
		my $worktime = $date -> {date}; #obtain localtime yyyymmdd
		my $resultfile = "result".$worktime;
		system "rm -rf ./$resultfile";
		system "mkdir ./$resultfile";
		
		system "cp ./data/primer_overall.txt ./$resultfile/SSR2Marker.txt";
		system "cp ./data/*.stat ./$resultfile/";
		
		open (RD, "./$resultfile/SSR2Marker.stat") || die ("\nError in $0: $! !\n\n");
		
		while (<RD>) {
			if (m/^\#/) {
				#nothing to do
			} else {
				print $_;
			}
		}
		
		close (RD);
		
		print "\nAll SSR Primers results as well as their corresponding statistics have already been generated by SSR2Marker\. \n"; #label
		print "You can view them in the folder of $resultfile now\. \n"; #label
	}
# ----------------------------- This is demarcation line -----------------------------
# ----------------------------- This is demarcation line -----------------------------
# ----------------------------- This is demarcation line -----------------------------
} else {
	$message = "\nError: Please enter correct arguments to run this program !\nYou may type '-help' or '-h' for further information\n\n";
	die $message;
}


sub getTime {
	my $time = shift || time();
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($time);
	$mon ++;
	$sec  = ($sec<10) ? "0".$sec:$sec;
	$min  = ($min<10) ? "0".$min:$min;
	$hour = ($hour<10) ? "0".$hour:$hour;
	$mday = ($mday<10) ? "0".$mday:$mday;
	$mon  = ($mon<10) ? "0".$mon:$mon;
	$year+=1900;
	
	my $weekday = ('Sun','Mon','Tue','Wed','Thu','Fri','Sat')[$wday];
	return {
		'second' => $sec,
		'minute' => $min,
		'hour'   => $hour,
		'day'    => $mday,
		'month'  => $mon,
		'year'   => $year,
		'weekNo' => $wday,
		'wday'   => $weekday,
		'yday'   => $yday,
		'date'   => "$year$mon$mday"
	};
}
