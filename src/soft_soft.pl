# Author: Junyang Yue
# Program name: soft_soft.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl soft_soft.pl
## 
## 
## USAGE: perl soft_soft.pl
## ________________________________________________________________________________
##


# check installed softwares #

`perl -version > SSR2Marker.log 2>&1`;

open (LOG, "./SSR2Marker.log");

my @array_perl;

while (<LOG>) {
	chomp;
	if ($_) {
		push @array_perl, $_;
	}
}

$count_perl = @array_perl;

if (($count_perl <= 2) && ("@array_perl" =~ m/command not found/i)) {
	$message = "\nError: The Perl software is not found. You should specify the path or install it firstly !\n\n";
	die $message;
} else {
	print "The Perl software could be correctly called\. \n";
}

close (LOG);

system "rm -rf ./SSR2Marker.log";


`blastn -version > SSR2Marker.log 2>&1`;

open (LOG, "./SSR2Marker.log");

my @array_blast;

while (<LOG>) {
	chomp;
	if ($_) {
		push @array_blast, $_;
	}
}

$count_blast = @array_blast;

if (($count_blast <= 2) && ("@array_blast" =~ m/command not found/i)) {
	$message = "\nError: The BLAST+ software is not found. You should specify the path or install it firstly !\n\n";
	die $message;
} else {
	print "The BLAST+ software could be correctly called\. \n";
}

close (LOG);

system "rm -rf ./SSR2Marker.log";


`mafft --version > SSR2Marker.log 2>&1`;

open (LOG, "./SSR2Marker.log");

my @array_mafft;

while (<LOG>) {
	chomp;
	if ($_) {
		push @array_mafft, $_;
	}
}

$count_mafft = @array_mafft;

if (($count_mafft <= 2) && ("@array_mafft" =~ m/command not found/i)) {
	$message = "\nError: The MAFFT software is not found. You should specify the path or install it firstly !\n\n";
	die $message;
} else {
	print "The MAFFT software could be correctly called\. \n";
}

close (LOG);

system "rm -rf ./SSR2Marker.log";


`primer3_core -version > SSR2Marker.log 2>&1`;

open (LOG, "./SSR2Marker.log");

my @array_primer;

while (<LOG>) {
	chomp;
	if ($_) {
		push @array_primer, $_;
	}
}

$count_primer = @array_primer;

if (($count_primer <= 2) && ("@array_primer" =~ m/command not found/i)) {
	$message = "\nError: The Primer3 software is not found. You should specify the path or install it firstly !\n\n";
	die $message;
} else {
	print "The Primer3 software could be correctly called\. \n";
}

close (LOG);

system "rm -rf ./SSR2Marker.log";


`re-PCR -version > SSR2Marker.log 2>&1`;

open (LOG, "./SSR2Marker.log");

my @array_epcr;

while (<LOG>) {
	chomp;
	if ($_) {
		push @array_epcr, $_;
	}
}

$count_epcr = @array_epcr;

if (($count_epcr <= 2) && ("@array_epcr" =~ m/command not found/i)) {
	$message = "\nError: The NCBI-ePCR software is not found. You should specify the path or install it firstly !\n\n";
	die $message;
} else {
	print "The NCBI-ePCR software could be correctly called\. \n";
}

close (LOG);

system "rm -rf ./SSR2Marker.log";

