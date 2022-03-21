#!/usr/bin/perl

# Program name: SSR2Marker.pl
# Release date: 28/01/2022
# Current version: Version 1.0

# Author: Junyang Yue
# Email: yuejy@ahau.edu.cn
# Thank you for using SSR2Marker. If you encounter any problems or find some bugs, 
# please contact us through the Email address.


#Load dependent modules

use strict;
use File::Copy;
use FileHandle;
use File::Basename qw/basename/;
use Cwd 'abs_path';
use threads;
use Getopt::Long;
use POSIX qw/mktime/;
no strict 'refs';

#Clear buffer
$| = 1;

#Current version information

my $current_version = 'SSR2Marker v1.0';

#Set default options

my ($fasta1, $fasta2, $misa1, $misa2, $motif, $flanking_sequence_length, $output_name, $version, $help);

our @suffix_list = qw(.gz .fa .fq .fasta .fastq .genome .tab .txt);

#Fetch options from command line

GetOptions (
	'fasta1|f1=s'        =>   \$fasta1,
	'fasta2|f2=s'        =>   \$fasta2,
	'misa1|m1=s'         =>   \$misa1,
	'misa2|m2=s'         =>   \$misa2,
	'motif|m!'           =>   \$motif,
	'flanking|f=i'       =>   \$flanking_sequence_length,
    'out|o=s'            =>   \$output_name,
	'version|v!'         =>   \$version,
	'help|h!'            =>   \$help,
) or die;

#Describe program information

my $usage = <<__USAGE__;
###############################################################################
Name:
  SSR2Marker - Identify candidate SSRs that are well suited for use as 
  molecular markers.
Description:
  An integrated pipeline for identification of SSR markers between two
  genome-scale sequences.
Usage:
  perl SSR2Marker.pl option1 <value1> option2 <value2> ... optionN <valueN>
Options:
  -f1|-fasta1      <str> : Single file in FASTA format containing the
                           sequence(s) from species A (must be provided).
  -f2|-fasta2      <str> : Single file in FASTA format containing the
                           sequence(s) from species B (must be provided).
  -m1|-misa1       <str> : The results including SSR information genrated from
                           'fasta1' file by the MISA program.
  -m2|-misa2       <str> : The results including SSR information genrated from
                           'fasta2' file by the MISA program.
  -m|-motif        <str> : Setting the motifs for SSR identification through
                           navigational operations.
  -f|-flanking     <int> : Length of the flanking sequences between each SSR
                           locus (must be an integer and larger than 100 (bp),
                           default: 200).
  -v|-version            : Showing the version information.
  -h|-help               : Showing the help information.
###############################################################################
__USAGE__

#Check the options

my $message = <<__GUIDE__;
\nYou may type '-help' or '-h' for further information !
__GUIDE__

die "$current_version\n" if $version;
die $usage if $help;
die "Error: options '-f1' and '-f2' must be provided !$message" if ! $fasta1 or ! $fasta2;
die "Error: options '-f1' and '-f2' must be different !$message" if $fasta1 eq $fasta2;
die "Error: option '-m2' must be provided !\n" if $misa1 and ! $misa2;
die "Error: option '-m1' must be provided !\n" if $misa2 and ! $misa1;
die "Error: options '-m1' and '-m2' must be different !$message" if $misa1 and $misa1 eq $misa2;
die "Error: option '-f' must be an integer larger than 100 !$message" if ($flanking_sequence_length && ($flanking_sequence_length <= 99)) or $flanking_sequence_length =~ /\.\d*[1-9]+/;

#Obtain folder basename
my $abs_fasta1 = abs_path ($fasta1);
my $base_fasta1 = basename ($fasta1);
my $prefix_fasta1 = basename ($fasta1, @suffix_list);
my $abs_fasta2 = abs_path ($fasta2);
my $base_fasta2 = basename ($fasta2);
my $prefix_fasta2 = basename ($fasta2, @suffix_list);

#Check the dependent software
print "Beginning to check the dependent softwares needed in this program \.\.\. \n"; #label
&checkInstalledSoft ("perl");
&checkInstalledSoft ("blastn");
&checkInstalledSoft ("mafft");
&checkInstalledSoft ("primer3_core");
&checkInstalledSoft ("re-PCR");

my $file_log = "./soft_check.log";
if (-e $file_log) {
	die ("\nThe dependent software is missing. Please install it firstly $! !\n\n");
} else {
	print "All dependent softwares could be correctly called here\. \n"; #label
}

#Obtain thread setting
our $thread = `cat /proc/cpuinfo| grep "processor"| wc -l`;
chomp ($thread);

#Obtain SSR motif
my $motifs;
if ($motif) {
	$motifs = &setSSRMotif;
} else {
	$motifs = '1=10,2=6,3=5,4=5,5=5,6=5';
}
my $flanking_sequence_lengths;
if ($flanking_sequence_length) {
	$flanking_sequence_lengths = $flanking_sequence_length;
} else {
	$flanking_sequence_lengths = 200;
}


#========================================#
#        Start main program ...          #
#========================================#

my $date_start = &getLocalTime();
my $worktime_detail = $date_start -> {date_detail}; #obtain localtime yyyymmddhhmmss
our $workfile = "Running_SSR2Marker_".$worktime_detail;

system "rm -rf ./$workfile/";
print "\nPreparing the working directory \.\.\. \n"; #label
mkdir "./$workfile/";
print "The working directory is successfully built\. \n"; #label
	
if (! -e -f "$base_fasta1") {
	system "ln -s $abs_fasta1 $base_fasta1";
}
if (! -e -f "$base_fasta2") {
	system "ln -s $abs_fasta2 $base_fasta2";
}

&copyFastaFile ("./$fasta1", "./$workfile/$prefix_fasta1.fasta", "./$workfile/$prefix_fasta1.rename", 6);
&copyFastaFile ("./$fasta2", "./$workfile/$prefix_fasta2.fasta", "./$workfile/$prefix_fasta2.rename", 6);

if ($misa1) {
	&changeMISAFormat ($misa1, "./$workfile/$prefix_fasta1.rename", "./$workfile/$prefix_fasta1.misa");
	&changeMISAFormat ($misa2, "./$workfile/$prefix_fasta2.rename", "./$workfile/$prefix_fasta2.misa");
} else {
	print "Beginning to identify the SSRs in the input files \.\.\. \n"; #label
	&runMISAProgram ("./$workfile/$prefix_fasta1.fasta", "./$workfile/$prefix_fasta1.misa", $motifs);
	&runMISAProgram ("./$workfile/$prefix_fasta2.fasta", "./$workfile/$prefix_fasta2.misa", $motifs);
	print "Finishing the identification of SSRs\. \n"; #label
}

print "Grabbing the SSR sequences as well as their upstream and downstream sequences \.\.\. \n"; #label
&handleMISAResult ("./$workfile/$prefix_fasta1.fasta", $flanking_sequence_lengths, 6);
&handleMISAResult ("./$workfile/$prefix_fasta2.fasta", $flanking_sequence_lengths, 6);

system "cat ./$workfile/*.up ./$workfile/*.down > ./$workfile/merge.flanking";
system "cat ./$workfile/$prefix_fasta1.unit ./$workfile/$prefix_fasta2.unit > ./$workfile/merge.unit";
if (-s "./$workfile/merge.flanking") {
	$/ = "\n";
	open (IN, "./$workfile/merge.flanking") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT1, ">./$workfile/$prefix_fasta1.flanking.fasta");
	open (OUT2, ">./$workfile/$prefix_fasta2.flanking.fasta");
	while (<IN>) {
		chomp;
		if ($_) {
			my @line = split m/\t/, $_;
			chomp ($line[0]);
			chomp ($line[5]);
			if ($line[5] !~ m/^N+$/i) {
				$line[5] =~ s/N//gi;
				$line[5] =~ s/(B|D|E|F|H|I|J|K|L|M|O|P|Q|R|S|U|V|W|X|Y|Z)//gi;
				if ($line[0] =~ m/$prefix_fasta1/) {
					print OUT1 ">".$line[0]."\n".$line[5]."\n";
				} elsif ($line[0] =~ m/$prefix_fasta2/) {
					print OUT2 ">".$line[0]."\n".$line[5]."\n";
				}
			}
		}
	}
	close (IN);
	close (OUT);
}

system "rm -rf ./temp/";
mkdir "./temp/";
mkdir "./temp/db/";
mkdir "./temp/input/";
mkdir "./temp/output/";

copy ("./$workfile/$prefix_fasta1.flanking.fasta", "./temp/db/reference.fasta");
copy ("./$workfile/$prefix_fasta2.flanking.fasta", "./temp/input/query.fasta");

&runBlastProgram ("./temp/db/reference.fasta", "./temp/input/query.fasta", "./temp/output/both_flanking_seq.blast", "1e-5", $thread);
copy ("./temp/output/both_flanking_seq.blast", "./$workfile/");

print "Beginning to identify single-copy SSRs \.\.\. \n"; #label
if (-s "./$workfile/both_flanking_seq.blast") {
	$/ = "\n";
	open (IN, "./$workfile/both_flanking_seq.blast") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">./$workfile/both_flanking_seq.blast_1_2.tab");
	open (OUT1, ">./$workfile/both_flanking_seq.blast_1.tab");
	open (OUT2, ">./$workfile/both_flanking_seq.blast_2.tab");
	while (<IN>) {
		chomp;
		my @line = split m/\t/, $_;
		chomp ($line[0]);
		chomp ($line[1]);
		chomp ($line[3]);
		chomp ($line[8]);
		chomp ($line[9]);
		if ($line[3] > $flanking_sequence_lengths/2) { #half of the select cut length
			if ($line[8] < $line[9]) {
				print OUT $line[0]."\t".$line[1]."\t"."F"."\n"; #Forward
				print OUT1 $line[0]."\n";
				print OUT2 $line[1]."\n";
			} elsif ($line[8] > $line[9]) {
				print OUT $line[0]."\t".$line[1]."\t"."R"."\n"; #Reverse
				print OUT1 $line[0]."\n";
				print OUT2 $line[1]."\n";
			}
		}
	}
	close (IN);
	close (OUT);
	close (OUT1);
	close (OUT2);
}

&countOneRow ("./$workfile/both_flanking_seq.blast_1.tab", "./$workfile/both_flanking_seq.blast_1.count");
&countOneRow ("./$workfile/both_flanking_seq.blast_2.tab", "./$workfile/both_flanking_seq.blast_2.count");

if (-s "./$workfile/both_flanking_seq.blast_1_2.tab") {
	$/ = "\n";
	open (CT1, "./$workfile/both_flanking_seq.blast_1.count") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (CT2, "./$workfile/both_flanking_seq.blast_2.count") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	my @ct1 = <CT1>;
	my %hash1;
	foreach (@ct1) {
		chomp ($_);
		if ($_) {
			my @line1 = split m/\t/, $_;
			chomp ($line1[0]);
			chomp ($line1[1]);
			$hash1{$line1[0]} = $line1[1];
		}
	}
	close (CT1);
	my @ct2 = <CT2>;
	my %hash2;
	foreach (@ct2) {
		chomp ($_);
		if ($_) {
			my @line2 = split m/\t/, $_;
			chomp ($line2[0]);
			chomp ($line2[1]);
			$hash2{$line2[0]} = $line2[1];
		}
	}
	close (CT2);
	open (IN, "./$workfile/both_flanking_seq.blast_1_2.tab") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">./$workfile/both_flanking_seq.blast_1_2.tab.match");
	my @in = <IN>;
	foreach (@in) {
		chomp ($_);
		if ($_) {
			my $this_line = $_;
			my @list = split m/\t/, $this_line;
			chomp ($list[0]);
			chomp ($list[1]);
			print OUT $this_line."\t".$hash1{$list[0]}."\t".$hash2{$list[1]}."\n";
		}
	}
	close (IN);
	close (OUT);
}

if (-s "./$workfile/both_flanking_seq.blast_1_2.tab.match") {
	$/ = "\n";
	open (IN, "./$workfile/both_flanking_seq.blast_1_2.tab.match") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">./$workfile/both_flanking_seq.blast_1_2.tab.singlecopy");
	while (<IN>) {
		chomp;
		my @line = split m/\t/, $_;
		chomp ($line[0]);
		chomp ($line[1]);
		chomp ($line[2]);
		chomp ($line[3]);
		chomp ($line[4]);
		if (($line[3] == 1) && ($line[4] == 1)) {
			print OUT $line[0]."\t".$line[1]."\t".$line[2]."\n";
		}
	}
	close (IN);
	close (OUT);
}
print "Finishing the identification of single-copy SSRs\. \n"; #label

if (-s "./$workfile/both_flanking_seq.blast_1_2.tab.singlecopy") {
	$/ = "\n";
	open (MG, "./$workfile/merge.flanking") || die ("\nError in $0: The merged file doesn't exist: $! !\n\n");
	my @mg = <MG>;
	my %hash_merge;
	foreach (@mg) {
		chomp ($_);
		if ($_) {
			my @line = split m/\t/, $_;
			chomp ($line[0]);
			chomp ($line[4]);
			$line[4] =~ s/\(([ATCGN]+)\)\d+/$1/;
			$hash_merge{$line[0]} = $line[4];
		}
	}
	close (MG);
	sub storeSSRClass;
	my $sub_hash = &storeSSRClass;
	
	open (IN, "./$workfile/both_flanking_seq.blast_1_2.tab.singlecopy") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">./$workfile/both_flanking_seq.blast_1_2.tab.singlecopy.match");
	my @in = <IN>;
	foreach (@in) {
		chomp;
		if ($_) {
			my $this_line = $_;
			my @list = split m/\t/, $this_line;
			chomp($list[0]);
			chomp($list[1]);
			print OUT $this_line."\t".$$sub_hash{$hash_merge{$list[0]}}."\t".$$sub_hash{$hash_merge{$list[1]}}."\n";
		}
	}
	close (IN);
	close (OUT);
}

if (-s "./$workfile/both_flanking_seq.blast_1_2.tab.singlecopy.match") {
	$/ = "\n";
	open (IN, "./$workfile/both_flanking_seq.blast_1_2.tab.singlecopy.match") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">./$workfile/both_flanking_seq.blast_1_2.tab.singlecopy.sameclass");
	while (<IN>) {
		chomp;
		my @line = split m/\t/, $_;
		chomp ($line[0]);
		chomp ($line[1]);
		chomp ($line[2]);
		chomp ($line[3]);
		chomp ($line[4]);
		my ($id1, $id2, $stream1, $stream2);
		if ($line[3] eq $line[4]) { #only reserve the same class
			$id1 = substr ($line[0],0,-1);
			$stream1 = substr ($line[0],-1,);
			$id2 = substr ($line[1],0,-1);
			$stream2 = substr ($line[1],-1,);
			print OUT $id1."-".$id2."_".$line[2]."-".$line[3]."\t".$stream1."-".$stream2."\n";
		}
	}
	close (IN);
	close (OUT);
}

if (-s "./$workfile/both_flanking_seq.blast_1_2.tab.singlecopy.sameclass") {
	$/ = "\n";
	open (IN, "./$workfile/both_flanking_seq.blast_1_2.tab.singlecopy.sameclass") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">./$workfile/both_flanking_seq.blast_1_2.tab.singlecopy.sameclass.flankinglabel");
	my @in = <IN>;
	my @name;
	my @value;
	foreach (@in) {
		chomp;
		if ($_) {
			my @line = split m/\t/, $_;
			chomp ($line[0]);
			chomp ($line[1]);
			push @name, $line[0];
			push @value, $line[1];
		}
	}
	my %hash;
	foreach (0..$#value) {
		$hash{$name[$_]} = $hash{$name[$_]}."_".$value[$_];
	}
	my @name_new;
	my @value_new;

	foreach (sort keys %hash) {
		push @name_new, $_;
		push @value_new, $hash{$_};
	}
	for (my $i=0; $i<=$#value_new; $i++) {
		$value_new[$i] =~ s/^_//;
		if ($value_new[$i] =~ m/_/) {
			print OUT $name_new[$i]."\t".$value_new[$i]."\n";
		}
	}
	close (IN);
	close (OUT);
}

mkdir "./temp/fasta/";
mkdir "./temp/mafft/";

if (-s "./$workfile/both_flanking_seq.blast_1_2.tab.singlecopy.sameclass.flankinglabel") {
	$/ = "\n";
	open (MG, "./$workfile/merge.flanking") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	my @mg = <MG>;
	my %hash_merge;
	foreach (@mg) {
		chomp ($_);
		if ($_) {
			my @line = split m/\t/, $_;
			chomp ($line[0]);
			chomp ($line[5]);
			$hash_merge{$line[0]} = $line[5];
		}
	}
	close (MG);
	
	open (IN, "./$workfile/both_flanking_seq.blast_1_2.tab.singlecopy.sameclass.flankinglabel") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">./temp/mafft.id");
	while (<IN>) {
		chomp;
		my @line = split m/\t/, $_;
		chomp ($line[0]);
		chomp ($line[1]);
		my @list0 = split m/_/, $line[0];
		my @list1 = split m/_/, $line[1];
		chomp ($list0[0]);
		chomp ($list0[1]);
		chomp ($list1[0]);
		chomp ($list1[1]);
		my @id = split m/-/, $list0[0]; #id
		chomp ($id[0]);
		chomp ($id[1]);
		my @di = split m/-/, $list0[1]; #direction
		chomp ($di[0]);
		chomp ($di[1]); #SSR class
		my @stream0 = split m/-/, $list1[0];
		chomp ($stream0[0]);
		chomp ($stream0[1]);
		my @stream1 = split m/-/, $list1[1];
		chomp ($stream1[0]);
		chomp ($stream1[1]);
		my $id1_up = $id[0].$stream0[0];
		my $id2_up = $id[1].$stream0[1];
		my $id1_down = $id[0].$stream1[0];
		my $id2_down = $id[1].$stream1[1];
		
		if ($di[0] eq "F") { #blast forward
			my $id1_up_seq = $hash_merge{$id1_up};
			my $id2_up_seq = $hash_merge{$id2_up};
			my $id1_down_seq = $hash_merge{$id1_down};
			my $id2_down_seq = $hash_merge{$id2_down};
			my $filenameFU = $id[0]."_".$id[1]."_"."FU"; #upstream
			my $filenameFD = $id[0]."_".$id[1]."_"."FD"; #downstream
			print OUT $filenameFU."\n";
			print OUT $filenameFD."\n";
			open (OUTFU, ">./temp/fasta/$filenameFU");
			open (OUTFD, ">./temp/fasta/$filenameFD");
			print OUTFU ">".$id1_up."\n".$id1_up_seq."\n".">".$id2_up."\n".$id2_up_seq."\n";
			print OUTFD ">".$id1_down."\n".$id1_down_seq."\n".">".$id2_down."\n".$id2_down_seq."\n";
		} elsif ($di[0] eq "R") { #blast reverse
			my $id1_up_seq = $hash_merge{$id1_up};
			my $id2_up_seq = $hash_merge{$id2_up};
			my $id2_up_seq_R = reverse ($id2_up_seq); #reverse_complementary_sequences
			$id2_up_seq_R =~ tr/ACGTUacgtt/TGCAAtgcaa/;
			my $id1_down_seq = $hash_merge{$id1_down};
			my $id2_down_seq = $hash_merge{$id2_down};
			my $id2_down_seq_R = reverse ($id2_down_seq); #reverse_complementary_sequences
			$id2_down_seq_R =~ tr/ACGTUacgtt/TGCAAtgcaa/;
			my $filenameRU = $id[0]."_".$id[1]."_"."RU"; #upstream
			my $filenameRD = $id[0]."_".$id[1]."_"."RD"; #downstream
			print OUT $filenameRU."\n";
			print OUT $filenameRD."\n";
			open (OUTRU, ">./temp/fasta/$filenameRU");
			open (OUTRD, ">./temp/fasta/$filenameRD");
			print OUTRU ">".$id1_up."\n".$id1_up_seq."\n".">".$id2_up."\n".$id2_up_seq_R."\n";
			print OUTRD ">".$id1_down."\n".$id1_down_seq."\n".">".$id2_down."\n".$id2_down_seq_R."\n";
		}
	}
	close (IN);
	close (OUT);
	close (OUTFU);
	close (OUTFD);
	close (OUTRU);
	close (OUTRD);
}
print "SSR, upstream and downstream sequences are obtained\. \n"; #label

print "Beginning to identify the identical regions within upstream and downstream sequences, respectively \.\.\. \n"; #label
&runMafftProgram ("./temp/fasta/", "./temp/mafft/", "./temp/mafft.id");
&handleMafftResult ("./temp/mafft/", "./$workfile/identified_consistent_seq_after_mafft.tab");
&selectLongestSeq ("./$workfile/identified_consistent_seq_after_mafft.tab", "./$workfile/longest_consistent_seq_after_mafft.tab", "./$workfile/primer_prepare.id");
print "Finishing the identification of identical regions\. \n"; #label

&prepareFilePrimer ("./$workfile/longest_consistent_seq_after_mafft.tab", "./$workfile/primer_prepare.id", "./$workfile/primer_prepare.txt");
system "primer3_core ./$workfile/primer_prepare.txt > ./$workfile/primer_bank.txt" || die ("\nError in $0: Primer3 running is broken: $! !\n\n");
&handlePrimerResult ("./$workfile/primer_bank.txt", "./$workfile/primer_extract.txt");

mkdir "./temp/epcr/";
mkdir "./temp/epcr/output/";

system "mv ./$workfile/$prefix_fasta1.fasta ./temp/epcr/";
system "mv ./$workfile/$prefix_fasta2.fasta ./temp/epcr/";

system "famap -t n -b ./temp/epcr/$prefix_fasta1.famap ./temp/epcr/$prefix_fasta1.fasta" || die ("\nError in $0: ePCR running is broken: $! !\n\n");
system "famap -t n -b ./temp/epcr/$prefix_fasta2.famap ./temp/epcr/$prefix_fasta2.fasta" || die ("\nError in $0: ePCR running is broken: $! !\n\n");
system "fahash -b ./temp/epcr/$prefix_fasta1.hash -w 12 -f 3 ./temp/epcr/$prefix_fasta1.famap" || die ("\nError in $0: ePCR running is broken: $! !\n\n");
system "fahash -b ./temp/epcr/$prefix_fasta2.hash -w 12 -f 3 ./temp/epcr/$prefix_fasta2.famap" || die ("\nError in $0: ePCR running is broken: $! !\n\n");
&primerToEPCR ("./$workfile/primer_extract.txt", "$prefix_fasta1.fasta");
&primerToEPCR ("./$workfile/primer_extract.txt", "$prefix_fasta2.fasta");
&handleEPCRResult ("./$workfile/primer_eliminate.txt", "./$workfile/primer_check.txt", "./$workfile/primer_shortlong.temp", "./$workfile/primer_multimatch.temp");

if (-s "./$workfile/primer_check.txt") {
	$/ = "\n";
	open (IN, "./$workfile/primer_check.txt") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">./$workfile/primer_number.txt");
	my %result;
	while (<IN>) {
		chomp;
		my @line = split m/\t/, $_;
		chomp ($line[0]);
		if ($line[0] =~ m/(.+\d{6}_.+\d{6})_[FR]\.p\d/) {
			my $id = $1;
			if (not exists $result{$id}){
				$result{$id} = 1;
			} else {
				$result{$id}++;
			}	
		}
	}
	foreach (sort keys %result) {
		print OUT $_."\t".$result{$_}."\n";
	}
	close (IN);
	close (OUT);
	
	open (NO, "./$workfile/primer_number.txt");
	open (ORD, ">./$workfile/primer_order.txt");
	my $i = 0;
	while (<NO>) {
		chomp;
		my $order = sprintf "%05d", $i+1;
		my $ssr2marker_id = "S2M".$order;
		print ORD $_."\t".$ssr2marker_id."\n";
		$i++;
	}
	close (NO);
	close (ORD);
}

system "mv ./temp/epcr/$prefix_fasta1.fasta ./$workfile/";
system "mv ./temp/epcr/$prefix_fasta2.fasta ./$workfile/";

if (-s "./$workfile/primer_order.txt") {
	$/ = "\n";
	open (HA0, "./$workfile/primer_order.txt") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	my @ha0 = <HA0>;
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
	open (HA1, "./$workfile/merge.unit") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	my @ha1 = <HA1>;
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
	
	open (HA6, "./$workfile/primer_extract.txt") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	my @ha6 = <HA6>;
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
			if (($line6[1] =~ m/^[A|C|G|T]+$/) && ($line6[3] =~ m/^[A|C|G|T]+$/)) {
				$hash6{$line6[0]} = $line6[1];
				$hash7{$line6[0]} = $line6[2];
				$hash8{$line6[0]} = $line6[3];
				$hash9{$line6[0]} = $line6[4];
			}
		}
	}
	close (HA6);
	sub storeSSRClass;
	my $sub_hash= &storeSSRClass;
	$/ = ">";
	open (HAA, "./$workfile/$prefix_fasta1.fasta") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	my @haa = <HAA>;
	my %hasha; #species 1 or A
	foreach (@haa) {
		chomp;
		my @linea = split m/\n/, $_;
		chomp ($linea[0]);
		$linea[0] =~ s/\r//;
		chomp ($linea[1]);
		$linea[1] =~ s/\r//;
		$hasha{$linea[0]} = $linea[1];
	}
	close (HAA);
	open (HAB, "./$workfile/$prefix_fasta2.fasta") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	my @hab = <HAB>;
	my %hashb; #species 2 or B
	foreach (@hab) {
		chomp;
		my @lineb = split m/\n/, $_;
		chomp ($lineb[0]);
		$lineb[0] =~ s/\r//;
		chomp ($lineb[1]);
		$lineb[1] =~ s/\r//;
		$hashb{$lineb[0]} = $lineb[1];
	}
	close (HAB);
	
	$/ = "\n";
	open (IN, "./$workfile/primer_check.txt") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">./$workfile/primer_overall.txt"); # totally 36 rows

	print OUT "\#"."Primer_ID"."\t"."SSR2Marker_ID"."\t"."SSR_Class"."\t"."Primer_Pair_Count"."\t"."Primer_Pair_Number"."\t"."Primer_Forward"."\t"."Primer_Forward_Length"."\t"."Primer_Forward_Tm"."\t"."Primer_Reverse"."\t"."Primer_Reverse_Length"."\t"."Primer_Reverse_Tm"."\t"."Target_Sequence_Difference"."\t"."SSR_ID_".$prefix_fasta1."\t"."SSR_Unit_".$prefix_fasta1."\t"."SSR_Repeat_".$prefix_fasta1."\t"."SSR_Length_".$prefix_fasta1."\t"."SSR_Chromosome_".$prefix_fasta1."\t"."SSR_Chromosome_Direction_".$prefix_fasta1."\t"."SSR_Chromosome_Start_".$prefix_fasta1."\t"."SSR_Chromosome_End_".$prefix_fasta1."\t"."Target_Chromosome_Start_".$prefix_fasta1."\t"."Target_Chromosome_End_".$prefix_fasta1."\t"."Target_Seqence_".$prefix_fasta1."\t"."Target_Sequence_Length_".$prefix_fasta1."\t"."SSR_ID_".$prefix_fasta2."\t"."SSR_Unit_".$prefix_fasta2."\t"."SSR_Repeat_".$prefix_fasta2."\t"."SSR_Length_".$prefix_fasta2."\t"."SSR_Chromosome_".$prefix_fasta2."\t"."SSR_Chromosome_Direction_".$prefix_fasta2."\t"."SSR_Chromosome_Start_".$prefix_fasta2."\t"."SSR_Chromosome_End_".$prefix_fasta2."\t"."Target_Chromosome_Start_".$prefix_fasta2."\t"."Target_Chromosome_End_".$prefix_fasta2."\t"."Target_Seqence_".$prefix_fasta2."\t"."Target_Sequence_Length_".$prefix_fasta2."\n"; #title
	
	my ($primer_id, $ssr2marker_id, $ssr_class, $primer_pair_count, $primer_pair_number, $primer_forward, $primer_forward_length, $primer_forward_tm, $primer_reverse, $primer_reverse_length, $primer_reverse_tm, $target_sequence_difference, $ssr_id_species1, $ssr_unit_species1, $ssr_repeat_species1, $ssr_length_species1, $ssr_chromosome_species1, $ssr_chromosome_direction_species1, $ssr_chromosome_start_species1, $ssr_chromosome_end_species1, $target_chromosome_start_species1, $target_chromosome_end_species1, $target_seqence_species1, $target_sequence_length_species1, $ssr_id_species2, $ssr_unit_species2, $ssr_repeat_species2, $ssr_length_species2, $ssr_chromosome_species2, $ssr_chromosome_direction_species2, $ssr_chromosome_start_species2, $ssr_chromosome_end_species2, $target_chromosome_start_species2, $target_chromosome_end_species2, $target_seqence_species2, $target_sequence_length_species2);
	my $i = 0;
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
		
		if ($line[0] =~ m/(.+\d{6})_(.+\d{6})_([FR])\.p(\d)/) {
			
			my $id = $1."_".$2;
			
			$ssr2marker_id = $hash0id{$id}; # row 2
			$primer_pair_count = $hash0{$id}; # row 4
			
			$ssr_id_species1 = $2; # row 13
			$ssr_id_species2 = $1; # row 25
			
			$primer_pair_number = $4; # row 5
			
			my @list = split m/\#\#/, $line[1];
			chomp ($list[0]);
			chomp ($list[1]);
			
			my @a = split m/:/, $list[0];
			my @b = split m/:/, $list[1];
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
			
			if ($prefix_fasta1 eq $a[0]) {
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
				
				$ssr_class = $$sub_hash{$ssr_unit_species1}; # row 3, only the same selected, the same of species1 and species2, only assign once
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
		
		if (($primer_forward =~ m/^[A|C|G|T]/) && ($primer_reverse =~ m/^[A|C|G|T]/)) {
			print OUT $primer_id."\t".$ssr2marker_id."\t".$ssr_class."\t".$primer_pair_count."\t".$primer_pair_number."\t".$primer_forward."\t".$primer_forward_length."\t".$primer_forward_tm."\t".$primer_reverse."\t".$primer_reverse_length."\t".$primer_reverse_tm."\t".$target_sequence_difference."\t".$ssr_id_species1."\t".$ssr_unit_species1."\t".$ssr_repeat_species1."\t".$ssr_length_species1."\t".$ssr_chromosome_species1."\t".$ssr_chromosome_direction_species1."\t".$ssr_chromosome_start_species1."\t".$ssr_chromosome_end_species1."\t".$target_chromosome_start_species1."\t".$target_chromosome_end_species1."\t".$target_seqence_species1."\t".$target_sequence_length_species1."\t".$ssr_id_species2."\t".$ssr_unit_species2."\t".$ssr_repeat_species2."\t".$ssr_length_species2."\t".$ssr_chromosome_species2."\t".$ssr_chromosome_direction_species2."\t".$ssr_chromosome_start_species2."\t".$ssr_chromosome_end_species2."\t".$target_chromosome_start_species2."\t".$target_chromosome_end_species2."\t".$target_seqence_species2."\t".$target_sequence_length_species2."\n"; #result
		}
		
		$i++;
		if ($i > 999999) {
			die ("\nError in $0: The system broke down under excessive loads due to the huge SSR sequences !\nYou could split the input file and run this program again !\n\n");
		}
	}

	close (IN);
	close (OUT);
}

if (-s "./$workfile/primer_overall.txt") {
	my $lastline = `tail -n 1 ./$workfile/primer_overall.txt`;

	my @lastline = split m/\t/, $lastline;
	chomp ($lastline[0]);
	chomp ($lastline[1]);

	my $primer_count = $lastline[0];
	my $ssr2marker_count = $lastline[1];

	$primer_count =~ s/Primer0+//; #line 1 #title a
	$primer_count =~ s/Primer//; #line 1 #title a
	$ssr2marker_count =~ s/S2M0+//; #line 2 #title b
	$ssr2marker_count =~ s/S2M//; #line 2 #title b

	open (IN, "./$workfile/primer_overall.txt") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (IN1, "./$workfile/primer_multimatch.temp") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (IN2, "./$workfile/primer_shortlong.temp") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">./$workfile/SSR2Marker.stat");
	open (OUT1, ">./$workfile/SSR_Class.stat");
	open (OUT2, ">./$workfile/Target_Sequence_Difference.stat");
	open (OUT3, ">./$workfile/SSR_Length_$prefix_fasta1.stat");
	open (OUT4, ">./$workfile/SSR_Chromosome_$prefix_fasta1.stat");
	open (OUT5, ">./$workfile/Target_Sequence_Length_$prefix_fasta1.stat");
	open (OUT6, ">./$workfile/SSR_Length_$prefix_fasta2.stat");
	open (OUT7, ">./$workfile/SSR_Chromosome_$prefix_fasta2.stat");
	open (OUT8, ">./$workfile/Target_Sequence_Length_$prefix_fasta2.stat");

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

	print OUT1 "\#"."SSR classes shared by $base_fasta1 and $base_fasta2"."\t"."Count"."\n";
	print OUT2 "\#"."Length difference of target sequences between $base_fasta1 and $base_fasta2"."\t"."Count"."\n";
	print OUT3 "\#"."SSR length in $base_fasta1"."\t"."Count"."\n";
	print OUT4 "\#"."Chromosome distribution of SSRs in $base_fasta1"."\t"."Count"."\n";
	print OUT5 "\#"."Target sequence length in $base_fasta1"."\t"."Count"."\n";
	print OUT6 "\#"."SSR length in $base_fasta2"."\t"."Count"."\n";
	print OUT7 "\#"."Chromosome distribution of SSRs in $base_fasta2"."\t"."Count"."\n";
	print OUT8 "\#"."Target sequence length in $base_fasta2"."\t"."Count"."\n";

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
			my @line = split m/\t/, $_;
			
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
		
		foreach (my $i=0;$i<=$#keys_target_sequence_difference;$i++) {
			$count_target_sequence_difference += $keys_target_sequence_difference[$i];
		}
		
		print OUT "The total number of target sequences with different length between $base_fasta1 and $base_fasta2 is $count_target_sequence_difference\.\n"; #title d
	}

	# only min and max #
	foreach (sort {$a <=> $b} keys %ssr_length_species1) {
		print OUT3 $_."\t".$ssr_length_species1{$_}."\n";
	}

	my $min_ssr_length_species1 = (sort {$a <=> $b} keys %ssr_length_species1)[0]; #title e
	my $max_ssr_length_species1 = (sort {$b <=> $a} keys %ssr_length_species1)[0]; #title e

	my @keys_ssr_length_species1 = keys %ssr_length_species1;
	my $size_ssr_length_species1 = @keys_ssr_length_species1;

	# only count #
	my %rename1;
	open (R1, "./$workfile/$prefix_fasta1.rename") || die $!;
	while (<R1>) {
		chomp;
		if ($_) {
			my @line = split m/\t/, $_;
			chomp ($line[0]);
			chomp ($line[1]);
			$rename1{$line[0]} = $line[1];
		}
	}
	close (R1);
	my %rename2;
	open (R2, "./$workfile/$prefix_fasta2.rename") || die $!;
	while (<R2>) {
		chomp;
		if ($_) {
			my @line = split m/\t/, $_;
			chomp ($line[0]);
			chomp ($line[1]);
			$rename2{$line[0]} = $line[1];
		}
	}
	close (R2);
	
	foreach (sort keys %ssr_chromosome_species1) {
		print OUT4 $rename1{$_}."\t".$ssr_chromosome_species1{$_}."\n";
	}

	my @keys_ssr_chromosome_species1 = keys %ssr_chromosome_species1;
	my $size_ssr_chromosome_species1 = @keys_ssr_chromosome_species1; #title f

	# only min and max #
	foreach (sort {$a <=> $b} keys %target_sequence_length_species1) {
		print OUT5 $_."\t".$target_sequence_length_species1{$_}."\n";
	}

	my $min_target_sequence_length_species1 = (sort {$a <=> $b} keys %target_sequence_length_species1)[0]; #title g
	my $max_target_sequence_length_species1 = (sort {$b <=> $a} keys %target_sequence_length_species1)[0]; #title g

	my @keys_target_sequence_length_species1 = keys %target_sequence_length_species1;
	my $size_target_sequence_length_species1 = @keys_target_sequence_length_species1;

	# only min and max #
	foreach (sort {$a <=> $b} keys %ssr_length_species2) {
		print OUT6 $_."\t".$ssr_length_species2{$_}."\n";
	}

	my $min_ssr_length_species2 = (sort {$a <=> $b} keys %ssr_length_species2)[0]; #title h
	my $max_ssr_length_species2 = (sort {$b <=> $a} keys %ssr_length_species2)[0]; #title h

	my @keys_ssr_length_species2 = keys %ssr_length_species2;
	my $size_ssr_length_species2 = @keys_ssr_length_species2;

	# only count #
	foreach (sort keys %ssr_chromosome_species2) {
		print OUT7 $rename2{$_}."\t".$ssr_chromosome_species2{$_}."\n";
	}

	my @keys_ssr_chromosome_species2 = keys %ssr_chromosome_species2;
	my $size_ssr_chromosome_species2 = @keys_ssr_chromosome_species2; #title i

	# only min and max #
	foreach (sort {$a <=> $b} keys %target_sequence_length_species2) {
		print OUT8 $_."\t".$target_sequence_length_species2{$_}."\n";
	}

	my $min_target_sequence_length_species2 = (sort {$a <=> $b} keys %target_sequence_length_species2)[0]; #title j
	my $max_target_sequence_length_species2 = (sort {$b <=> $a} keys %target_sequence_length_species2)[0]; #title g

	my @keys_target_sequence_length_species2 = keys %target_sequence_length_species2;
	my $size_target_sequence_length_species2 = @keys_target_sequence_length_species2;

	print OUT "The length of SSRs in $base_fasta1 ranging from $min_ssr_length_species1 to $max_ssr_length_species1\.\n"; #title e
	print OUT "The total number of chromosomes or sequences with identified SSRs in $base_fasta1 is $size_ssr_chromosome_species1\.\n"; #title f
	print OUT "The length of target sequences in $base_fasta1 ranging from $min_target_sequence_length_species1 to $max_target_sequence_length_species1\.\n"; #title g
	print OUT "The length of SSRs in $base_fasta2 ranging from $min_ssr_length_species2 to $max_ssr_length_species2\.\n"; #title h
	print OUT "The total number of chromosomes with identified SSRs in $base_fasta2 is $size_ssr_chromosome_species2\.\n"; #title i
	print OUT "The length of target sequences in $base_fasta2 ranging from $min_target_sequence_length_species2 to $max_target_sequence_length_species2\.\n"; #title j

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
}

my $date_end = &getLocalTime();
my $workouttime = $date_end -> {date}; #obtain localtime yyyymmdd
my $workouttime_detail = $date_end -> {date_detail}; #obtain localtime yyyymmddhhmmss
my $running_second = &convertTimeSecond ($workouttime_detail) - &convertTimeSecond ($worktime_detail);
my $running_time = &convertSecondTime ($running_second);

my $resultfile;
if ($output_name) {
	$resultfile = $output_name;
} else {
	$resultfile = "Result_SSR2Marker_".$workouttime;
}

system "rm -rf ./$resultfile/";
system "mkdir ./$resultfile/";

system "cp ./$workfile/primer_overall.txt ./$resultfile/SSR2Marker.txt";
system "cp ./$workfile/*.stat ./$resultfile/";

open (RD, "./$resultfile/SSR2Marker.stat") || die ("\nError in $0: $! !\n\n");
while (<RD>) {
	if ($_ !~ /^\#/) {
		print $_;
	}
}
close (RD);

system "rm -rf ./$workfile/*.temp";
system "rm -rf ./*.temp";
system "rm -rf ./temp/";

print "\nAll identified molecular markers as well as their corresponding statistics have already been generated by SSR2Marker\. \n"; #label
print "You can view them in the folder of $resultfile now\. \n"; #label
print "The job takes about $running_time minutes\. \n"; #label


#========================================#
#          end main program !!!          #
#========================================#

#All subroutines used by the eMarker program

#Subroutine 1: Check the dependent softwares or programs
sub checkInstalledSoft { ##(softname["name"]<string>)
	$/ = "\n";
	
	`$_[0] -version > soft_check.log 2>&1`;
	
	open (LOG, "./soft_check.log");
	
	my @array;
	
	while (<LOG>) {
		chomp;
		if ($_) {
			push @array, $_;
		}
	}
	
	my $count = @array;
	my $softname = uc($_[0]);
	
	if (($count <= 2) && ("@array" =~ m/command not found/i)) {
		my $message = "\nError: The $softname software is not found. You should specify the path or install it firstly !\n\n";
		die $message;
	}
	else {
		print "The $softname software could be correctly called\. \n";
	}
	
	close (LOG);
	
	system "rm -rf ./soft_check.log";
}

#Subroutine 2: Get the local time
sub getLocalTime { ##()
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
		'date'   => "$year$mon$mday",
		'date_detail'   => "$year$mon$mday$hour$min$sec"
	};
}

#Subroutine 3: Convert the time to seconds for counting their difference
sub convertTimeSecond { ##(time[yyyymmddhhmmss]<int>)
	$_[0] =~ /(\d{4})\D?(\d{2})\D?(\d{2})(\d{2})\D?(\d{2})\D?(\d{2})/;
	my ($year, $mon, $day, $hour, $min, $sec) = ($1, $2, $3, $4, $5, $6);
	return mktime($sec, $min, $hour, $day, $mon - 1, $year - 1900, 0, 0);
}

#Subroutine 4: Convert the seconds to time for seeing the time
sub convertSecondTime { ##(time[number]<int>)
	my ($time, $day, $hour, $min, $sec, $day_remainder, $hour_remainder, $min_remainder);
	$day = int ($_[0]/86400);
	$day_remainder = ($_[0]%86400);
	$hour = int ($day_remainder/3600);
	$hour_remainder = ($day_remainder%3600);
	$min = int ($hour_remainder/60);
	$min_remainder = ($hour_remainder%60);
	$sec = $min_remainder;
	if ($day > 0) {
		$time = $day." day\(s\) ".$hour." hour\(s\) ".$min." minute\(s\) ".$sec." second\(s\)";
	} elsif ($hour > 0) {
		$time = $hour." hour\(s\) ".$min." minute\(s\) ".$sec." second\(s\)";
	} elsif ($min > 0) {
		$time = $min." minute\(s\) ".$sec." second\(s\)";
	} else {
		$time = $sec." second\(s\)";
	}
	return $time;
}

#Subroutine 5: Set the SSR motif
sub setSSRMotif { ##(motif[<STDIN>]<int>)
	print "Please type the minimum repeat number for each motif.\nmono-nucleotide (only positive integer number):";
	my $mono_type = <STDIN>;
	print "di-nucleotide (only positive integer number):";
	my $di_type = <STDIN>;
	print "tri-nucleotide (only positive integer number):";
	my $tri_type = <STDIN>;
	print "tetra-nucleotide (only positive integer number):";
	my $tetra_type = <STDIN>;
	print "penta-nucleotide (only positive integer number):";
	my $penta_type = <STDIN>;
	print "hexa-nucleotide (only positive integer number):";
	my $hexa_type = <STDIN>;
	my $motifs_default = '1=10,2=7,3=6,4=5,5=4,6=4';
	if ($mono_type =~ m/^\d+/) {
		$mono_type =~ s/(\r|\n)//;
	} else {
		die "Please type the correct parameter !\n;"
	}
	if ($di_type =~ m/^\d+/) {
		$di_type =~ s/(\r|\n)//;
	} else {
		die "Please type the correct parameter !\n;"
	}
	if ($tri_type =~ m/^\d+/) {
		$tri_type =~ s/(\r|\n)//;
	} else {
		die "Please type the correct parameter !\n;"
	}
	if ($tetra_type =~ m/^\d+/) {
		$tetra_type =~ s/(\r|\n)//;
	} else {
		die "Please type the correct parameter !\n;"
	}
	if ($penta_type =~ m/^\d+/) {
		$penta_type =~ s/(\r|\n)//;
	} else {
		die "Please type the correct parameter !\n;"
	}
	if ($hexa_type =~ m/^\d+/) {
		$hexa_type =~ s/(\r|\n)//;
	} else {
		die "Please type the correct parameter !\n;"
	}
	$motifs_default =~ s/1\=10/1\=$mono_type/;
	$motifs_default =~ s/2\=7/2\=$di_type/;
	$motifs_default =~ s/3\=6/3\=$tri_type/;
	$motifs_default =~ s/4\=5/4\=$tetra_type/;
	$motifs_default =~ s/5\=4/5\=$penta_type/;
	$motifs_default =~ s/6\=4/6\=$hexa_type/;
	return $motifs_default;
}

#Subroutine 6: Copy the fasta file and rename each sequence with a serial identifiers
sub copyFastaFile { ##(inputFile[>id\nseq]<path>, outputFile[>id\nseq]<path>, renameSeqFile[id\tid]<path>, floatingNumber[number]<int>)
	$/ = ">";
	
	my $base_thisfile = basename ($_[1]);	
	my $prefixthisfile = basename ($_[1], @suffix_list);	
	my $i = 0;
	
	open (IN, "$_[0]") || die ("\nError in $0: The file doesn't exist: $! !\n\n"); #input file
	open (OUT1, ">$_[1]"); #output file
	open (OUT2, ">>$_[2]"); #name convert file
	while (<IN>) {
		chomp;
		if ($_ =~ m/(.+?)\n(.+)/sx) {
			my $id = $1;
			my $seq = $2;
			$id =~ s/(\r|\n)//g;
			$seq =~ s/(\r|\n|\*)//g;
			my $sequence = uc ($seq);
			$sequence =~ s/U/T/g;
			$sequence =~ s/[BDEFHIJKLMOPQRSVWXYZ]/N/g;
			my $id_serial = sprintf "%0$_[3]d", $i+1; #serial number
			print OUT1 ">".$prefixthisfile.$id_serial."\n".$sequence."\n";
			print OUT2 $prefixthisfile.$id_serial."\t".$id."\n";
			
			$i++;
		}
		
		my $max_value = 9 x $_[3];
		
		if ($i >= $max_value) {
			die ("\nError in $0: The system broke down under excessive loads due to the huge input sequences !\nYou could split the input file and run this program again !\n\n");
		}
	}
	close (IN);
	close (OUT1);
	close (OUT2);
}

#Subroutine 7: Run MISA tools
sub runMISAProgram { ##(inputFile[>id\nseq]<path>, outputFile[\t]<path>, setSSRMotif[%]<string>)
	$/ = ">";
	
	open (IN, "$_[0]") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">$_[1]");
	my $motifs_info = $_[2];
	my @motif_digit = $motifs_info =~ /(\d+)=(\d+)/g;
	my %typrep = @motif_digit;
	my @typ = sort {$a <=> $b} keys %typrep;

	my $max_repeats = 1;
	my $min_repeats = 1000;
	my (%count_motif,%count_class);
	my ($number_sequences,$size_sequences,%ssr_containing_seqs);
	my $ssr_in_compound = 0;
	my ($id,$seq);
	
	while (<IN>) {
		next unless (($id,$seq) = /(.*?)\n(.*)/s);
		my ($nr,%start,@order,%end,%motif,%repeats);
		
		$seq =~ s/[\d\s>]//g;
		$id =~ s/^\s*//g;
		$id =~ s/\s*$//g;
		$id =~ s/\s/_/g;
		
		$number_sequences++;
		$size_sequences += length $seq;
		
		for (my $i=0;$i<scalar(@typ);$i++)	{
			my $motiflen = $typ[$i];
			my $minreps = $typrep{$typ[$i]} - 1;
			if ($min_repeats > $typrep{$typ[$i]}) {
				$min_repeats = $typrep{$typ[$i]}
			};
			my $search = "(([acgt]{$motiflen})\\2{$minreps,})";
			my %count_motifs;
			while ($seq =~ /$search/ig) {
				my $motif = uc $2;
				my $redundant;
				for (my $j = $motiflen-1;$j>0;$j--) {
					my $redmotif = "([ACGT]{$j})\\1{".($motiflen/$j-1)."}";
					$redundant = 1 if ( $motif =~ /$redmotif/ )
				};
				next if $redundant;
				$motif{++$nr} = $motif;
				my $ssr = uc $1;
				$repeats{$nr} = length($ssr) / $motiflen;
				$end{$nr} = pos($seq);
				$start{$nr} = $end{$nr} - length($ssr) + 1;
				$count_motifs{$motif{$nr}}++;
				$motif{$nr}->{$repeats{$nr}}++;
				$count_class{$typ[$i]}++;
				if ($max_repeats < $repeats{$nr}) {
					$max_repeats = $repeats{$nr}
				};
			};
		};
		next if (!$nr);
		$ssr_containing_seqs{$nr}++;
		@order = sort {$start{$a} <=> $start{$b}} keys %start;
		my $i = 0;
		my $count_seq;
		my ($start,$end,$ssrseq,$ssrtype,$size);
		while ($i < $nr) {
			my $space = 1;
			if (!$order[$i+1]) {
				$count_seq++;
				my $motiflen = length ($motif{$order[$i]});
				$ssrtype = "p".$motiflen;
				$ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}";
				$start = $start{$order[$i]};
				$end = $end{$order[$i++]};
				next
			};
			if (($start{$order[$i+1]} - $end{$order[$i]}) > $space) {
				$count_seq++;
				my $motiflen = length ($motif{$order[$i]});
				$ssrtype = "p".$motiflen;
				$ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}";
				$start = $start{$order[$i]}; $end = $end{$order[$i++]};
				next
			};
			my ($interssr);
			if (($start{$order[$i+1]} - $end{$order[$i]}) < 1) {
				$count_seq++; $ssr_in_compound++;
				$ssrtype = 'c*';
				$ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}($motif{$order[$i+1]})$repeats{$order[$i+1]}*";
				$start = $start{$order[$i]}; $end = $end{$order[$i+1]}
			} else {
				$count_seq++; $ssr_in_compound++;
				$interssr = lc substr($seq,$end{$order[$i]},($start{$order[$i+1]} - $end{$order[$i]}) - 1);
				$ssrtype = 'c';
				$ssrseq = "($motif{$order[$i]})$repeats{$order[$i]}$interssr($motif{$order[$i+1]})$repeats{$order[$i+1]}";
				$start = $start{$order[$i]};  $end = $end{$order[$i+1]};
			};
			while ($order[++$i + 1] and (($start{$order[$i+1]} - $end{$order[$i]}) <= $space)) {
				if (($start{$order[$i+1]} - $end{$order[$i]}) < 1) {
					$ssr_in_compound++;
					$ssrseq .= "($motif{$order[$i+1]})$repeats{$order[$i+1]}*";
					$ssrtype = 'c*';
					$end = $end{$order[$i+1]}
				} else {
					$ssr_in_compound++;
					$interssr = lc substr($seq,$end{$order[$i]},($start{$order[$i+1]} - $end{$order[$i]}) - 1);
					$ssrseq .= "$interssr($motif{$order[$i+1]})$repeats{$order[$i+1]}";
					$end = $end{$order[$i+1]};
				}
			};
			$i++;
		}
		continue {
			if ($ssrtype =~ m/c/) {
				print OUT "$id\t$count_seq\t$ssrtype\t\(N\)",($end - $start + 1),"\t",($end - $start + 1),"\t$start\t$end\t$ssrseq\n"
			} else {
				print OUT "$id\t$count_seq\t$ssrtype\t$ssrseq\t",($end - $start + 1),"\t$start\t$end\t$ssrseq\n"
			}
		};
	};
	close (IN);
	close (OUT);
}

#Subroutine 8: Handle the MISA result
sub handleMISAResult { ##(inputFastaFile[>id\nseq]<path>, flankingSeqLength[number]<int>, floatingNumber[number]<int>)
	$/ = ">";
	undef my %hash;
	my $base_thisfile = basename ($_[0]);
	my $prefix_thisfile = basename ($_[0], @suffix_list);
	open (FA, "./$workfile/$prefix_thisfile.fasta") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	my @fa = <FA>;
	foreach (@fa) {
		chomp;
		if ($_) {
			my @line = split m/\n/, $_;
			chomp ($line[0]);
			chomp ($line[1]);
			$hash{$line[0]} = $line[1];
		}
	}
	close (FA);
	
	$/ = "\n";
	open (MI, "./$workfile/$prefix_thisfile.misa") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (UP, ">./$workfile/$prefix_thisfile.up");
	open (FULL, ">./$workfile/$prefix_thisfile.full");
	open (DOWN, ">./$workfile/$prefix_thisfile.down");
	open (UNIT, ">./$workfile/$prefix_thisfile.unit");
	
	my $get_length = $_[1];
	my $max_value = 9 x $_[2];
	my @mi = <MI>;
	my $i = 0;
	foreach (@mi) {
		if ($_) {
			chomp;
			my $id = sprintf "%0$_[2]d", $i+1;
			my @list = split m/\t/, $_;
			chomp ($list[0]);
			chomp ($list[3]);
			chomp ($list[5]);
			chomp ($list[6]);
			chomp ($list[7]);
			if ($list[5] > $get_length) {
				print UP $prefix_thisfile.$id."U"."\t".$list[0]."\t".($list[5]-$get_length)."\t".($list[5]-1)."\t".$list[3]."\t".substr($hash{$list[0]},$list[5]-$get_length-1,$get_length)."\n";
				print FULL $prefix_thisfile.$id."\t".$list[0]."\t".$list[5]."\t".$list[6]."\t".$list[3]."\t".($get_length+1)."\t".($list[6]-$list[5]+$get_length+1)."\t".substr($hash{$list[0]},$list[5]-$get_length-1,$list[6]-$list[5]+$get_length*2+1)."\n";
				print DOWN $prefix_thisfile.$id."D"."\t".$list[0]."\t".($list[6]+1)."\t".($list[6]+$get_length)."\t".$list[3]."\t".substr($hash{$list[0]},$list[6],$get_length)."\n";
			} elsif ($list[5] == 1) {
				print UP $prefix_thisfile.$id."U"."\t".$list[0]."\t"."0"."\t"."0"."\t".$list[3]."\t"."N"."\n";
				print FULL $prefix_thisfile.$id."\t".$list[0]."\t".$list[5]."\t".$list[6]."\t".$list[3]."\t".$list[5]."\t".$list[6]."\t".substr($hash{$list[0]},0,$list[6]+$get_length)."\n";
				print DOWN $prefix_thisfile.$id."D"."\t".$list[0]."\t".($list[6]+1)."\t".($list[6]+$get_length+1)."\t".$list[3]."\t".substr($hash{$list[0]},$list[6],$get_length)."\n";
			} else {
				print UP $prefix_thisfile.$id."U"."\t".$list[0]."\t"."1"."\t".($list[5]-1)."\t".$list[3]."\t".substr($hash{$list[0]},0,$list[5]-1)."\n";
				print FULL $prefix_thisfile.$id."\t".$list[0]."\t".$list[5]."\t".$list[6]."\t".$list[3]."\t".$list[5]."\t".$list[6]."\t".substr($hash{$list[0]},0,$list[6]+$get_length)."\n";
				print DOWN $prefix_thisfile.$id."D"."\t".$list[0]."\t".($list[6]+1)."\t".($list[6]+$get_length)."\t".$list[3]."\t".substr($hash{$list[0]},$list[6],$get_length)."\n";
			}
			print UNIT $prefix_thisfile.$id."\t".$list[5]."\t".$list[6]."\t".$list[3]."\t".$list[7]."\n";
			$i++;
			if ($i >= $max_value) {
				die ("\nError in $0: The system broke down under excessive loads due to the huge sequences !\nYou could split the input file and run this program again !\n\n");
			}
		}
	}

	close (MI);
	close (UP);
	close (FULL);
	close (DOWN);
	close (UNIT);
}

#Subroutine 9: Run BLAST program
sub runBlastProgram { ##(referenceSeqFile[>id\nseq]<path>, querySeqFile[>id\nseq]<path>, outputTabFile[\t]<path>, evalue[number]<sprintf>, identity[number]<int>)
	system "makeblastdb -in $_[0] -parse_seqids -hash_index -dbtype nucl";
	system "blastn -db $_[0] -query $_[1] -out $_[2] -evalue $_[3] -outfmt 6 -perc_identity $_[4] -max_target_seqs 5 -num_threads $thread";
}

#subroutine 10: Count the ID in one row
sub countOneRow { ##(inputFile[id]<path>, outputFile[id\tseq]<path>)
	$/ = "\n";
	
	open (IN, "$_[0]") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">$_[1]");
	my %result;
	while (<IN>) {
		chomp;
		my $line = $_;
		if (not exists $result{$line}){
			$result{$line} = 1;
		}
		else {
			$result{$line}++;
		}
	}
	foreach (keys %result) {
		if ($result{$_} == 1){
			print OUT $_."\t"."1"."\n";
		}
		else {
			print OUT $_."\t".$result{$_}."\n";
		}
	}
	close (IN);
	close (OUT);
}

#subroutine 11: Store the SSR class to hash value
sub storeSSRClass { ##(%SSR_motif_class)
	my %SSR_motif_class = ("A","A","AA","A","AAA","A","AAAA","A","AAAAA","A","AAAAAA","A","T","A","TT","A","TTT","A","TTTT","A","TTTTT","A","TTTTTT","A","C","C","CC","C","CCC","C","CCCC","C","CCCCC","C","CCCCCC","C","G","C","GG","C","GGG","C","GGGG","C","GGGGG","C","GGGGGG","C","AT","AT","ATAT","AT","ATATAT","AT","TA","AT","TATA","AT","TATATA","AT","AC","AC","ACAC","AC","ACACAC","AC","CA","AC","CACA","AC","CACACA","AC","GT","AC","GTGT","AC","GTGTGT","AC","TG","AC","TGTG","AC","TGTGTG","AC","AG","AG","AGAG","AG","AGAGAG","AG","CT","AG","CTCT","AG","CTCTCT","AG","GA","AG","GAGA","AG","GAGAGA","AG","TC","AG","TCTC","AG","TCTCTC","AG","CG","CG","CGCG","CG","CGCGCG","CG","GC","CG","GCGC","CG","GCGCGC","CG","AAT","AAT","AATAAT","AAT","ATA","AAT","ATAATA","AAT","ATT","AAT","ATTATT","AAT","TAA","AAT","TAATAA","AAT","TAT","AAT","TATTAT","AAT","TTA","AAT","TTATTA","AAT","AAC","AAC","AACAAC","AAC","ACA","AAC","ACAACA","AAC","CAA","AAC","CAACAA","AAC","GTT","AAC","GTTGTT","AAC","TGT","AAC","TGTTGT","AAC","TTG","AAC","TTGTTG","AAC","AAG","AAG","AAGAAG","AAG","AGA","AAG","AGAAGA","AAG","CTT","AAG","CTTCTT","AAG","GAA","AAG","GAAGAA","AAG","TCT","AAG","TCTTCT","AAG","TTC","AAG","TTCTTC","AAG","AGT","AGT","AGTAGT","AGT","ATC","AGT","ATCATC","AGT","CAT","AGT","CATCAT","AGT","GTA","AGT","GTAGTA","AGT","TAG","AGT","TAGTAG","AGT","TCA","AGT","TCATCA","AGT","ACT","ACT","ACTACT","ACT","ATG","ACT","ATGATG","ACT","CTA","ACT","CTACTA","ACT","GAT","ACT","GATGAT","ACT","TAC","ACT","TACTAC","ACT","TGA","ACT","TGATGA","ACT","ACC","ACC","ACCACC","ACC","CAC","ACC","CACCAC","ACC","CCA","ACC","CCACCA","ACC","GGT","ACC","GGTGGT","ACC","GTG","ACC","GTGGTG","ACC","TGG","ACC","TGGTGG","ACC","ACG","ACG","ACGACG","ACG","CGA","ACG","CGACGA","ACG","CTG","ACG","CTGCTG","ACG","GAC","ACG","GACGAC","ACG","GCT","ACG","GCTGCT","ACG","TGC","ACG","TGCTGC","ACG","AGC","AGC","AGCAGC","AGC","CAG","AGC","CAGCAG","AGC","CGT","AGC","CGTCGT","AGC","GCA","AGC","GCAGCA","AGC","GTC","AGC","GTCGTC","AGC","TCG","AGC","TCGTCG","AGC","AGG","AGG","AGGAGG","AGG","CCT","AGG","CCTCCT","AGG","CTC","AGG","CTCCTC","AGG","GAG","AGG","GAGGAG","AGG","GGA","AGG","GGAGGA","AGG","TCC","AGG","TCCTCC","AGG","CCG","CCG","CCGCCG","CCG","CGC","CCG","CGCCGC","CCG","CGG","CCG","CGGCGG","CCG","GCC","CCG","GCCGCC","CCG","GCG","CCG","GCGGCG","CCG","GGC","CCG","GGCGGC","CCG","AAAT","AAAT","AATA","AAAT","ATAA","AAAT","ATTT","AAAT","TAAA","AAAT","TATT","AAAT","TTAT","AAAT","TTTA","AAAT","AAAC","AAAC","AACA","AAAC","ACAA","AAAC","CAAA","AAAC","GTTT","AAAC","TGTT","AAAC","TTGT","AAAC","TTTG","AAAC","AAAG","AAAG","AAGA","AAAG","AGAA","AAAG","CTTT","AAAG","GAAA","AAAG","TCTT","AAAG","TTCT","AAAG","TTTC","AAAG","AATT","AATT","ATTA","AATT","TAAT","AATT","TTAA","AATT","AATC","AATC","AGTT","AATC","ATCA","AATC","CAAT","AATC","GTTA","AATC","TAGT","AATC","TCAA","AATC","TTAG","AATC","AATG","AATG","ACTT","AATG","ATGA","AATG","CTTA","AATG","GAAT","AATG","TACT","AATG","TGAA","AATG","TTAC","AATG","AACT","AACT","ACTA","AACT","ATTG","AACT","CTAA","AACT","GATT","AACT","TAAC","AACT","TGAT","AACT","TTGA","AACT","AACC","AACC","ACCA","AACC","CAAC","AACC","CCAA","AACC","GGTT","AACC","GTTG","AACC","TGGT","AACC","TTGG","AACC","AACG","AACG","ACGA","AACG","CGAA","AACG","CTTG","AACG","GAAC","AACG","GCTT","AACG","TGCT","AACG","TTGC","AACG","AAGT","AAGT","AGTA","AAGT","ATTC","AAGT","CATT","AAGT","GTAA","AAGT","TAAG","AAGT","TCAT","AAGT","TTCA","AAGT","AAGC","AAGC","AGCA","AAGC","CAAG","AAGC","CGTT","AAGC","GCAA","AAGC","GTTC","AAGC","TCGT","AAGC","TTCG","AAGC","AAGG","AAGG","AGGA","AAGG","CCTT","AAGG","CTTC","AAGG","GAAG","AAGG","GGAA","AAGG","TCCT","AAGG","TTCC","AAGG","ACAT","ACAT","ATAC","ACAT","ATGT","ACAT","CATA","ACAT","GTAT","ACAT","TACA","ACAT","TATG","ACAT","TGTA","ACAT","AGAT","AGAT","ATAG","AGAT","ATCT","AGAT","CTAT","AGAT","GATA","AGAT","TAGA","AGAT","TATC","AGAT","TCTA","AGAT","AGGT","AGGT","ATCC","AGGT","CATC","AGGT","CCAT","AGGT","GGTA","AGGT","GTAG","AGGT","TAGG","AGGT","TCCA","AGGT","AGCT","AGCT","ATCG","AGCT","CGAT","AGCT","CTAG","AGCT","GATC","AGCT","GCTA","AGCT","TAGC","AGCT","TCGA","AGCT","ACGT","ACGT","ATGC","ACGT","CATG","ACGT","CGTA","ACGT","GCAT","ACGT","GTAC","ACGT","TACG","ACGT","TGCA","ACGT","ACCT","ACCT","ATGG","ACCT","CCTA","ACCT","CTAC","ACCT","GATG","ACCT","GGAT","ACCT","TACC","ACCT","TGGA","ACCT","ACAG","ACAG","AGAC","ACAG","CAGA","ACAG","CTGT","ACAG","GACA","ACAG","GTCT","ACAG","TCTG","ACAG","TGTC","ACAG","ACTC","ACTC","AGTG","ACTC","CACT","ACTC","CTCA","ACTC","GAGT","ACTC","GTGA","ACTC","TCAC","ACTC","TGAG","ACTC","ACTG","ACTG","CTGA","ACTG","GACT","ACTG","TGAC","ACTG","ACCC","ACCC","CACC","ACCC","CCAC","ACCC","CCCA","ACCC","GGGT","ACCC","GGTG","ACCC","GTGG","ACCC","TGGG","ACCC","ACCG","ACCG","CCGA","ACCG","CGAC","ACCG","CTGG","ACCG","GACC","ACCG","GCTG","ACCG","GGCT","ACCG","TGGC","ACCG","ACGC","ACGC","CACG","ACGC","CGCA","ACGC","CGTG","ACGC","GCAC","ACGC","GCGT","ACGC","GTGC","ACGC","TGCG","ACGC","ACGG","ACGG","CCTG","ACGG","CGGA","ACGG","CTGC","ACGG","GACG","ACGG","GCCT","ACGG","GGAC","ACGG","TGCC","ACGG","AGTC","AGTC","CAGT","AGTC","GTCA","AGTC","TCAG","AGTC","AGCC","AGCC","CAGC","AGCC","CCAG","AGCC","CGGT","AGCC","GCCA","AGCC","GGTC","AGCC","GTCG","AGCC","TCGG","AGCC","AGCG","AGCG","CGAG","AGCG","CGCT","AGCG","CTCG","AGCG","GAGC","AGCG","GCGA","AGCG","GCTC","AGCG","TCGC","AGCG","AGGC","AGGC","CAGG","AGGC","CCGT","AGGC","CGTC","AGGC","GCAG","AGGC","GGCA","AGGC","GTCC","AGGC","TCCG","AGGC","AGGG","AGGG","CCCT","AGGG","CCTC","AGGG","CTCC","AGGG","GAGG","AGGG","GGAG","AGGG","GGGA","AGGG","TCCC","AGGG","CCCG","CCCG","CCGC","CCCG","CGCC","CCCG","CGGG","CCCG","GCCC","CCCG","GCGG","CCCG","GGCG","CCCG","GGGC","CCCG","CCGG","CCGG","CGGC","CCGG","GCCG","CCGG","GGCC","CCGG","AAAAT","AAAAT","AAATA","AAAAT","AATAA","AAAAT","ATAAA","AAAAT","ATTTT","AAAAT","TAAAA","AAAAT","TATTT","AAAAT","TTATT","AAAAT","TTTAT","AAAAT","TTTTA","AAAAT","AAAAC","AAAAC","AAACA","AAAAC","AACAA","AAAAC","ACAAA","AAAAC","CAAAA","AAAAC","GTTTT","AAAAC","TGTTT","AAAAC","TTGTT","AAAAC","TTTGT","AAAAC","TTTTG","AAAAC","AAAAG","AAAAG","AAAGA","AAAAG","AAGAA","AAAAG","AGAAA","AAAAG","CTTTT","AAAAG","GAAAA","AAAAG","TCTTT","AAAAG","TTCTT","AAAAG","TTTCT","AAAAG","TTTTC","AAAAG","AAATT","AAATT","AATTA","AAATT","AATTT","AAATT","ATTAA","AAATT","ATTTA","AAATT","TAAAT","AAATT","TAATT","AAATT","TTAAA","AAATT","TTAAT","AAATT","TTTAA","AAATT","AAATC","AAATC","AATCA","AAATC","AGTTT","AAATC","ATCAA","AAATC","CAAAT","AAATC","GTTTA","AAATC","TAGTT","AAATC","TCAAA","AAATC","TTAGT","AAATC","TTTAG","AAATC","AAATG","AAATG","AATGA","AAATG","ACTTT","AAATG","ATGAA","AAATG","CTTTA","AAATG","GAAAT","AAATG","TACTT","AAATG","TGAAA","AAATG","TTACT","AAATG","TTTAC","AAATG","AAACT","AAACT","AACTA","AAACT","ACTAA","AAACT","ATTTG","AAACT","CTAAA","AAACT","GATTT","AAACT","TAAAC","AAACT","TGATT","AAACT","TTGAT","AAACT","TTTGA","AAACT","AAACC","AAACC","AACCA","AAACC","ACCAA","AAACC","CAAAC","AAACC","CCAAA","AAACC","GGTTT","AAACC","GTTTG","AAACC","TGGTT","AAACC","TTGGT","AAACC","TTTGG","AAACC","AAACG","AAACG","AACGA","AAACG","ACGAA","AAACG","CGAAA","AAACG","CTTTG","AAACG","GAAAC","AAACG","GCTTT","AAACG","TGCTT","AAACG","TTGCT","AAACG","TTTGC","AAACG","AAAGT","AAAGT","AAGTA","AAAGT","AGTAA","AAAGT","ATTTC","AAAGT","CATTT","AAAGT","GTAAA","AAAGT","TAAAG","AAAGT","TCATT","AAAGT","TTCAT","AAAGT","TTTCA","AAAGT","AAAGC","AAAGC","AAGCA","AAAGC","AGCAA","AAAGC","CAAAG","AAAGC","CGTTT","AAAGC","GCAAA","AAAGC","GTTTC","AAAGC","TCGTT","AAAGC","TTCGT","AAAGC","TTTCG","AAAGC","AAAGG","AAAGG","AAGGA","AAAGG","AGGAA","AAAGG","CCTTT","AAAGG","CTTTC","AAAGG","GAAAG","AAAGG","GGAAA","AAAGG","TCCTT","AAAGG","TTCCT","AAAGG","TTTCC","AAAGG","AATAT","AATAT","ATAAT","AATAT","ATATA","AATAT","ATATT","AATAT","ATTAT","AATAT","TAATA","AATAT","TATAA","AATAT","TATAT","AATAT","TATTA","AATAT","TTATA","AATAT","AATAC","AATAC","ACAAT","AATAC","ATACA","AATAC","ATGTT","AATAC","CAATA","AATAC","GTTAT","AATAC","TACAA","AATAC","TATGT","AATAC","TGTTA","AATAC","TTATG","AATAC","AATAG","AATAG","AGAAT","AATAG","ATAGA","AATAG","ATCTT","AATAG","CTTAT","AATAG","GAATA","AATAG","TAGAA","AATAG","TATCT","AATAG","TCTTA","AATAG","TTATC","AATAG","AAGTT","AAGTT","AATTC","AAGTT","AGTTA","AAGTT","ATTCA","AAGTT","CAATT","AAGTT","GTTAA","AAGTT","TAAGT","AAGTT","TCAAT","AAGTT","TTAAG","AAGTT","TTCAA","AAGTT","AACTT","AACTT","AATTG","AACTT","ACTTA","AACTT","ATTGA","AACTT","CTTAA","AACTT","GAATT","AACTT","TAACT","AACTT","TGAAT","AACTT","TTAAC","AACTT","TTGAA","AACTT","AATCT","AATCT","AGATT","AATCT","ATCTA","AATCT","ATTAG","AATCT","CTAAT","AATCT","GATTA","AATCT","TAATC","AATCT","TAGAT","AATCT","TCTAA","AATCT","TTAGA","AATCT","AATCC","AATCC","AGGTT","AATCC","ATCCA","AATCC","CAATC","AATCC","CCAAT","AATCC","GGTTA","AATCC","GTTAG","AATCC","TAGGT","AATCC","TCCAA","AATCC","TTAGG","AATCC","AATCG","AATCG","AGCTT","AATCG","ATCGA","AATCG","CGAAT","AATCG","CTTAG","AATCG","GAATC","AATCG","GCTTA","AATCG","TAGCT","AATCG","TCGAA","AATCG","TTAGC","AATCG","AATGT","AATGT","ACATT","AATGT","ATGTA","AATGT","ATTAC","AATGT","CATTA","AATGT","GTAAT","AATGT","TAATG","AATGT","TACAT","AATGT","TGTAA","AATGT","TTACA","AATGT","AATGC","AATGC","ACGTT","AATGC","ATGCA","AATGC","CAATG","AATGC","CGTTA","AATGC","GCAAT","AATGC","GTTAC","AATGC","TACGT","AATGC","TGCAA","AATGC","TTACG","AATGC","AATGG","AATGG","ACCTT","AATGG","ATGGA","AATGG","CCTTA","AATGG","CTTAC","AATGG","GAATG","AATGG","GGAAT","AATGG","TACCT","AATGG","TGGAA","AATGG","TTACC","AATGG","AACAT","AACAT","ACATA","AACAT","ATAAC","AACAT","ATTGT","AACAT","CATAA","AACAT","GTATT","AACAT","TAACA","AACAT","TATTG","AACAT","TGTAT","AACAT","TTGTA","AACAT","AACAC","AACAC","ACAAC","AACAC","ACACA","AACAC","CAACA","AACAC","CACAA","AACAC","GTGTT","AACAC","GTTGT","AACAC","TGTGT","AACAC","TGTTG","AACAC","TTGTG","AACAC","AACAG","AACAG","ACAGA","AACAG","AGAAC","AACAG","CAGAA","AACAG","CTTGT","AACAG","GAACA","AACAG","GTCTT","AACAG","TCTTG","AACAG","TGTCT","AACAG","TTGTC","AACAG","AACTC","AACTC","ACTCA","AACTC","AGTTG","AACTC","CAACT","AACTC","CTCAA","AACTC","GAGTT","AACTC","GTTGA","AACTC","TCAAC","AACTC","TGAGT","AACTC","TTGAG","AACTC","AACTG","AACTG","ACTGA","AACTG","ACTTG","AACTG","CTGAA","AACTG","CTTGA","AACTG","GAACT","AACTG","GACTT","AACTG","TGAAC","AACTG","TGACT","AACTG","TTGAC","AACTG","AACCT","AACCT","ACCTA","AACCT","ATTGG","AACCT","CCTAA","AACCT","CTAAC","AACCT","GATTG","AACCT","GGATT","AACCT","TAACC","AACCT","TGGAT","AACCT","TTGGA","AACCT","AACCC","AACCC","ACCCA","AACCC","CAACC","AACCC","CCAAC","AACCC","CCCAA","AACCC","GGGTT","AACCC","GGTTG","AACCC","GTTGG","AACCC","TGGGT","AACCC","TTGGG","AACCC","AACCG","AACCG","ACCGA","AACCG","CCGAA","AACCG","CGAAC","AACCG","CTTGG","AACCG","GAACC","AACCG","GCTTG","AACCG","GGCTT","AACCG","TGGCT","AACCG","TTGGC","AACCG","AACGT","AACGT","ACGTA","AACGT","ATTGC","AACGT","CATTG","AACGT","CGTAA","AACGT","GCATT","AACGT","GTAAC","AACGT","TAACG","AACGT","TGCAT","AACGT","TTGCA","AACGT","AACGC","AACGC","ACGCA","AACGC","CAACG","AACGC","CGCAA","AACGC","CGTTG","AACGC","GCAAC","AACGC","GCGTT","AACGC","GTTGC","AACGC","TGCGT","AACGC","TTGCG","AACGC","AACGG","AACGG","ACGGA","AACGG","CCTTG","AACGG","CGGAA","AACGG","CTTGC","AACGG","GAACG","AACGG","GCCTT","AACGG","GGAAC","AACGG","TGCCT","AACGG","TTGCC","AACGG","AAGAT","AAGAT","AGATA","AAGAT","ATAAG","AAGAT","ATTCT","AAGAT","CTATT","AAGAT","GATAA","AAGAT","TAAGA","AAGAT","TATTC","AAGAT","TCTAT","AAGAT","TTCTA","AAGAT","AAGAC","AAGAC","ACAAG","AAGAC","AGACA","AAGAC","CAAGA","AAGAC","CTGTT","AAGAC","GACAA","AAGAC","GTTCT","AAGAC","TCTGT","AAGAC","TGTTC","AAGAC","TTCTG","AAGAC","AAGAG","AAGAG","AGAAG","AAGAG","AGAGA","AAGAG","CTCTT","AAGAG","CTTCT","AAGAG","GAAGA","AAGAG","GAGAA","AAGAG","TCTCT","AAGAG","TCTTC","AAGAG","TTCTC","AAGAG","AAGTC","AAGTC","AGTCA","AAGTC","AGTTC","AAGTC","CAAGT","AAGTC","CAGTT","AAGTC","GTCAA","AAGTC","GTTCA","AAGTC","TCAAG","AAGTC","TCAGT","AAGTC","TTCAG","AAGTC","AAGTG","AAGTG","ACTTC","AAGTG","AGTGA","AAGTG","CACTT","AAGTG","CTTCA","AAGTG","GAAGT","AAGTG","GTGAA","AAGTG","TCACT","AAGTG","TGAAG","AAGTG","TTCAC","AAGTG","AAGCT","AAGCT","AGCTA","AAGCT","ATTCG","AAGCT","CGATT","AAGCT","CTAAG","AAGCT","GATTC","AAGCT","GCTAA","AAGCT","TAAGC","AAGCT","TCGAT","AAGCT","TTCGA","AAGCT","AAGCC","AAGCC","AGCCA","AAGCC","CAAGC","AAGCC","CCAAG","AAGCC","CGGTT","AAGCC","GCCAA","AAGCC","GGTTC","AAGCC","GTTCG","AAGCC","TCGGT","AAGCC","TTCGG","AAGCC","AAGCG","AAGCG","AGCGA","AAGCG","CGAAG","AAGCG","CGCTT","AAGCG","CTTCG","AAGCG","GAAGC","AAGCG","GCGAA","AAGCG","GCTTC","AAGCG","TCGCT","AAGCG","TTCGC","AAGCG","AAGGT","AAGGT","AGGTA","AAGGT","ATTCC","AAGGT","CATTC","AAGGT","CCATT","AAGGT","GGTAA","AAGGT","GTAAG","AAGGT","TAAGG","AAGGT","TCCAT","AAGGT","TTCCA","AAGGT","AAGGC","AAGGC","AGGCA","AAGGC","CAAGG","AAGGC","CCGTT","AAGGC","CGTTC","AAGGC","GCAAG","AAGGC","GGCAA","AAGGC","GTTCC","AAGGC","TCCGT","AAGGC","TTCCG","AAGGC","AAGGG","AAGGG","AGGGA","AAGGG","CCCTT","AAGGG","CCTTC","AAGGG","CTTCC","AAGGG","GAAGG","AAGGG","GGAAG","AAGGG","GGGAA","AAGGG","TCCCT","AAGGG","TTCCC","AAGGG","AGTAT","AGTAT","ATAGT","AGTAT","ATATC","AGTAT","ATCAT","AGTAT","CATAT","AGTAT","GTATA","AGTAT","TAGTA","AGTAT","TATAG","AGTAT","TATCA","AGTAT","TCATA","AGTAT","ACTAT","ACTAT","ATACT","ACTAT","ATATG","ACTAT","ATGAT","ACTAT","CTATA","ACTAT","GATAT","ACTAT","TACTA","ACTAT","TATAC","ACTAT","TATGA","ACTAT","TGATA","ACTAT","ACCAT","ACCAT","ATACC","ACCAT","ATGGT","ACCAT","CATAC","ACCAT","CCATA","ACCAT","GGTAT","ACCAT","GTATG","ACCAT","TACCA","ACCAT","TATGG","ACCAT","TGGTA","ACCAT","ACGAT","ACGAT","ATACG","ACGAT","ATGCT","ACGAT","CGATA","ACGAT","CTATG","ACGAT","GATAC","ACGAT","GCTAT","ACGAT","TACGA","ACGAT","TATGC","ACGAT","TGCTA","ACGAT","AGCAT","AGCAT","ATAGC","AGCAT","ATCGT","AGCAT","CATAG","AGCAT","CGTAT","AGCAT","GCATA","AGCAT","GTATC","AGCAT","TAGCA","AGCAT","TATCG","AGCAT","TCGTA","AGCAT","AGGAT","AGGAT","ATAGG","AGGAT","ATCCT","AGGAT","CCTAT","AGGAT","CTATC","AGGAT","GATAG","AGGAT","GGATA","AGGAT","TAGGA","AGGAT","TATCC","AGGAT","TCCTA","AGGAT","ACATC","ACATC","AGTGT","ACATC","ATCAC","ACATC","CACAT","ACATC","CATCA","ACATC","GTAGT","ACATC","GTGTA","ACATC","TAGTG","ACATC","TCACA","ACATC","TGTAG","ACATC","AGATC","AGATC","AGTCT","AGATC","ATCAG","AGATC","CAGAT","AGATC","CTAGT","AGATC","GATCA","AGATC","GTCTA","AGATC","TAGTC","AGATC","TCAGA","AGATC","TCTAG","AGATC","AGAGT","AGAGT","AGTAG","AGAGT","ATCTC","AGAGT","CATCT","AGAGT","CTCAT","AGAGT","GAGTA","AGAGT","GTAGA","AGAGT","TAGAG","AGAGT","TCATC","AGAGT","TCTCA","AGAGT","ACTAG","ACTAG","AGACT","ACTAG","ATCTG","ACTAG","CTAGA","ACTAG","CTGAT","ACTAG","GACTA","ACTAG","GATCT","ACTAG","TAGAC","ACTAG","TCTGA","ACTAG","TGATC","ACTAG","AGGGT","AGGGT","ATCCC","AGGGT","CATCC","AGGGT","CCATC","AGGGT","CCCAT","AGGGT","GGGTA","AGGGT","GGTAG","AGGGT","GTAGG","AGGGT","TAGGG","AGGGT","TCCCA","AGGGT","AGGCT","AGGCT","ATCCG","AGGCT","CCGAT","AGGCT","CGATC","AGGCT","CTAGG","AGGCT","GATCC","AGGCT","GCTAG","AGGCT","GGCTA","AGGCT","TAGGC","AGGCT","TCCGA","AGGCT","AGCGT","AGCGT","ATCGC","AGCGT","CATCG","AGCGT","CGCAT","AGCGT","CGTAG","AGCGT","GCATC","AGCGT","GCGTA","AGCGT","GTAGC","AGCGT","TAGCG","AGCGT","TCGCA","AGCGT","AGCCT","AGCCT","ATCGG","AGCCT","CCTAG","AGCCT","CGGAT","AGCCT","CTAGC","AGCCT","GATCG","AGCCT","GCCTA","AGCCT","GGATC","AGCCT","TAGCC","AGCCT","TCGGA","AGCCT","ACATG","ACATG","ACTGT","ACATG","ATGAC","ACATG","CATGA","ACATG","CTGTA","ACATG","GACAT","ACATG","GTACT","ACATG","TACTG","ACATG","TGACA","ACATG","TGTAC","ACATG","ACTCT","ACTCT","AGATG","ACTCT","ATGAG","ACTCT","CTACT","ACTCT","CTCTA","ACTCT","GAGAT","ACTCT","GATGA","ACTCT","TACTC","ACTCT","TCTAC","ACTCT","TGAGA","ACTCT","ACAGT","ACAGT","AGTAC","ACAGT","ATGTC","ACAGT","CAGTA","ACAGT","CATGT","ACAGT","GTACA","ACAGT","GTCAT","ACAGT","TACAG","ACAGT","TCATG","ACAGT","TGTCA","ACAGT","ACACT","ACACT","ACTAC","ACACT","ATGTG","ACACT","CACTA","ACACT","CTACA","ACACT","GATGT","ACACT","GTGAT","ACACT","TACAC","ACACT","TGATG","ACACT","TGTGA","ACACT","ACGGT","ACGGT","ATGCC","ACGGT","CATGC","ACGGT","CCATG","ACGGT","CGGTA","ACGGT","GCCAT","ACGGT","GGTAC","ACGGT","GTACG","ACGGT","TACGG","ACGGT","TGCCA","ACGGT","ACGCT","ACGCT","ATGCG","ACGCT","CGATG","ACGCT","CGCTA","ACGCT","CTACG","ACGCT","GATGC","ACGCT","GCGAT","ACGCT","GCTAC","ACGCT","TACGC","ACGCT","TGCGA","ACGCT","ACCGT","ACCGT","ATGGC","ACCGT","CATGG","ACCGT","CCGTA","ACCGT","CGTAC","ACCGT","GCATG","ACCGT","GGCAT","ACCGT","GTACC","ACCGT","TACCG","ACCGT","TGGCA","ACCGT","ACCCT","ACCCT","ATGGG","ACCCT","CCCTA","ACCCT","CCTAC","ACCCT","CTACC","ACCCT","GATGG","ACCCT","GGATG","ACCCT","GGGAT","ACCCT","TACCC","ACCCT","TGGGA","ACCCT","ACACC","ACACC","ACCAC","ACACC","CACAC","ACACC","CACCA","ACACC","CCACA","ACACC","GGTGT","ACACC","GTGGT","ACACC","GTGTG","ACACC","TGGTG","ACACC","TGTGG","ACACC","ACACG","ACACG","ACGAC","ACACG","CACGA","ACACG","CGACA","ACACG","CTGTG","ACACG","GACAC","ACACG","GCTGT","ACACG","GTGCT","ACACG","TGCTG","ACACG","TGTGC","ACACG","ACAGC","ACAGC","AGCAC","ACAGC","CACAG","ACAGC","CAGCA","ACAGC","CGTGT","ACAGC","GCACA","ACAGC","GTCGT","ACAGC","GTGTC","ACAGC","TCGTG","ACAGC","TGTCG","ACAGC","ACAGG","ACAGG","AGGAC","ACAGG","CAGGA","ACAGG","CCTGT","ACAGG","CTGTC","ACAGG","GACAG","ACAGG","GGACA","ACAGG","GTCCT","ACAGG","TCCTG","ACAGG","TGTCC","ACAGG","ACTCC","ACTCC","AGGTG","ACTCC","CACTC","ACTCC","CCACT","ACTCC","CTCCA","ACTCC","GAGGT","ACTCC","GGTGA","ACTCC","GTGAG","ACTCC","TCCAC","ACTCC","TGAGG","ACTCC","ACTCG","ACTCG","AGCTG","ACTCG","CGACT","ACTCG","CTCGA","ACTCG","CTGAG","ACTCG","GACTC","ACTCG","GAGCT","ACTCG","GCTGA","ACTCG","TCGAC","ACTCG","TGAGC","ACTCG","ACGTG","ACGTG","ACTGC","ACGTG","CACTG","ACGTG","CGTGA","ACGTG","CTGCA","ACGTG","GACGT","ACGTG","GCACT","ACGTG","GTGAC","ACGTG","TGACG","ACGTG","TGCAC","ACGTG","ACCTG","ACCTG","ACTGG","ACCTG","CCTGA","ACCTG","CTGAC","ACCTG","CTGGA","ACCTG","GACCT","ACCTG","GACTG","ACCTG","GGACT","ACCTG","TGACC","ACCTG","TGGAC","ACCTG","ACCAG","ACCAG","AGACC","ACCAG","CAGAC","ACCAG","CCAGA","ACCAG","CTGGT","ACCAG","GACCA","ACCAG","GGTCT","ACCAG","GTCTG","ACCAG","TCTGG","ACCAG","TGGTC","ACCAG","ACCTC","ACCTC","AGTGG","ACCTC","CACCT","ACCTC","CCTCA","ACCTC","CTCAC","ACCTC","GAGTG","ACCTC","GGAGT","ACCTC","GTGGA","ACCTC","TCACC","ACCTC","TGGAG","ACCTC","ACCCC","ACCCC","CACCC","ACCCC","CCACC","ACCCC","CCCAC","ACCCC","CCCCA","ACCCC","GGGGT","ACCCC","GGGTG","ACCCC","GGTGG","ACCCC","GTGGG","ACCCC","TGGGG","ACCCC","ACCCG","ACCCG","CCCGA","ACCCG","CCGAC","ACCCG","CGACC","ACCCG","CTGGG","ACCCG","GACCC","ACCCG","GCTGG","ACCCG","GGCTG","ACCCG","GGGCT","ACCCG","TGGGC","ACCCG","ACCGC","ACCGC","CACCG","ACCGC","CCGCA","ACCGC","CGCAC","ACCGC","CGTGG","ACCGC","GCACC","ACCGC","GCGTG","ACCGC","GGCGT","ACCGC","GTGGC","ACCGC","TGGCG","ACCGC","ACCGG","ACCGG","CCGGA","ACCGG","CCTGG","ACCGG","CGGAC","ACCGG","CTGGC","ACCGG","GACCG","ACCGG","GCCTG","ACCGG","GGACC","ACCGG","GGCCT","ACCGG","TGGCC","ACCGG","ACGAG","ACGAG","AGACG","ACGAG","CGAGA","ACGAG","CTCTG","ACGAG","CTGCT","ACGAG","GACGA","ACGAG","GAGAC","ACGAG","GCTCT","ACGAG","TCTGC","ACGAG","TGCTC","ACGAG","ACGTC","ACGTC","AGTGC","ACGTC","CACGT","ACGTC","CAGTG","ACGTC","CGTCA","ACGTC","GCAGT","ACGTC","GTCAC","ACGTC","GTGCA","ACGTC","TCACG","ACGTC","TGCAG","ACGTC","ACGCC","ACGCC","CACGC","ACGCC","CCACG","ACGCC","CGCCA","ACGCC","CGGTG","ACGCC","GCCAC","ACGCC","GCGGT","ACGCC","GGTGC","ACGCC","GTGCG","ACGCC","TGCGG","ACGCC","ACGCG","ACGCG","CGACG","ACGCG","CGCGA","ACGCG","CGCTG","ACGCG","CTGCG","ACGCG","GACGC","ACGCG","GCGAC","ACGCG","GCGCT","ACGCG","GCTGC","ACGCG","TGCGC","ACGCG","ACGGC","ACGGC","CACGG","ACGGC","CCGTG","ACGGC","CGGCA","ACGGC","CGTGC","ACGGC","GCACG","ACGGC","GCCGT","ACGGC","GGCAC","ACGGC","GTGCC","ACGGC","TGCCG","ACGGC","ACGGG","ACGGG","CCCTG","ACGGG","CCTGC","ACGGG","CGGGA","ACGGG","CTGCC","ACGGG","GACGG","ACGGG","GCCCT","ACGGG","GGACG","ACGGG","GGGAC","ACGGG","TGCCC","ACGGG","AGAGC","AGAGC","AGCAG","AGAGC","CAGAG","AGAGC","CGTCT","AGAGC","CTCGT","AGAGC","GAGCA","AGAGC","GCAGA","AGAGC","GTCTC","AGAGC","TCGTC","AGAGC","TCTCG","AGAGC","AGAGG","AGAGG","AGGAG","AGAGG","CCTCT","AGAGG","CTCCT","AGAGG","CTCTC","AGAGG","GAGAG","AGAGG","GAGGA","AGAGG","GGAGA","AGAGG","TCCTC","AGAGG","TCTCC","AGAGG","AGGTC","AGGTC","AGTCC","AGGTC","CAGGT","AGGTC","CAGTC","AGGTC","CCAGT","AGGTC","GGTCA","AGGTC","GTCAG","AGGTC","GTCCA","AGGTC","TCAGG","AGGTC","TCCAG","AGGTC","AGCTC","AGCTC","AGTCG","AGCTC","CAGCT","AGCTC","CGAGT","AGCTC","CTCAG","AGCTC","GAGTC","AGCTC","GCTCA","AGCTC","GTCGA","AGCTC","TCAGC","AGCTC","TCGAG","AGCTC","AGCCC","AGCCC","CAGCC","AGCCC","CCAGC","AGCCC","CCCAG","AGCCC","CGGGT","AGCCC","GCCCA","AGCCC","GGGTC","AGCCC","GGTCG","AGCCC","GTCGG","AGCCC","TCGGG","AGCCC","AGCCG","AGCCG","CCGAG","AGCCG","CGAGC","AGCCG","CGGCT","AGCCG","CTCGG","AGCCG","GAGCC","AGCCG","GCCGA","AGCCG","GCTCG","AGCCG","GGCTC","AGCCG","TCGGC","AGCCG","AGCGC","AGCGC","CAGCG","AGCGC","CGCAG","AGCGC","CGCGT","AGCGC","CGTCG","AGCGC","GCAGC","AGCGC","GCGCA","AGCGC","GCGTC","AGCGC","GTCGC","AGCGC","TCGCG","AGCGC","AGCGG","AGCGG","CCTCG","AGCGG","CGCCT","AGCGG","CGGAG","AGCGG","CTCGC","AGCGG","GAGCG","AGCGG","GCCTC","AGCGG","GCGGA","AGCGG","GGAGC","AGCGG","TCGCC","AGCGG","AGGCC","AGGCC","CAGGC","AGGCC","CCAGG","AGGCC","CCGGT","AGGCC","CGGTC","AGGCC","GCCAG","AGGCC","GGCCA","AGGCC","GGTCC","AGGCC","GTCCG","AGGCC","TCCGG","AGGCC","AGGCG","AGGCG","CCGCT","AGGCG","CGAGG","AGGCG","CGCTC","AGGCG","CTCCG","AGGCG","GAGGC","AGGCG","GCGAG","AGGCG","GCTCC","AGGCG","GGCGA","AGGCG","TCCGC","AGGCG","AGGGC","AGGGC","CAGGG","AGGGC","CCCGT","AGGGC","CCGTC","AGGGC","CGTCC","AGGGC","GCAGG","AGGGC","GGCAG","AGGGC","GGGCA","AGGGC","GTCCC","AGGGC","TCCCG","AGGGC","AGGGG","AGGGG","CCCCT","AGGGG","CCCTC","AGGGG","CCTCC","AGGGG","CTCCC","AGGGG","GAGGG","AGGGG","GGAGG","AGGGG","GGGAG","AGGGG","GGGGA","AGGGG","TCCCC","AGGGG","CCCCG","CCCCG","CCCGC","CCCCG","CCGCC","CCCCG","CGCCC","CCCCG","CGGGG","CCCCG","GCCCC","CCCCG","GCGGG","CCCCG","GGCGG","CCCCG","GGGCG","CCCCG","GGGGC","CCCCG","CCCGG","CCCGG","CCGGC","CCCGG","CCGGG","CCCGG","CGGCC","CCCGG","CGGGC","CCCGG","GCCCG","CCCGG","GCCGG","CCCGG","GGCCC","CCCGG","GGCCG","CCCGG","GGGCC","CCCGG","CCGCG","CCGCG","CGCCG","CCGCG","CGCGC","CCGCG","CGCGG","CCGCG","CGGCG","CCGCG","GCCGC","CCGCG","GCGCC","CCGCG","GCGCG","CCGCG","GCGGC","CCGCG","GGCGC","CCGCG","AAAAAT","AAAAAT","AAAATA","AAAAAT","AAATAA","AAAAAT","AATAAA","AAAAAT","ATAAAA","AAAAAT","ATTTTT","AAAAAT","TAAAAA","AAAAAT","TATTTT","AAAAAT","TTATTT","AAAAAT","TTTATT","AAAAAT","TTTTAT","AAAAAT","TTTTTA","AAAAAT","AAAAAC","AAAAAC","AAAACA","AAAAAC","AAACAA","AAAAAC","AACAAA","AAAAAC","ACAAAA","AAAAAC","CAAAAA","AAAAAC","GTTTTT","AAAAAC","TGTTTT","AAAAAC","TTGTTT","AAAAAC","TTTGTT","AAAAAC","TTTTGT","AAAAAC","TTTTTG","AAAAAC","AAAAAG","AAAAAG","AAAAGA","AAAAAG","AAAGAA","AAAAAG","AAGAAA","AAAAAG","AGAAAA","AAAAAG","CTTTTT","AAAAAG","GAAAAA","AAAAAG","TCTTTT","AAAAAG","TTCTTT","AAAAAG","TTTCTT","AAAAAG","TTTTCT","AAAAAG","TTTTTC","AAAAAG","AAAATT","AAAATT","AAATTA","AAAATT","AATTAA","AAAATT","AATTTT","AAAATT","ATTAAA","AAAATT","ATTTTA","AAAATT","TAAAAT","AAAATT","TAATTT","AAAATT","TTAAAA","AAAATT","TTAATT","AAAATT","TTTAAT","AAAATT","TTTTAA","AAAATT","AAAATC","AAAATC","AAATCA","AAAATC","AATCAA","AAAATC","AGTTTT","AAAATC","ATCAAA","AAAATC","CAAAAT","AAAATC","GTTTTA","AAAATC","TAGTTT","AAAATC","TCAAAA","AAAATC","TTAGTT","AAAATC","TTTAGT","AAAATC","TTTTAG","AAAATC","AAAATG","AAAATG","AAATGA","AAAATG","AATGAA","AAAATG","ACTTTT","AAAATG","ATGAAA","AAAATG","CTTTTA","AAAATG","GAAAAT","AAAATG","TACTTT","AAAATG","TGAAAA","AAAATG","TTACTT","AAAATG","TTTACT","AAAATG","TTTTAC","AAAATG","AAAACT","AAAACT","AAACTA","AAAACT","AACTAA","AAAACT","ACTAAA","AAAACT","ATTTTG","AAAACT","CTAAAA","AAAACT","GATTTT","AAAACT","TAAAAC","AAAACT","TGATTT","AAAACT","TTGATT","AAAACT","TTTGAT","AAAACT","TTTTGA","AAAACT","AAAACC","AAAACC","AAACCA","AAAACC","AACCAA","AAAACC","ACCAAA","AAAACC","CAAAAC","AAAACC","CCAAAA","AAAACC","GGTTTT","AAAACC","GTTTTG","AAAACC","TGGTTT","AAAACC","TTGGTT","AAAACC","TTTGGT","AAAACC","TTTTGG","AAAACC","AAAACG","AAAACG","AAACGA","AAAACG","AACGAA","AAAACG","ACGAAA","AAAACG","CGAAAA","AAAACG","CTTTTG","AAAACG","GAAAAC","AAAACG","GCTTTT","AAAACG","TGCTTT","AAAACG","TTGCTT","AAAACG","TTTGCT","AAAACG","TTTTGC","AAAACG","AAAAGT","AAAAGT","AAAGTA","AAAAGT","AAGTAA","AAAAGT","AGTAAA","AAAAGT","ATTTTC","AAAAGT","CATTTT","AAAAGT","GTAAAA","AAAAGT","TAAAAG","AAAAGT","TCATTT","AAAAGT","TTCATT","AAAAGT","TTTCAT","AAAAGT","TTTTCA","AAAAGT","AAAAGC","AAAAGC","AAAGCA","AAAAGC","AAGCAA","AAAAGC","AGCAAA","AAAAGC","CAAAAG","AAAAGC","CGTTTT","AAAAGC","GCAAAA","AAAAGC","GTTTTC","AAAAGC","TCGTTT","AAAAGC","TTCGTT","AAAAGC","TTTCGT","AAAAGC","TTTTCG","AAAAGC","AAAAGG","AAAAGG","AAAGGA","AAAAGG","AAGGAA","AAAAGG","AGGAAA","AAAAGG","CCTTTT","AAAAGG","CTTTTC","AAAAGG","GAAAAG","AAAAGG","GGAAAA","AAAAGG","TCCTTT","AAAAGG","TTCCTT","AAAAGG","TTTCCT","AAAAGG","TTTTCC","AAAAGG","AAATAT","AAATAT","AATATA","AAATAT","ATAAAT","AAATAT","ATATAA","AAATAT","ATATTT","AAATAT","ATTTAT","AAATAT","TAAATA","AAATAT","TATAAA","AAATAT","TATATT","AAATAT","TATTTA","AAATAT","TTATAT","AAATAT","TTTATA","AAATAT","AAATAC","AAATAC","AATACA","AAATAC","ACAAAT","AAATAC","ATACAA","AAATAC","ATGTTT","AAATAC","CAAATA","AAATAC","GTTTAT","AAATAC","TACAAA","AAATAC","TATGTT","AAATAC","TGTTTA","AAATAC","TTATGT","AAATAC","TTTATG","AAATAC","AAATAG","AAATAG","AATAGA","AAATAG","AGAAAT","AAATAG","ATAGAA","AAATAG","ATCTTT","AAATAG","CTTTAT","AAATAG","GAAATA","AAATAG","TAGAAA","AAATAG","TATCTT","AAATAG","TCTTTA","AAATAG","TTATCT","AAATAG","TTTATC","AAATAG","AAATTT","AAATTT","AATTTA","AAATTT","ATTTAA","AAATTT","TAAATT","AAATTT","TTAAAT","AAATTT","TTTAAA","AAATTT","AAATTC","AAATTC","AAGTTT","AAATTC","AATTCA","AAATTC","AGTTTA","AAATTC","ATTCAA","AAATTC","CAAATT","AAATTC","GTTTAA","AAATTC","TAAGTT","AAATTC","TCAAAT","AAATTC","TTAAGT","AAATTC","TTCAAA","AAATTC","TTTAAG","AAATTC","AAATTG","AAATTG","AACTTT","AAATTG","AATTGA","AAATTG","ACTTTA","AAATTG","ATTGAA","AAATTG","CTTTAA","AAATTG","GAAATT","AAATTG","TAACTT","AAATTG","TGAAAT","AAATTG","TTAACT","AAATTG","TTGAAA","AAATTG","TTTAAC","AAATTG","AAATCT","AAATCT","AATCTA","AAATCT","AGATTT","AAATCT","ATCTAA","AAATCT","ATTTAG","AAATCT","CTAAAT","AAATCT","GATTTA","AAATCT","TAAATC","AAATCT","TAGATT","AAATCT","TCTAAA","AAATCT","TTAGAT","AAATCT","TTTAGA","AAATCT","AAATCC","AAATCC","AATCCA","AAATCC","AGGTTT","AAATCC","ATCCAA","AAATCC","CAAATC","AAATCC","CCAAAT","AAATCC","GGTTTA","AAATCC","GTTTAG","AAATCC","TAGGTT","AAATCC","TCCAAA","AAATCC","TTAGGT","AAATCC","TTTAGG","AAATCC","AAATCG","AAATCG","AATCGA","AAATCG","AGCTTT","AAATCG","ATCGAA","AAATCG","CGAAAT","AAATCG","CTTTAG","AAATCG","GAAATC","AAATCG","GCTTTA","AAATCG","TAGCTT","AAATCG","TCGAAA","AAATCG","TTAGCT","AAATCG","TTTAGC","AAATCG","AAATGT","AAATGT","AATGTA","AAATGT","ACATTT","AAATGT","ATGTAA","AAATGT","ATTTAC","AAATGT","CATTTA","AAATGT","GTAAAT","AAATGT","TAAATG","AAATGT","TACATT","AAATGT","TGTAAA","AAATGT","TTACAT","AAATGT","TTTACA","AAATGT","AAATGC","AAATGC","AATGCA","AAATGC","ACGTTT","AAATGC","ATGCAA","AAATGC","CAAATG","AAATGC","CGTTTA","AAATGC","GCAAAT","AAATGC","GTTTAC","AAATGC","TACGTT","AAATGC","TGCAAA","AAATGC","TTACGT","AAATGC","TTTACG","AAATGC","AAATGG","AAATGG","AATGGA","AAATGG","ACCTTT","AAATGG","ATGGAA","AAATGG","CCTTTA","AAATGG","CTTTAC","AAATGG","GAAATG","AAATGG","GGAAAT","AAATGG","TACCTT","AAATGG","TGGAAA","AAATGG","TTACCT","AAATGG","TTTACC","AAATGG","AAACAT","AAACAT","AACATA","AAACAT","ACATAA","AAACAT","ATAAAC","AAACAT","ATTTGT","AAACAT","CATAAA","AAACAT","GTATTT","AAACAT","TAAACA","AAACAT","TATTTG","AAACAT","TGTATT","AAACAT","TTGTAT","AAACAT","TTTGTA","AAACAT","AAACAC","AAACAC","AACACA","AAACAC","ACAAAC","AAACAC","ACACAA","AAACAC","CAAACA","AAACAC","CACAAA","AAACAC","GTGTTT","AAACAC","GTTTGT","AAACAC","TGTGTT","AAACAC","TGTTTG","AAACAC","TTGTGT","AAACAC","TTTGTG","AAACAC","AAACAG","AAACAG","AACAGA","AAACAG","ACAGAA","AAACAG","AGAAAC","AAACAG","CAGAAA","AAACAG","CTTTGT","AAACAG","GAAACA","AAACAG","GTCTTT","AAACAG","TCTTTG","AAACAG","TGTCTT","AAACAG","TTGTCT","AAACAG","TTTGTC","AAACAG","AAACTT","AAACTT","AACTTA","AAACTT","AATTTG","AAACTT","ACTTAA","AAACTT","ATTTGA","AAACTT","CTTAAA","AAACTT","GAATTT","AAACTT","TAAACT","AAACTT","TGAATT","AAACTT","TTAAAC","AAACTT","TTGAAT","AAACTT","TTTGAA","AAACTT","AAACTC","AAACTC","AACTCA","AAACTC","ACTCAA","AAACTC","AGTTTG","AAACTC","CAAACT","AAACTC","CTCAAA","AAACTC","GAGTTT","AAACTC","GTTTGA","AAACTC","TCAAAC","AAACTC","TGAGTT","AAACTC","TTGAGT","AAACTC","TTTGAG","AAACTC","AAACTG","AAACTG","AACTGA","AAACTG","ACTGAA","AAACTG","ACTTTG","AAACTG","CTGAAA","AAACTG","CTTTGA","AAACTG","GAAACT","AAACTG","GACTTT","AAACTG","TGAAAC","AAACTG","TGACTT","AAACTG","TTGACT","AAACTG","TTTGAC","AAACTG","AAACCT","AAACCT","AACCTA","AAACCT","ACCTAA","AAACCT","ATTTGG","AAACCT","CCTAAA","AAACCT","CTAAAC","AAACCT","GATTTG","AAACCT","GGATTT","AAACCT","TAAACC","AAACCT","TGGATT","AAACCT","TTGGAT","AAACCT","TTTGGA","AAACCT","AAACCC","AAACCC","AACCCA","AAACCC","ACCCAA","AAACCC","CAAACC","AAACCC","CCAAAC","AAACCC","CCCAAA","AAACCC","GGGTTT","AAACCC","GGTTTG","AAACCC","GTTTGG","AAACCC","TGGGTT","AAACCC","TTGGGT","AAACCC","TTTGGG","AAACCC","AAACCG","AAACCG","AACCGA","AAACCG","ACCGAA","AAACCG","CCGAAA","AAACCG","CGAAAC","AAACCG","CTTTGG","AAACCG","GAAACC","AAACCG","GCTTTG","AAACCG","GGCTTT","AAACCG","TGGCTT","AAACCG","TTGGCT","AAACCG","TTTGGC","AAACCG","AAACGT","AAACGT","AACGTA","AAACGT","ACGTAA","AAACGT","ATTTGC","AAACGT","CATTTG","AAACGT","CGTAAA","AAACGT","GCATTT","AAACGT","GTAAAC","AAACGT","TAAACG","AAACGT","TGCATT","AAACGT","TTGCAT","AAACGT","TTTGCA","AAACGT","AAACGC","AAACGC","AACGCA","AAACGC","ACGCAA","AAACGC","CAAACG","AAACGC","CGCAAA","AAACGC","CGTTTG","AAACGC","GCAAAC","AAACGC","GCGTTT","AAACGC","GTTTGC","AAACGC","TGCGTT","AAACGC","TTGCGT","AAACGC","TTTGCG","AAACGC","AAACGG","AAACGG","AACGGA","AAACGG","ACGGAA","AAACGG","CCTTTG","AAACGG","CGGAAA","AAACGG","CTTTGC","AAACGG","GAAACG","AAACGG","GCCTTT","AAACGG","GGAAAC","AAACGG","TGCCTT","AAACGG","TTGCCT","AAACGG","TTTGCC","AAACGG","AAAGAT","AAAGAT","AAGATA","AAAGAT","AGATAA","AAAGAT","ATAAAG","AAAGAT","ATTTCT","AAAGAT","CTATTT","AAAGAT","GATAAA","AAAGAT","TAAAGA","AAAGAT","TATTTC","AAAGAT","TCTATT","AAAGAT","TTCTAT","AAAGAT","TTTCTA","AAAGAT","AAAGAC","AAAGAC","AAGACA","AAAGAC","ACAAAG","AAAGAC","AGACAA","AAAGAC","CAAAGA","AAAGAC","CTGTTT","AAAGAC","GACAAA","AAAGAC","GTTTCT","AAAGAC","TCTGTT","AAAGAC","TGTTTC","AAAGAC","TTCTGT","AAAGAC","TTTCTG","AAAGAC","AAAGAG","AAAGAG","AAGAGA","AAAGAG","AGAAAG","AAAGAG","AGAGAA","AAAGAG","CTCTTT","AAAGAG","CTTTCT","AAAGAG","GAAAGA","AAAGAG","GAGAAA","AAAGAG","TCTCTT","AAAGAG","TCTTTC","AAAGAG","TTCTCT","AAAGAG","TTTCTC","AAAGAG","AAAGTT","AAAGTT","AAGTTA","AAAGTT","AATTTC","AAAGTT","AGTTAA","AAAGTT","ATTTCA","AAAGTT","CAATTT","AAAGTT","GTTAAA","AAAGTT","TAAAGT","AAAGTT","TCAATT","AAAGTT","TTAAAG","AAAGTT","TTCAAT","AAAGTT","TTTCAA","AAAGTT","AAAGTC","AAAGTC","AAGTCA","AAAGTC","AGTCAA","AAAGTC","AGTTTC","AAAGTC","CAAAGT","AAAGTC","CAGTTT","AAAGTC","GTCAAA","AAAGTC","GTTTCA","AAAGTC","TCAAAG","AAAGTC","TCAGTT","AAAGTC","TTCAGT","AAAGTC","TTTCAG","AAAGTC","AAAGTG","AAAGTG","AAGTGA","AAAGTG","ACTTTC","AAAGTG","AGTGAA","AAAGTG","CACTTT","AAAGTG","CTTTCA","AAAGTG","GAAAGT","AAAGTG","GTGAAA","AAAGTG","TCACTT","AAAGTG","TGAAAG","AAAGTG","TTCACT","AAAGTG","TTTCAC","AAAGTG","AAAGCT","AAAGCT","AAGCTA","AAAGCT","AGCTAA","AAAGCT","ATTTCG","AAAGCT","CGATTT","AAAGCT","CTAAAG","AAAGCT","GATTTC","AAAGCT","GCTAAA","AAAGCT","TAAAGC","AAAGCT","TCGATT","AAAGCT","TTCGAT","AAAGCT","TTTCGA","AAAGCT","AAAGCC","AAAGCC","AAGCCA","AAAGCC","AGCCAA","AAAGCC","CAAAGC","AAAGCC","CCAAAG","AAAGCC","CGGTTT","AAAGCC","GCCAAA","AAAGCC","GGTTTC","AAAGCC","GTTTCG","AAAGCC","TCGGTT","AAAGCC","TTCGGT","AAAGCC","TTTCGG","AAAGCC","AAAGCG","AAAGCG","AAGCGA","AAAGCG","AGCGAA","AAAGCG","CGAAAG","AAAGCG","CGCTTT","AAAGCG","CTTTCG","AAAGCG","GAAAGC","AAAGCG","GCGAAA","AAAGCG","GCTTTC","AAAGCG","TCGCTT","AAAGCG","TTCGCT","AAAGCG","TTTCGC","AAAGCG","AAAGGT","AAAGGT","AAGGTA","AAAGGT","AGGTAA","AAAGGT","ATTTCC","AAAGGT","CATTTC","AAAGGT","CCATTT","AAAGGT","GGTAAA","AAAGGT","GTAAAG","AAAGGT","TAAAGG","AAAGGT","TCCATT","AAAGGT","TTCCAT","AAAGGT","TTTCCA","AAAGGT","AAAGGC","AAAGGC","AAGGCA","AAAGGC","AGGCAA","AAAGGC","CAAAGG","AAAGGC","CCGTTT","AAAGGC","CGTTTC","AAAGGC","GCAAAG","AAAGGC","GGCAAA","AAAGGC","GTTTCC","AAAGGC","TCCGTT","AAAGGC","TTCCGT","AAAGGC","TTTCCG","AAAGGC","AAAGGG","AAAGGG","AAGGGA","AAAGGG","AGGGAA","AAAGGG","CCCTTT","AAAGGG","CCTTTC","AAAGGG","CTTTCC","AAAGGG","GAAAGG","AAAGGG","GGAAAG","AAAGGG","GGGAAA","AAAGGG","TCCCTT","AAAGGG","TTCCCT","AAAGGG","TTTCCC","AAAGGG","AACAAT","AACAAT","AATAAC","AACAAT","ACAATA","AACAAT","ATAACA","AACAAT","ATTGTT","AACAAT","CAATAA","AACAAT","GTTATT","AACAAT","TAACAA","AACAAT","TATTGT","AACAAT","TGTTAT","AACAAT","TTATTG","AACAAT","TTGTTA","AACAAT","AAGAAT","AAGAAT","AATAAG","AAGAAT","AGAATA","AAGAAT","ATAAGA","AAGAAT","ATTCTT","AAGAAT","CTTATT","AAGAAT","GAATAA","AAGAAT","TAAGAA","AAGAAT","TATTCT","AAGAAT","TCTTAT","AAGAAT","TTATTC","AAGAAT","TTCTTA","AAGAAT","AATATT","AATATT","AATTAT","AATATT","ATAATT","AATATT","ATATTA","AATATT","ATTAAT","AATATT","ATTATA","AATATT","TAATAT","AATATT","TAATTA","AATATT","TATAAT","AATATT","TATTAA","AATATT","TTAATA","AATATT","TTATAA","AATATT","AATATC","AATATC","AGTTAT","AATATC","ATAGTT","AATATC","ATATCA","AATATC","ATCAAT","AATATC","CAATAT","AATATC","GTTATA","AATATC","TAGTTA","AATATC","TATAGT","AATATC","TATCAA","AATATC","TCAATA","AATATC","TTATAG","AATATC","AATATG","AATATG","ACTTAT","AATATG","ATACTT","AATATG","ATATGA","AATATG","ATGAAT","AATATG","CTTATA","AATATG","GAATAT","AATATG","TACTTA","AATATG","TATACT","AATATG","TATGAA","AATATG","TGAATA","AATATG","TTATAC","AATATG","AATACT","AATACT","ACTAAT","AATACT","ATACTA","AATACT","ATGATT","AATACT","ATTATG","AATACT","CTAATA","AATACT","GATTAT","AATACT","TAATAC","AATACT","TACTAA","AATACT","TATGAT","AATACT","TGATTA","AATACT","TTATGA","AATACT","AATACC","AATACC","ACCAAT","AATACC","ATACCA","AATACC","ATGGTT","AATACC","CAATAC","AATACC","CCAATA","AATACC","GGTTAT","AATACC","GTTATG","AATACC","TACCAA","AATACC","TATGGT","AATACC","TGGTTA","AATACC","TTATGG","AATACC","AATACG","AATACG","ACGAAT","AATACG","ATACGA","AATACG","ATGCTT","AATACG","CGAATA","AATACG","CTTATG","AATACG","GAATAC","AATACG","GCTTAT","AATACG","TACGAA","AATACG","TATGCT","AATACG","TGCTTA","AATACG","TTATGC","AATACG","AATAGT","AATAGT","AGTAAT","AATAGT","ATAGTA","AATAGT","ATCATT","AATAGT","ATTATC","AATAGT","CATTAT","AATAGT","GTAATA","AATAGT","TAATAG","AATAGT","TAGTAA","AATAGT","TATCAT","AATAGT","TCATTA","AATAGT","TTATCA","AATAGT","AATAGC","AATAGC","AGCAAT","AATAGC","ATAGCA","AATAGC","ATCGTT","AATAGC","CAATAG","AATAGC","CGTTAT","AATAGC","GCAATA","AATAGC","GTTATC","AATAGC","TAGCAA","AATAGC","TATCGT","AATAGC","TCGTTA","AATAGC","TTATCG","AATAGC","AATAGG","AATAGG","AGGAAT","AATAGG","ATAGGA","AATAGG","ATCCTT","AATAGG","CCTTAT","AATAGG","CTTATC","AATAGG","GAATAG","AATAGG","GGAATA","AATAGG","TAGGAA","AATAGG","TATCCT","AATAGG","TCCTTA","AATAGG","TTATCC","AATAGG","AATGTT","AATGTT","AATTAC","AATGTT","ACAATT","AATGTT","ATGTTA","AATGTT","ATTACA","AATGTT","CAATTA","AATGTT","GTTAAT","AATGTT","TAATGT","AATGTT","TACAAT","AATGTT","TGTTAA","AATGTT","TTAATG","AATGTT","TTACAA","AATGTT","AATCTT","AATCTT","AATTAG","AATCTT","AGAATT","AATCTT","ATCTTA","AATCTT","ATTAGA","AATCTT","CTTAAT","AATCTT","GAATTA","AATCTT","TAATCT","AATCTT","TAGAAT","AATCTT","TCTTAA","AATCTT","TTAATC","AATCTT","TTAGAA","AATCTT","AAGATT","AAGATT","AATTCT","AAGATT","AGATTA","AAGATT","ATTAAG","AAGATT","ATTCTA","AAGATT","CTAATT","AAGATT","GATTAA","AAGATT","TAAGAT","AAGATT","TAATTC","AAGATT","TCTAAT","AAGATT","TTAAGA","AAGATT","TTCTAA","AAGATT","AAGGTT","AAGGTT","AATTCC","AAGGTT","AGGTTA","AAGGTT","ATTCCA","AAGGTT","CAATTC","AAGGTT","CCAATT","AAGGTT","GGTTAA","AAGGTT","GTTAAG","AAGGTT","TAAGGT","AAGGTT","TCCAAT","AAGGTT","TTAAGG","AAGGTT","TTCCAA","AAGGTT","AAGCTT","AAGCTT","AATTCG","AAGCTT","AGCTTA","AAGCTT","ATTCGA","AAGCTT","CGAATT","AAGCTT","CTTAAG","AAGCTT","GAATTC","AAGCTT","GCTTAA","AAGCTT","TAAGCT","AAGCTT","TCGAAT","AAGCTT","TTAAGC","AAGCTT","TTCGAA","AAGCTT","AACATT","AACATT","AATTGT","AACATT","ACATTA","AACATT","ATTAAC","AACATT","ATTGTA","AACATT","CATTAA","AACATT","GTAATT","AACATT","TAACAT","AACATT","TAATTG","AACATT","TGTAAT","AACATT","TTAACA","AACATT","TTGTAA","AACATT","AACGTT","AACGTT","AATTGC","AACGTT","ACGTTA","AACGTT","ATTGCA","AACGTT","CAATTG","AACGTT","CGTTAA","AACGTT","GCAATT","AACGTT","GTTAAC","AACGTT","TAACGT","AACGTT","TGCAAT","AACGTT","TTAACG","AACGTT","TTGCAA","AACGTT","AACCTT","AACCTT","AATTGG","AACCTT","ACCTTA","AACCTT","ATTGGA","AACCTT","CCTTAA","AACCTT","CTTAAC","AACCTT","GAATTG","AACCTT","GGAATT","AACCTT","TAACCT","AACCTT","TGGAAT","AACCTT","TTAACC","AACCTT","TTGGAA","AACCTT","AATCAT","AATCAT","AGTATT","AATCAT","ATAATC","AATCAT","ATCATA","AATCAT","ATTAGT","AATCAT","CATAAT","AATCAT","GTATTA","AATCAT","TAATCA","AATCAT","TAGTAT","AATCAT","TATTAG","AATCAT","TCATAA","AATCAT","TTAGTA","AATCAT","AATCAC","AATCAC","ACAATC","AATCAC","AGTGTT","AATCAC","ATCACA","AATCAC","CAATCA","AATCAC","CACAAT","AATCAC","GTGTTA","AATCAC","GTTAGT","AATCAC","TAGTGT","AATCAC","TCACAA","AATCAC","TGTTAG","AATCAC","TTAGTG","AATCAC","AATCAG","AATCAG","AGAATC","AATCAG","AGTCTT","AATCAG","ATCAGA","AATCAG","CAGAAT","AATCAG","CTTAGT","AATCAG","GAATCA","AATCAG","GTCTTA","AATCAG","TAGTCT","AATCAG","TCAGAA","AATCAG","TCTTAG","AATCAG","TTAGTC","AATCAG","AATCTC","AATCTC","AGAGTT","AATCTC","AGTTAG","AATCTC","ATCTCA","AATCTC","CAATCT","AATCTC","CTCAAT","AATCTC","GAGTTA","AATCTC","GTTAGA","AATCTC","TAGAGT","AATCTC","TCAATC","AATCTC","TCTCAA","AATCTC","TTAGAG","AATCTC","AATCTG","AATCTG","ACTTAG","AATCTG","AGACTT","AATCTG","ATCTGA","AATCTG","CTGAAT","AATCTG","CTTAGA","AATCTG","GAATCT","AATCTG","GACTTA","AATCTG","TAGACT","AATCTG","TCTGAA","AATCTG","TGAATC","AATCTG","TTAGAC","AATCTG","AATCCT","AATCCT","AGGATT","AATCCT","ATCCTA","AATCCT","ATTAGG","AATCCT","CCTAAT","AATCCT","CTAATC","AATCCT","GATTAG","AATCCT","GGATTA","AATCCT","TAATCC","AATCCT","TAGGAT","AATCCT","TCCTAA","AATCCT","TTAGGA","AATCCT","AATCCC","AATCCC","AGGGTT","AATCCC","ATCCCA","AATCCC","CAATCC","AATCCC","CCAATC","AATCCC","CCCAAT","AATCCC","GGGTTA","AATCCC","GGTTAG","AATCCC","GTTAGG","AATCCC","TAGGGT","AATCCC","TCCCAA","AATCCC","TTAGGG","AATCCC","AATCCG","AATCCG","AGGCTT","AATCCG","ATCCGA","AATCCG","CCGAAT","AATCCG","CGAATC","AATCCG","CTTAGG","AATCCG","GAATCC","AATCCG","GCTTAG","AATCCG","GGCTTA","AATCCG","TAGGCT","AATCCG","TCCGAA","AATCCG","TTAGGC","AATCCG","AATCGT","AATCGT","AGCATT","AATCGT","ATCGTA","AATCGT","ATTAGC","AATCGT","CATTAG","AATCGT","CGTAAT","AATCGT","GCATTA","AATCGT","GTAATC","AATCGT","TAATCG","AATCGT","TAGCAT","AATCGT","TCGTAA","AATCGT","TTAGCA","AATCGT","AATCGC","AATCGC","AGCGTT","AATCGC","ATCGCA","AATCGC","CAATCG","AATCGC","CGCAAT","AATCGC","CGTTAG","AATCGC","GCAATC","AATCGC","GCGTTA","AATCGC","GTTAGC","AATCGC","TAGCGT","AATCGC","TCGCAA","AATCGC","TTAGCG","AATCGC","AATCGG","AATCGG","AGCCTT","AATCGG","ATCGGA","AATCGG","CCTTAG","AATCGG","CGGAAT","AATCGG","CTTAGC","AATCGG","GAATCG","AATCGG","GCCTTA","AATCGG","GGAATC","AATCGG","TAGCCT","AATCGG","TCGGAA","AATCGG","TTAGCC","AATCGG","AATGAT","AATGAT","ACTATT","AATGAT","ATAATG","AATGAT","ATGATA","AATGAT","ATTACT","AATGAT","CTATTA","AATGAT","GATAAT","AATGAT","TAATGA","AATGAT","TACTAT","AATGAT","TATTAC","AATGAT","TGATAA","AATGAT","TTACTA","AATGAT","AATGAC","AATGAC","ACAATG","AATGAC","ACTGTT","AATGAC","ATGACA","AATGAC","CAATGA","AATGAC","CTGTTA","AATGAC","GACAAT","AATGAC","GTTACT","AATGAC","TACTGT","AATGAC","TGACAA","AATGAC","TGTTAC","AATGAC","TTACTG","AATGAC","AATGAG","AATGAG","ACTCTT","AATGAG","AGAATG","AATGAG","ATGAGA","AATGAG","CTCTTA","AATGAG","CTTACT","AATGAG","GAATGA","AATGAG","GAGAAT","AATGAG","TACTCT","AATGAG","TCTTAC","AATGAG","TGAGAA","AATGAG","TTACTC","AATGAG","AATGTC","AATGTC","ACAGTT","AATGTC","AGTTAC","AATGTC","ATGTCA","AATGTC","CAATGT","AATGTC","CAGTTA","AATGTC","GTCAAT","AATGTC","GTTACA","AATGTC","TACAGT","AATGTC","TCAATG","AATGTC","TGTCAA","AATGTC","TTACAG","AATGTC","AATGTG","AATGTG","ACACTT","AATGTG","ACTTAC","AATGTG","ATGTGA","AATGTG","CACTTA","AATGTG","CTTACA","AATGTG","GAATGT","AATGTG","GTGAAT","AATGTG","TACACT","AATGTG","TGAATG","AATGTG","TGTGAA","AATGTG","TTACAC","AATGTG","AATGCT","AATGCT","ACGATT","AATGCT","ATGCTA","AATGCT","ATTACG","AATGCT","CGATTA","AATGCT","CTAATG","AATGCT","GATTAC","AATGCT","GCTAAT","AATGCT","TAATGC","AATGCT","TACGAT","AATGCT","TGCTAA","AATGCT","TTACGA","AATGCT","AATGCC","AATGCC","ACGGTT","AATGCC","ATGCCA","AATGCC","CAATGC","AATGCC","CCAATG","AATGCC","CGGTTA","AATGCC","GCCAAT","AATGCC","GGTTAC","AATGCC","GTTACG","AATGCC","TACGGT","AATGCC","TGCCAA","AATGCC","TTACGG","AATGCC","AATGCG","AATGCG","ACGCTT","AATGCG","ATGCGA","AATGCG","CGAATG","AATGCG","CGCTTA","AATGCG","CTTACG","AATGCG","GAATGC","AATGCG","GCGAAT","AATGCG","GCTTAC","AATGCG","TACGCT","AATGCG","TGCGAA","AATGCG","TTACGC","AATGCG","AATGGT","AATGGT","ACCATT","AATGGT","ATGGTA","AATGGT","ATTACC","AATGGT","CATTAC","AATGGT","CCATTA","AATGGT","GGTAAT","AATGGT","GTAATG","AATGGT","TAATGG","AATGGT","TACCAT","AATGGT","TGGTAA","AATGGT","TTACCA","AATGGT","AATGGC","AATGGC","ACCGTT","AATGGC","ATGGCA","AATGGC","CAATGG","AATGGC","CCGTTA","AATGGC","CGTTAC","AATGGC","GCAATG","AATGGC","GGCAAT","AATGGC","GTTACC","AATGGC","TACCGT","AATGGC","TGGCAA","AATGGC","TTACCG","AATGGC","AATGGG","AATGGG","ACCCTT","AATGGG","ATGGGA","AATGGG","CCCTTA","AATGGG","CCTTAC","AATGGG","CTTACC","AATGGG","GAATGG","AATGGG","GGAATG","AATGGG","GGGAAT","AATGGG","TACCCT","AATGGG","TGGGAA","AATGGG","TTACCC","AATGGG","AACAAG","AACAAG","AAGAAC","AACAAG","ACAAGA","AACAAG","AGAACA","AACAAG","CAAGAA","AACAAG","CTTGTT","AACAAG","GAACAA","AACAAG","GTTCTT","AACAAG","TCTTGT","AACAAG","TGTTCT","AACAAG","TTCTTG","AACAAG","TTGTTC","AACAAG","AACATC","AACATC","ACATCA","AACATC","AGTTGT","AACATC","ATCAAC","AACATC","CAACAT","AACATC","CATCAA","AACATC","GTAGTT","AACATC","GTTGTA","AACATC","TAGTTG","AACATC","TCAACA","AACATC","TGTAGT","AACATC","TTGTAG","AACATC","AACATG","AACATG","ACATGA","AACATG","ACTTGT","AACATG","ATGAAC","AACATG","CATGAA","AACATG","CTTGTA","AACATG","GAACAT","AACATG","GTACTT","AACATG","TACTTG","AACATG","TGAACA","AACATG","TGTACT","AACATG","TTGTAC","AACATG","AACACT","AACACT","ACACTA","AACACT","ACTAAC","AACACT","ATTGTG","AACACT","CACTAA","AACACT","CTAACA","AACACT","GATTGT","AACACT","GTGATT","AACACT","TAACAC","AACACT","TGATTG","AACACT","TGTGAT","AACACT","TTGTGA","AACACT","AACACC","AACACC","ACACCA","AACACC","ACCAAC","AACACC","CAACAC","AACACC","CACCAA","AACACC","CCAACA","AACACC","GGTTGT","AACACC","GTGGTT","AACACC","GTTGTG","AACACC","TGGTTG","AACACC","TGTGGT","AACACC","TTGTGG","AACACC","AACACG","AACACG","ACACGA","AACACG","ACGAAC","AACACG","CACGAA","AACACG","CGAACA","AACACG","CTTGTG","AACACG","GAACAC","AACACG","GCTTGT","AACACG","GTGCTT","AACACG","TGCTTG","AACACG","TGTGCT","AACACG","TTGTGC","AACACG","AACAGT","AACAGT","ACAGTA","AACAGT","AGTAAC","AACAGT","ATTGTC","AACAGT","CAGTAA","AACAGT","CATTGT","AACAGT","GTAACA","AACAGT","GTCATT","AACAGT","TAACAG","AACAGT","TCATTG","AACAGT","TGTCAT","AACAGT","TTGTCA","AACAGT","AACAGC","AACAGC","ACAGCA","AACAGC","AGCAAC","AACAGC","CAACAG","AACAGC","CAGCAA","AACAGC","CGTTGT","AACAGC","GCAACA","AACAGC","GTCGTT","AACAGC","GTTGTC","AACAGC","TCGTTG","AACAGC","TGTCGT","AACAGC","TTGTCG","AACAGC","AACAGG","AACAGG","ACAGGA","AACAGG","AGGAAC","AACAGG","CAGGAA","AACAGG","CCTTGT","AACAGG","CTTGTC","AACAGG","GAACAG","AACAGG","GGAACA","AACAGG","GTCCTT","AACAGG","TCCTTG","AACAGG","TGTCCT","AACAGG","TTGTCC","AACAGG","AACTAT","AACTAT","ACTATA","AACTAT","ATAACT","AACTAT","ATATTG","AACTAT","ATTGAT","AACTAT","CTATAA","AACTAT","GATATT","AACTAT","TAACTA","AACTAT","TATAAC","AACTAT","TATTGA","AACTAT","TGATAT","AACTAT","TTGATA","AACTAT","AACTAC","AACTAC","ACAACT","AACTAC","ACTACA","AACTAC","ATGTTG","AACTAC","CAACTA","AACTAC","CTACAA","AACTAC","GATGTT","AACTAC","GTTGAT","AACTAC","TACAAC","AACTAC","TGATGT","AACTAC","TGTTGA","AACTAC","TTGATG","AACTAC","AACTAG","AACTAG","ACTAGA","AACTAG","AGAACT","AACTAG","ATCTTG","AACTAG","CTAGAA","AACTAG","CTTGAT","AACTAG","GAACTA","AACTAG","GATCTT","AACTAG","TAGAAC","AACTAG","TCTTGA","AACTAG","TGATCT","AACTAG","TTGATC","AACTAG","AACTTC","AACTTC","AAGTTG","AACTTC","ACTTCA","AACTTC","AGTTGA","AACTTC","CAACTT","AACTTC","CTTCAA","AACTTC","GAAGTT","AACTTC","GTTGAA","AACTTC","TCAACT","AACTTC","TGAAGT","AACTTC","TTCAAC","AACTTC","TTGAAG","AACTTC","AACTTG","AACTTG","ACTTGA","AACTTG","CTTGAA","AACTTG","GAACTT","AACTTG","TGAACT","AACTTG","TTGAAC","AACTTG","AACTCT","AACTCT","ACTCTA","AACTCT","AGATTG","AACTCT","ATTGAG","AACTCT","CTAACT","AACTCT","CTCTAA","AACTCT","GAGATT","AACTCT","GATTGA","AACTCT","TAACTC","AACTCT","TCTAAC","AACTCT","TGAGAT","AACTCT","TTGAGA","AACTCT","AACTCC","AACTCC","ACTCCA","AACTCC","AGGTTG","AACTCC","CAACTC","AACTCC","CCAACT","AACTCC","CTCCAA","AACTCC","GAGGTT","AACTCC","GGTTGA","AACTCC","GTTGAG","AACTCC","TCCAAC","AACTCC","TGAGGT","AACTCC","TTGAGG","AACTCC","AACTCG","AACTCG","ACTCGA","AACTCG","AGCTTG","AACTCG","CGAACT","AACTCG","CTCGAA","AACTCG","CTTGAG","AACTCG","GAACTC","AACTCG","GAGCTT","AACTCG","GCTTGA","AACTCG","TCGAAC","AACTCG","TGAGCT","AACTCG","TTGAGC","AACTCG","AACTGT","AACTGT","ACATTG","AACTGT","ACTGTA","AACTGT","ATTGAC","AACTGT","CATTGA","AACTGT","CTGTAA","AACTGT","GACATT","AACTGT","GTAACT","AACTGT","TAACTG","AACTGT","TGACAT","AACTGT","TGTAAC","AACTGT","TTGACA","AACTGT","AACTGC","AACTGC","ACGTTG","AACTGC","ACTGCA","AACTGC","CAACTG","AACTGC","CGTTGA","AACTGC","CTGCAA","AACTGC","GACGTT","AACTGC","GCAACT","AACTGC","GTTGAC","AACTGC","TGACGT","AACTGC","TGCAAC","AACTGC","TTGACG","AACTGC","AACTGG","AACTGG","ACCTTG","AACTGG","ACTGGA","AACTGG","CCTTGA","AACTGG","CTGGAA","AACTGG","CTTGAC","AACTGG","GAACTG","AACTGG","GACCTT","AACTGG","GGAACT","AACTGG","TGACCT","AACTGG","TGGAAC","AACTGG","TTGACC","AACTGG","AACCAT","AACCAT","ACCATA","AACCAT","ATAACC","AACCAT","ATTGGT","AACCAT","CATAAC","AACCAT","CCATAA","AACCAT","GGTATT","AACCAT","GTATTG","AACCAT","TAACCA","AACCAT","TATTGG","AACCAT","TGGTAT","AACCAT","TTGGTA","AACCAT","AACCAC","AACCAC","ACAACC","AACCAC","ACCACA","AACCAC","CAACCA","AACCAC","CACAAC","AACCAC","CCACAA","AACCAC","GGTGTT","AACCAC","GTGTTG","AACCAC","GTTGGT","AACCAC","TGGTGT","AACCAC","TGTTGG","AACCAC","TTGGTG","AACCAC","AACCAG","AACCAG","ACCAGA","AACCAG","AGAACC","AACCAG","CAGAAC","AACCAG","CCAGAA","AACCAG","CTTGGT","AACCAG","GAACCA","AACCAG","GGTCTT","AACCAG","GTCTTG","AACCAG","TCTTGG","AACCAG","TGGTCT","AACCAG","TTGGTC","AACCAG","AACCTC","AACCTC","ACCTCA","AACCTC","AGTTGG","AACCTC","CAACCT","AACCTC","CCTCAA","AACCTC","CTCAAC","AACCTC","GAGTTG","AACCTC","GGAGTT","AACCTC","GTTGGA","AACCTC","TCAACC","AACCTC","TGGAGT","AACCTC","TTGGAG","AACCTC","AACCTG","AACCTG","ACCTGA","AACCTG","ACTTGG","AACCTG","CCTGAA","AACCTG","CTGAAC","AACCTG","CTTGGA","AACCTG","GAACCT","AACCTG","GACTTG","AACCTG","GGACTT","AACCTG","TGAACC","AACCTG","TGGACT","AACCTG","TTGGAC","AACCTG","AACCCT","AACCCT","ACCCTA","AACCCT","ATTGGG","AACCCT","CCCTAA","AACCCT","CCTAAC","AACCCT","CTAACC","AACCCT","GATTGG","AACCCT","GGATTG","AACCCT","GGGATT","AACCCT","TAACCC","AACCCT","TGGGAT","AACCCT","TTGGGA","AACCCT","AACCCC","AACCCC","ACCCCA","AACCCC","CAACCC","AACCCC","CCAACC","AACCCC","CCCAAC","AACCCC","CCCCAA","AACCCC","GGGGTT","AACCCC","GGGTTG","AACCCC","GGTTGG","AACCCC","GTTGGG","AACCCC","TGGGGT","AACCCC","TTGGGG","AACCCC","AACCCG","AACCCG","ACCCGA","AACCCG","CCCGAA","AACCCG","CCGAAC","AACCCG","CGAACC","AACCCG","CTTGGG","AACCCG","GAACCC","AACCCG","GCTTGG","AACCCG","GGCTTG","AACCCG","GGGCTT","AACCCG","TGGGCT","AACCCG","TTGGGC","AACCCG","AACCGT","AACCGT","ACCGTA","AACCGT","ATTGGC","AACCGT","CATTGG","AACCGT","CCGTAA","AACCGT","CGTAAC","AACCGT","GCATTG","AACCGT","GGCATT","AACCGT","GTAACC","AACCGT","TAACCG","AACCGT","TGGCAT","AACCGT","TTGGCA","AACCGT","AACCGC","AACCGC","ACCGCA","AACCGC","CAACCG","AACCGC","CCGCAA","AACCGC","CGCAAC","AACCGC","CGTTGG","AACCGC","GCAACC","AACCGC","GCGTTG","AACCGC","GGCGTT","AACCGC","GTTGGC","AACCGC","TGGCGT","AACCGC","TTGGCG","AACCGC","AACCGG","AACCGG","ACCGGA","AACCGG","CCGGAA","AACCGG","CCTTGG","AACCGG","CGGAAC","AACCGG","CTTGGC","AACCGG","GAACCG","AACCGG","GCCTTG","AACCGG","GGAACC","AACCGG","GGCCTT","AACCGG","TGGCCT","AACCGG","TTGGCC","AACCGG","AACGAT","AACGAT","ACGATA","AACGAT","ATAACG","AACGAT","ATTGCT","AACGAT","CGATAA","AACGAT","CTATTG","AACGAT","GATAAC","AACGAT","GCTATT","AACGAT","TAACGA","AACGAT","TATTGC","AACGAT","TGCTAT","AACGAT","TTGCTA","AACGAT","AACGAC","AACGAC","ACAACG","AACGAC","ACGACA","AACGAC","CAACGA","AACGAC","CGACAA","AACGAC","CTGTTG","AACGAC","GACAAC","AACGAC","GCTGTT","AACGAC","GTTGCT","AACGAC","TGCTGT","AACGAC","TGTTGC","AACGAC","TTGCTG","AACGAC","AACGAG","AACGAG","ACGAGA","AACGAG","AGAACG","AACGAG","CGAGAA","AACGAG","CTCTTG","AACGAG","CTTGCT","AACGAG","GAACGA","AACGAG","GAGAAC","AACGAG","GCTCTT","AACGAG","TCTTGC","AACGAG","TGCTCT","AACGAG","TTGCTC","AACGAG","AACGTC","AACGTC","ACGTCA","AACGTC","AGTTGC","AACGTC","CAACGT","AACGTC","CAGTTG","AACGTC","CGTCAA","AACGTC","GCAGTT","AACGTC","GTCAAC","AACGTC","GTTGCA","AACGTC","TCAACG","AACGTC","TGCAGT","AACGTC","TTGCAG","AACGTC","AACGTG","AACGTG","ACGTGA","AACGTG","ACTTGC","AACGTG","CACTTG","AACGTG","CGTGAA","AACGTG","CTTGCA","AACGTG","GAACGT","AACGTG","GCACTT","AACGTG","GTGAAC","AACGTG","TGAACG","AACGTG","TGCACT","AACGTG","TTGCAC","AACGTG","AACGCT","AACGCT","ACGCTA","AACGCT","ATTGCG","AACGCT","CGATTG","AACGCT","CGCTAA","AACGCT","CTAACG","AACGCT","GATTGC","AACGCT","GCGATT","AACGCT","GCTAAC","AACGCT","TAACGC","AACGCT","TGCGAT","AACGCT","TTGCGA","AACGCT","AACGCC","AACGCC","ACGCCA","AACGCC","CAACGC","AACGCC","CCAACG","AACGCC","CGCCAA","AACGCC","CGGTTG","AACGCC","GCCAAC","AACGCC","GCGGTT","AACGCC","GGTTGC","AACGCC","GTTGCG","AACGCC","TGCGGT","AACGCC","TTGCGG","AACGCC","AACGCG","AACGCG","ACGCGA","AACGCG","CGAACG","AACGCG","CGCGAA","AACGCG","CGCTTG","AACGCG","CTTGCG","AACGCG","GAACGC","AACGCG","GCGAAC","AACGCG","GCGCTT","AACGCG","GCTTGC","AACGCG","TGCGCT","AACGCG","TTGCGC","AACGCG","AACGGT","AACGGT","ACGGTA","AACGGT","ATTGCC","AACGGT","CATTGC","AACGGT","CCATTG","AACGGT","CGGTAA","AACGGT","GCCATT","AACGGT","GGTAAC","AACGGT","GTAACG","AACGGT","TAACGG","AACGGT","TGCCAT","AACGGT","TTGCCA","AACGGT","AACGGC","AACGGC","ACGGCA","AACGGC","CAACGG","AACGGC","CCGTTG","AACGGC","CGGCAA","AACGGC","CGTTGC","AACGGC","GCAACG","AACGGC","GCCGTT","AACGGC","GGCAAC","AACGGC","GTTGCC","AACGGC","TGCCGT","AACGGC","TTGCCG","AACGGC","AACGGG","AACGGG","ACGGGA","AACGGG","CCCTTG","AACGGG","CCTTGC","AACGGG","CGGGAA","AACGGG","CTTGCC","AACGGG","GAACGG","AACGGG","GCCCTT","AACGGG","GGAACG","AACGGG","GGGAAC","AACGGG","TGCCCT","AACGGG","TTGCCC","AACGGG","AAGATC","AAGATC","AGATCA","AAGATC","AGTTCT","AAGATC","ATCAAG","AAGATC","CAAGAT","AAGATC","CTAGTT","AAGATC","GATCAA","AAGATC","GTTCTA","AAGATC","TAGTTC","AAGATC","TCAAGA","AAGATC","TCTAGT","AAGATC","TTCTAG","AAGATC","AAGATG","AAGATG","ACTTCT","AAGATG","AGATGA","AAGATG","ATGAAG","AAGATG","CTACTT","AAGATG","CTTCTA","AAGATG","GAAGAT","AAGATG","GATGAA","AAGATG","TACTTC","AAGATG","TCTACT","AAGATG","TGAAGA","AAGATG","TTCTAC","AAGATG","AAGACT","AAGACT","ACTAAG","AAGACT","AGACTA","AAGACT","ATTCTG","AAGACT","CTAAGA","AAGACT","CTGATT","AAGACT","GACTAA","AAGACT","GATTCT","AAGACT","TAAGAC","AAGACT","TCTGAT","AAGACT","TGATTC","AAGACT","TTCTGA","AAGACT","AAGACC","AAGACC","ACCAAG","AAGACC","AGACCA","AAGACC","CAAGAC","AAGACC","CCAAGA","AAGACC","CTGGTT","AAGACC","GACCAA","AAGACC","GGTTCT","AAGACC","GTTCTG","AAGACC","TCTGGT","AAGACC","TGGTTC","AAGACC","TTCTGG","AAGACC","AAGACG","AAGACG","ACGAAG","AAGACG","AGACGA","AAGACG","CGAAGA","AAGACG","CTGCTT","AAGACG","CTTCTG","AAGACG","GAAGAC","AAGACG","GACGAA","AAGACG","GCTTCT","AAGACG","TCTGCT","AAGACG","TGCTTC","AAGACG","TTCTGC","AAGACG","AAGAGT","AAGAGT","AGAGTA","AAGAGT","AGTAAG","AAGAGT","ATTCTC","AAGAGT","CATTCT","AAGAGT","CTCATT","AAGAGT","GAGTAA","AAGAGT","GTAAGA","AAGAGT","TAAGAG","AAGAGT","TCATTC","AAGAGT","TCTCAT","AAGAGT","TTCTCA","AAGAGT","AAGAGC","AAGAGC","AGAGCA","AAGAGC","AGCAAG","AAGAGC","CAAGAG","AAGAGC","CGTTCT","AAGAGC","CTCGTT","AAGAGC","GAGCAA","AAGAGC","GCAAGA","AAGAGC","GTTCTC","AAGAGC","TCGTTC","AAGAGC","TCTCGT","AAGAGC","TTCTCG","AAGAGC","AAGAGG","AAGAGG","AGAGGA","AAGAGG","AGGAAG","AAGAGG","CCTTCT","AAGAGG","CTCCTT","AAGAGG","CTTCTC","AAGAGG","GAAGAG","AAGAGG","GAGGAA","AAGAGG","GGAAGA","AAGAGG","TCCTTC","AAGAGG","TCTCCT","AAGAGG","TTCTCC","AAGAGG","AAGTAT","AAGTAT","AGTATA","AAGTAT","ATAAGT","AAGTAT","ATATTC","AAGTAT","ATTCAT","AAGTAT","CATATT","AAGTAT","GTATAA","AAGTAT","TAAGTA","AAGTAT","TATAAG","AAGTAT","TATTCA","AAGTAT","TCATAT","AAGTAT","TTCATA","AAGTAT","AAGTAC","AAGTAC","ACAAGT","AAGTAC","AGTACA","AAGTAC","ATGTTC","AAGTAC","CAAGTA","AAGTAC","CATGTT","AAGTAC","GTACAA","AAGTAC","GTTCAT","AAGTAC","TACAAG","AAGTAC","TCATGT","AAGTAC","TGTTCA","AAGTAC","TTCATG","AAGTAC","AAGTAG","AAGTAG","AGAAGT","AAGTAG","AGTAGA","AAGTAG","ATCTTC","AAGTAG","CATCTT","AAGTAG","CTTCAT","AAGTAG","GAAGTA","AAGTAG","GTAGAA","AAGTAG","TAGAAG","AAGTAG","TCATCT","AAGTAG","TCTTCA","AAGTAG","TTCATC","AAGTAG","AAGTTC","AAGTTC","AGTTCA","AAGTTC","CAAGTT","AAGTTC","GTTCAA","AAGTTC","TCAAGT","AAGTTC","TTCAAG","AAGTTC","AAGTCT","AAGTCT","AGATTC","AAGTCT","AGTCTA","AAGTCT","ATTCAG","AAGTCT","CAGATT","AAGTCT","CTAAGT","AAGTCT","GATTCA","AAGTCT","GTCTAA","AAGTCT","TAAGTC","AAGTCT","TCAGAT","AAGTCT","TCTAAG","AAGTCT","TTCAGA","AAGTCT","AAGTCC","AAGTCC","AGGTTC","AAGTCC","AGTCCA","AAGTCC","CAAGTC","AAGTCC","CAGGTT","AAGTCC","CCAAGT","AAGTCC","GGTTCA","AAGTCC","GTCCAA","AAGTCC","GTTCAG","AAGTCC","TCAGGT","AAGTCC","TCCAAG","AAGTCC","TTCAGG","AAGTCC","AAGTCG","AAGTCG","AGCTTC","AAGTCG","AGTCGA","AAGTCG","CAGCTT","AAGTCG","CGAAGT","AAGTCG","CTTCAG","AAGTCG","GAAGTC","AAGTCG","GCTTCA","AAGTCG","GTCGAA","AAGTCG","TCAGCT","AAGTCG","TCGAAG","AAGTCG","TTCAGC","AAGTCG","AAGTGT","AAGTGT","ACATTC","AAGTGT","AGTGTA","AAGTGT","ATTCAC","AAGTGT","CACATT","AAGTGT","CATTCA","AAGTGT","GTAAGT","AAGTGT","GTGTAA","AAGTGT","TAAGTG","AAGTGT","TCACAT","AAGTGT","TGTAAG","AAGTGT","TTCACA","AAGTGT","AAGTGC","AAGTGC","ACGTTC","AAGTGC","AGTGCA","AAGTGC","CAAGTG","AAGTGC","CACGTT","AAGTGC","CGTTCA","AAGTGC","GCAAGT","AAGTGC","GTGCAA","AAGTGC","GTTCAC","AAGTGC","TCACGT","AAGTGC","TGCAAG","AAGTGC","TTCACG","AAGTGC","AAGTGG","AAGTGG","ACCTTC","AAGTGG","AGTGGA","AAGTGG","CACCTT","AAGTGG","CCTTCA","AAGTGG","CTTCAC","AAGTGG","GAAGTG","AAGTGG","GGAAGT","AAGTGG","GTGGAA","AAGTGG","TCACCT","AAGTGG","TGGAAG","AAGTGG","TTCACC","AAGTGG","AAGCAT","AAGCAT","AGCATA","AAGCAT","ATAAGC","AAGCAT","ATTCGT","AAGCAT","CATAAG","AAGCAT","CGTATT","AAGCAT","GCATAA","AAGCAT","GTATTC","AAGCAT","TAAGCA","AAGCAT","TATTCG","AAGCAT","TCGTAT","AAGCAT","TTCGTA","AAGCAT","AAGCAC","AAGCAC","ACAAGC","AAGCAC","AGCACA","AAGCAC","CAAGCA","AAGCAC","CACAAG","AAGCAC","CGTGTT","AAGCAC","GCACAA","AAGCAC","GTGTTC","AAGCAC","GTTCGT","AAGCAC","TCGTGT","AAGCAC","TGTTCG","AAGCAC","TTCGTG","AAGCAC","AAGCAG","AAGCAG","AGAAGC","AAGCAG","AGCAGA","AAGCAG","CAGAAG","AAGCAG","CGTCTT","AAGCAG","CTTCGT","AAGCAG","GAAGCA","AAGCAG","GCAGAA","AAGCAG","GTCTTC","AAGCAG","TCGTCT","AAGCAG","TCTTCG","AAGCAG","TTCGTC","AAGCAG","AAGCTC","AAGCTC","AGCTCA","AAGCTC","AGTTCG","AAGCTC","CAAGCT","AAGCTC","CGAGTT","AAGCTC","CTCAAG","AAGCTC","GAGTTC","AAGCTC","GCTCAA","AAGCTC","GTTCGA","AAGCTC","TCAAGC","AAGCTC","TCGAGT","AAGCTC","TTCGAG","AAGCTC","AAGCTG","AAGCTG","ACTTCG","AAGCTG","AGCTGA","AAGCTG","CGACTT","AAGCTG","CTGAAG","AAGCTG","CTTCGA","AAGCTG","GAAGCT","AAGCTG","GACTTC","AAGCTG","GCTGAA","AAGCTG","TCGACT","AAGCTG","TGAAGC","AAGCTG","TTCGAC","AAGCTG","AAGCCT","AAGCCT","AGCCTA","AAGCCT","ATTCGG","AAGCCT","CCTAAG","AAGCCT","CGGATT","AAGCCT","CTAAGC","AAGCCT","GATTCG","AAGCCT","GCCTAA","AAGCCT","GGATTC","AAGCCT","TAAGCC","AAGCCT","TCGGAT","AAGCCT","TTCGGA","AAGCCT","AAGCCC","AAGCCC","AGCCCA","AAGCCC","CAAGCC","AAGCCC","CCAAGC","AAGCCC","CCCAAG","AAGCCC","CGGGTT","AAGCCC","GCCCAA","AAGCCC","GGGTTC","AAGCCC","GGTTCG","AAGCCC","GTTCGG","AAGCCC","TCGGGT","AAGCCC","TTCGGG","AAGCCC","AAGCCG","AAGCCG","AGCCGA","AAGCCG","CCGAAG","AAGCCG","CGAAGC","AAGCCG","CGGCTT","AAGCCG","CTTCGG","AAGCCG","GAAGCC","AAGCCG","GCCGAA","AAGCCG","GCTTCG","AAGCCG","GGCTTC","AAGCCG","TCGGCT","AAGCCG","TTCGGC","AAGCCG","AAGCGT","AAGCGT","AGCGTA","AAGCGT","ATTCGC","AAGCGT","CATTCG","AAGCGT","CGCATT","AAGCGT","CGTAAG","AAGCGT","GCATTC","AAGCGT","GCGTAA","AAGCGT","GTAAGC","AAGCGT","TAAGCG","AAGCGT","TCGCAT","AAGCGT","TTCGCA","AAGCGT","AAGCGC","AAGCGC","AGCGCA","AAGCGC","CAAGCG","AAGCGC","CGCAAG","AAGCGC","CGCGTT","AAGCGC","CGTTCG","AAGCGC","GCAAGC","AAGCGC","GCGCAA","AAGCGC","GCGTTC","AAGCGC","GTTCGC","AAGCGC","TCGCGT","AAGCGC","TTCGCG","AAGCGC","AAGCGG","AAGCGG","AGCGGA","AAGCGG","CCTTCG","AAGCGG","CGCCTT","AAGCGG","CGGAAG","AAGCGG","CTTCGC","AAGCGG","GAAGCG","AAGCGG","GCCTTC","AAGCGG","GCGGAA","AAGCGG","GGAAGC","AAGCGG","TCGCCT","AAGCGG","TTCGCC","AAGCGG","AAGGAT","AAGGAT","AGGATA","AAGGAT","ATAAGG","AAGGAT","ATTCCT","AAGGAT","CCTATT","AAGGAT","CTATTC","AAGGAT","GATAAG","AAGGAT","GGATAA","AAGGAT","TAAGGA","AAGGAT","TATTCC","AAGGAT","TCCTAT","AAGGAT","TTCCTA","AAGGAT","AAGGAC","AAGGAC","ACAAGG","AAGGAC","AGGACA","AAGGAC","CAAGGA","AAGGAC","CCTGTT","AAGGAC","CTGTTC","AAGGAC","GACAAG","AAGGAC","GGACAA","AAGGAC","GTTCCT","AAGGAC","TCCTGT","AAGGAC","TGTTCC","AAGGAC","TTCCTG","AAGGAC","AAGGAG","AAGGAG","AGAAGG","AAGGAG","AGGAGA","AAGGAG","CCTCTT","AAGGAG","CTCTTC","AAGGAG","CTTCCT","AAGGAG","GAAGGA","AAGGAG","GAGAAG","AAGGAG","GGAGAA","AAGGAG","TCCTCT","AAGGAG","TCTTCC","AAGGAG","TTCCTC","AAGGAG","AAGGTC","AAGGTC","AGGTCA","AAGGTC","AGTTCC","AAGGTC","CAAGGT","AAGGTC","CAGTTC","AAGGTC","CCAGTT","AAGGTC","GGTCAA","AAGGTC","GTCAAG","AAGGTC","GTTCCA","AAGGTC","TCAAGG","AAGGTC","TCCAGT","AAGGTC","TTCCAG","AAGGTC","AAGGTG","AAGGTG","ACTTCC","AAGGTG","AGGTGA","AAGGTG","CACTTC","AAGGTG","CCACTT","AAGGTG","CTTCCA","AAGGTG","GAAGGT","AAGGTG","GGTGAA","AAGGTG","GTGAAG","AAGGTG","TCCACT","AAGGTG","TGAAGG","AAGGTG","TTCCAC","AAGGTG","AAGGCT","AAGGCT","AGGCTA","AAGGCT","ATTCCG","AAGGCT","CCGATT","AAGGCT","CGATTC","AAGGCT","CTAAGG","AAGGCT","GATTCC","AAGGCT","GCTAAG","AAGGCT","GGCTAA","AAGGCT","TAAGGC","AAGGCT","TCCGAT","AAGGCT","TTCCGA","AAGGCT","AAGGCC","AAGGCC","AGGCCA","AAGGCC","CAAGGC","AAGGCC","CCAAGG","AAGGCC","CCGGTT","AAGGCC","CGGTTC","AAGGCC","GCCAAG","AAGGCC","GGCCAA","AAGGCC","GGTTCC","AAGGCC","GTTCCG","AAGGCC","TCCGGT","AAGGCC","TTCCGG","AAGGCC","AAGGCG","AAGGCG","AGGCGA","AAGGCG","CCGCTT","AAGGCG","CGAAGG","AAGGCG","CGCTTC","AAGGCG","CTTCCG","AAGGCG","GAAGGC","AAGGCG","GCGAAG","AAGGCG","GCTTCC","AAGGCG","GGCGAA","AAGGCG","TCCGCT","AAGGCG","TTCCGC","AAGGCG","AAGGGT","AAGGGT","AGGGTA","AAGGGT","ATTCCC","AAGGGT","CATTCC","AAGGGT","CCATTC","AAGGGT","CCCATT","AAGGGT","GGGTAA","AAGGGT","GGTAAG","AAGGGT","GTAAGG","AAGGGT","TAAGGG","AAGGGT","TCCCAT","AAGGGT","TTCCCA","AAGGGT","AAGGGC","AAGGGC","AGGGCA","AAGGGC","CAAGGG","AAGGGC","CCCGTT","AAGGGC","CCGTTC","AAGGGC","CGTTCC","AAGGGC","GCAAGG","AAGGGC","GGCAAG","AAGGGC","GGGCAA","AAGGGC","GTTCCC","AAGGGC","TCCCGT","AAGGGC","TTCCCG","AAGGGC","AAGGGG","AAGGGG","AGGGGA","AAGGGG","CCCCTT","AAGGGG","CCCTTC","AAGGGG","CCTTCC","AAGGGG","CTTCCC","AAGGGG","GAAGGG","AAGGGG","GGAAGG","AAGGGG","GGGAAG","AAGGGG","GGGGAA","AAGGGG","TCCCCT","AAGGGG","TTCCCC","AAGGGG","ACATAT","ACATAT","ATACAT","ACATAT","ATATAC","ACATAT","ATATGT","ACATAT","ATGTAT","ACATAT","CATATA","ACATAT","GTATAT","ACATAT","TACATA","ACATAT","TATACA","ACATAT","TATATG","ACATAT","TATGTA","ACATAT","TGTATA","ACATAT","AGATAT","AGATAT","ATAGAT","AGATAT","ATATAG","AGATAT","ATATCT","AGATAT","ATCTAT","AGATAT","CTATAT","AGATAT","GATATA","AGATAT","TAGATA","AGATAT","TATAGA","AGATAT","TATATC","AGATAT","TATCTA","AGATAT","TCTATA","AGATAT","AGGTAT","AGGTAT","ATAGGT","AGGTAT","ATATCC","AGGTAT","ATCCAT","AGGTAT","CATATC","AGGTAT","CCATAT","AGGTAT","GGTATA","AGGTAT","GTATAG","AGGTAT","TAGGTA","AGGTAT","TATAGG","AGGTAT","TATCCA","AGGTAT","TCCATA","AGGTAT","AGCTAT","AGCTAT","ATAGCT","AGCTAT","ATATCG","AGCTAT","ATCGAT","AGCTAT","CGATAT","AGCTAT","CTATAG","AGCTAT","GATATC","AGCTAT","GCTATA","AGCTAT","TAGCTA","AGCTAT","TATAGC","AGCTAT","TATCGA","AGCTAT","TCGATA","AGCTAT","ACGTAT","ACGTAT","ATACGT","ACGTAT","ATATGC","ACGTAT","ATGCAT","ACGTAT","CATATG","ACGTAT","CGTATA","ACGTAT","GCATAT","ACGTAT","GTATAC","ACGTAT","TACGTA","ACGTAT","TATACG","ACGTAT","TATGCA","ACGTAT","TGCATA","ACGTAT","ACCTAT","ACCTAT","ATACCT","ACCTAT","ATATGG","ACCTAT","ATGGAT","ACCTAT","CCTATA","ACCTAT","CTATAC","ACCTAT","GATATG","ACCTAT","GGATAT","ACCTAT","TACCTA","ACCTAT","TATACC","ACCTAT","TATGGA","ACCTAT","TGGATA","ACCTAT","ACACAT","ACACAT","ACATAC","ACACAT","ATACAC","ACACAT","ATGTGT","ACACAT","CACATA","ACACAT","CATACA","ACACAT","GTATGT","ACACAT","GTGTAT","ACACAT","TACACA","ACACAT","TATGTG","ACACAT","TGTATG","ACACAT","TGTGTA","ACACAT","ACAGAT","ACAGAT","AGATAC","ACAGAT","ATACAG","ACAGAT","ATGTCT","ACAGAT","CAGATA","ACAGAT","CTATGT","ACAGAT","GATACA","ACAGAT","GTCTAT","ACAGAT","TACAGA","ACAGAT","TATGTC","ACAGAT","TCTATG","ACAGAT","TGTCTA","ACAGAT","ACTCAT","ACTCAT","AGTATG","ACTCAT","ATACTC","ACTCAT","ATGAGT","ACTCAT","CATACT","ACTCAT","CTCATA","ACTCAT","GAGTAT","ACTCAT","GTATGA","ACTCAT","TACTCA","ACTCAT","TATGAG","ACTCAT","TCATAC","ACTCAT","TGAGTA","ACTCAT","ACTATG","ACTATG","ACTGAT","ACTATG","ATACTG","ACTATG","ATGACT","ACTATG","CTATGA","ACTATG","CTGATA","ACTATG","GACTAT","ACTATG","GATACT","ACTATG","TACTGA","ACTATG","TATGAC","ACTATG","TGACTA","ACTATG","TGATAC","ACTATG","ACCCAT","ACCCAT","ATACCC","ACCCAT","ATGGGT","ACCCAT","CATACC","ACCCAT","CCATAC","ACCCAT","CCCATA","ACCCAT","GGGTAT","ACCCAT","GGTATG","ACCCAT","GTATGG","ACCCAT","TACCCA","ACCCAT","TATGGG","ACCCAT","TGGGTA","ACCCAT","ACCGAT","ACCGAT","ATACCG","ACCGAT","ATGGCT","ACCGAT","CCGATA","ACCGAT","CGATAC","ACCGAT","CTATGG","ACCGAT","GATACC","ACCGAT","GCTATG","ACCGAT","GGCTAT","ACCGAT","TACCGA","ACCGAT","TATGGC","ACCGAT","TGGCTA","ACCGAT","ACGCAT","ACGCAT","ATACGC","ACGCAT","ATGCGT","ACGCAT","CATACG","ACGCAT","CGCATA","ACGCAT","CGTATG","ACGCAT","GCATAC","ACGCAT","GCGTAT","ACGCAT","GTATGC","ACGCAT","TACGCA","ACGCAT","TATGCG","ACGCAT","TGCGTA","ACGCAT","ACGGAT","ACGGAT","ATACGG","ACGGAT","ATGCCT","ACGGAT","CCTATG","ACGGAT","CGGATA","ACGGAT","CTATGC","ACGGAT","GATACG","ACGGAT","GCCTAT","ACGGAT","GGATAC","ACGGAT","TACGGA","ACGGAT","TATGCC","ACGGAT","TGCCTA","ACGGAT","ACATAG","ACATAG","AGACAT","ACATAG","ATAGAC","ACATAG","ATCTGT","ACATAG","CATAGA","ACATAG","CTGTAT","ACATAG","GACATA","ACATAG","GTATCT","ACATAG","TAGACA","ACATAG","TATCTG","ACATAG","TCTGTA","ACATAG","TGTATC","ACATAG","AGAGAT","AGAGAT","AGATAG","AGAGAT","ATAGAG","AGAGAT","ATCTCT","AGAGAT","CTATCT","AGAGAT","CTCTAT","AGAGAT","GAGATA","AGAGAT","GATAGA","AGAGAT","TAGAGA","AGAGAT","TATCTC","AGAGAT","TCTATC","AGAGAT","TCTCTA","AGAGAT","AGTATC","AGTATC","AGTCAT","AGTATC","ATAGTC","AGTATC","ATCAGT","AGTATC","CAGTAT","AGTATC","CATAGT","AGTATC","GTATCA","AGTATC","GTCATA","AGTATC","TAGTCA","AGTATC","TATCAG","AGTATC","TCAGTA","AGTATC","TCATAG","AGTATC","ACTATC","ACTATC","AGTGAT","ACTATC","ATAGTG","ACTATC","ATCACT","ACTATC","CACTAT","ACTATC","CTATCA","ACTATC","GATAGT","ACTATC","GTGATA","ACTATC","TAGTGA","ACTATC","TATCAC","ACTATC","TCACTA","ACTATC","TGATAG","ACTATC","AGCCAT","AGCCAT","ATAGCC","AGCCAT","ATCGGT","AGCCAT","CATAGC","AGCCAT","CCATAG","AGCCAT","CGGTAT","AGCCAT","GCCATA","AGCCAT","GGTATC","AGCCAT","GTATCG","AGCCAT","TAGCCA","AGCCAT","TATCGG","AGCCAT","TCGGTA","AGCCAT","AGCGAT","AGCGAT","ATAGCG","AGCGAT","ATCGCT","AGCGAT","CGATAG","AGCGAT","CGCTAT","AGCGAT","CTATCG","AGCGAT","GATAGC","AGCGAT","GCGATA","AGCGAT","GCTATC","AGCGAT","TAGCGA","AGCGAT","TATCGC","AGCGAT","TCGCTA","AGCGAT","AGGCAT","AGGCAT","ATAGGC","AGGCAT","ATCCGT","AGGCAT","CATAGG","AGGCAT","CCGTAT","AGGCAT","CGTATC","AGGCAT","GCATAG","AGGCAT","GGCATA","AGGCAT","GTATCC","AGGCAT","TAGGCA","AGGCAT","TATCCG","AGGCAT","TCCGTA","AGGCAT","AGGGAT","AGGGAT","ATAGGG","AGGGAT","ATCCCT","AGGGAT","CCCTAT","AGGGAT","CCTATC","AGGGAT","CTATCC","AGGGAT","GATAGG","AGGGAT","GGATAG","AGGGAT","GGGATA","AGGGAT","TAGGGA","AGGGAT","TATCCC","AGGGAT","TCCCTA","AGGGAT","ACTAGT","ACTAGT","AGTACT","ACTAGT","ATCATG","ACTAGT","ATGATC","ACTAGT","CATGAT","ACTAGT","CTAGTA","ACTAGT","GATCAT","ACTAGT","GTACTA","ACTAGT","TACTAG","ACTAGT","TAGTAC","ACTAGT","TCATGA","ACTAGT","TGATCA","ACTAGT","ACCATC","ACCATC","AGTGGT","ACCATC","ATCACC","ACCATC","CACCAT","ACCATC","CATCAC","ACCATC","CCATCA","ACCATC","GGTAGT","ACCATC","GTAGTG","ACCATC","GTGGTA","ACCATC","TAGTGG","ACCATC","TCACCA","ACCATC","TGGTAG","ACCATC","ACGATC","ACGATC","AGTGCT","ACGATC","ATCACG","ACGATC","CACGAT","ACGATC","CGATCA","ACGATC","CTAGTG","ACGATC","GATCAC","ACGATC","GCTAGT","ACGATC","GTGCTA","ACGATC","TAGTGC","ACGATC","TCACGA","ACGATC","TGCTAG","ACGATC","AGCATC","AGCATC","AGTCGT","AGCATC","ATCAGC","AGCATC","CAGCAT","AGCATC","CATCAG","AGCATC","CGTAGT","AGCATC","GCATCA","AGCATC","GTAGTC","AGCATC","GTCGTA","AGCATC","TAGTCG","AGCATC","TCAGCA","AGCATC","TCGTAG","AGCATC","AGGATC","AGGATC","AGTCCT","AGGATC","ATCAGG","AGGATC","CAGGAT","AGGATC","CCTAGT","AGGATC","CTAGTC","AGGATC","GATCAG","AGGATC","GGATCA","AGGATC","GTCCTA","AGGATC","TAGTCC","AGGATC","TCAGGA","AGGATC","TCCTAG","AGGATC","ACATCT","ACATCT","AGATGT","ACATCT","ATCTAC","ACATCT","ATGTAG","ACATCT","CATCTA","ACATCT","CTACAT","ACATCT","GATGTA","ACATCT","GTAGAT","ACATCT","TACATC","ACATCT","TAGATG","ACATCT","TCTACA","ACATCT","TGTAGA","ACATCT","AGATCT","AGATCT","ATCTAG","AGATCT","CTAGAT","AGATCT","GATCTA","AGATCT","TAGATC","AGATCT","TCTAGA","AGATCT","AGAGGT","AGAGGT","AGGTAG","AGAGGT","ATCTCC","AGAGGT","CATCTC","AGAGGT","CCATCT","AGAGGT","CTCCAT","AGAGGT","GAGGTA","AGAGGT","GGTAGA","AGAGGT","GTAGAG","AGAGGT","TAGAGG","AGAGGT","TCCATC","AGAGGT","TCTCCA","AGAGGT","AGAGCT","AGAGCT","AGCTAG","AGAGCT","ATCTCG","AGAGCT","CGATCT","AGAGCT","CTAGAG","AGAGCT","CTCGAT","AGAGCT","GAGCTA","AGAGCT","GATCTC","AGAGCT","GCTAGA","AGAGCT","TAGAGC","AGAGCT","TCGATC","AGAGCT","TCTCGA","AGAGCT","ACGTAG","ACGTAG","AGACGT","ACGTAG","ATCTGC","ACGTAG","CATCTG","ACGTAG","CGTAGA","ACGTAG","CTGCAT","ACGTAG","GACGTA","ACGTAG","GCATCT","ACGTAG","GTAGAC","ACGTAG","TAGACG","ACGTAG","TCTGCA","ACGTAG","TGCATC","ACGTAG","ACCTAG","ACCTAG","AGACCT","ACCTAG","ATCTGG","ACCTAG","CCTAGA","ACCTAG","CTAGAC","ACCTAG","CTGGAT","ACCTAG","GACCTA","ACCTAG","GATCTG","ACCTAG","GGATCT","ACCTAG","TAGACC","ACCTAG","TCTGGA","ACCTAG","TGGATC","ACCTAG","ACATCC","ACATCC","AGGTGT","ACATCC","ATCCAC","ACATCC","CACATC","ACATCC","CATCCA","ACATCC","CCACAT","ACATCC","GGTGTA","ACATCC","GTAGGT","ACATCC","GTGTAG","ACATCC","TAGGTG","ACATCC","TCCACA","ACATCC","TGTAGG","ACATCC","AGATCC","AGATCC","AGGTCT","AGATCC","ATCCAG","AGATCC","CAGATC","AGATCC","CCAGAT","AGATCC","CTAGGT","AGATCC","GATCCA","AGATCC","GGTCTA","AGATCC","GTCTAG","AGATCC","TAGGTC","AGATCC","TCCAGA","AGATCC","TCTAGG","AGATCC","AGGAGT","AGGAGT","AGTAGG","AGGAGT","ATCCTC","AGGAGT","CATCCT","AGGAGT","CCTCAT","AGGAGT","CTCATC","AGGAGT","GAGTAG","AGGAGT","GGAGTA","AGGAGT","GTAGGA","AGGAGT","TAGGAG","AGGAGT","TCATCC","AGGAGT","TCCTCA","AGGAGT","ACTAGG","ACTAGG","AGGACT","ACTAGG","ATCCTG","ACTAGG","CCTGAT","ACTAGG","CTAGGA","ACTAGG","CTGATC","ACTAGG","GACTAG","ACTAGG","GATCCT","ACTAGG","GGACTA","ACTAGG","TAGGAC","ACTAGG","TCCTGA","ACTAGG","TGATCC","ACTAGG","AGGGGT","AGGGGT","ATCCCC","AGGGGT","CATCCC","AGGGGT","CCATCC","AGGGGT","CCCATC","AGGGGT","CCCCAT","AGGGGT","GGGGTA","AGGGGT","GGGTAG","AGGGGT","GGTAGG","AGGGGT","GTAGGG","AGGGGT","TAGGGG","AGGGGT","TCCCCA","AGGGGT","AGGGCT","AGGGCT","ATCCCG","AGGGCT","CCCGAT","AGGGCT","CCGATC","AGGGCT","CGATCC","AGGGCT","CTAGGG","AGGGCT","GATCCC","AGGGCT","GCTAGG","AGGGCT","GGCTAG","AGGGCT","GGGCTA","AGGGCT","TAGGGC","AGGGCT","TCCCGA","AGGGCT","AGGCGT","AGGCGT","ATCCGC","AGGCGT","CATCCG","AGGCGT","CCGCAT","AGGCGT","CGCATC","AGGCGT","CGTAGG","AGGCGT","GCATCC","AGGCGT","GCGTAG","AGGCGT","GGCGTA","AGGCGT","GTAGGC","AGGCGT","TAGGCG","AGGCGT","TCCGCA","AGGCGT","AGGCCT","AGGCCT","ATCCGG","AGGCCT","CCGGAT","AGGCCT","CCTAGG","AGGCCT","CGGATC","AGGCCT","CTAGGC","AGGCCT","GATCCG","AGGCCT","GCCTAG","AGGCCT","GGATCC","AGGCCT","GGCCTA","AGGCCT","TAGGCC","AGGCCT","TCCGGA","AGGCCT","ACATCG","ACATCG","AGCTGT","ACATCG","ATCGAC","ACATCG","CATCGA","ACATCG","CGACAT","ACATCG","CTGTAG","ACATCG","GACATC","ACATCG","GCTGTA","ACATCG","GTAGCT","ACATCG","TAGCTG","ACATCG","TCGACA","ACATCG","TGTAGC","ACATCG","AGATCG","AGATCG","AGCTCT","AGATCG","ATCGAG","AGATCG","CGAGAT","AGATCG","CTAGCT","AGATCG","CTCTAG","AGATCG","GAGATC","AGATCG","GATCGA","AGATCG","GCTCTA","AGATCG","TAGCTC","AGATCG","TCGAGA","AGATCG","TCTAGC","AGATCG","AGCAGT","AGCAGT","AGTAGC","AGCAGT","ATCGTC","AGCAGT","CAGTAG","AGCAGT","CATCGT","AGCAGT","CGTCAT","AGCAGT","GCAGTA","AGCAGT","GTAGCA","AGCAGT","GTCATC","AGCAGT","TAGCAG","AGCAGT","TCATCG","AGCAGT","TCGTCA","AGCAGT","ACTAGC","ACTAGC","AGCACT","ACTAGC","ATCGTG","ACTAGC","CACTAG","ACTAGC","CGTGAT","ACTAGC","CTAGCA","ACTAGC","GATCGT","ACTAGC","GCACTA","ACTAGC","GTGATC","ACTAGC","TAGCAC","ACTAGC","TCGTGA","ACTAGC","TGATCG","ACTAGC","AGCGGT","AGCGGT","ATCGCC","AGCGGT","CATCGC","AGCGGT","CCATCG","AGCGGT","CGCCAT","AGCGGT","CGGTAG","AGCGGT","GCCATC","AGCGGT","GCGGTA","AGCGGT","GGTAGC","AGCGGT","GTAGCG","AGCGGT","TAGCGG","AGCGGT","TCGCCA","AGCGGT","AGCGCT","AGCGCT","ATCGCG","AGCGCT","CGATCG","AGCGCT","CGCGAT","AGCGCT","CGCTAG","AGCGCT","CTAGCG","AGCGCT","GATCGC","AGCGCT","GCGATC","AGCGCT","GCGCTA","AGCGCT","GCTAGC","AGCGCT","TAGCGC","AGCGCT","TCGCGA","AGCGCT","AGCCGT","AGCCGT","ATCGGC","AGCCGT","CATCGG","AGCCGT","CCGTAG","AGCCGT","CGGCAT","AGCCGT","CGTAGC","AGCCGT","GCATCG","AGCCGT","GCCGTA","AGCCGT","GGCATC","AGCCGT","GTAGCC","AGCCGT","TAGCCG","AGCCGT","TCGGCA","AGCCGT","AGCCCT","AGCCCT","ATCGGG","AGCCCT","CCCTAG","AGCCCT","CCTAGC","AGCCCT","CGGGAT","AGCCCT","CTAGCC","AGCCCT","GATCGG","AGCCCT","GCCCTA","AGCCCT","GGATCG","AGCCCT","GGGATC","AGCCCT","TAGCCC","AGCCCT","TCGGGA","AGCCCT","ACCATG","ACCATG","ACTGGT","ACCATG","ATGACC","ACCATG","CATGAC","ACCATG","CCATGA","ACCATG","CTGGTA","ACCATG","GACCAT","ACCATG","GGTACT","ACCATG","GTACTG","ACCATG","TACTGG","ACCATG","TGACCA","ACCATG","TGGTAC","ACCATG","ACGATG","ACGATG","ACTGCT","ACGATG","ATGACG","ACGATG","CGATGA","ACGATG","CTACTG","ACGATG","CTGCTA","ACGATG","GACGAT","ACGATG","GATGAC","ACGATG","GCTACT","ACGATG","TACTGC","ACGATG","TGACGA","ACGATG","TGCTAC","ACGATG","ACTCGT","ACTCGT","AGCATG","ACTCGT","ATGAGC","ACTCGT","CATGAG","ACTCGT","CGTACT","ACTCGT","CTCGTA","ACTCGT","GAGCAT","ACTCGT","GCATGA","ACTCGT","GTACTC","ACTCGT","TACTCG","ACTCGT","TCGTAC","ACTCGT","TGAGCA","ACTCGT","ACTCCT","ACTCCT","AGGATG","ACTCCT","ATGAGG","ACTCCT","CCTACT","ACTCCT","CTACTC","ACTCCT","CTCCTA","ACTCCT","GAGGAT","ACTCCT","GATGAG","ACTCCT","GGATGA","ACTCCT","TACTCC","ACTCCT","TCCTAC","ACTCCT","TGAGGA","ACTCCT","ACATGT","ACATGT","ATGTAC","ACATGT","CATGTA","ACATGT","GTACAT","ACATGT","TACATG","ACATGT","TGTACA","ACATGT","ACAGGT","ACAGGT","AGGTAC","ACAGGT","ATGTCC","ACAGGT","CAGGTA","ACAGGT","CATGTC","ACAGGT","CCATGT","ACAGGT","GGTACA","ACAGGT","GTACAG","ACAGGT","GTCCAT","ACAGGT","TACAGG","ACAGGT","TCCATG","ACAGGT","TGTCCA","ACAGGT","ACAGCT","ACAGCT","AGCTAC","ACAGCT","ATGTCG","ACAGCT","CAGCTA","ACAGCT","CGATGT","ACAGCT","CTACAG","ACAGCT","GATGTC","ACAGCT","GCTACA","ACAGCT","GTCGAT","ACAGCT","TACAGC","ACAGCT","TCGATG","ACAGCT","TGTCGA","ACAGCT","ACACGT","ACACGT","ACGTAC","ACACGT","ATGTGC","ACACGT","CACGTA","ACACGT","CATGTG","ACACGT","CGTACA","ACACGT","GCATGT","ACACGT","GTACAC","ACACGT","GTGCAT","ACACGT","TACACG","ACACGT","TGCATG","ACACGT","TGTGCA","ACACGT","ACACCT","ACACCT","ACCTAC","ACACCT","ATGTGG","ACACCT","CACCTA","ACACCT","CCTACA","ACACCT","CTACAC","ACACCT","GATGTG","ACACCT","GGATGT","ACACCT","GTGGAT","ACACCT","TACACC","ACACCT","TGGATG","ACACCT","TGTGGA","ACACCT","ACATGC","ACATGC","ACGTGT","ACATGC","ATGCAC","ACATGC","CACATG","ACATGC","CATGCA","ACATGC","CGTGTA","ACATGC","GCACAT","ACATGC","GTACGT","ACATGC","GTGTAC","ACATGC","TACGTG","ACATGC","TGCACA","ACATGC","TGTACG","ACATGC","ACGTCT","ACGTCT","AGATGC","ACGTCT","ATGCAG","ACGTCT","CAGATG","ACGTCT","CGTCTA","ACGTCT","CTACGT","ACGTCT","GATGCA","ACGTCT","GCAGAT","ACGTCT","GTCTAC","ACGTCT","TACGTC","ACGTCT","TCTACG","ACGTCT","TGCAGA","ACGTCT","ACGAGT","ACGAGT","AGTACG","ACGAGT","ATGCTC","ACGAGT","CATGCT","ACGAGT","CGAGTA","ACGAGT","CTCATG","ACGAGT","GAGTAC","ACGAGT","GCTCAT","ACGAGT","GTACGA","ACGAGT","TACGAG","ACGAGT","TCATGC","ACGAGT","TGCTCA","ACGAGT","ACGACT","ACGACT","ACTACG","ACGACT","ATGCTG","ACGACT","CGACTA","ACGACT","CTACGA","ACGACT","CTGATG","ACGACT","GACTAC","ACGACT","GATGCT","ACGACT","GCTGAT","ACGACT","TACGAC","ACGACT","TGATGC","ACGACT","TGCTGA","ACGACT","ACGGGT","ACGGGT","ATGCCC","ACGGGT","CATGCC","ACGGGT","CCATGC","ACGGGT","CCCATG","ACGGGT","CGGGTA","ACGGGT","GCCCAT","ACGGGT","GGGTAC","ACGGGT","GGTACG","ACGGGT","GTACGG","ACGGGT","TACGGG","ACGGGT","TGCCCA","ACGGGT","ACGGCT","ACGGCT","ATGCCG","ACGGCT","CCGATG","ACGGCT","CGATGC","ACGGCT","CGGCTA","ACGGCT","CTACGG","ACGGCT","GATGCC","ACGGCT","GCCGAT","ACGGCT","GCTACG","ACGGCT","GGCTAC","ACGGCT","TACGGC","ACGGCT","TGCCGA","ACGGCT","ACGCGT","ACGCGT","ATGCGC","ACGCGT","CATGCG","ACGCGT","CGCATG","ACGCGT","CGCGTA","ACGCGT","CGTACG","ACGCGT","GCATGC","ACGCGT","GCGCAT","ACGCGT","GCGTAC","ACGCGT","GTACGC","ACGCGT","TACGCG","ACGCGT","TGCGCA","ACGCGT","ACGCCT","ACGCCT","ATGCGG","ACGCCT","CCTACG","ACGCCT","CGCCTA","ACGCCT","CGGATG","ACGCCT","CTACGC","ACGCCT","GATGCG","ACGCCT","GCCTAC","ACGCCT","GCGGAT","ACGCCT","GGATGC","ACGCCT","TACGCC","ACGCCT","TGCGGA","ACGCCT","ACATGG","ACATGG","ACCTGT","ACATGG","ATGGAC","ACATGG","CATGGA","ACATGG","CCTGTA","ACATGG","CTGTAC","ACATGG","GACATG","ACATGG","GGACAT","ACATGG","GTACCT","ACATGG","TACCTG","ACATGG","TGGACA","ACATGG","TGTACC","ACATGG","ACCTCT","ACCTCT","AGATGG","ACCTCT","ATGGAG","ACCTCT","CCTCTA","ACCTCT","CTACCT","ACCTCT","CTCTAC","ACCTCT","GAGATG","ACCTCT","GATGGA","ACCTCT","GGAGAT","ACCTCT","TACCTC","ACCTCT","TCTACC","ACCTCT","TGGAGA","ACCTCT","ACCAGT","ACCAGT","AGTACC","ACCAGT","ATGGTC","ACCAGT","CAGTAC","ACCAGT","CATGGT","ACCAGT","CCAGTA","ACCAGT","GGTCAT","ACCAGT","GTACCA","ACCAGT","GTCATG","ACCAGT","TACCAG","ACCAGT","TCATGG","ACCAGT","TGGTCA","ACCAGT","ACCACT","ACCACT","ACTACC","ACCACT","ATGGTG","ACCACT","CACTAC","ACCACT","CCACTA","ACCACT","CTACCA","ACCACT","GATGGT","ACCACT","GGTGAT","ACCACT","GTGATG","ACCACT","TACCAC","ACCACT","TGATGG","ACCACT","TGGTGA","ACCACT","ACCGGT","ACCGGT","ATGGCC","ACCGGT","CATGGC","ACCGGT","CCATGG","ACCGGT","CCGGTA","ACCGGT","CGGTAC","ACCGGT","GCCATG","ACCGGT","GGCCAT","ACCGGT","GGTACC","ACCGGT","GTACCG","ACCGGT","TACCGG","ACCGGT","TGGCCA","ACCGGT","ACCGCT","ACCGCT","ATGGCG","ACCGCT","CCGCTA","ACCGCT","CGATGG","ACCGCT","CGCTAC","ACCGCT","CTACCG","ACCGCT","GATGGC","ACCGCT","GCGATG","ACCGCT","GCTACC","ACCGCT","GGCGAT","ACCGCT","TACCGC","ACCGCT","TGGCGA","ACCGCT","ACCCGT","ACCCGT","ATGGGC","ACCCGT","CATGGG","ACCCGT","CCCGTA","ACCCGT","CCGTAC","ACCCGT","CGTACC","ACCCGT","GCATGG","ACCCGT","GGCATG","ACCCGT","GGGCAT","ACCCGT","GTACCC","ACCCGT","TACCCG","ACCCGT","TGGGCA","ACCCGT","ACCCCT","ACCCCT","ATGGGG","ACCCCT","CCCCTA","ACCCCT","CCCTAC","ACCCCT","CCTACC","ACCCCT","CTACCC","ACCCCT","GATGGG","ACCCCT","GGATGG","ACCCCT","GGGATG","ACCCCT","GGGGAT","ACCCCT","TACCCC","ACCCCT","TGGGGA","ACCCCT","ACACAG","ACACAG","ACAGAC","ACACAG","AGACAC","ACACAG","CACAGA","ACACAG","CAGACA","ACACAG","CTGTGT","ACACAG","GACACA","ACACAG","GTCTGT","ACACAG","GTGTCT","ACACAG","TCTGTG","ACACAG","TGTCTG","ACACAG","TGTGTC","ACACAG","ACACTC","ACACTC","ACTCAC","ACACTC","AGTGTG","ACACTC","CACACT","ACACTC","CACTCA","ACACTC","CTCACA","ACACTC","GAGTGT","ACACTC","GTGAGT","ACACTC","GTGTGA","ACACTC","TCACAC","ACACTC","TGAGTG","ACACTC","TGTGAG","ACACTC","ACACTG","ACACTG","ACTGAC","ACACTG","ACTGTG","ACACTG","CACTGA","ACACTG","CTGACA","ACACTG","CTGTGA","ACACTG","GACACT","ACACTG","GACTGT","ACACTG","GTGACT","ACACTG","TGACAC","ACACTG","TGACTG","ACACTG","TGTGAC","ACACTG","ACACCC","ACACCC","ACCCAC","ACACCC","CACACC","ACACCC","CACCCA","ACACCC","CCACAC","ACACCC","CCCACA","ACACCC","GGGTGT","ACACCC","GGTGTG","ACACCC","GTGGGT","ACACCC","GTGTGG","ACACCC","TGGGTG","ACACCC","TGTGGG","ACACCC","ACACCG","ACACCG","ACCGAC","ACACCG","CACCGA","ACACCG","CCGACA","ACACCG","CGACAC","ACACCG","CTGTGG","ACACCG","GACACC","ACACCG","GCTGTG","ACACCG","GGCTGT","ACACCG","GTGGCT","ACACCG","TGGCTG","ACACCG","TGTGGC","ACACCG","ACACGC","ACACGC","ACGCAC","ACACGC","CACACG","ACACGC","CACGCA","ACACGC","CGCACA","ACACGC","CGTGTG","ACACGC","GCACAC","ACACGC","GCGTGT","ACACGC","GTGCGT","ACACGC","GTGTGC","ACACGC","TGCGTG","ACACGC","TGTGCG","ACACGC","ACACGG","ACACGG","ACGGAC","ACACGG","CACGGA","ACACGG","CCTGTG","ACACGG","CGGACA","ACACGG","CTGTGC","ACACGG","GACACG","ACACGG","GCCTGT","ACACGG","GGACAC","ACACGG","GTGCCT","ACACGG","TGCCTG","ACACGG","TGTGCC","ACACGG","ACAGAG","ACAGAG","AGACAG","ACAGAG","AGAGAC","ACAGAG","CAGAGA","ACAGAG","CTCTGT","ACAGAG","CTGTCT","ACAGAG","GACAGA","ACAGAG","GAGACA","ACAGAG","GTCTCT","ACAGAG","TCTCTG","ACAGAG","TCTGTC","ACAGAG","TGTCTC","ACAGAG","ACAGTC","ACAGTC","AGTCAC","ACAGTC","AGTGTC","ACAGTC","CACAGT","ACAGTC","CAGTCA","ACAGTC","CAGTGT","ACAGTC","GTCACA","ACAGTC","GTCAGT","ACAGTC","GTGTCA","ACAGTC","TCACAG","ACAGTC","TCAGTG","ACAGTC","TGTCAG","ACAGTC","ACAGTG","ACAGTG","ACTGTC","ACAGTG","AGTGAC","ACAGTG","CACTGT","ACAGTG","CAGTGA","ACAGTG","CTGTCA","ACAGTG","GACAGT","ACAGTG","GTCACT","ACAGTG","GTGACA","ACAGTG","TCACTG","ACAGTG","TGACAG","ACAGTG","TGTCAC","ACAGTG","ACAGCC","ACAGCC","AGCCAC","ACAGCC","CACAGC","ACAGCC","CAGCCA","ACAGCC","CCACAG","ACAGCC","CGGTGT","ACAGCC","GCCACA","ACAGCC","GGTGTC","ACAGCC","GTCGGT","ACAGCC","GTGTCG","ACAGCC","TCGGTG","ACAGCC","TGTCGG","ACAGCC","ACAGCG","ACAGCG","AGCGAC","ACAGCG","CAGCGA","ACAGCG","CGACAG","ACAGCG","CGCTGT","ACAGCG","CTGTCG","ACAGCG","GACAGC","ACAGCG","GCGACA","ACAGCG","GCTGTC","ACAGCG","GTCGCT","ACAGCG","TCGCTG","ACAGCG","TGTCGC","ACAGCG","ACAGGC","ACAGGC","AGGCAC","ACAGGC","CACAGG","ACAGGC","CAGGCA","ACAGGC","CCGTGT","ACAGGC","CGTGTC","ACAGGC","GCACAG","ACAGGC","GGCACA","ACAGGC","GTCCGT","ACAGGC","GTGTCC","ACAGGC","TCCGTG","ACAGGC","TGTCCG","ACAGGC","ACAGGG","ACAGGG","AGGGAC","ACAGGG","CAGGGA","ACAGGG","CCCTGT","ACAGGG","CCTGTC","ACAGGG","CTGTCC","ACAGGG","GACAGG","ACAGGG","GGACAG","ACAGGG","GGGACA","ACAGGG","GTCCCT","ACAGGG","TCCCTG","ACAGGG","TGTCCC","ACAGGG","ACTCAG","ACTCAG","AGACTC","ACTCAG","AGTCTG","ACTCAG","CAGACT","ACTCAG","CTCAGA","ACTCAG","CTGAGT","ACTCAG","GACTCA","ACTCAG","GAGTCT","ACTCAG","GTCTGA","ACTCAG","TCAGAC","ACTCAG","TCTGAG","ACTCAG","TGAGTC","ACTCAG","ACTCTC","ACTCTC","AGAGTG","ACTCTC","AGTGAG","ACTCTC","CACTCT","ACTCTC","CTCACT","ACTCTC","CTCTCA","ACTCTC","GAGAGT","ACTCTC","GAGTGA","ACTCTC","GTGAGA","ACTCTC","TCACTC","ACTCTC","TCTCAC","ACTCTC","TGAGAG","ACTCTC","ACTCTG","ACTCTG","ACTGAG","ACTCTG","AGACTG","ACTCTG","CTCTGA","ACTCTG","CTGACT","ACTCTG","CTGAGA","ACTCTG","GACTCT","ACTCTG","GACTGA","ACTCTG","GAGACT","ACTCTG","TCTGAC","ACTCTG","TGACTC","ACTCTG","TGAGAC","ACTCTG","ACTCCC","ACTCCC","AGGGTG","ACTCCC","CACTCC","ACTCCC","CCACTC","ACTCCC","CCCACT","ACTCCC","CTCCCA","ACTCCC","GAGGGT","ACTCCC","GGGTGA","ACTCCC","GGTGAG","ACTCCC","GTGAGG","ACTCCC","TCCCAC","ACTCCC","TGAGGG","ACTCCC","ACTCCG","ACTCCG","AGGCTG","ACTCCG","CCGACT","ACTCCG","CGACTC","ACTCCG","CTCCGA","ACTCCG","CTGAGG","ACTCCG","GACTCC","ACTCCG","GAGGCT","ACTCCG","GCTGAG","ACTCCG","GGCTGA","ACTCCG","TCCGAC","ACTCCG","TGAGGC","ACTCCG","ACTCGC","ACTCGC","AGCGTG","ACTCGC","CACTCG","ACTCGC","CGCACT","ACTCGC","CGTGAG","ACTCGC","CTCGCA","ACTCGC","GAGCGT","ACTCGC","GCACTC","ACTCGC","GCGTGA","ACTCGC","GTGAGC","ACTCGC","TCGCAC","ACTCGC","TGAGCG","ACTCGC","ACTCGG","ACTCGG","AGCCTG","ACTCGG","CCTGAG","ACTCGG","CGGACT","ACTCGG","CTCGGA","ACTCGG","CTGAGC","ACTCGG","GACTCG","ACTCGG","GAGCCT","ACTCGG","GCCTGA","ACTCGG","GGACTC","ACTCGG","TCGGAC","ACTCGG","TGAGCC","ACTCGG","ACGGTG","ACGGTG","ACTGCC","ACGGTG","CACTGC","ACGGTG","CCACTG","ACGGTG","CGGTGA","ACGGTG","CTGCCA","ACGGTG","GACGGT","ACGGTG","GCCACT","ACGGTG","GGTGAC","ACGGTG","GTGACG","ACGGTG","TGACGG","ACGGTG","TGCCAC","ACGGTG","ACGCTG","ACGCTG","ACTGCG","ACGCTG","CGACTG","ACGCTG","CGCTGA","ACGCTG","CTGACG","ACGCTG","CTGCGA","ACGCTG","GACGCT","ACGCTG","GACTGC","ACGCTG","GCGACT","ACGCTG","GCTGAC","ACGCTG","TGACGC","ACGCTG","TGCGAC","ACGCTG","ACCGTG","ACCGTG","ACTGGC","ACCGTG","CACTGG","ACCGTG","CCGTGA","ACCGTG","CGTGAC","ACCGTG","CTGGCA","ACCGTG","GACCGT","ACCGTG","GCACTG","ACCGTG","GGCACT","ACCGTG","GTGACC","ACCGTG","TGACCG","ACCGTG","TGGCAC","ACCGTG","ACCCTG","ACCCTG","ACTGGG","ACCCTG","CCCTGA","ACCCTG","CCTGAC","ACCCTG","CTGACC","ACCCTG","CTGGGA","ACCCTG","GACCCT","ACCCTG","GACTGG","ACCCTG","GGACTG","ACCCTG","GGGACT","ACCCTG","TGACCC","ACCCTG","TGGGAC","ACCCTG","ACCACG","ACCACG","ACGACC","ACCACG","CACGAC","ACCACG","CCACGA","ACCACG","CGACCA","ACCACG","CTGGTG","ACCACG","GACCAC","ACCACG","GCTGGT","ACCACG","GGTGCT","ACCACG","GTGCTG","ACCACG","TGCTGG","ACCACG","TGGTGC","ACCACG","ACCAGC","ACCAGC","AGCACC","ACCAGC","CACCAG","ACCAGC","CAGCAC","ACCAGC","CCAGCA","ACCAGC","CGTGGT","ACCAGC","GCACCA","ACCAGC","GGTCGT","ACCAGC","GTCGTG","ACCAGC","GTGGTC","ACCAGC","TCGTGG","ACCAGC","TGGTCG","ACCAGC","ACCAGG","ACCAGG","AGGACC","ACCAGG","CAGGAC","ACCAGG","CCAGGA","ACCAGG","CCTGGT","ACCAGG","CTGGTC","ACCAGG","GACCAG","ACCAGG","GGACCA","ACCAGG","GGTCCT","ACCAGG","GTCCTG","ACCAGG","TCCTGG","ACCAGG","TGGTCC","ACCAGG","ACCTCC","ACCTCC","AGGTGG","ACCTCC","CACCTC","ACCTCC","CCACCT","ACCTCC","CCTCCA","ACCTCC","CTCCAC","ACCTCC","GAGGTG","ACCTCC","GGAGGT","ACCTCC","GGTGGA","ACCTCC","GTGGAG","ACCTCC","TCCACC","ACCTCC","TGGAGG","ACCTCC","ACCTCG","ACCTCG","AGCTGG","ACCTCG","CCTCGA","ACCTCG","CGACCT","ACCTCG","CTCGAC","ACCTCG","CTGGAG","ACCTCG","GACCTC","ACCTCG","GAGCTG","ACCTCG","GCTGGA","ACCTCG","GGAGCT","ACCTCG","TCGACC","ACCTCG","TGGAGC","ACCTCG","ACCTGC","ACCTGC","ACGTGG","ACCTGC","CACCTG","ACCTGC","CCTGCA","ACCTGC","CGTGGA","ACCTGC","CTGCAC","ACCTGC","GACGTG","ACCTGC","GCACCT","ACCTGC","GGACGT","ACCTGC","GTGGAC","ACCTGC","TGCACC","ACCTGC","TGGACG","ACCTGC","ACCTGG","ACCTGG","CCTGGA","ACCTGG","CTGGAC","ACCTGG","GACCTG","ACCTGG","GGACCT","ACCTGG","TGGACC","ACCTGG","ACCCAG","ACCCAG","AGACCC","ACCCAG","CAGACC","ACCCAG","CCAGAC","ACCCAG","CCCAGA","ACCCAG","CTGGGT","ACCCAG","GACCCA","ACCCAG","GGGTCT","ACCCAG","GGTCTG","ACCCAG","GTCTGG","ACCCAG","TCTGGG","ACCCAG","TGGGTC","ACCCAG","ACCCTC","ACCCTC","AGTGGG","ACCCTC","CACCCT","ACCCTC","CCCTCA","ACCCTC","CCTCAC","ACCCTC","CTCACC","ACCCTC","GAGTGG","ACCCTC","GGAGTG","ACCCTC","GGGAGT","ACCCTC","GTGGGA","ACCCTC","TCACCC","ACCCTC","TGGGAG","ACCCTC","ACCCCC","ACCCCC","CACCCC","ACCCCC","CCACCC","ACCCCC","CCCACC","ACCCCC","CCCCAC","ACCCCC","CCCCCA","ACCCCC","GGGGGT","ACCCCC","GGGGTG","ACCCCC","GGGTGG","ACCCCC","GGTGGG","ACCCCC","GTGGGG","ACCCCC","TGGGGG","ACCCCC","ACCCCG","ACCCCG","CCCCGA","ACCCCG","CCCGAC","ACCCCG","CCGACC","ACCCCG","CGACCC","ACCCCG","CTGGGG","ACCCCG","GACCCC","ACCCCG","GCTGGG","ACCCCG","GGCTGG","ACCCCG","GGGCTG","ACCCCG","GGGGCT","ACCCCG","TGGGGC","ACCCCG","ACCCGC","ACCCGC","CACCCG","ACCCGC","CCCGCA","ACCCGC","CCGCAC","ACCCGC","CGCACC","ACCCGC","CGTGGG","ACCCGC","GCACCC","ACCCGC","GCGTGG","ACCCGC","GGCGTG","ACCCGC","GGGCGT","ACCCGC","GTGGGC","ACCCGC","TGGGCG","ACCCGC","ACCCGG","ACCCGG","CCCGGA","ACCCGG","CCGGAC","ACCCGG","CCTGGG","ACCCGG","CGGACC","ACCCGG","CTGGGC","ACCCGG","GACCCG","ACCCGG","GCCTGG","ACCCGG","GGACCC","ACCCGG","GGCCTG","ACCCGG","GGGCCT","ACCCGG","TGGGCC","ACCCGG","ACCGAG","ACCGAG","AGACCG","ACCGAG","CCGAGA","ACCGAG","CGAGAC","ACCGAG","CTCTGG","ACCGAG","CTGGCT","ACCGAG","GACCGA","ACCGAG","GAGACC","ACCGAG","GCTCTG","ACCGAG","GGCTCT","ACCGAG","TCTGGC","ACCGAG","TGGCTC","ACCGAG","ACCGTC","ACCGTC","AGTGGC","ACCGTC","CACCGT","ACCGTC","CAGTGG","ACCGTC","CCGTCA","ACCGTC","CGTCAC","ACCGTC","GCAGTG","ACCGTC","GGCAGT","ACCGTC","GTCACC","ACCGTC","GTGGCA","ACCGTC","TCACCG","ACCGTC","TGGCAG","ACCGTC","ACCGCC","ACCGCC","CACCGC","ACCGCC","CCACCG","ACCGCC","CCGCCA","ACCGCC","CGCCAC","ACCGCC","CGGTGG","ACCGCC","GCCACC","ACCGCC","GCGGTG","ACCGCC","GGCGGT","ACCGCC","GGTGGC","ACCGCC","GTGGCG","ACCGCC","TGGCGG","ACCGCC","ACCGCG","ACCGCG","CCGCGA","ACCGCG","CGACCG","ACCGCG","CGCGAC","ACCGCG","CGCTGG","ACCGCG","CTGGCG","ACCGCG","GACCGC","ACCGCG","GCGACC","ACCGCG","GCGCTG","ACCGCG","GCTGGC","ACCGCG","GGCGCT","ACCGCG","TGGCGC","ACCGCG","ACCGGC","ACCGGC","CACCGG","ACCGGC","CCGGCA","ACCGGC","CCGTGG","ACCGGC","CGGCAC","ACCGGC","CGTGGC","ACCGGC","GCACCG","ACCGGC","GCCGTG","ACCGGC","GGCACC","ACCGGC","GGCCGT","ACCGGC","GTGGCC","ACCGGC","TGGCCG","ACCGGC","ACCGGG","ACCGGG","CCCTGG","ACCGGG","CCGGGA","ACCGGG","CCTGGC","ACCGGG","CGGGAC","ACCGGG","CTGGCC","ACCGGG","GACCGG","ACCGGG","GCCCTG","ACCGGG","GGACCG","ACCGGG","GGCCCT","ACCGGG","GGGACC","ACCGGG","TGGCCC","ACCGGG","ACGAGC","ACGAGC","AGCACG","ACGAGC","CACGAG","ACGAGC","CGAGCA","ACGAGC","CGTGCT","ACGAGC","CTCGTG","ACGAGC","GAGCAC","ACGAGC","GCACGA","ACGAGC","GCTCGT","ACGAGC","GTGCTC","ACGAGC","TCGTGC","ACGAGC","TGCTCG","ACGAGC","ACGAGG","ACGAGG","AGGACG","ACGAGG","CCTGCT","ACGAGG","CGAGGA","ACGAGG","CTCCTG","ACGAGG","CTGCTC","ACGAGG","GACGAG","ACGAGG","GAGGAC","ACGAGG","GCTCCT","ACGAGG","GGACGA","ACGAGG","TCCTGC","ACGAGG","TGCTCC","ACGAGG","ACGTCC","ACGTCC","AGGTGC","ACGTCC","CACGTC","ACGTCC","CAGGTG","ACGTCC","CCACGT","ACGTCC","CGTCCA","ACGTCC","GCAGGT","ACGTCC","GGTGCA","ACGTCC","GTCCAC","ACGTCC","GTGCAG","ACGTCC","TCCACG","ACGTCC","TGCAGG","ACGTCC","ACGTCG","ACGTCG","AGCTGC","ACGTCG","CAGCTG","ACGTCG","CGACGT","ACGTCG","CGTCGA","ACGTCG","CTGCAG","ACGTCG","GACGTC","ACGTCG","GCAGCT","ACGTCG","GCTGCA","ACGTCG","GTCGAC","ACGTCG","TCGACG","ACGTCG","TGCAGC","ACGTCG","ACGTGC","ACGTGC","CACGTG","ACGTGC","CGTGCA","ACGTGC","GCACGT","ACGTGC","GTGCAC","ACGTGC","TGCACG","ACGTGC","ACGCAG","ACGCAG","AGACGC","ACGCAG","CAGACG","ACGCAG","CGCAGA","ACGCAG","CGTCTG","ACGCAG","CTGCGT","ACGCAG","GACGCA","ACGCAG","GCAGAC","ACGCAG","GCGTCT","ACGCAG","GTCTGC","ACGCAG","TCTGCG","ACGCAG","TGCGTC","ACGCAG","ACGCTC","ACGCTC","AGTGCG","ACGCTC","CACGCT","ACGCTC","CGAGTG","ACGCTC","CGCTCA","ACGCTC","CTCACG","ACGCTC","GAGTGC","ACGCTC","GCGAGT","ACGCTC","GCTCAC","ACGCTC","GTGCGA","ACGCTC","TCACGC","ACGCTC","TGCGAG","ACGCTC","ACGCCC","ACGCCC","CACGCC","ACGCCC","CCACGC","ACGCCC","CCCACG","ACGCCC","CGCCCA","ACGCCC","CGGGTG","ACGCCC","GCCCAC","ACGCCC","GCGGGT","ACGCCC","GGGTGC","ACGCCC","GGTGCG","ACGCCC","GTGCGG","ACGCCC","TGCGGG","ACGCCC","ACGCCG","ACGCCG","CCGACG","ACGCCG","CGACGC","ACGCCG","CGCCGA","ACGCCG","CGGCTG","ACGCCG","CTGCGG","ACGCCG","GACGCC","ACGCCG","GCCGAC","ACGCCG","GCGGCT","ACGCCG","GCTGCG","ACGCCG","GGCTGC","ACGCCG","TGCGGC","ACGCCG","ACGCGC","ACGCGC","CACGCG","ACGCGC","CGCACG","ACGCGC","CGCGCA","ACGCGC","CGCGTG","ACGCGC","CGTGCG","ACGCGC","GCACGC","ACGCGC","GCGCAC","ACGCGC","GCGCGT","ACGCGC","GCGTGC","ACGCGC","GTGCGC","ACGCGC","TGCGCG","ACGCGC","ACGCGG","ACGCGG","CCTGCG","ACGCGG","CGCCTG","ACGCGG","CGCGGA","ACGCGG","CGGACG","ACGCGG","CTGCGC","ACGCGG","GACGCG","ACGCGG","GCCTGC","ACGCGG","GCGCCT","ACGCGG","GCGGAC","ACGCGG","GGACGC","ACGCGG","TGCGCC","ACGCGG","ACGGAG","ACGGAG","AGACGG","ACGGAG","CCTCTG","ACGGAG","CGGAGA","ACGGAG","CTCTGC","ACGGAG","CTGCCT","ACGGAG","GACGGA","ACGGAG","GAGACG","ACGGAG","GCCTCT","ACGGAG","GGAGAC","ACGGAG","TCTGCC","ACGGAG","TGCCTC","ACGGAG","ACGGTC","ACGGTC","AGTGCC","ACGGTC","CACGGT","ACGGTC","CAGTGC","ACGGTC","CCAGTG","ACGGTC","CGGTCA","ACGGTC","GCCAGT","ACGGTC","GGTCAC","ACGGTC","GTCACG","ACGGTC","GTGCCA","ACGGTC","TCACGG","ACGGTC","TGCCAG","ACGGTC","ACGGCC","ACGGCC","CACGGC","ACGGCC","CCACGG","ACGGCC","CCGGTG","ACGGCC","CGGCCA","ACGGCC","CGGTGC","ACGGCC","GCCACG","ACGGCC","GCCGGT","ACGGCC","GGCCAC","ACGGCC","GGTGCC","ACGGCC","GTGCCG","ACGGCC","TGCCGG","ACGGCC","ACGGCG","ACGGCG","CCGCTG","ACGGCG","CGACGG","ACGGCG","CGCTGC","ACGGCG","CGGCGA","ACGGCG","CTGCCG","ACGGCG","GACGGC","ACGGCG","GCCGCT","ACGGCG","GCGACG","ACGGCG","GCTGCC","ACGGCG","GGCGAC","ACGGCG","TGCCGC","ACGGCG","ACGGGC","ACGGGC","CACGGG","ACGGGC","CCCGTG","ACGGGC","CCGTGC","ACGGGC","CGGGCA","ACGGGC","CGTGCC","ACGGGC","GCACGG","ACGGGC","GCCCGT","ACGGGC","GGCACG","ACGGGC","GGGCAC","ACGGGC","GTGCCC","ACGGGC","TGCCCG","ACGGGC","ACGGGG","ACGGGG","CCCCTG","ACGGGG","CCCTGC","ACGGGG","CCTGCC","ACGGGG","CGGGGA","ACGGGG","CTGCCC","ACGGGG","GACGGG","ACGGGG","GCCCCT","ACGGGG","GGACGG","ACGGGG","GGGACG","ACGGGG","GGGGAC","ACGGGG","TGCCCC","ACGGGG","AGAGTC","AGAGTC","AGTCAG","AGAGTC","AGTCTC","AGAGTC","CAGAGT","AGAGTC","CAGTCT","AGAGTC","CTCAGT","AGAGTC","GAGTCA","AGAGTC","GTCAGA","AGAGTC","GTCTCA","AGAGTC","TCAGAG","AGAGTC","TCAGTC","AGAGTC","TCTCAG","AGAGTC","AGAGCC","AGAGCC","AGCCAG","AGAGCC","CAGAGC","AGAGCC","CCAGAG","AGAGCC","CGGTCT","AGAGCC","CTCGGT","AGAGCC","GAGCCA","AGAGCC","GCCAGA","AGAGCC","GGTCTC","AGAGCC","GTCTCG","AGAGCC","TCGGTC","AGAGCC","TCTCGG","AGAGCC","AGAGCG","AGAGCG","AGCGAG","AGAGCG","CGAGAG","AGAGCG","CGCTCT","AGAGCG","CTCGCT","AGAGCG","CTCTCG","AGAGCG","GAGAGC","AGAGCG","GAGCGA","AGAGCG","GCGAGA","AGAGCG","GCTCTC","AGAGCG","TCGCTC","AGAGCG","TCTCGC","AGAGCG","AGAGGC","AGAGGC","AGGCAG","AGAGGC","CAGAGG","AGAGGC","CCGTCT","AGAGGC","CGTCTC","AGAGGC","CTCCGT","AGAGGC","GAGGCA","AGAGGC","GCAGAG","AGAGGC","GGCAGA","AGAGGC","GTCTCC","AGAGGC","TCCGTC","AGAGGC","TCTCCG","AGAGGC","AGAGGG","AGAGGG","AGGGAG","AGAGGG","CCCTCT","AGAGGG","CCTCTC","AGAGGG","CTCCCT","AGAGGG","CTCTCC","AGAGGG","GAGAGG","AGAGGG","GAGGGA","AGAGGG","GGAGAG","AGAGGG","GGGAGA","AGAGGG","TCCCTC","AGAGGG","TCTCCC","AGAGGG","AGGGTC","AGGGTC","AGTCCC","AGGGTC","CAGGGT","AGGGTC","CAGTCC","AGGGTC","CCAGTC","AGGGTC","CCCAGT","AGGGTC","GGGTCA","AGGGTC","GGTCAG","AGGGTC","GTCAGG","AGGGTC","GTCCCA","AGGGTC","TCAGGG","AGGGTC","TCCCAG","AGGGTC","AGGCTC","AGGCTC","AGTCCG","AGGCTC","CAGGCT","AGGCTC","CCGAGT","AGGCTC","CGAGTC","AGGCTC","CTCAGG","AGGCTC","GAGTCC","AGGCTC","GCTCAG","AGGCTC","GGCTCA","AGGCTC","GTCCGA","AGGCTC","TCAGGC","AGGCTC","TCCGAG","AGGCTC","AGCGTC","AGCGTC","AGTCGC","AGCGTC","CAGCGT","AGCGTC","CAGTCG","AGCGTC","CGCAGT","AGCGTC","CGTCAG","AGCGTC","GCAGTC","AGCGTC","GCGTCA","AGCGTC","GTCAGC","AGCGTC","GTCGCA","AGCGTC","TCAGCG","AGCGTC","TCGCAG","AGCGTC","AGCCTC","AGCCTC","AGTCGG","AGCCTC","CAGCCT","AGCCTC","CCTCAG","AGCCTC","CGGAGT","AGCCTC","CTCAGC","AGCCTC","GAGTCG","AGCCTC","GCCTCA","AGCCTC","GGAGTC","AGCCTC","GTCGGA","AGCCTC","TCAGCC","AGCCTC","TCGGAG","AGCCTC","AGCAGG","AGCAGG","AGGAGC","AGCAGG","CAGGAG","AGCAGG","CCTCGT","AGCAGG","CGTCCT","AGCAGG","CTCGTC","AGCAGG","GAGCAG","AGCAGG","GCAGGA","AGCAGG","GGAGCA","AGCAGG","GTCCTC","AGCAGG","TCCTCG","AGCAGG","TCGTCC","AGCAGG","AGCTCC","AGCTCC","AGGTCG","AGCTCC","CAGCTC","AGCTCC","CCAGCT","AGCTCC","CGAGGT","AGCTCC","CTCCAG","AGCTCC","GAGGTC","AGCTCC","GCTCCA","AGCTCC","GGTCGA","AGCTCC","GTCGAG","AGCTCC","TCCAGC","AGCTCC","TCGAGG","AGCTCC","AGCTCG","AGCTCG","CGAGCT","AGCTCG","CTCGAG","AGCTCG","GAGCTC","AGCTCG","GCTCGA","AGCTCG","TCGAGC","AGCTCG","AGCCCC","AGCCCC","CAGCCC","AGCCCC","CCAGCC","AGCCCC","CCCAGC","AGCCCC","CCCCAG","AGCCCC","CGGGGT","AGCCCC","GCCCCA","AGCCCC","GGGGTC","AGCCCC","GGGTCG","AGCCCC","GGTCGG","AGCCCC","GTCGGG","AGCCCC","TCGGGG","AGCCCC","AGCCCG","AGCCCG","CCCGAG","AGCCCG","CCGAGC","AGCCCG","CGAGCC","AGCCCG","CGGGCT","AGCCCG","CTCGGG","AGCCCG","GAGCCC","AGCCCG","GCCCGA","AGCCCG","GCTCGG","AGCCCG","GGCTCG","AGCCCG","GGGCTC","AGCCCG","TCGGGC","AGCCCG","AGCCGC","AGCCGC","CAGCCG","AGCCGC","CCGCAG","AGCCGC","CGCAGC","AGCCGC","CGGCGT","AGCCGC","CGTCGG","AGCCGC","GCAGCC","AGCCGC","GCCGCA","AGCCGC","GCGTCG","AGCCGC","GGCGTC","AGCCGC","GTCGGC","AGCCGC","TCGGCG","AGCCGC","AGCCGG","AGCCGG","CCGGAG","AGCCGG","CCTCGG","AGCCGG","CGGAGC","AGCCGG","CGGCCT","AGCCGG","CTCGGC","AGCCGG","GAGCCG","AGCCGG","GCCGGA","AGCCGG","GCCTCG","AGCCGG","GGAGCC","AGCCGG","GGCCTC","AGCCGG","TCGGCC","AGCCGG","AGCGCC","AGCGCC","CAGCGC","AGCGCC","CCAGCG","AGCGCC","CGCCAG","AGCGCC","CGCGGT","AGCGCC","CGGTCG","AGCGCC","GCCAGC","AGCGCC","GCGCCA","AGCGCC","GCGGTC","AGCGCC","GGTCGC","AGCGCC","GTCGCG","AGCGCC","TCGCGG","AGCGCC","AGCGCG","AGCGCG","CGAGCG","AGCGCG","CGCGAG","AGCGCG","CGCGCT","AGCGCG","CGCTCG","AGCGCG","CTCGCG","AGCGCG","GAGCGC","AGCGCG","GCGAGC","AGCGCG","GCGCGA","AGCGCG","GCGCTC","AGCGCG","GCTCGC","AGCGCG","TCGCGC","AGCGCG","AGCGGC","AGCGGC","CAGCGG","AGCGGC","CCGTCG","AGCGGC","CGCCGT","AGCGGC","CGGCAG","AGCGGC","CGTCGC","AGCGGC","GCAGCG","AGCGGC","GCCGTC","AGCGGC","GCGGCA","AGCGGC","GGCAGC","AGCGGC","GTCGCC","AGCGGC","TCGCCG","AGCGGC","AGCGGG","AGCGGG","CCCTCG","AGCGGG","CCTCGC","AGCGGG","CGCCCT","AGCGGG","CGGGAG","AGCGGG","CTCGCC","AGCGGG","GAGCGG","AGCGGG","GCCCTC","AGCGGG","GCGGGA","AGCGGG","GGAGCG","AGCGGG","GGGAGC","AGCGGG","TCGCCC","AGCGGG","AGGTCC","AGGTCC","CAGGTC","AGGTCC","CCAGGT","AGGTCC","GGTCCA","AGGTCC","GTCCAG","AGGTCC","TCCAGG","AGGTCC","AGGCCC","AGGCCC","CAGGCC","AGGCCC","CCAGGC","AGGCCC","CCCAGG","AGGCCC","CCGGGT","AGGCCC","CGGGTC","AGGCCC","GCCCAG","AGGCCC","GGCCCA","AGGCCC","GGGTCC","AGGCCC","GGTCCG","AGGCCC","GTCCGG","AGGCCC","TCCGGG","AGGCCC","AGGCCG","AGGCCG","CCGAGG","AGGCCG","CCGGCT","AGGCCG","CGAGGC","AGGCCG","CGGCTC","AGGCCG","CTCCGG","AGGCCG","GAGGCC","AGGCCG","GCCGAG","AGGCCG","GCTCCG","AGGCCG","GGCCGA","AGGCCG","GGCTCC","AGGCCG","TCCGGC","AGGCCG","AGGCGC","AGGCGC","CAGGCG","AGGCGC","CCGCGT","AGGCGC","CGCAGG","AGGCGC","CGCGTC","AGGCGC","CGTCCG","AGGCGC","GCAGGC","AGGCGC","GCGCAG","AGGCGC","GCGTCC","AGGCGC","GGCGCA","AGGCGC","GTCCGC","AGGCGC","TCCGCG","AGGCGC","AGGCGG","AGGCGG","CCGCCT","AGGCGG","CCTCCG","AGGCGG","CGCCTC","AGGCGG","CGGAGG","AGGCGG","CTCCGC","AGGCGG","GAGGCG","AGGCGG","GCCTCC","AGGCGG","GCGGAG","AGGCGG","GGAGGC","AGGCGG","GGCGGA","AGGCGG","TCCGCC","AGGCGG","AGGGCC","AGGGCC","CAGGGC","AGGGCC","CCAGGG","AGGGCC","CCCGGT","AGGGCC","CCGGTC","AGGGCC","CGGTCC","AGGGCC","GCCAGG","AGGGCC","GGCCAG","AGGGCC","GGGCCA","AGGGCC","GGTCCC","AGGGCC","GTCCCG","AGGGCC","TCCCGG","AGGGCC","AGGGCG","AGGGCG","CCCGCT","AGGGCG","CCGCTC","AGGGCG","CGAGGG","AGGGCG","CGCTCC","AGGGCG","CTCCCG","AGGGCG","GAGGGC","AGGGCG","GCGAGG","AGGGCG","GCTCCC","AGGGCG","GGCGAG","AGGGCG","GGGCGA","AGGGCG","TCCCGC","AGGGCG","AGGGGC","AGGGGC","CAGGGG","AGGGGC","CCCCGT","AGGGGC","CCCGTC","AGGGGC","CCGTCC","AGGGGC","CGTCCC","AGGGGC","GCAGGG","AGGGGC","GGCAGG","AGGGGC","GGGCAG","AGGGGC","GGGGCA","AGGGGC","GTCCCC","AGGGGC","TCCCCG","AGGGGC","AGGGGG","AGGGGG","CCCCCT","AGGGGG","CCCCTC","AGGGGG","CCCTCC","AGGGGG","CCTCCC","AGGGGG","CTCCCC","AGGGGG","GAGGGG","AGGGGG","GGAGGG","AGGGGG","GGGAGG","AGGGGG","GGGGAG","AGGGGG","GGGGGA","AGGGGG","TCCCCC","AGGGGG","CCCCCG","CCCCCG","CCCCGC","CCCCCG","CCCGCC","CCCCCG","CCGCCC","CCCCCG","CGCCCC","CCCCCG","CGGGGG","CCCCCG","GCCCCC","CCCCCG","GCGGGG","CCCCCG","GGCGGG","CCCCCG","GGGCGG","CCCCCG","GGGGCG","CCCCCG","GGGGGC","CCCCCG","CCCCGG","CCCCGG","CCCGGC","CCCCGG","CCGGCC","CCCCGG","CCGGGG","CCCCGG","CGGCCC","CCCCGG","CGGGGC","CCCCGG","GCCCCG","CCCCGG","GCCGGG","CCCCGG","GGCCCC","CCCCGG","GGCCGG","CCCCGG","GGGCCG","CCCCGG","GGGGCC","CCCCGG","CCCGCG","CCCGCG","CCGCGC","CCCGCG","CGCCCG","CCCGCG","CGCGCC","CCCGCG","CGCGGG","CCCGCG","CGGGCG","CCCGCG","GCCCGC","CCCGCG","GCGCCC","CCCGCG","GCGCGG","CCCGCG","GCGGGC","CCCGCG","GGCGCG","CCCGCG","GGGCGC","CCCGCG","CCCGGG","CCCGGG","CCGGGC","CCCGGG","CGGGCC","CCCGGG","GCCCGG","CCCGGG","GGCCCG","CCCGGG","GGGCCC","CCCGGG","CCGCGG","CCGCGG","CCGGCG","CCGCGG","CGCCGG","CCGCGG","CGCGGC","CCGCGG","CGGCCG","CCGCGG","CGGCGC","CCGCGG","GCCGCG","CCGCGG","GCCGGC","CCGCGG","GCGCCG","CCGCGG","GCGGCC","CCGCGG","GGCCGC","CCGCGG","GGCGCC","CCGCGG","N","N");
	return \%SSR_motif_class;
}

#Subroutine 12: Run MAFFT program
sub runMafftProgram { ##(fastaFolder[>id\nseq]<path>, mafftFolder[>id\nseq]<path>, mafftIDFile[id]<path>)
	$/ = "\n";
	open (MA, "$_[2]") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	my @mafft_id = <MA>;
	for (my $i=0;$i<=$#mafft_id;$i++) {
		chomp ($mafft_id[$i]);
		my $input = "$_[0]/$mafft_id[$i]";
		my $output = "$_[1]/$mafft_id[$i].mafft";
		system "mafft --thread $thread --quiet $input > $output";
	}
	close (MA);
}

#Subroutine 13: Handle MAFFT results
sub handleMafftResult { ##(mafftFolder[>id\nseq]<path>, identifySeqFile[id\tseq]<path>)
	$/ = ">";
	
	mkdir "./temp/handle/";
	my $dirname1 = "$_[0]";
	opendir (DIR1, $dirname1);
	while ((my $filename1 = readdir(DIR1))) {
		open (IN, "$_[0]/$filename1") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
		open ($filename1, ">>./temp/handle/$filename1");
		
		my $j=0;
		my @data;
		while (<IN>) {
			chomp;
			my $sequence;
			if ($_) {
				my @line = split m/\n/, $_;
				for (my $i=1;$i<=$#line;$i++) {
					$sequence .= $line[$i];
				}
				my @list = split m//, $sequence;
				for my $k (0 .. $#list) {
					$data[$k][$j]=$list[$k];
				}
				$sequence = "";
				$j++;
			}
		}
		for (@data) {
			print $filename1 join("",@{$_}),"\n";
		}
		close ($filename1);
	}
	close (IN);
	closedir(DIR1);
	
	$/ = "\n";
	open (OUT, ">$_[1]");
	my $dirname2 = "./temp/handle/";
	opendir (DIR2, $dirname2);
	while ((my $filename2 = readdir(DIR2))) {
		if (-f "./temp/handle/$filename2") {
			open (TMP, "./temp/handle/$filename2") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
			my $id = $filename2;
			$id =~ s/\.mafft//;
			print OUT $id."\t";
			while (<TMP>) {
				chomp;
				if ($_) {
					my %hash;
					my @line = split m//, $_;
					for (my $i=0;$i<=$#line;$i++) {
						chomp ($line[$i]);
						if (not exists $hash{$line[$i]}) {
							$hash{$line[$i]} = 1;
						}
						else {
							$hash{$line[$i]}++;
						}
					}
					foreach (keys %hash) {
						if ($hash{$_} == 1) {
							print OUT $_;
						}
						else {
							print OUT $_;
						}
					}
				}
				print OUT "##";
			}
			print OUT "\n";
		}
	}
	close (TMP);
	close (OUT);
	closedir(DIR2);
}

#Subroutine 14: Select the consistent and longest sequences
sub selectLongestSeq { ##(identifySeqFile[id\tseq]<path>, longestSeqFile[>id\nseq]<path>, primerIDFile[\n]<path>)
	$/ = "\n";
	
	open (IN, "$_[0]") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">$_[1]");
	open (OUTID, ">$_[2]");

	while (<IN>) {
		chomp;
		my @line = split m/\t/, $_;
		chomp ($line[0]);
		chomp ($line[1]);
		my $seq = uc($line[1]);
		my @list = split m/\#\#/, $seq;
		my $sequence;
		my $base;
		for (my $i=0;$i<=$#list;$i++) {
			chomp ($list[$i]);
			if ($list[$i] =~ /^[ATCG]$/) {
				$base = $list[$i];
			} else {
				$base = "-";
			}
			$sequence .= $base;
		}
		
		my @seq = split m/-/, $sequence;
		my %hash;
		my @value;
		for (my $j=0;$j<=$#seq;$j++) {
			chomp ($seq[$j]);
			if ($seq[$j]) {
				$hash{$seq[$j]} = length($seq[$j]);
			}
		}
		foreach my $key (sort {$hash{$b} <=> $hash{$a}} keys %hash) {
			push @value, $key;
		}
		chomp ($value[0]);
		my $primer = uc($value[0]);
		print OUT $line[0]."\t".$primer."\n";
		
		my $id = $line[0];
		if ($id =~ m/U$/) {
			$id =~ s/U$//;
			print OUTID $id."\n";
		}
	}
	close (IN);
	close (OUT);
	close (OUTID);
}

#Subroutine 15: Prepare files for design primers
sub prepareFilePrimer { ##(longestSeqFile[id\tseq]<path>, primerIDFile[\n]<path>, outputFile[=]<path>)
	$/ = "\n";
	
	open (SEQ, "$_[0]") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	my @seq = <SEQ>;
	my %hash;
	foreach (@seq) {
		chomp;
		if ($_) {
			my @line = split m/\t/, $_;
			chomp ($line[0]);
			chomp ($line[1]);
			$hash{$line[0]} = $line[1];
		}
	}
	close (SEQ);
	
	open (IN, "$_[1]") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">$_[2]");
	while (<IN>) {
		chomp;
		my $id = $_;
		my $up_stream = $_."U";
		my $target = "N" x 100; #mimic the target sequences
		my $down_stream = $_."D";
		my $seq = $hash{$up_stream}.$target.$hash{$down_stream};
		my $start = length($hash{$up_stream});
		print OUT "SEQUENCE_ID=".$id."\n"; #needed
		print OUT "SEQUENCE_TEMPLATE=".$seq."\n"; #needed
		print OUT "PRIMER_TASK=generic"."\n";
		print OUT "PRIMER_PICK_LEFT_PRIMER=1"."\n";
		print OUT "PRIMER_PICK_INTERNAL_OLIGO=0"."\n";
		print OUT "PRIMER_PICK_RIGHT_PRIMER=1"."\n";
		print OUT "PRIMER_OPT_SIZE=20"."\n";
		print OUT "PRIMER_MIN_SIZE=16"."\n";
		print OUT "PRIMER_MAX_SIZE=28"."\n";
		print OUT "PRIMER_PRODUCT_SIZE_RANGE=100-500"."\n";
		print OUT "PRIMER_EXPLAIN_NLAG=1"."\n";
		print OUT "SEQUENCE_TARGET=".$start.","."100"."\n"; #needed
		print OUT "="."\n"; #needed
	}
	close (IN);
	close (OUT);
}

#Subroutine 16: Handle Primer-3 results
sub handlePrimerResult { ##(inputFile[=]<path>, outputFile[\t]<path>)
	$/ = "=\n";
	
	open (IN, "$_[0]") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">$_[1]");
	
	while (<IN>) {
		my ($id) = (/SEQUENCE_ID=(.+?\n)/);
		$id =~ s/\r//;
		$id =~ s/\n//;
		
		# ------------------------- Primer 1 ------------------------- #
		/PRIMER_LEFT_0_SEQUENCE=(.*)/ || do {next};
		my $info = $id."\.p1\t$1\t";
		/PRIMER_LEFT_0_TM=(.*)/; $info .= "$1\t";

		/PRIMER_RIGHT_0_SEQUENCE=(.*)/;	$info .= "$1\t";
		/PRIMER_RIGHT_0_TM=(.*)/; $info .= "$1\n";
		
		# ------------------------- Primer 2 ------------------------- #
		/PRIMER_LEFT_1_SEQUENCE=(.*)/; $info .= $id."\.p2\t$1\t";
		/PRIMER_LEFT_1_TM=(.*)/; $info .= "$1\t";
			
		/PRIMER_RIGHT_1_SEQUENCE=(.*)/;	$info .= "$1\t";
		/PRIMER_RIGHT_1_TM=(.*)/; $info .= "$1\n";
		
		# ------------------------- Primer 3 ------------------------- #
		/PRIMER_LEFT_2_SEQUENCE=(.*)/; $info .= $id."\.p3\t$1\t";
		/PRIMER_LEFT_2_TM=(.*)/; $info .= "$1\t";
			
		/PRIMER_RIGHT_2_SEQUENCE=(.*)/;	$info .= "$1\t";
		/PRIMER_RIGHT_2_TM=(.*)/; $info .= "$1\n";
		
		# ------------------------- Primer 4 ------------------------- #
		/PRIMER_LEFT_3_SEQUENCE=(.*)/; $info .= $id."\.p4\t$1\t";
		/PRIMER_LEFT_3_TM=(.*)/; $info .= "$1\t";
			
		/PRIMER_RIGHT_3_SEQUENCE=(.*)/;	$info .= "$1\t";
		/PRIMER_RIGHT_3_TM=(.*)/; $info .= "$1\n";
		
		# ------------------------- Primer 5 ------------------------- #
		/PRIMER_LEFT_4_SEQUENCE=(.*)/; $info .= $id."\.p5\t$1\t";
		/PRIMER_LEFT_4_TM=(.*)/; $info .= "$1\t";
			
		/PRIMER_RIGHT_4_SEQUENCE=(.*)/;	$info .= "$1\t";
		/PRIMER_RIGHT_4_TM=(.*)/; $info .= "$1\n";
		
		print OUT $info;
	};
	
	close (IN);
	close (OUT);
}

#Subroutine 17: Primer to e-PCR
sub primerToEPCR { ##(inputFile[\t]<path>, seqFileName[$]<string>)
	$/ = "\n";
	
	open (IN, "$_[0]") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	while (<IN>) {
		chomp;
		my @line = split m/\t/, $_;
		chomp ($line[0]);
		chomp ($line[1]);
		chomp ($line[3]);
		my $prefix_thisfile = basename ($_[1], @suffix_list);
		my $out_filename = $line[0].".".$prefix_thisfile;
		system "re-PCR -s ./temp/epcr/$prefix_thisfile.hash -n 1 -g 1 $line[1] $line[3] 50-1000 > ./temp/epcr/output/$out_filename";
	}
	close (IN);
}

#Subroutine 18: Handle e-PCR results
sub handleEPCRResult { ##(inputFile[\t]<path>, outputFile1[\t]<path>, outputFile2[\t]<path>, outputFile3[\t]<path>)
	$/ = "\n";
	
	my $dirname = "./temp/epcr/output/";
	opendir (DIR, $dirname);
	open (OUT1, ">$_[0]"); #reserved
	open (TEMP1, ">$_[2]"); #eliminated, just for statistics
	while ((my $filename = readdir(DIR))) {
		open (IN, "./temp/epcr/output/$filename") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
		if ($filename =~ m/(.+\d{6}_.+\d{6}_[FR]\.p\d)\.(.+)/) {
			my $id = $1;
			my $file = $2;
			while (<IN>) {
				chomp;
				if ($_ !~ /^\#/) {
					print OUT1 $id."\t".$file."\t".$_."\n";
				}
				if (m/Error/) {
					print TEMP1 $id."\t".$file."\t".$_."\n";
				}
			}
		}
		close (IN);
	}
	close (OUT1);
	close (TEMP1);
	closedir(DIR);
	
	open (OUT1, "$_[0]") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT2, ">$_[1]"); #reserved
	open (TEMP2, ">$_[3]"); #eliminated, just for statistics
	my @file = <OUT1>;
	my @name;
	my @value;
	my $long;
	my $list;
	foreach (@file) {
		chomp;
		if ($_) {
			my @line = split m/\t/, $_;
			chomp ($line[0]);
			chomp ($line[1]);
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
	my %hash;
	foreach (0..$#value) {
		$hash{$name[$_]} = $hash{$name[$_]}."##".$value[$_];
	}
	my @name_new;
	my @value_new;
	foreach (sort keys %hash) {
		push @name_new, $_;
		push @value_new, $hash{$_};
	}
	for (my $i=0; $i<=$#value_new; $i++) {
		$value_new[$i] =~ s/^\#\#//;
		my $count = () = ($value_new[$i] =~ m/\#\#/g);
		if ($count == 1) {
			print OUT2 $name_new[$i]."\t".$value_new[$i]."\n";
		} else {
			print TEMP2 $name_new[$i]."\t".$value_new[$i]."\n";
		}
	}
	close (OUT1);
	close (OUT2);
	close (TEMP2);
}

#Subroutine 19: Change the format of input MISA results
sub changeMISAFormat { ##(inputMISAFile[\t]<path>, outputMISAFile[\t]<path>)
	$/ = "\n";
	
	undef my %hash1;
	undef my %hash2;
	open (HA, "$_[1]") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	while (<HA>) {
		chomp;
		if ($_) {
			my @id = split m/\t/, $_;
			chomp ($id[0]); #novel ID
			chomp ($id[1]); #original ID
			$hash1{$id[1]} = $id[0];
			$hash2{$id[0]} = $id[1];
		}
	}
	close (HA);
	
	open (IN, "$_[0]") || die ("\nError in $0: The file doesn't exist: $! !\n\n");
	open (OUT, ">$_[2].temp");
	while (<IN>) {
		chomp;
		my @line = split m/\t/, $_;
		chomp ($line[0]); #chromosome id
		chomp ($line[1]); #serial number
		chomp ($line[2]); #class
		chomp ($line[3]); #unit
		chomp ($line[4]); #length
		chomp ($line[5]); #start
		chomp ($line[6]); #end
		if ($hash1{$line[0]}) {
			if ($line[2] =~ m/c/) {
				print OUT $hash1{$line[0]}."\t".$line[1]."\t".$line[2]."\t"."\(N\)$line[4]"."\t".$line[4]."\t".$line[5]."\t".$line[6]."\t".$line[3]."\n";
			} else {
				print OUT $hash1{$line[0]}."\t".$line[1]."\t".$line[2]."\t".$line[3]."\t".$line[4]."\t".$line[5]."\t".$line[6]."\t".$line[3]."\n";
			}
		} elsif ($hash2{$line[0]}) {
			if ($line[2] =~ m/c/) {
				print OUT $line[0]."\t".$line[1]."\t".$line[2]."\t"."\(N\)$line[4]"."\t".$line[4]."\t".$line[5]."\t".$line[6]."\t".$line[3]."\n";
			} else {
				print OUT $_."\t".$line[3]."\n";
			}
		} else {
			die "Error: The IDs of MISA results do not match the sequence names of their corresponding FASTA files !$message";
		}
	}
	close (IN);
	close (OUT);
	copy ("$_[2].temp","$_[2]");
}
