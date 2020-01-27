# Author: Junyang Yue
# Program name: flank2fasta.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl flank2fasta.pl
## 
## 
## USAGE: perl flank2fasta.pl
## ________________________________________________________________________________
##


# prepare FASTA file for mafft #

open (FA, "./data/merge.tab") || die ("\nError in $0: The merged file doesn't exist: $! !\n\n");

@fa = <FA>;
my %hash;
foreach (@fa) {
	chomp ($_);
	if ($_) {
		my @line = split m/\t/, $_;
		chomp ($line[0]);
		chomp ($line[5]);
		$hash{$line[0]} = $line[5];
	}
}
close (FA);


open (IN, "./data/output/merge_blast_match_select_class_compare_flank.tab") || die ("\nError in $0: The BLAST result file with flanking information doesn't exist: $! !\n\n");
open (OUT, ">./data/mafft.id"); #with no postfix, such as fasta

while (<IN>) {
	chomp;
	@line = split m/\t/, $_;
	chomp ($line[0]);
	chomp ($line[1]);
	@list0 = split m/_/, $line[0];
	@list1 = split m/_/, $line[1];
	chomp ($list0[0]);
	chomp ($list0[1]);
	chomp ($list1[0]);
	chomp ($list1[1]);
	@id = split m/-/, $list0[0]; #id
	chomp ($id[0]);
	chomp ($id[1]);
	@di = split m/-/, $list0[1]; #direction
	chomp ($di[0]);
	chomp ($di[1]); #SSR class
	@stream0 = split m/-/, $list1[0];
	chomp ($stream0[0]);
	chomp ($stream0[1]);
	@stream1 = split m/-/, $list1[1];
	chomp ($stream1[0]);
	chomp ($stream1[1]);
	$id1_up = $id[0].$stream0[0];
	$id2_up = $id[1].$stream0[1];
	$id1_down = $id[0].$stream1[0];
	$id2_down = $id[1].$stream1[1];
	
	if ($di[0] eq "F") { #blast forward
		$id1_up_seq = $hash{$id1_up};
		$id2_up_seq = $hash{$id2_up};
		$id1_down_seq = $hash{$id1_down};
		$id2_down_seq = $hash{$id2_down};
		$filenameFU = $id[0]."_".$id[1]."_"."FU"; #upstream
		$filenameFD = $id[0]."_".$id[1]."_"."FD"; #downstream
		print OUT $filenameFU."\n";
		print OUT $filenameFD."\n";
		open (OUTFU, ">./data/fasta/$filenameFU\.fasta");
		open (OUTFD, ">./data/fasta/$filenameFD\.fasta");
		print OUTFU ">".$id1_up."\n".$id1_up_seq."\n".">".$id2_up."\n".$id2_up_seq."\n";
		print OUTFD ">".$id1_down."\n".$id1_down_seq."\n".">".$id2_down."\n".$id2_down_seq."\n";
	} elsif ($di[0] eq "R") { #blast reverse
		$id1_up_seq = $hash{$id1_up};
		$id2_up_seq = $hash{$id2_up};
		$id2_up_seq_R = reverse ($id2_up_seq); #reverse_complementary_sequences
		$id2_up_seq_R =~ tr/ACGTUacgtt/TGCAAtgcaa/;
		$id1_down_seq = $hash{$id1_down};
		$id2_down_seq = $hash{$id2_down};
		$id2_down_seq_R = reverse ($id2_down_seq); #reverse_complementary_sequences
		$id2_down_seq_R =~ tr/ACGTUacgtt/TGCAAtgcaa/;
		$filenameRU = $id[0]."_".$id[1]."_"."RU"; #upstream
		$filenameRD = $id[0]."_".$id[1]."_"."RD"; #downstream
		print OUT $filenameRU."\n";
		print OUT $filenameRD."\n";
		open (OUTRU, ">./data/fasta/$filenameRU\.fasta");
		open (OUTRD, ">./data/fasta/$filenameRD\.fasta");
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
