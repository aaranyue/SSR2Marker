# Author: Junyang Yue
# Program name: identical2longest.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl identical2longest.pl
## 
## 
## USAGE: perl identical2longest.pl
## ________________________________________________________________________________
##


# identify the longest region from mafft results #

open (IN, "./data/identical_region.tab") || die ("\nError in $0: The result file with information of identical regions doesn't exist: $! !\n\n");
open (OUT, ">./data/longest_identical_region.tab");
open (OUTID, ">./data/primer_prepare.id");



while (<IN>) {
	chomp;
	my %hash;
	my @seq;
	my $primer;
	
	@line = split m/\t/, $_;
	chomp ($line[0]);
	chomp ($line[1]);
	
	$title = $line[0];
	$title =~ s/\.mafft//;
	
	@list = split m/-/, $line[1];
	
	for ($i=0;$i<=$#list;$i++) {
		if ($list[$i]) {
			$hash{$list[$i]} = length($list[$i]);
		}
	}

	foreach my $key (sort {$hash{$b} <=> $hash{$a}} keys %hash) {
		my $value = $hash{$key};
		#print OUT $key,"=>",$hash{$key},"\n";
		push @seq, $key;
	}
	
	chomp ($seq[0]);
	$primer = uc($seq[0]);
	
	print OUT $title."\t".$primer."\n";
	
	$id = $title;
	
	if ($id =~ m/U$/) {
		$id =~ s/U$//;
		print OUTID $id."\n";
	}
}

close (IN);
close (OUT);
close (OUTID);
