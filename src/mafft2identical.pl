# Author: Junyang Yue
# Program name: mafft2identical.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl mafft2identical.pl
## 
## 
## USAGE: perl mafft2identical.pl
## ________________________________________________________________________________
##


# find identical region from mafft results #

$dirname = "./data/mafft/";
opendir (DIR, $dirname);

open (OUT, ">./data/identical_region.tab");

while (($filename = readdir(DIR))) {
	open (IN, "./data/mafft/$filename") || die ("\nError in $0: The MAFFT result file doesn't exist: $! !\n\n");
	
	if ($filename =~ m/\d+/) {
		print OUT $filename."\t";
		
		$/ = ">";
		my @seq;
		
		while (<IN>) {
			chomp;
			if ($_) {
				$_ =~ s/\r//;
				$_ =~ s/U\n/Y/;
				$_ =~ s/D\n/Y/;
				$_ =~ s/\n//g;
				@line = split m/Y/, $_;
				if ($line[1]) {
					push @seq, $line[1];
				}
			}
		}
		
		@list1 = split m//, $seq[0];
		@list2 = split m//, $seq[1];
		
		for ($i=0;$i<=$#list1;$i++) {
			if (($list1[$i] eq '-')|($list2[$i] eq '-')) {
				print OUT '-';
			} elsif ($list1[$i] ne $list2[$i]) {
				print OUT '-';
			} else {
				print OUT $list1[$i];
			}
		}
		
		print OUT "\n";
	}
}

close (IN);
close (OUT);
closedir(DIR);
