# Author: Junyang Yue
# Program name: epcr2eliminate.pl

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   perl epcr2eliminate.pl
## 
## 
## USAGE: perl epcr2eliminate.pl
## ________________________________________________________________________________
##


# eliminate primers with multiple matches #

$dirname = "./data/epcr/output";
opendir (DIR, $dirname);

open (OUT, ">./data/primer_eliminate.txt"); #reserved
open (OUT1, ">./data/primer_shortlong.txt"); #eliminated, just for statistics

while (($filename = readdir(DIR))) {
	open (IN, "./data/epcr/output/$filename") || die ("\nError in $0: The ePCR result file doesn't exist: $! !\n\n");
	
	if ($filename =~ m/(.+\d{6}_.+\d{6}_[A-Z]\d)\.(.+)/) {
		
		$id = $1;
		$species = $2;
		
		while (<IN>) {
			chomp;
			if (m/^\#/) {
				#nothing to do
			} elsif ($_) {
				print OUT $id."\t".$species."\t".$_."\n";
			} 
			if (m/Error/) {
				print OUT1 $id."\t".$species."\t".$_."\n";
			}
		}
	}
}

close (IN);
close (OUT);
close (OUT1);
closedir(DIR);
