# SSR2Marker
An integrated pipeline for identification of SSR markers

The test data could be downloaded from the website: ftp://www.atcgn.com/SSR2Marker, including Hongyang.fasta and White.fasta.
Otherwise, users could also use their own sequence data only required in the FASTA format.

After obtaining the data (e.g., Hongyang.fasta, White.fasta), users simply need to type one command to run this pipeline:

[user1@localhost]$ cd <your directory path>
[user1@localhost]$ perl SSR2Marker.pl Hongyang.fasta White.fasta

After running, a total of 10 result files, including detailed information of SSR motifs, primer pairs, amplified fragments, sequence sizes, length polymorphisms and statistics calculations, are obtained.

More information is provided in the manual of SSR2Marker pipeline.

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
