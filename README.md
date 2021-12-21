# User Guide of SSR2Marker

Author
======
Junyang Yue
########################################################################################

Version
=======
SSR2Marker v1.0
########################################################################################

Date
====
2021-11-28
########################################################################################

Introduction
============
SSR2Marker is an integrated pipeline for identification of SSR markers, classification
of SSR categories and design of primer pairs between any two given genome-scale 
sequences. It can find out both the monomorphic and dimorphic molecular markers with a 
definite comparison value of sequence length. This program is written in Perl and 
integrated with BLAST, MAFFT, Primer3 and e-PCR.

PERL5 is freely available at https://www.perl.org/.
BLAST is freely available at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST.
MAFFT is freely available at https://mafft.cbrc.jp/alignment/software/source.html.
Primer3 is freely available at https://github.com/primer3-org/primer3.
e-PCR is freely available at http://ftp.debian.org/debian/pool/main/e/epcr/.
########################################################################################

Download
========
SSR2Marker is freely available at from https://github.com/aaranyue/SSR2Marker.
########################################################################################

Install
=======
Users don't need to do anything extra if the dependant softwares are already installed,
use it directly !
########################################################################################

Usage
=====
perl SSR2Marker.pl [option1] <value1> [option2] <value2> ... [optionN] <valueN>
########################################################################################

Options
=======
-f1|-fasta1     <str> : A single file in FASTA format containing the genomic or
                        transcriptomic sequence(s) from species A.
-f2|-fasta2     <str> : A single file in FASTA format containing the genomic or
                        transcriptomic sequence(s) from species B.
-m1|-misa1      <str> : The MISA result of species A produced by the MISA program.
-m2|-misa2      <str> : The MISA result of species B produced by the MISA program.
-m |-motif      <str> : Setting the motifs for SSR identification through navigational
                        operations (default: 1=10,2=6,3=5,4=5,5=5,6=5).
-f |-flanking   <int> : Length of the flanking sequences between each SSR locus (must
                        be an integer and larger than 100 (bp), default: 200).
-v |-version          : Showing the version information.
-h |-help             : Showing the help information.

**Note** :      The FASTA files are necessary. The MISA files are optional. If users
                have existing MISA results, they are suggested to provide them in the
                current step. This will save about 17% operation times to obtain the
                SSR primers. While the MISA results were not provided, SSR2Marker will
                begin with performing the MISA analysis until obtaining the SSR
                primers. What the users need to know is that the MISA results may
                present some difference due to the different parameter settings.
                Anyway, both the FASTA and MISA (if provided) files from two species
                should be provided together.
########################################################################################

Syntax 1 - Suppose users just have two FASTA files named Hongyang.fasta and White.fasta.
========
[user@localhost]$ cd <your directory path>
[user@localhost]$ perl SSR2Marker.pl -f1 Hongyang.fasta -f2 White.fasta
########################################################################################

Syntax 2 - Suppose users have two FASTA files as well as their respective MISA results
           named Hongyang.fasta, Hongyang.fasta.misa, White.fasta and White.fasta.misa.
========
[user@localhost]$ cd <your directory path>
[user@localhost]$ perl SSR2Marker.pl -f1 Hongyang.fasta -f2 White.fasta\
                  -m1 Hongyang.fasta.misa -m2 White.fasta.misa
########################################################################################

Syntax 3 - Suppose users want to set the SSR motifs.
========
[user@localhost]$ cd <your directory path>
[user@localhost]$ perl SSR2Marker.pl -f1 Hongyang.fasta -f2 White.fasta -m
########################################################################################

Syntax 4 - Suppose users want to set a different length for selection of flanking
           sequences.
========
[user@localhost]$ cd <your directory path>
[user@localhost]$ perl SSR2Marker.pl -f1 Hongyang.fasta -f2 White.fasta -f 150
########################################################################################

Test
====
The test data could be downloaded from the website: https://github.com/aaranyue/SSR2Marker,
including Hongyang.fasta and White.fasta. Otherwise, users could also use their own
sequence data only required in the FASTA format. After obtaining the data (e.g.,
Hongyang.fasta, White.fasta), users are simply needed to type one command to run this
pipeline. After running, a total of 10 result files in a new folder, including detailed
information of SSR motifs, primer pairs, amplified fragments, sequence sizes, length
polymorphisms and statistics calculations, are obtained.

Email
=====
Thank you for using SSR2Marker. If you have any question or suggestion, please contact
us via the Email address: yuejy@ahau.edu.cn.
