# SSR2Marker
An integrated pipeline for identification of SSR markers

The test data could be downloaded from the website: [ftp://www.atcgn.com/SSR2Marker](ftp://www.atcgn.com/SSR2Marker), including Hongyang.fasta and White.fasta.

Otherwise, users could also use their own sequence data only required in the FASTA format.

After obtaining the data (e.g., Hongyang.fasta, White.fasta), users simply need to type one command to run this pipeline:

  [user1@localhost]$ cd \<your directory path\>
  [user1@localhost]$ perl SSR2Marker.pl Hongyang.fasta White.fasta

After running, a total of 10 result files, including detailed information of SSR motifs, primer pairs, amplified fragments, sequence sizes, length polymorphisms and statistics calculations, are obtained.

More information is provided in the manual of SSR2Marker pipeline.
