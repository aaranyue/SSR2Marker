# Author: Junyang Yue
# Program name: fasta2mafft.sh

 
## ________________________________________________________________________________
## 
## DESCRIPTION: Tool for the identification and localization of conserved SSR among close species.
## 
## SYNTAX:   sh fasta2mafft.sh
## 
## 
## USAGE: sh fasta2mafft.sh
## ________________________________________________________________________________
##


# run mafft program #

#! /bin/bash

cpu_count=`cat /proc/cpuinfo| grep "processor"| wc -l`;

while read line

do
	input=${line}".fasta"
	mafft --thread $cpu_count --quiet ./data/fasta/$input > ./data/mafft/${line}.mafft
done < ./data/mafft.id
