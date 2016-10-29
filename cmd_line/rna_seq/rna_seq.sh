#!/bin/bash
# This script contains commands to do a simple RNA-seq experiment to determine genes that are differentially expressed
#
# required tools used:
# - samtools    v1.2
# - bowtie2     v2.2.2
# - tophat      v2.0.14
# - cufflinks   v2.2.1
# - cuffmerege  v2.2.1
# - cuffcompare v2.2.1
# - cuffdiff    v2.2.1

# create directory for athal index
# uncommented out as no need to create athal dir
#mkdir ./data/athal_index

# generate bowtie2 index of athal genome 
#bowtie2-build ./data/athal_chr.fa ./data/athal_index/athal

# include a copy of the genome with the name 'athal.fa' in the index directory
#cp ./data/athal_chr.fa ./data/athal_index/athal.fa

# dir should now look like this

# [guest@centos6 rna_seq]$ ls -al ./data/athal_index
# total 9440
# drwxrwx---. 1 root vboxsf   65536 Oct 29 17:18 .
# drwxrwx---. 1 root vboxsf       0 Oct 29 17:13 ..
# -rwxrwx---. 1 root vboxsf 4361256 Oct 29 17:15 athal.1.bt2
# -rwxrwx---. 1 root vboxsf  125008 Oct 29 17:15 athal.2.bt2
# -rwxrwx---. 1 root vboxsf      17 Oct 29 17:15 athal.3.bt2
# -rwxrwx---. 1 root vboxsf  125000 Oct 29 17:15 athal.4.bt2
# -rwxrwx---. 1 root vboxsf  500073 Oct 29 17:18 athal.fa
# -rwxrwx---. 1 root vboxsf 4361256 Oct 29 17:15 athal.rev.1.bt2
# -rwxrwx---. 1 root vboxsf  125008 Oct 29 17:15 athal.rev.2.bt2

