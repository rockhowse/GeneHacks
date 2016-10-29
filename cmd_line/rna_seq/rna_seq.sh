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

# make directories needed for this 
mkdir -p ./data/
mkdir -p ./data/tophat/
mkdir -p ./data/annotation/

# create directory for athal index
# uncommented out as no need to create athal dir
mkdir ./data/athal_index

# generate bowtie2 index of athal genome 
bowtie2-build ./data/athal_chr.fa ./data/athal_index/athal

# include a copy of the genome with the name 'athal.fa' in the index directory
cp ./data/athal_chr.fa ./data/athal_index/athal.fa

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

# copy the annotations to the annotation directory
cp ./data/athal_genes.gtf ./data/annotation/athal_genes.gtf

# not using these as the dirs are short enough for clean cmd line
# not used: DATA_DIR=./data/
# not used: WORK_DIR=$DATA_DIR/tophat
# not used: ANNOT=$DATA_DIR/annotation/athal_genes.gtf
# not used: ANNOT_IDX=
# not used: BWT2_IDX=$DATA_DIR/atal_index

# run tophat using supplied gene annotations and pre-generated bowtie2 index from above
tophat2 -o ./data/tophat/athal \               	# output directory
	-p 10 \                                	# 10 threads
	-G ./data/annotation/athal_genes.gtf \ 	# gene annotations
	./data/athal_index/athal \		# bowtie2 index
	./data/Day8.fastq ./data/Day16.fastq	# data to seq

# SUCCESS!
# [2016-10-29 18:05:34] Run complete: 00:00:21 elapsed
# directory structure after tophat completes
# [guest@centos6 rna_seq]$ ls -al ./data/tophat/athal/
# total 5493
# drwxrwx---. 1 root vboxsf   65536 Oct 29 18:05 .
# drwxrwx---. 1 root vboxsf       0 Oct 29 18:03 ..
# -rwxrwx---. 1 root vboxsf 5457908 Oct 29 18:05 accepted_hits.bam
# -rwxrwx---. 1 root vboxsf     732 Oct 29 18:05 align_summary.txt
# -rwxrwx---. 1 root vboxsf    2820 Oct 29 18:05 deletions.bed
# -rwxrwx---. 1 root vboxsf    2081 Oct 29 18:05 insertions.bed
# -rwxrwx---. 1 root vboxsf   22770 Oct 29 18:05 junctions.bed
# drwxrwx---. 1 root vboxsf   65536 Oct 29 18:05 logs
# -rwxrwx---. 1 root vboxsf     270 Oct 29 18:05 prep_reads.info
# -rwxrwx---. 1 root vboxsf    5124 Oct 29 18:05 unmapped.bam




