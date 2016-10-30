#!/bin/bash
# This script contains commands to do a simple RNA-seq experiment to determine genes that are differentially expressed
#
# Interesting addtional information on RNA-seq techniques used here:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2672628/
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

# run tophat2 using supplied gene annotations and pre-generated bowtie2 index from above on Day8
#tophat2 -o ./data/tophat/athal/Day8 \
#	-p 10 \
#	-G ./data/annotation/athal_genes.gtf \
#	./data/athal_index/athal \i
#	./data/Day8.fastq

# SUCCESS!
# [2016-10-29 18:42:24] Run complete: 00:00:11 elapsed

# directory structure after tophat2 completes
# guest@centos6 rna_seq]$ ls -al ./data/tophat/athal/Day8
# total 2399
# drwxrwx---. 1 root vboxsf   65536 Oct 29 18:42 .
# drwxrwx---. 1 root vboxsf       0 Oct 29 18:40 ..
# -rwxrwx---. 1 root vboxsf 2301786 Oct 29 18:42 accepted_hits.bam
# -rwxrwx---. 1 root vboxsf     199 Oct 29 18:42 align_summary.txt
# -rwxrwx---. 1 root vboxsf    1260 Oct 29 18:42 deletions.bed
# -rwxrwx---. 1 root vboxsf     583 Oct 29 18:42 insertions.bed
# -rwxrwx---. 1 root vboxsf   15148 Oct 29 18:42 junctions.bed
# drwxrwx---. 1 root vboxsf   65536 Oct 29 18:42 logs
# -rwxrwx---. 1 root vboxsf      64 Oct 29 18:42 prep_reads.info
# -rwxrwx---. 1 root vboxsf    3839 Oct 29 18:42 unmapped.bam

# run tophat2 using supplied gene annotations and pre-generated bowtie2 index from above on Day1
#tophat2 -o ./data/tophat/athal/Day16 \
#	-p 10 \
#	-G ./data/annotation/athal_genes.gtf \
#        ./data/athal_index/athal \
#        ./data/Day16.fastq

# SUCCESS!
# [2016-10-29 18:58:27] Run complete: 00:00:11 elapsed

# directory structure after tophat2 completes
# [guest@centos6 rna_seq]$ ls -al ./data/tophat/athal/Day16
# total 2207
# drwxrwx---. 1 root vboxsf   65536 Oct 29 18:58 .
# drwxrwx---. 1 root vboxsf       0 Oct 29 18:58 ..
# -rwxrwx---. 1 root vboxsf 2099352 Oct 29 18:58 accepted_hits.bam
# -rwxrwx---. 1 root vboxsf     199 Oct 29 18:58 align_summary.txt
# -rwxrwx---. 1 root vboxsf    1728 Oct 29 18:58 deletions.bed
# -rwxrwx---. 1 root vboxsf    1749 Oct 29 18:58 insertions.bed
# -rwxrwx---. 1 root vboxsf   21601 Oct 29 18:58 junctions.bed
# drwxrwx---. 1 root vboxsf   65536 Oct 29 18:58 logs
# -rwxrwx---. 1 root vboxsf      64 Oct 29 18:58 prep_reads.info
# -rwxrwx---. 1 root vboxsf    1996 Oct 29 18:58 unmapped.bam

# evidently using the gtf annotations at this juncture is WRONG as -G is NOT a default non-required parameter
tophat2 -o ./data/tophat/athal/Day8 ./data/athal_index/athal ./data/Day8.fastq
tophat2 -o ./data/tophat/athal/Day16 ./data/athal_index/athal ./data/Day16.fastq

# get number of reads/mapped/mutiple alignments
cat ./data/tophat/athal/Day8/align_summary.txt
cat ./data/tophat/athal/Day16/align_summary.txt

# get the number of spliced alignments, grep for 'N' in the CIGAR string
# https://www.biostars.org/p/165061/
samtools view  ./data/tophat/athal/Day8/accepted_hits.bam | awk '($6 ~ /N/)' | wc -l
samtools view  ./data/tophat/athal/Day16/accepted_hits.bam | awk '($6 ~ /N/)' | wc -l

# count the number of splice junctions
# this wasn't covered in either tophat video (?) 
# found demo here:
# https://wikis.utexas.edu/display/bioiteam/Mapping+with+Tophat+Exercises
samtools view ./data/tophat/athal/Day8/accepted_hits.bam | cut -f 1,6 | grep 'N'|head | wc -l
samtools view ./data/tophat/athal/Day16/accepted_hits.bam | cut -f 1,6 | grep 'N'|head | wc -l

