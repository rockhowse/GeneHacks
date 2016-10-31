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
# - cuffmerge  v2.2.1
# - cuffcompare v2.2.1
# - cuffdiff    v2.2.1

# make directories needed for this process
mkdir -p ./data/
mkdir -p ./data/tophat/
mkdir -p ./data/annotation/
mkdir -p ./data/cufflinks/
mkdir -p ./data/cuffcompare/
mkdir -p ./data/cuffmerge/

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

# as our data includes non-mapped alignments, let's exclude these when getting the number of alignments (NOT READS!)
samtools view ./data/tophat/athal/Day8/accepted_hits.bam | cut -f3 | grep -v '*' | wc -l
samtools view ./data/tophat/athal/Day16/accepted_hits.bam | cut -f3 | grep -v '*' | wc -l

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

# assemble the aligned RNA-seq reads into genes and transcripts using cufflinks
cufflinks -o ./data/cufflinks/athal/Day8 ./data/tophat/athal/Day8/accepted_hits.bam
cufflinks -o ./data/cufflinks/athal/Day16 ./data/tophat/athal/Day16/accepted_hits.bam

# output dirs
# [guest@centos6 rna_seq]$ ls -al ./data/cufflinks/athal/Day8/
# total 152
# drwxrwx---. 1 root vboxsf      0 Oct 30 10:08 .
# drwxrwx---. 1 root vboxsf      0 Oct 30 10:08 ..
# -rwxrwx---. 1 root vboxsf  13437 Oct 30 10:08 genes.fpkm_tracking
# -rwxrwx---. 1 root vboxsf  15768 Oct 30 10:08 isoforms.fpkm_tracking
# -rwxrwx---. 1 root vboxsf      0 Oct 30 10:08 skipped.gtf
# -rwxrwx---. 1 root vboxsf 124933 Oct 30 10:08 transcripts.gtf
# [guest@centos6 rna_seq]$ ls -al ./data/cufflinks/athal/Day16/
# total 119
# drwxrwx---. 1 root vboxsf      0 Oct 30 10:10 .
# drwxrwx---. 1 root vboxsf      0 Oct 30 10:10 ..
# -rwxrwx---. 1 root vboxsf   5684 Oct 30 10:10 genes.fpkm_tracking
# -rwxrwx---. 1 root vboxsf   7509 Oct 30 10:10 isoforms.fpkm_tracking
# -rwxrwx---. 1 root vboxsf      0 Oct 30 10:10 skipped.gtf
# -rwxrwx---. 1 root vboxsf 107303 Oct 30 10:10 transcripts.gtf

# get the number of genes, use field 9, split by space, and then field 2
cut -f9 ./data/cufflinks/athal/Day8/transcripts.gtf | cut -d ' ' -f2 | sort -u | wc -l
cut -f9 ./data/cufflinks/athal/Day16/transcripts.gtf | cut -d ' ' -f2 | sort -u | wc -l

# get the number of transcripts, use field 9, split by space, and then field 4
cut -f9 ./data/cufflinks/athal/Day8/transcripts.gtf | cut -d ' ' -f4 | sort -u | wc -l
cut -f9 ./data/cufflinks/athal/Day16/transcripts.gtf | cut -d ' ' -f4 | sort -u | wc -l

# get the number of single transcript genes, group by field 2 and get count of only unique values
cut -f3,9 ./data/cufflinks/athal/Day8/transcripts.gtf | grep '^transcript' | cut -f2 | cut -d ' ' -f4 | uniq -cu | wc -l
cut -f3,9 ./data/cufflinks/athal/Day16/transcripts.gtf | grep '^transcript' | cut -f2 | cut -d ' ' -f4 | uniq -cu | wc -l

# get the number of single-exon transcripts find exon records, group by field 4 and get count of only unique values (uniq -cu)
cut -f9 ./data/cufflinks/athal/Day8/transcripts.gtf | grep exon | cut -d ' ' -f4 | uniq -cu | wc -l
cut -f9 ./data/cufflinks/athal/Day16/transcripts.gtf | grep exon | cut -d ' ' -f4 | uniq -cu | wc -l

# get the number of multi-exon transcripts find exon records, group by field 4 and get count of only duplicate values (uniq -cd)
cut -f9 ./data/cufflinks/athal/Day8/transcripts.gtf | grep exon | cut -d ' ' -f4 | uniq -cd | wc -l
cut -f9 ./data/cufflinks/athal/Day16/transcripts.gtf | grep exon | cut -d ' ' -f4 | uniq -cd | wc -l

# make directories for cuffcompare as it doesn't appear to make the output directory structure when used with -o
mkdir -p ./data/cuffcompare/athal/Day8
mkdir -p ./data/cuffcompare/athal/Day16

# run cuffcompare using the reference gene annotation and use -R to only consider the reference transcripts that overlap input transfrags
cuffcompare -r ./data/athal_genes.gtf -R -o ./data/cuffcompare/athal/Day8/athal ./data/cufflinks/athal/Day8/transcripts.gtf
cuffcompare -r ./data/athal_genes.gtf -R -o ./data/cuffcompare/athal/Day16/athal ./data/cufflinks/athal/Day16/transcripts.gtf

# directory structure
# [guest@centos6 rna_seq]$ ls -al ./data/cuffcompare/athal/Day8/
# total 117
# drwxrwx---. 1 root vboxsf     0 Oct 30 13:56 .
# drwxrwx---. 1 root vboxsf     0 Oct 30 13:59 ..
# -rwxrwx---. 1 root vboxsf 86277 Oct 30 14:00 athal.combined.gtf
# -rwxrwx---. 1 root vboxsf  6266 Oct 30 14:00 athal.loci
# -rwxrwx---. 1 root vboxsf  1211 Oct 30 14:00 athal.stats
# -rwxrwx---. 1 root vboxsf 24464 Oct 30 14:00 athal.tracking
# [guest@centos6 rna_seq]$ ls -al ./data/cuffcompare/athal/Day16/
# total 106
# drwxrwx---. 1 root vboxsf     0 Oct 30 14:00 .
# drwxrwx---. 1 root vboxsf     0 Oct 30 13:59 ..
# -rwxrwx---. 1 root vboxsf 89155 Oct 30 14:00 athal.combined.gtf
# -rwxrwx---. 1 root vboxsf  4782 Oct 30 14:00 athal.loci
# -rwxrwx---. 1 root vboxsf  1214 Oct 30 14:00 athal.stats
# -rwxrwx---. 1 root vboxsf 11391 Oct 30 14:00 athal.tracking

# get the number of transcripts that fully the reference transcripts (class_code ~ '=')
cut -f3 ./data/cufflinks/athal/Day8/athal.transcripts.gtf.tmap | sort | uniq -c | grep =
cut -f3 ./data/cufflinks/athal/Day16/athal.transcripts.gtf.tmap | sort | uniq -c | grep =

# get the number of splice variants for a specific gene
cat ./data/cufflinks/athal/Day8/athal.transcripts.gtf.tmap | grep AT4G20240 | wc -l
cat ./data/cufflinks/athal/Day16/athal.transcripts.gtf.tmap | grep AT4G20240 | wc -l

# get the number of transcripts that are partial reconstructions of the reference transcrips (class_code ~ 'c')
cut -f3 ./data/cufflinks/athal/Day8/athal.transcripts.gtf.tmap | sort | uniq -c | grep 'c$'
cut -f3 ./data/cufflinks/athal/Day16/athal.transcripts.gtf.tmap | sort | uniq -c | grep 'c$'

# get the number of transcripts that are potentially novel splice variants (class_code ~ 'j')
cut -f3 ./data/cufflinks/athal/Day8/athal.transcripts.gtf.tmap | sort | uniq -c | grep 'j$'
cut -f3 ./data/cufflinks/athal/Day16/athal.transcripts.gtf.tmap | sort | uniq -c | grep 'j$'

# get the number of transcripts that fall entirely within a reference intron (class_code ~ 'i')
cut -f3 ./data/cufflinks/athal/Day8/athal.transcripts.gtf.tmap | sort | uniq -c | grep 'i$'
cut -f3 ./data/cufflinks/athal/Day16/athal.transcripts.gtf.tmap | sort | uniq -c | grep 'i$'

# create the list of gtf files to compare
mkdir -p ./data/cuffmerge/athal/
rm ./data/cuffmerge/athal/transcript_gtf_list.txt
echo "./data/cufflinks/athal/Day8/transcripts.gtf" >> ./data/cuffmerge/athal/transcript_gtf_list.txt
echo "./data/cufflinks/athal/Day16/transcripts.gtf" >> ./data/cuffmerge/athal/transcript_gtf_list.txt

# file should now include the following
# [guest@centos6 rna_seq]$ cat ./data/cuffmerge/athal/transcript_gtf_list.txt 
# ./data/cufflinks/athal/Day8/transcripts.gtf
# ./data/cufflinks/athal/Day16/transcripts.gtf

# use cuffmerge and provided annotation to merge and reconcile the two ests of cufflinks transcripts (Day8 and Day16)
cuffmerge -o ./data/cuffmerge/athal/ -g ./data/athal_genes.gtf ./data/cuffmerge/athal/transcript_gtf_list.txt

# get the number of genes (loci) in the merged.gtf
cut -f9 ./data/cuffmerge/athal/merged.gtf | cut -d ' ' -f2 | sort -u | wc -l

# get the number of transcripts in the merged.gtf
cut -f9 ./data/cuffmerge/athal/merged.gtf | cut -d ' ' -f4 | sort -u | wc -l

# use cuffdiff to perform the differential expression analysis



