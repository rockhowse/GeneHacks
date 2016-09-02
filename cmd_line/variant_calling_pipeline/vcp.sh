~!#/bin/bash
# command line files that implements a very simple variant calling pipline using the following:
# samtools v.1.2
# bowtie v.2.2.5
# bcftools v.1.2

# make the data dir
# ciommented out as we only need it once
# mkdir data

# copy in the data to the data directory
# wu_0.v7.fas		~ fasta formatted genome file
# wu_wu_0_A_wgs.fastq	~ fastq formatted reads

#how many sequences were in the genome (fasta file)
cat ./data/wu_0.v7.fas | grep ">" |  wc -l

#What was the name of the third sequence in the genome file?
# Give the name only, without the “>” sign.
# Read 3rd value
cat ./data/wu_0.v7.fas | grep ">" 

#What was the name of the last sequence in the genome file? 
#Give the name only, without the “>” sign.
# Read last value
cat ./data/wu_0.v7.fas | grep ">"

# create directory for wu_0
# uncommented out as no need to create wu_0 dir
#mkdir ./data/wu_0

# generate bowtie2 index of wu_0_A genome 
# uncommented out as it takes 1.15 min
#bowtie2-build ./data/wu_0.v7.fas ./data/wu_0/wu_0

#How many reads were in the original fastq file?
# divide this number by 4 (4 lines per read in fastq format)
cat ./data/wu_0_A_wgs.fastq | wc -l

#how many alignments  are reported for the original(full-match) setting?
#excle lines in the file containing unmapped reads
bowtie2 -p 4 -x ./data/wu_0/wu_0 ./data/wu_0_A_wgs.fastq -S ./data/wu_0_A_wgs.bt2.sam
cat ./data/wu_0_A_wgs.bt2.sam | grep -v "^@" | cut -f3 | grep -v "*" | wc -l

#How many matches (alignments) were reported with the local-match setting? Exclude lines in the file containing unmapped reads.
bowtie2 --local -p 4 -x ./data/wu_0/wu_0 ./data/wu_0_A_wgs.fastq -S ./data/wu_0_A_wgs.bt2.local.sam
cat ./data/wu_0_A_wgs.bt2.local.sam | grep -v "^@" | cut -f3 | grep -v "*" | wc -l
