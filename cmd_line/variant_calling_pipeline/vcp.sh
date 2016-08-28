~!#/bin/bash
# command line files that implements a very simple variant calling pipline using the following:
# samtools v.1.2
# bowtie v.2.2.5
# bcftools v.1.2

#how many sequences were in the genome
cat ./data/wu_0.v7.fas | grep ">" |  wc -l

#What was the name of the third sequence in the genome file?
# Give the name only, without the “>” sign.
cat ./data/wu_0.v7.fas | grep ">" 

#What was the name of the last sequence in the genome file? 
#Give the name only, without the “>” sign.
cat ./data/wu_0.v7.fas | grep ">"

# create directory for wu_0
# icd ./data/
# mkdir wu_0

# generate bowtie2 index of wu_0_A genome 
# bowtie2-build wu_0.v7.fas wu_0/wu_0


