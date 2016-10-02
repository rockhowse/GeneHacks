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
# could also use this and take "mapped"
samtools flagstat ./data/wu_0_A_wgs.bt2.sam

#how many alignments  are reported for the original(full-match) setting?
#excle lines in the file containing unmapped reads
bowtie2 -p 4 -x ./data/wu_0/wu_0 ./data/wu_0_A_wgs.fastq -S ./data/wu_0_A_wgs.bt2.sam
cat ./data/wu_0_A_wgs.bt2.sam | grep -v "^@" | cut -f3 | grep -v "*" | wc -l
# could also use this and tape "QC passed"
samtools flagstat ./data/wu_0_A_wgs.bt2.sam

#How many matches (alignments) were reported with the local-match setting? Exclude lines in the file containing unmapped reads.
bowtie2 --local -p 4 -x ./data/wu_0/wu_0 ./data/wu_0_A_wgs.fastq -S ./data/wu_0_A_wgs.bt2.local.sam
cat ./data/wu_0_A_wgs.bt2.local.sam | grep -v "^@" | cut -f3 | grep -v "*" | wc -l
# could also use this and take "mapped"
samtools flagstat ./data/wu_0_A_wgs.bt2.local.sam

#how many alignments contained insertions and/or deletions
samtools view ./data/wu_0_A_wgs.bt2.sam | cut -f6 | grep -cE 'I|D'
samtools view ./data/wu_0_A_wgs.bt2.local.sam | cut -f6 | grep -cE 'I|D'

# get a number of entries for specific chromosome
# convert the .sam file to .bam, need the genome
samtools view -bT ./data/wu_0.v7.fas ./data/wu_0_A_wgs.bt2.sam > ./data/wu_0_A_wgs.bt2.bam

#first we need to sort the data
samtools sort ./data/wu_0_A_wgs.bt2.bam ./data/wu_0_A_wgs.bt2.sorted

# generate the mpilup format
# use -u and -v for uncompressed VCF format
samtools mpileup -f ./data/wu_0.v7.fas -uv ./data/wu_0_A_wgs.bt2.sorted.bam > ./data/wu_0_A_wgs.mpileup.vcf

# count the number of entries in the VCF file for specific chromosome
cat ./data/wu_0_A_wgs.mpileup.vcf | grep -v "^#" | cut -f1 | grep -c "^Chr3"

# get the number of entries that contain the letter of the genome
cat ./data/wu_0_A_wgs.mpileup.vcf | grep -v "^#" | cut -f4 | grep -c "^A$"

#Get number of entries with specified number of supporting reads (depth)
cat ./data/wu_0_A_wgs.mpileup.vcf | grep -v "^#" | grep -c "DP=20"

#Get number of entries that represent indels
cat ./data/wu_0_A_wgs.mpileup.vcf | grep -v "^#" | grep -c "INDEL"

#Get number of entries for specific position
#mpilup has chr and position data in fields 1 and 2
#yse $ for end of string for position
cat ./data/wu_0_A_wgs.mpileup.vcf | grep -v "^#" | cut -f1,2 | grep -Pc "^Chr1\t175672$"

#how many variants are on a specific Chromosome
#get data in BCF format
samtools mpileup -f ./data/wu_0.v7.fas -g ./data/wu_0_A_wgs.bt2.sorted.bam > ./data/wu_0_A_wgs.mpileup.bcf

#cal varients using bcftools call, multi-allelic (-m), variants only (-v) and getting uncompressed vCF ("-O -v")
bcftools call -m -v -Ov ./data/wu_0_A_wgs.mpileup.bcf > ./data/wu_0_A_wgs.mpileup.2.vcf

# now we can get variants by filtering on the field 1
cat ./data/wu_0_A_wgs.mpileup.2.vcf | grep -v "^#" | cut -f1 | grep "Chr3" | wc -l

#number of variants as TSNP of speific type
cat ./data/wu_0_A_wgs.mpileup.2.vcf | grep -v "^#" | cut -f4,5 | grep -P "^A\tT$" | wc -l

#number of variants that are indels
cat ./data/wu_0_A_wgs.mpileup.2.vcf | grep -v "^#" | grep -c "INDEL"

#number of variant entries at a specific depth
cat ./data/wu_0_A_wgs.mpileup.2.vcf | grep -v "^#" | grep -c "DP=20"

#type of variant at specific position on a given chromosome
#use your mind-brain to grock this one
cat ./data/wu_0_A_wgs.mpileup.2.vcf | grep -v "^#" | grep "Chr3" | grep "11937923"

