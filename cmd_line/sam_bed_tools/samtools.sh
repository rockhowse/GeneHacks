!#/bin/bash
# get number of alignments ina a sam file
samtools flagstat athal_wu_0_A.bam | head -n 1
# Note that, if the file was created with a tool that includes unmapped reads into the BAM file, we would need to exclude the lines representing unmapped reads, i.e. with a ‘*’ in column 3 (chrom):
samtools view athal_wu_0_A.bam | cut -f3 | grep -v '*' | wc -l

# get the alignments showing the read mate's unmapped (FLAG 0x8)
#WRONG:samtools view -f 0x8 athal_wu_0_A.bam | wc -l
#An alignment with an unmapped mate is marked with a ‘*’ in column 7. Note that the question asks for alignments, not reads, so we simply count the number of lines in the SAM file with a ‘*’ in column 7:
samtools view athal_wu_0_A.bam | cut -f7 | grep -c '*'

# get the count of reads that contain deletions
#WRONG: samtools view athal_wu_0_A.bam | grep D | wc -l
#Deletions are be marked with the letter ‘D’ in the CIGAR string for the alignment, shown in column 6:
samtools view athal_wu_0_A.bam | cut -f6 | grep -c 'D'
