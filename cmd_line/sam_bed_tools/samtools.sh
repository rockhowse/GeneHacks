!#/bin/bash
# these commands were used on samtools v1.2 htslib 1.2.1

# get number of alignments ina a sam file
samtools flagstat ./data/athal_wu_0_A.bam | head -n 1
# Note that, if the file was created with a tool that includes unmapped reads into the BAM file, we would need to exclude the lines representing unmapped reads, i.e. with a ‘*’ in column 3 (chrom):
samtools view ./data/athal_wu_0_A.bam | cut -f3 | grep -v '*' | wc -l 

# get the alignments showing the read mate's unmapped (FLAG 0x8)
#WRONG:samtools view -f 0x8 athal_wu_0_A.bam | wc -l
#An alignment with an unmapped mate is marked with a ‘*’ in column 7. Note that the question asks for alignments, not reads, so we simply count the number of lines in the SAM file with a ‘*’ in column 7:
samtools view ./data/athal_wu_0_A.bam | cut -f7 | grep -c '*'

# get the count of alignments that contain deletions
#WRONG: samtools view athal_wu_0_A.bam | grep D | wc -l
#Deletions are be marked with the letter ‘D’ in the CIGAR string for the alignment, shown in column 6:
samtools view ./data/athal_wu_0_A.bam | cut -f6 | grep -c 'D'

# get the number of alignments who's mate is mapped to the same chromosome
#SAM docs column 7. RNEXT: Reference sequence name of the primary alignment of the NEXT read in the template. For
#the last read, the next read is the first read in the template. If @SQ header lines are present, RNEXT (if
#not ‘*’ or ‘=’) must be present in one of the SQ-SN tag. This field is set as ‘*’ when the information is
#unavailable, and set as ‘=’ if RNEXT is identical RNAME. If not ‘=’ and the next read in the template
#has one primary mapping (see also bit 0x100 in FLAG), this field is identical to RNAME at the primary
#line of the next read. If RNEXT is ‘*’, no assumptions can be made on PNEXT and bit 0x20
samtools view ./data/athal_wu_0_A.bam | cut -f7 | grep -c '='

# find the number of spliced alignments
# value of 'N' in CIGAR field 6
samtools view ./data/athal_wu_0_A.bam | cut -f6 | grep -c 'N'

# sort our bam file so we can index it
samtools sort ./data/athal_wu_0_A.bam ./data/athal_wu_0_A.sorted

#index sorted file (generates .bai file)
samtools index ./data/athal_wu_0_A.sorted.bam

#extract specific range using bed file
# make sure to include header
samtools view -h -L ./data/athal_wu_0_A.bed ./data/athal_wu_0_A.sorted.bam > ./data/athal_wu_0_A.extracted.bam 

#get number of alignments in extracted region
samtools view ./data/athal_wu_0_A.extracted.bam | cut -f3 | grep -v '*' | wc -l

#get number of alignments with mate unmapped in extracted region
samtools view ./data/athal_wu_0_A.extracted.bam | cut -f7 | grep -c '*'

# get number of alighnments with a deletion in extracted region
samtools view ./data/athal_wu_0_A.extracted.bam | cut -f6 | grep -c 'D'

# get number of alignments that show reads mate mapped to the same chromosome in extracted region
samtools view ./data/athal_wu_0_A.extracted.bam | cut -f7 | grep -c '='

# get the number of spliced alignments in the extracted region
samtools view ./data/athal_wu_0_A.extracted.bam | cut -f6 | grep -c 'N'
