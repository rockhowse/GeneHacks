!#/bin/bash
# useful bedtools commands
# bedtools v2.24.0

# get a count of overlaps in filtered data set
# http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
#We start by running BEDtools on the alignment set restricted to the specified region (Chr3:11777000-11794000) and the GTF annotation file listed above. To allow the input to be read directly from the BAM file, we use the option ‘-abam’; in this case we will need to also specify ‘-bed’ for the BAM alignment information to be shown in BED column format in the output:
bedtools intersect -abam ./data/athal_wu_0_A.region.bam -b ./data/athal_wu_0_A_annot.gtf -bed -wo > ./data/athal_wu_0_A.overlaps.bed
more ./data/athal_wu_0_A.overlaps.bed | wc -l

# how many have 10 bases or longer
#The size of the overlap is listed in column 22 of the ‘overlaps.bed’ file. To determine those longer than 10 bases, we extract the column, sort numerically in decreasing order, and simply determine by visual inspection of the file the number of such records. For instance, in ‘vim’ we search for the first line listing ‘9’ (‘:/9’), then determine its line number (Ctrl+g). Alternatively, one can use grep with option ‘-n’ to list the lines and corresponding line numbers:
# make sure to subtract one from this answer!
cut -f22 ./data/athal_wu_0_A.overlaps.bed | sort -nrk1 | grep "^9" | head -1

# get the number of annotations that overlap the alignments
# Columns 1-12 define the alignments:
# Potential pitfalls: Multiple reads may map at the same coordinates, so the information in columns 1-3 is insufficient. The minimum information needed to define the alignments is contained in columns 1-5, which include the read ID and the flag, specifying whether this is read 1 or read 2 in a pair with the same read ID).
cut -f1-12 ./data/athal_wu_0_A.overlaps.bed | sort -u | wc -l

# how many exons have reads mapped to them
#Columns 13-21 define the exons:
cut -f13-21 ./data/athal_wu_0_A.overlaps.bed | sort -u | wc -l


