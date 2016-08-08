#!/bin/bash

# how many chromosomes
grep ">" apple.genome | wc -l

# number of genes
cut -f1 apple.genes | uniq | wc -l

# number of varients (f2)
cut -f2 apple.genes | sort -u | wc -l
# number of single splice varients
cut -f1 apple.genes | uniq -c | grep " 1 " | wc -l
# number of multi splice varients
cut -f1 apple.genes | uniq -c | grep -v " 1 " | wc -l

# How many genes are there on the ‘+’ strand
cut -f1,4 apple.genes | sort -u | grep "+" | wc -l
# how many genes are there on the '-' strand
cut -f1,4 apple.genes | sort -u | grep "-" | wc -l

#How many genes are there on chromosome chr1
cut -f1,3 apple.genes | sort -u | grep "chr1" | wc -l
#How many genes are there on chromosome chr2
cut -f1,3 apple.genes | sort -u | grep "chr2" | wc -l
#How many genes are there on chromosome chr3
cut -f1,3 apple.genes | sort -u | grep "chr3" | wc -l

#How many transcripts are there on chromosome chr1
cut -f2,3 apple.genes | sort -u | grep "chr1" | wc -l
#How many transcripts are there on chromosome chr2
cut -f2,3 apple.genes | sort -u | grep "chr2" | wc -l
#How many transcripts are there on chromosome chr3
cut -f2,3 apple.genes | sort -u | grep "chr3" | wc -l

# pull out unique genes for A, B and C
cut -f1 apple.conditionA | sort -u > apple.cA.sorted
cut -f1 apple.conditionB | sort -u > apple.cB.sorted
cut -f1 apple.conditionC | sort -u > apple.cC.sorted

# number of genes in both A and B
comm -1 -2 apple.cA.sorted apple.cB.sorted | wc -l
# number of genes specific to A
comm -2 -3 apple.cA.sorted apple.cB.sorted | wc -l
# number of genes specific to B
comm -1 -3 apple.cA.sorted apple.cB.sorted | wc -l
# number of genes common in A, B and C
cat apple.c{A,B,C}.sorted | sort | uniq -c | grep " 3 " | wc -l


