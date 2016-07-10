# GeneHacks
Set of tools used for processing Geneomic and Bioinformatic data implemented in Python 2

# dna
Package containing utility data structures and functions for handling generic DNA data

## dna_utils
Useful data structures and functions for handling DNA data

It supports the following data structures:

* dna_standard_alphabet ~ standard DNA text characters: "ACGT"
* complement ~ dictionary of the complement of a DNA character, used for reverse compliment

It supports the following functions:

* generate_random_seq(num_in_seq) ~ function to get back a random seq made from the dna_standard_alphabet
* reverse_complement(dna_seq) ~ given a DNA sequence, return the reverse complement of that sequence
* get_random_reads(genome, num_reads, read_len) ~ returns a list containing num_reads of random "reads" from the genome of read_len

# fasta
Package containing utility data structures and functions for handling fasta formatted data

## fasta_utils
Useful functions for fast files with a single geneome

It supports the following functions:

read_genome(file_name) ~ reads an entire genome in a fasta file given the fasta file name

## fasta_utils_multi
Useful functions for handling fast files with multiple reads in them

It supports the following data structures:

* CodonInfo() containing the following:
  1. codons ~ list of all codons extracted from a given single sided DNA sequence
  2. open_reading_frames ~ list of seq staring with ATG and ending with TAA, TAG, or TGA

* FastaRecord() containing the following:
  1. id ~ unique identifier
  2. header
  3. sequence (combined)
  4. frame_1_codon_info list of all codons and ORFs starting at position 0 in DNA seq
  5. frame_2_codon_info list of all codons and ORFs starting at position 1 in DNA seq
  6. frame_3_codon_info list of all codons and ORFs starting at position 2 in DNA seq

It supports the following functions:

* get a count of the number of valid and potentially erroneous records
* get a list containing just the record headers
* get a dictionary containing a list of FastaRecord() objects pulled from the file keyed on the id in the FastaRecord()

# fastq
Package containing utility data structures and functions for handling fastq formatted data

## fastq_utils
Useful data structures and functions for handling fastq data

It supports the following functions:

* q_to_phred_33(Q) ~ Turn Q into Phred+33 ASCII-encoded quality
* phread_33_to_q(qual) ~ Turn Phred+33 ASCII-encoded quality into Q
* read_fastq(file_name) ~ given a fastq file, return a list of reads and their corresponding qualities

## Package Dependencies

* matplotlib ~ used for graphing

# read_alignment
Package containing utility data structures and functions and tests for aligning reads to sequences

## read_alignment_utils
Utility data structures and functions for aligning reads to sequences

It supports the following functions:

* naive_exact(read, sequence) ~ returns a list of offsets where the pattern occurs in the sequence as well as the number of matched and mismatched character reads