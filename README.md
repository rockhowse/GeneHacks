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

It supports the following data structures:

global_alignment_score ~ 2x2 penalty matrix with scores for alignment mis-matches using dynamic programming and global alignment. transitions (A<->G, C<->T) [2] vs tranversions [4] and gaps [8]
KMerIndex ~ class used for calculating and querying a Kmer index on a given sequence
BoyerMoore ~ class used to do Boyer-Moore pre-processing on a given read before doing alignment

It supports the following functions:

* naive_exact(read, sequence) ~ returns a list of offsets where the pattern occurs in the sequence as well as the number of matched and mismatched character reads
* naive_exact_with_rc(read, sequence) ~ returns a list of offset where the read occurs on either the forward OR reverse strand
* naive_exact_with_counts(read, sequence) ~ returns the same information as the standard naive_exact() function, but with two additional counts: num_alignments, and num_characters_compared as a measure of work done
* naive_mm_allowed(read, sequence, num_mm_allowed=2) ~  Matches exact read in DNA sequence, returning a list of the occurrences (offsets from start of sequence), allows up to num_mm_allowed mismatches, with a default of 2

* boyer_moore(read, p_bm, sequence) ~ uses the boyer moore algorithm to return matches, requires that you pre-process the read with the BoyerMoore() object which is passed in as the second parameter
* boyer_moore_with_counts(read, p_pm, sequence) ~ returns the same information as the standard boyer_moore() function, but with two additional counts: num_alignments, and num_characters_compared as a measure of work done
* approximate_match_boyer_moore(read, sequence, num_allowed_edits) ~ approximate matching function that uses num_allowed_edits+1 segments and pigeon-hole matching with boyer-moore for exact matching per segment

* query_k_mer_index(read, sequence, index) ~ searches a pre-indexed sequence stored in a KMerIndex object for the given read
* approximate_match_kmer_index(read, sequence, num_allowed_edits, kmer_index) ~ approximate matching function that uses num_allowed_edits+1 segments and pigeon-hole matching with the kmer-index for exact matching per segment

* approximate_match_subsequence_index(read, sequence, num_allowed_edits, subsequence_index) ~ approximate matching function that uses num_allowed_edits+1 segments and pigeon-hole matching with the subsequence-index for exact matching per segment

* get_hamming_distance(str_1, str_2) ~ returns the hamming distance between two strings of equal length, if the strings are not equal length, return -1
* get_edit_distance_recursive(str_1, str_2) ~ returns the edit distance between two strings implemented using a recursive technique (SLOW!!!)
* get_edit_distance_dynamic_programming(str_1, str_2) ~ returns the edit distance between two strings implemented using a dynamic programming technique
* get_edit_distance_dynamic_programming_global_alignment(str_1, str_2) ~ returns the global edit distance between two strings using a scoring matrix
* get_edit_distance_dynamic_programming_approximate(str_1, str_2) ~ returns the approximate edit distance implemented using dynamic programming 

#assembly
Package containing utility data structures and functions for dealing with genomic assembly

It supports the following functions:

* overlap(str_1, str_2, min_overlap_length=3) ~ returns the number of overlapping characters between the suffix of str_1 overlaps and the prefix of str_2 with at LEAST min_overlap_len characters matching. Default of 3.
* naive_overlap_map(reads, min_overlap_length) ~ returns a map of overlapping reads using a simple naive suffix -> prefix overlapping structure with the nodes being the reads and the edges containing the number of characters that overlapped
* overlap_all_pairs(reads, k) ~ given a list of reads and a kmer value of k, it returns a list of tuples representing each overlapping read that matches exact suffix to prefix 
* shortest_common_super_string(set_of_strings) ~ SLOW!!!! N! run time =*( <--- saddest panda
* shortest_common_super_string_list(set_of_strings) ~ SLOW!! returns a list of all the possible shortest common super strings
* pick_maximal_overlap(reads, k) ~ Finds the two reads with the maximal overlap, combines them and returns both reads and the overlap length
* greedy_shortest_common_super_string(reads, k) ~ relies on the pick_maximal_overlap() function, returns the shortest common super string in a "greedy fasion", will need to find optimal value of k


## Package Dependencies

* bisect ~ used for bisect left (binary search), used in the k-mer and subsequence index implementations
* permutations ~ used for naive_overlap_map and shortest_common_super_string 
