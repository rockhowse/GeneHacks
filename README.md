# GeneHacks
Set of tools used for processing Geneomic and Bioinformatic data

# dna
Package containing utility data structures and functions for handling generic dna data

It supports the following data structures:

* dna_standard_alphabet ~ standard DNA text characters: "ACGT"

It supports the following functions:

* generate_random_seq(num_in_seq) ~ function to get back a random seq made from the dna_standard_alphabet

# fasta
Package containing utility data structures and functions for handling fasta formatted data

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
