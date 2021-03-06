"""
File that contains utility data structures and functions for handling DNA data
"""
import random
import collections

# standard DNA characters
dna_standard_alphabet = "ACGT"

# map used to get the reverse compliment of a single strand
complement = {'A': 'T',
              'C': 'G',
              'G': 'C',
              'T': 'A',
              'N': 'N'}


def generate_random_seq(num_in_seq):
    """
    function to generate a random DNA sequence of dna of length num_in_seq

    :param num_in_seq:
    :return:
    """

    seq = ''

    for _ in range(num_in_seq):
        seq += random.choice(dna_standard_alphabet)

    return seq


def reverse_complement(dna_seq):
    """
    given a DNA seq, return the reverse complement

    :param dna_seq:
    :return:
    """

    rc = ""

    for base in dna_seq:
        # make sure to pre-pend
        rc = complement[base] + rc

    return rc


def get_frequency_counts(sequence):
    """
    returns a dict containing the frequency of each base for a given sequence

    :param sequence:
    :return:
    """

    ''' doesn't handle reads with 'N' in fastq files
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    for base in genome:
        base_counts[base] += 1
    '''

    count = collections.Counter()

    count.update(sequence)

    return count


def get_random_reads(genome, num_reads, read_len):
    """
    Gets a list of randomly generated "reads" from a given genome
    :param genome:
    :param num_reads:
    :param read_len:
    :return reads:
    """

    reads = []
    for _ in range(num_reads):
        start = random.randint(0, len(genome)-read_len) - 1
        reads.append(genome[start:start+read_len])

    return reads
