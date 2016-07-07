"""
File that contains utility data structures and functions for handling DNA data
"""
import random

# standard DNA characters
dna_standard_alphabet = "ACGT"

# map used to get the reverse compliment of a single strand
complement = {'A': 'T',
              'C': 'G',
              'G': 'C',
              'T': 'A'}


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
