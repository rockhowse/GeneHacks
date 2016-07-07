"""
File that contains utility data structures and functions for handling DNA data
"""
import random

dna_standard_alphabet = "ACGT"


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