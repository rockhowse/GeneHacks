"""
Simple tests for use with the dna_utils.py
"""

# http://stackoverflow.com/questions/1054271/how-to-import-a-python-class-that-is-in-a-directory-above
import sys
# Adds higher directory to python modules path.
sys.path.append("../fasta/")

import unittest
import dna_utils as du
import fasta.fasta_utils as fau

debug = False


class TestDNAUtils(unittest.TestCase):

    def test_generate_random_seq(self):
        """
        Simple test to our random seq generator function
        """
        seq = du.generate_random_seq(10)

        self.assertEquals(seq.__len__(), 10)

    def test_reverse_complement(self):
        """
        Simple test to get the reverse complement of a DNA sequence
        :return:
        """
        seq = 'AGGGTCACCGTTACCTGAACCCGGGCG'

        reverse_complement = du.reverse_complement(seq)

        self.assertEquals(reverse_complement, 'CGCCCGGGTTCAGGTAACGGTGACCCT')

    def test_get_frequency_counts(self):
        """
        Simple test that gets the frequency of each base in a DNA sequence
        """

        dir_name = "../fasta/data/"
        file_name = "lambda_virus.fa"
        full_file_name = dir_name + file_name

        genome = fau.read_genome(full_file_name)

        counts = du.get_frequency_counts(genome)

        self.assertEqual(counts['A'], 12334)
        self.assertEqual(counts['C'], 11362)
        self.assertEqual(counts['G'], 12820)
        self.assertEqual(counts['T'], 11986)

"""
    Test all dna util functions
"""
if __name__ == '__main__':
    unittest.main()