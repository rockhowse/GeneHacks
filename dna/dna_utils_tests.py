"""
Simple tests for use with the dna_utils.py
"""

# http://stackoverflow.com/questions/1054271/how-to-import-a-python-class-that-is-in-a-directory-above
import sys
# Adds higher directory to python modules path.
sys.path.append("../fasta/")

import unittest
import dna_utils as dnau
import fasta.fasta_utils as fau

dir_name = "../fasta/data/"
file_name = "lambda_virus.fa"
full_file_name = dir_name + file_name

debug = False

class TestDNAUtils(unittest.TestCase):

    def test_generate_random_seq(self):
        """
        Simple test to our random seq generator function
        """
        seq = dnau.generate_random_seq(10)

        self.assertEquals(seq.__len__(), 10)

    def test_reverse_complement(self):
        """
        Simple test to get the reverse complement of a DNA sequence
        :return:
        """
        seq = 'AGGGTCACCGTTACCTGAACCCGGGCG'

        reverse_complement = dnau.reverse_complement(seq)

        self.assertEquals(reverse_complement, 'CGCCCGGGTTCAGGTAACGGTGACCCT')

    def test_get_frequency_counts(self):
        """
        Simple test that gets the frequency of each base in a DNA sequence
        """

        genome = fau.read_genome(full_file_name)

        counts = dnau.get_frequency_counts(genome)

        self.assertEqual(counts['A'], 12334)
        self.assertEqual(counts['C'], 11362)
        self.assertEqual(counts['G'], 12820)
        self.assertEqual(counts['T'], 11986)

    def test_get_random_reads(self):
        """
        Simple test that gets random reads from a given genome
        """

        genome = fau.read_genome(full_file_name)

        random_reads = dnau.get_random_reads(genome, 5, 100)

        self.assertEqual(len(random_reads), 5)

        for read in random_reads:
            self.assertEqual(len(read), 100)


"""
    Test all dna util functions
"""
if __name__ == '__main__':
    unittest.main()