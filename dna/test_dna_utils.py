"""
Simple tests for use with the dna_utils.py
"""

import unittest
import dna_utils as du

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

"""
    Test all dna util functions
"""
if __name__ == '__main__':
    unittest.main()