"""
Simple tests for use with the dna_utils.py
"""

import unittest
import dna_utils as dnau

debug = False

class TestDNAUtils(unittest.TestCase):
    """
    Simple test to our random seq generator function
    """

    def test_generate_random_seq(self):

        seq = dnau.generate_random_seq(10)

        self.assertEquals(seq.__len__(), 10)

"""
    Test all dna util functions
"""
if __name__ == '__main__':
    unittest.main()