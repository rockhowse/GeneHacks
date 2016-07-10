import unittest
# http://stackoverflow.com/questions/1054271/how-to-import-a-python-class-that-is-in-a-directory-above
import sys
# Adds higher directory to python modules path.
sys.path.append("../fasta/")
sys.path.append("../dna/")

import read_alilgnment_utils as rau
import fasta.fasta_utils as fau
import dna.dna_utils as dnau

data_dir = "../fasta/data/"
test_file_name = "phix.fa"
full_file_name = data_dir + test_file_name

debug = False


class TestReadAlignmentUtils(unittest.TestCase):

    def test_naive_exact(self):
        """
        Tests naive exact matching algorithm
        :return:
        """

        p = "word"
        t = "There would have been a time for such a word"

        occurrences, matches, mismatches = rau.naive_exact(p, t)

        self.assertEqual(len(occurrences), 1)
        self.assertEqual(occurrences[0], 40)
        self.assertEqual(matches, 6)
        self.assertEqual(mismatches, 40)

        p = "AAA"
        t = "AAATAA"

        occurrences, matches, mismatches = rau.naive_exact(p, t)

        self.assertEqual(len(occurrences), 1)
        self.assertEqual(occurrences[0], 0)
        self.assertEqual(matches, 6)
        self.assertEqual(mismatches, 3)

        genome = fau.read_genome(full_file_name)

        random_reads = dnau.get_random_reads(genome, 100, 100)

        self.assertEqual(len(random_reads), 100)

        num_matched = 0

        for random_read in random_reads:
            self.assertEqual(len(random_read), 100)

            occurrences, matches, mismatches = rau.naive_exact(random_read, genome)

            # since we KNOW there is an alignment there has to be one or more
            self.assertGreaterEqual(len(occurrences), 1)

            num_matched += 1

        self.assertEqual(num_matched, 100)

"""
    Test all read/alignment functions
"""
if __name__ == '__main__':
    unittest.main()
