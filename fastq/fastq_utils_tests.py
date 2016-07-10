"""
tests for fastq_utils
"""

import unittest
import fastq_utils as fqu
import dna.dna_utils as dnau

import matplotlib.pyplot as plt

data_dir = "./data/"
test_file_name = "SRR835775_1.first1000.fastq"
full_file_name = data_dir + test_file_name

debug = False


class TestFastqUtils(unittest.TestCase):

    def test_q_to_phred_33(self):
        """
        Tests converting a quality character to Phred+33 format
        :return:
        """

        phread_33 = fqu.q_to_phred_33(40)

        self.assertEqual(phread_33, 'I')

    def test_phread_33_to_q(self):
        """
        Tests converting a quality character to Phred+33 format
        :return:
        """

        q = fqu.phread_33_to_q('I')

        self.assertEqual(q, 40)

    def test_read_qualities(self):
        """
        Tests the ability to get a list of bases and their corresponding qualties from a fastq file
        :return:
        """

        sequences, qualities = fqu.read_fastq(full_file_name)

        self.assertGreater(len(sequences), 0)
        self.assertGreater(len(qualities), 0)

    def test_create_hist(self):
        """
        Tests the building of a histigram of quality scores from fastq file
        :return:
        """

        sequences, qualities = fqu.read_fastq(full_file_name)

        hist = fqu.create_hist(qualities)

        self.assertEqual(len(hist), 50)

        #hist of read qualities
        plt.bar(range(len(hist)), hist)
        plt.show()

    def test_find_gc_by_pos(self):
        """
        Tests the GC by position function
        :return:
        """

        sequences, qualities = fqu.read_fastq(full_file_name)

        gc_by_pos = fqu.find_gc_by_pos(sequences)

        self.assertEqual(len(gc_by_pos), 100)

        # line plot of GC ratio for these reads
        plt.plot(range(len(gc_by_pos)), gc_by_pos)
        plt.show()

    def test_fastq_base_dist(self):
        """
        Quick test to get the base distribution of our sequences
        :return:
        """

        sequences, qualities = fqu.read_fastq(full_file_name)

        sequence = "".join(sequences)

        base_dist = dnau.get_frequency_counts(sequence)

        self.assertEqual(base_dist['A'], 21132)
        self.assertEqual(base_dist['C'], 28272)
        self.assertEqual(base_dist['G'], 28742)
        self.assertEqual(base_dist['T'], 21836)
        # 'N' means not confident
        self.assertEqual(base_dist['N'], 18)

"""
    Test all fastq util functions
"""
if __name__ == '__main__':
    unittest.main()
