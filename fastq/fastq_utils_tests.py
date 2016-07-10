"""
tests for fastq_utils
"""

import unittest
import fastq_utils as fqu

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

        ''' hist of read qualities
        import matplotlib.pyplot as plt
        plt.bar(range(len(hist)), hist)
        plt.show()
        '''

        self.assertEqual(len(hist), 50)

"""
    Test all fastq util functions
"""
if __name__ == '__main__':
    unittest.main()