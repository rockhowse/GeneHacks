"""
tests for fastq_utils
"""

import unittest
import fasq_utils as fqu

data_dir = "./data/"
test_file_name = ""
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

"""
    Test all fastq util functions
"""
if __name__ == '__main__':
    unittest.main()