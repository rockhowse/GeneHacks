"""
Testing script for calling fasta_utils functions
"""
import unittest
import fasta_utils as fau

data_dir = "./data/"
test_file_name = "lambda_virus.fa"
full_file_name = data_dir + test_file_name

debug = False

class TestFastaUtils(unittest.TestCase):

    def test_read_genome(self):
        """
        Tests a fasta file genome read
        :return:
        """

        genome = fau.read_genome(full_file_name)

        self.assertEqual(len(genome), 48502)

"""
    Test all fasta util functions
"""
if __name__ == '__main__':
    unittest.main()