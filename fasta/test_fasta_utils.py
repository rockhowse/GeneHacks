"""
Testing script for calling fasta_utils functiosn
"""
import unittest
import fasta_utils as fau

data_dir = "./data/"
test_file_name = "herp.dna.example.fasta"
full_file_name = data_dir + test_file_name

class TestFastaUtils(unittest.TestCase):

    def test_fasta_num_records(self):
        """
        Tests a fasta file to get the number of records in it
        :return:
        """

        num_records,num_errors = fau.get_num_records(full_file_name)

        print(full_file_name + " has: " + str(num_records) + " records.")

        self.assertGreaterEqual(num_records, 0)
        self.assertEqual(num_errors, 0)

    def test_fasta_get_headers(self):
        """
        Tests the fasta file to get the list containing the full header line for each record
        :return:
        """

        header_records = fau.get_record_headers(full_file_name)

        for header_record in header_records:
            print header_record.strip()

        self.assertGreaterEqual(len(header_records), 0)

"""
    Test all fasta util functions
"""
if __name__ == '__main__':
    unittest.main()
