"""
Testing script for calling fasta_utils functiosn
"""
import unittest
import fasta_utils as fau

data_dir = "./data/"
test_file_name = "herp.dna.example.fasta"
full_file_name = data_dir + test_file_name

debug = False

class TestFastaUtils(unittest.TestCase):

    def test_is_fasta_header(self):
        """
        Tests a fasta file to see if the first line is a fasta header line
        :return:
        """

        is_fasta_header = False

        with open(full_file_name, "r") as in_file:
            for line in in_file:
                is_fasta_header = fau.is_header_line(line)

                # only testing the first line
                break

        self.assertEqual(is_fasta_header, True)

    def test_fasta_num_records(self):
        """
        Tests a fasta file to get the number of records in it
        :return:
        """

        num_records,num_errors = fau.get_num_records(full_file_name)

        if debug:
            print(full_file_name + " has: " + str(num_records) + " records.")

        self.assertGreaterEqual(num_records, 0)
        self.assertEqual(num_errors, 0)

    def test_fasta_get_headers(self):
        """
        Tests the fasta file to get the list containing the full header line for each record
        :return:
        """

        header_records = fau.get_record_headers(full_file_name)

        if debug:
            for header_record in header_records:
                print header_record.strip()

        self.assertGreaterEqual(len(header_records), 0)

    def test_fasta_get_id(self):
        """
        Tests the fasta file allowing you to get the unique id given a header record
        :return:
        """

        header_records = fau.get_record_headers(full_file_name)

        self.assertGreaterEqual(len(header_records), 0)

        unique_id = fau.get_record_id(header_records[0])

        # checks agains the first id in the first record in the supplied data file
        self.assertEqual(unique_id, "g~142022655~gb~EQ086233.1")

    def test_fasta_get_sequences(self):
        """
        Tests the fasta file to get the dictionary containing the unpacked id, header, sequence data
        :return:
        """

        sequences = fau.get_sequences(full_file_name)

        self.assertGreaterEqual(len(sequences), 0)

"""
    Test all fasta util functions
"""
if __name__ == '__main__':
    unittest.main()
