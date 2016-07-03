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
        self.assertEqual(unique_id, "gi|142022655|gb|EQ086233.1|43")

    def test_fasta_get_sequences(self):
        """
        Tests the fasta file to get the dictionary containing the unpacked id, header, sequence data
        :return:
        """

        sequences = fau.get_sequences(full_file_name)

        self.assertGreaterEqual(len(sequences), 0)

    def test_fasta_get_sequence_lengths(self):
        """
        Tests the lengths of the sequences in a fasta file
        :return:
        """

        sequence_lengths = fau.get_sequence_lengths(full_file_name)

        self.assertGreaterEqual(len(sequence_lengths), 0)

        for length_val in sequence_lengths:

            if debug:
                print "length: " + str(length_val)

            self.assertGreaterEqual(length_val, 0)

        if debug:
            print "min: " + str(min(sequence_lengths))
            print "max: " + str(max(sequence_lengths))

    def test_fasta_get_codons_from_seq(self):
        """
        Tests the splitting of a single stranded DNA seq into codons
        :return:
        """

        codon_list = fau.get_codons_from_sequence("AGGTGACACCGCAAGCCTTATATTAGC")

        if debug:
            for codon in codon_list:
                print codon

        self.assertGreaterEqual(len(codon_list), 0)


    def test_fasta_get_codons_from_seq_all_frames(self):
        """
        Tests the splitting of a single stranded DNA seq into codons, test frames 1, 2, 3
        :return:
        """

        codon_list_f1, codon_list_f2,codon_list_f3  = fau.get_codons_from_sequence_three_framed("AGGTGACACCGCAAGCCTTATATTAGC")

        if debug:
            print "condons in frame1:"
            for codon in codon_list_f1:
                print "\t" + codon

        self.assertGreaterEqual(len(codon_list_f1), 0)

        if debug:
            print "condons in frame2:"
            for codon in codon_list_f2:
                print "\t" + codon

        self.assertGreaterEqual(len(codon_list_f2), 0)

        if debug:
            print "condons in frame3:"
            for codon in codon_list_f3:
                print "\t" + codon

        self.assertGreaterEqual(len(codon_list_f3), 0)

"""
    Test all fasta util functions
"""
if __name__ == '__main__':
    unittest.main()
