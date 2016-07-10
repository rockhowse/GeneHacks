import unittest
# http://stackoverflow.com/questions/1054271/how-to-import-a-python-class-that-is-in-a-directory-above
import sys
# Adds higher directory to python modules path.
sys.path.append("../fasta/")
sys.path.append("../dna/")

import read_alilgnment_utils as rau
import fasta.fasta_utils as fau
import dna.dna_utils as dnau
import fastq.fastq_utils as fqu

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

        # ----- Test artificial random reads ----- #

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

        # ----- Test real-world random reads ----- #

        data_dir2 = "../fastq/data/"
        test_file_name2 = "ERR266411_1.first1000.fastq"
        full_file_name2 = data_dir2 + test_file_name2

        phix_reads, _ = fqu.read_fastq(full_file_name2)

        num_matched = 0
        num_read = 0
        for read in phix_reads:
            occurrences, _, _ = rau.naive_exact(read, genome)
            num_read += 1

            if len(occurrences) > 0:
                num_matched += 1

        # only 7 reads matched exactly out of 1000 in the single side
        self.assertEqual(num_matched, 7)
        self.assertEqual(num_read, 1000)

        print ('%d / %d reads matched the genome', num_matched, num_read)

        # ----- Test real-world random reads first 30 bases only ----- #

        data_dir2 = "../fastq/data/"
        test_file_name2 = "ERR266411_1.first1000.fastq"
        full_file_name2 = data_dir2 + test_file_name2

        phix_reads, _ = fqu.read_fastq(full_file_name2)

        num_matched = 0
        num_read = 0
        for read in phix_reads:
            occurrences, _, _ = rau.naive_exact(read[:30], genome)
            num_read += 1

            if len(occurrences) > 0:
                num_matched += 1

        # only 439 reads matched exactly out of 1000 in the single side
        self.assertEqual(num_matched, 459)
        self.assertEqual(num_read, 1000)

        print ('%d / %d reads matched the genome', num_matched, num_read)
        
        # ----- Test real-world random reads on both sides of the DNA ----- #

        data_dir2 = "../fastq/data/"
        test_file_name2 = "ERR266411_1.first1000.fastq"
        full_file_name2 = data_dir2 + test_file_name2

        phix_reads, _ = fqu.read_fastq(full_file_name2)

        num_matched = 0
        num_read = 0
        for read in phix_reads:
            r = read[:30]
            occurrences, _, _ = rau.naive_exact(r, genome)

            occurrences_reverse_compliment, _, _ = rau.naive_exact(dnau.reverse_complement(r), genome)

            if len(occurrences_reverse_compliment) < 0:
                derp = 27

            occurrences.extend(occurrences_reverse_compliment)

            num_read += 1

            if len(occurrences) > 0:
                num_matched += 1

        # only 439 reads matched exactly out of 1000 in the single side
        self.assertEqual(num_matched, 932)
        self.assertEqual(num_read, 1000)

        print ('%d / %d reads matched the genome', num_matched, num_read)

"""
    Test all read/alignment functions
"""
if __name__ == '__main__':
    unittest.main()
