import unittest
# http://stackoverflow.com/questions/1054271/how-to-import-a-python-class-that-is-in-a-directory-above
import sys
# Adds higher directory to python modules path.
sys.path.append("../fasta/")
sys.path.append("../dna/")

import read_alilgnment_utils as rau
import boyer_moore as bm
import kmer_index as kmi
import subsequent_index as ssi
import datetime as dt

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

        ''' take a long time, un-comment out for full testing
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
        '''

    def test_naive_exact_with_rc(self):
        """
        Tests the naive exact match with reverse compliement awareness
        :return:
        """

        p = 'CCC'
        ten_as = 'AAAAAAAAAA'
        t = ten_as + 'CCC' + ten_as + 'GGG' + ten_as
        occurrences, _, _  = rau.naive_exact_with_rc(p, t)
        print(occurrences)

        self.assertEqual(occurrences[0], 10)
        self.assertEqual(occurrences[1], 23)

        p = 'CGCG'
        t = ten_as + 'CGCG' + ten_as + 'CGCG' + ten_as
        occurrences, _, _ = rau.naive_exact_with_rc(p, t)
        print(occurrences)

        self.assertEqual(occurrences[0], 24)
        self.assertEqual(occurrences[1], 10)

        # ----- Real world alignment ----- #
        data_dir = "../fasta/data/"
        test_file_name = "phix.fa"
        full_file_name = data_dir + test_file_name

        phix_genome = fau.read_genome(full_file_name)

        occurrences, _, _ = rau.naive_exact_with_rc('ATTA', phix_genome)

        print('offset of leftmost occurrence: %d' % min(occurrences))
        print('# occurrences: %d' % len(occurrences))

        self.assertEqual(min(occurrences), 62)
        self.assertEqual(len(occurrences), 60)

    def test_naive_mm_allowed(self):
        """
        Tests the naive exact match with 2 mismatches allowed
        :return:
        """

        p = 'CTGT'
        ten_as = 'AAAAAAAAAA'
        t = ten_as + 'CTGT' + ten_as + 'CTTT' + ten_as + 'CGGG' + ten_as

        occurrences, _, _ = rau.naive_mm_allowed(p, t, 2)
        print(occurrences)

        self.assertEqual(len(occurrences), 3)
        self.assertEqual(occurrences[0], 10)
        self.assertEqual(occurrences[1], 24)
        self.assertEqual(occurrences[2], 38)

        # ----- Real world alignment ----- #
        data_dir = "../fasta/data/"
        test_file_name = "phix.fa"
        full_file_name = data_dir + test_file_name

        phix_genome = fau.read_genome(full_file_name)

        occurrences, _, _ = rau.naive_mm_allowed('GATTACA', phix_genome, 2)

        print('offset of leftmost occurrence: %d' % min(occurrences))
        print('# occurrences: %d' % len(occurrences))

        self.assertEqual(min(occurrences), 10)
        self.assertEqual(len(occurrences), 79)

    def test_boyer_moore(self):
        """
        tests the list of occurrences found when using boyer_moore matching
        :return:
        """

        t = 'GCTACGATCTAGAATCTA'
        p = 'TCTA'

        # pre-process pattern
        p_bm = bm.BoyerMoore(p)
        occurrences, matches, mismatches = rau.boyer_moore(p, p_bm, t)

        self.assertEquals(matches, 2)
        self.assertEquals(mismatches, 5)
        self.assertEquals(occurrences[0], 7)
        self.assertEquals(occurrences[1], 14)

        # check to make sure all occurrences in our returned list match up with the offset
        for occurrence in occurrences:
            self.assertEquals(t[occurrence:occurrence+4], p)

    def test_k_mer_index(self):
        """
        tests k_mer index implementation
        :return:
        """

        t = 'GCTACGATCTAGAATCTA'
        p = 'TCTA'

        index = kmi.KMerIndex(t, 2)
        occurrences = rau.query_kmer_index(p, t, index)

        self.assertEquals(len(occurrences), 2)
        self.assertEquals(occurrences[0], 7)
        self.assertEquals(occurrences[1], 14)

    def test_approximate_match(self):
        """
        tests our approximate matching function that uses num_segments+1 segments and pigeon-hole matching with boyer-moore
        :return:
        """

        t = 'CACTTAATTTG'
        p = 'AACTTG'

        occurrences = rau.approximate_match_boyer_moore(p, t, 2)

        self.assertEquals(len(occurrences), 2)
        self.assertEquals(occurrences[0], 0)
        self.assertEquals(occurrences[1], 5)

    def test_naive_exact_with_counts(self):
        """
        Tests out modified naive exact match that now includes num_alignments and num_characters as a measure of "work"
        :return:
        """

        # example 1
        p = 'word'
        t = 'there would have been a time for such a word'
        occurrences, _, _, num_alignments, num_character_comparisons = rau.naive_exact_with_counts(p, t)
        print(occurrences, num_alignments, num_character_comparisons)

        self.assertEqual(len(occurrences), 1)
        self.assertEqual(occurrences[0], 40)
        self.assertEqual(num_alignments, 41)
        self.assertEqual(num_character_comparisons, 46)

        # example 2
        p = 'needle'
        t = 'needle need noodle needle'
        occurrences, _, _, num_alignments, num_character_comparisons = rau.naive_exact_with_counts(p, t)
        print(occurrences, num_alignments, num_character_comparisons)

        self.assertEqual(len(occurrences), 2)
        self.assertEqual(occurrences[0], 0)
        self.assertEqual(occurrences[1], 19)
        self.assertEqual(num_alignments, 20)
        self.assertEqual(num_character_comparisons, 35)

    def test_boyer_moore_with_counts(self):
        """
        Tests out modified boyer_moore that now includes num_alignments and num_characters_compared
        :return:
        """

        # example 1
        p = 'word'
        t = 'there would have been a time for such a word'
        lowercase_alphabet = 'abcdefghijklmnopqrstuvwxyz '
        p_bm = bm.BoyerMoore(p, lowercase_alphabet)
        occurrences, _, _, num_alignments, num_character_comparisons = rau.boyer_moore_with_counts(p, p_bm, t)
        print(occurrences, num_alignments, num_character_comparisons)

        self.assertEqual(len(occurrences), 1)
        self.assertEqual(occurrences[0], 40)
        self.assertEqual(num_alignments, 12)
        self.assertEqual(num_character_comparisons, 15)

        #example 2
        p = 'needle'
        t = 'needle need noodle needle'
        p_bm = bm.BoyerMoore(p, lowercase_alphabet)
        occurrences, _, _, num_alignments, num_character_comparisons = rau.boyer_moore_with_counts(p, p_bm, t)
        print(occurrences, num_alignments, num_character_comparisons)

        self.assertEqual(len(occurrences), 2)
        self.assertEqual(occurrences[0], 0)
        self.assertEqual(occurrences[1], 19)
        self.assertEqual(num_alignments, 5)
        self.assertEqual(num_character_comparisons, 18)

    def test_approximate_match_subsequence_index(self):

        #example 1
        p = 'to-morrow and to-morrow '
        t = 'to-morrow and to-morrow and to-morrow creeps in this petty pace'
        num_allowed_edits = 2
        i = 3
        k = 8
        subsequence_index = ssi.SubseqIndex(t, k, i)
        occurrences, num_index_hits = rau.approximate_match_subsequence_index(p, t, num_allowed_edits, subsequence_index)

        print(occurrences)
        self.assertEqual(len(occurrences), 2)
        self.assertEqual(occurrences[0], 0)
        self.assertEqual(occurrences[1], 14)

        print(num_index_hits)
        self.assertEqual(num_index_hits, 6)

        #example 2
        # downloaded test data from: http://www.gutenberg.org/cache/epub/1110/pg1110.txt
        data_dir = "./data/"
        file_name = 'pg1110.txt'
        full_file_name = data_dir + file_name

        t = open(full_file_name).read()
        p = 'English measure backward'
        num_allowed_edits = 2
        i = 3
        k = 8
        subsequence_index = ssi.SubseqIndex(t, k, i)
        occurrences, num_index_hits = rau.approximate_match_subsequence_index(p, t, num_allowed_edits, subsequence_index)

        print(occurrences)
        self.assertEqual(len(occurrences), 1)
        # not sure if I am getting the correct download but close enough to 135249
        self.assertEqual(occurrences[0], 132188)

        print(num_index_hits)
        self.assertEqual(num_index_hits, 3)

    def test_get_hamming_distance(self):
        """
        Tests hamming distance function

        :return:
        """

        str_1 = "GATCAACCGGTAC"
        str_2 = "GACTAAGGGGTAC"

        hamming_distance = rau.get_hamming_distance(str_1, str_2)

        self.assertEqual(hamming_distance, 4)

        str_1 = "GATCAACCGGTA"
        str_2 = "GATCAACCGGTAC"

        hamming_distance = rau.get_hamming_distance(str_1, str_2)

        self.assertEqual(hamming_distance, -1)

    def test_get_edit_distance_recursive(self):
        """
        Tests edit distance using the recursive solution (SLOW!!!)

        :return:
        """

        str_1 = "GATCAACC"
        str_2 = "GATTAAGG"

        ''' un-comment out for SLOW tests
        edit_distance = rau.get_edit_distance_recursive(str_1, str_2)

        self.assertEqual(edit_distance, 3)

        str_1 = "GATCAACC"
        str_2 = "GATCAACCA"

        edit_distance = rau.get_edit_distance_recursive(str_1, str_2)

        self.assertEqual(edit_distance, 1)


        str_1 = "Shakespear"
        str_2 = "shake spear"

        # test running time
        start_time = dt.datetime.now()

        edit_distance = rau.get_edit_distance_recursive(str_1, str_2)

        end_time = dt.datetime.now()
        total_time_sec = (end_time-start_time).total_seconds()

        # we assert this takes longer than 1 second,
        # NOTE: this is based on current generation hardware, pentium i7-4700 HQ, 2.4Ghz
        self.assertGreaterEqual(total_time_sec, 1)
        '''

    def test_get_edit_distance_dynamic_programming(self):
        """
        Tests edit distance using the dynamic programming solution

        :return:
        """

        str_1 = "GATCAACC"
        str_2 = "GATTAAGG"

        edit_distance = rau.get_edit_distance_dynamic_programming(str_1, str_2)

        self.assertEqual(edit_distance, 3)

        str_1 = "GATCAACC"
        str_2 = "GATCAACCA"

        edit_distance = rau.get_edit_distance_dynamic_programming(str_1, str_2)

        self.assertEqual(edit_distance, 1)

        str_1 = "Shakespear"
        str_2 = "shake spear"

        # test running time
        start_time = dt.datetime.now()

        edit_distance = rau.get_edit_distance_dynamic_programming(str_1, str_2)

        end_time = dt.datetime.now()
        total_time_sec = (end_time - start_time).total_seconds()

        # we assert this takes less than 1 second,
        # NOTE: this is based on current generation hardware, pentium i7-4700 HQ, 2.4Ghz
        self.assertLessEqual(total_time_sec, 1)

"""
    Test all read/alignment functions
"""
if __name__ == '__main__':
    unittest.main()
