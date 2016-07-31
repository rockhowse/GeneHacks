"""
Simple tests for use with the assembly_utils.py
"""

import unittest
import assembly_utils as au

debug = False


class TestAssemblyUtils(unittest.TestCase):

    def test_overlap(self):
        str_1 = "TTACGT"
        str_2 = "CGTACCGT"

        num_overlap = au.overlap(str_1, str_2, 3)

        self.assertEqual(3, num_overlap)

        str_1 = "TTACGT"
        str_2 = "GTACCGT"

        num_overlap = au.overlap(str_1, str_2, 3)

        # we expect an overlap of 3 but this only overlaps 2 so it returns 0
        self.assertEqual(0, num_overlap)

    def test_naive_overlap_map(self):
        """
        Tests the naive_overlap_map function, which reurns a map of overlapping suffix/prefixes given a set of reads

        :return:
        """

        reads = ["ACGTTATC", "GATCAAGT", "TTCACGGA"]
        min_overlap_length = 3

        overlaps = au.naive_overlap_map(reads, min_overlap_length)

        # these reads have NO overlaps so the map is empty
        self.assertEqual(len(overlaps), 0)

        reads = ["ACGGATGATC", "GATCAAGT", "TTCACGGA"]
        min_overlap_length = 3

        overlaps = au.naive_overlap_map(reads, min_overlap_length)

        # there should now be two overlaps
        self.assertEqual(len(overlaps), 2)

        self.assertEqual(overlaps[(reads[0], reads[1])], 4)
        self.assertEqual(overlaps[(reads[2], reads[0])], 5)

    def test_overlap_all_pairs(self):
        """
        Tests the abiliyt to get all overlaping paris given a set of reads and a minimum match length of k
        :return:
        """

        # Test 1
        reads = ['ABCDEFG', 'EFGHIJ', 'HIJABC']

        all_pairs = au.overlap_all_pairs(reads, 3)

        self.assertEqual(len(all_pairs), 3)

        # Test 2
        reads = ["CGTACG", "TACGTA", "GTACGT", "ACGTAC", "GTACGA", "TACGAT"]

        all_pairs = au.overlap_all_pairs(reads, 4)

        self.assertEqual(len(all_pairs), 12)

        # Test 3
        all_pairs = au.overlap_all_pairs(reads, 5)

        self.assertEqual(len(all_pairs), 6)

    def test_shortest_common_super_string(self):

        string_set = ['ACGGTACGAGC', 'GAGCTTCGGA', 'GACACGG']

        shortest_common_super_string = au.shortest_common_super_string(string_set)

        self.assertEquals(shortest_common_super_string, 'GACACGGTACGAGCTTCGGA')

    def test_shortest_common_super_string_list(self):
        """
        Tests the function that returns the list of ALL possible shortest common super strings for a given set of strings

        :return:
        """

        # Example 1
        strings = ['ABC', 'BCA', 'CAB']

        shortest_common_super_string = au.shortest_common_super_string(strings)

        self.assertEqual('ABCAB', shortest_common_super_string)

        # Returns list of all superstrings that are tied for shortest
        shortest_common_super_strings = au.shortest_common_super_string_list(strings)

        self.assertEquals(len(shortest_common_super_strings), 3)
        self.assertEqual(shortest_common_super_strings[0], 'ABCAB')
        self.assertEqual(shortest_common_super_strings[1], 'BCABC')
        self.assertEqual(shortest_common_super_strings[2], 'CABCA')

        # Example 2
        strings = ['GAT', 'TAG', 'TCG', 'TGC', 'AAT', 'ATA']

        # Returns just one shortest superstring
        shortest_common_super_string = au.shortest_common_super_string(strings)

        self.assertEqual('TCGATGCAATAG', shortest_common_super_string)

        # Returns list of all superstrings that are tied for shortest
        shortest_common_super_strings = au.shortest_common_super_string_list(strings)

        self.assertEquals(len(shortest_common_super_strings), 10)
        self.assertEqual(shortest_common_super_strings[0], 'AATAGATCGTGC')
        self.assertEqual(shortest_common_super_strings[5], 'TCGAATAGATGC')
        self.assertEqual(shortest_common_super_strings[9], 'TGCAATCGATAG')

"""
    Test all assembly util functions
"""
if __name__ == '__main__':
    unittest.main()
