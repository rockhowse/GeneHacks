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

"""
    Test all assembly util functions
"""
if __name__ == '__main__':
    unittest.main()
