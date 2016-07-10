import unittest
import read_alilgnment_utils as rau

data_dir = "./data/"
test_file_name = ""
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

"""
    Test all read/alignment functions
"""
if __name__ == '__main__':
    unittest.main()
