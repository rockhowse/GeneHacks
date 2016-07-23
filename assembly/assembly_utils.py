"""
Useful data structures and functions to help solve the DNA assembly problem
"""


def overlap(str_1, str_2, min_overlap_length=3):
    """
    function that returns the number of overlapping characters between the suffix of str_1 overlaps and
    the prefix of str_2 with at LEAST min_overlap_len characters matching. Default of 3.

    :param str_1:
    :param str_2:
    :param min_overlap_length:
    :return:
    """

    start = 0

    while True:
        start = str_1.find(str_2[:min_overlap_length], start)

        # if we didn't find the suffix of str_1 in str_2, return 0
        if start == -1:
            return 0

        # if str_2 string starts with the suffix of str_1 at offset start
        # return the number of overlapping characters with str_1
        if str_2.startswith(str_1[start:]):
            return len(str_1)-start

        start += 1
