"""
Useful data structures and functions to help solve the DNA assembly problem
"""
from itertools import permutations


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


def naive_overlap_map(reads, min_overlap_length):
    """
    creates a dictionary of overlapped reads using naieve exact matching of suffxes from read_1 to prefixes of read_2

    :param reads:
    :param min_overlap_length:
    :return:
    """

    overlaps = {}

    for read_1, read_2 in permutations(reads, 2):
        overlap_length = overlap(read_1, read_2, min_overlap_length=min_overlap_length)

        if overlap_length > 0:
            overlaps[(read_1, read_2)] = overlap_length

    return overlaps