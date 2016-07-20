"""
Functions for calculating distances between two strings
"""


def hamming_distance(str_1, str_2):
    """
    returns the hamming distance between two given strings

    :param str_1:
    :param str_2:
    :return:
    """

    distance = 0

    # if the two strings are NOT the same length, return -1
    if len(str_1) != len(str_2):
        return -1

    for i in range(0, len(str_1)):
        if str_1[i] != str_2[i]:
            distance += 1

    return distance
