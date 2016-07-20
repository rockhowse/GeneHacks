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


def edit_distance_recursive(str_1, str_2):
    """
    recursive implementation of the edit distance algorithm

    :param str_1:
    :param str_2:
    :return:
    """

    if len(str_1) == 0:
        return len(str_2)
    if len(str_2) == 0:
        return len(str_1)

    delta = 1 if str_1[-1] != str_2[-1] else 0

    return min(edit_distance_recursive(str_1[:-1], str_2[:-1]) + delta,
               edit_distance_recursive(str_1[:-1], str_2) + 1,
               edit_distance_recursive(str_1, str_2[:-1]) + 1)


def edit_distance_dynamic_programming(str_1, str_2):
    """
    dynamic programming implementation of the edit distance

    :param str_1:
    :param str_2:
    :return:
    """

    distance_matrix = []

    # fill the matrix with str_1 + 1 x str2_1 + 1 of zeroes (+1 to account for empty string
    matrix_height = len(str_1) + 1
    matrix_width = len(str_2) + 1

    for i in range(matrix_height):
        distance_matrix.append([0] * matrix_width)

    # initialize first row and column to default empty string values
    for i in range(matrix_height):
        distance_matrix[i][0] = i
    for i in range(matrix_width):
        distance_matrix[0][i] = i

    # start at 1 due to empty string
    for i in range(1, matrix_height):
        for j in range(1, matrix_width):
            horizontal_distance = distance_matrix[i][j-1] + 1
            vertical_distance = distance_matrix[i-1][j] + 1

            delta = 0 if str_1[i-1] == str_2[j-1] else 1

            diagonal_distance = distance_matrix[i-1][j-1] + delta

            distance_matrix[i][j] = min(horizontal_distance,
                                        vertical_distance,
                                        diagonal_distance)

    # at this point, the distance should be the lower right hand corner of the matrix
    return distance_matrix[-1][-1]

