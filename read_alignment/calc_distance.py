"""
Functions for calculating distances between two strings
"""
import sys
sys.path.append("../dna/")

# used to get the standard dna alphabet
import dna.dna_utils as dnau


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


def edit_distance_dynamic_programming_approximate(str_1, str_2):
    """
    dynamic programming implementation of the edit distance with approximate matching

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

    # initialize the first column to default empty string values
    for i in range(matrix_height):
        distance_matrix[i][0] = i
    # initialize first row to all zeros
    for i in range(matrix_width):
        distance_matrix[0][i] = 0

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

    min_distance = -1

    # minimum distances on the bottom row is the smallest edit distance
    for i in range(matrix_width):

        check_distance = distance_matrix[-1][i]

        if min_distance == -1 or check_distance < min_distance:
            min_distance = check_distance

    # at this point, the distance should be the lower right hand corner of the matrix
    return min_distance


# penalty for transitions (A<->G, C<->T) [2] vs tranversions [4] and gaps [8]
global_alignment_score = [[0, 4, 2, 4, 8],
                          [4, 0, 4, 2, 8],
                          [2, 4, 0, 4, 8],
                          [4, 2, 4, 0, 8],
                          [8, 8, 8, 8, 8]]

alphabet = dnau.dna_standard_alphabet

def edit_distance_dynamic_programming_global_alignment(str_1, str_2):
    """
    dynamic programming implementation of the edit distance, using global alignment score matrix

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

    # skip first row and use the global_alignment_score for penalty of skipped
    for i in range(1, matrix_height):
        distance_matrix[i][0] = distance_matrix[i-1][0] + \
                                global_alignment_score[alphabet.index(str_1[i-1])][-1]

    # skip first column and use the global_alignment_score for penalty of skipped
    for i in range(1, matrix_width):
        distance_matrix[0][i] = distance_matrix[0][i-1] + \
                                global_alignment_score[-1][alphabet.index(str_2[i-1])]

    # start at 1 due to empty string
    for i in range(1, matrix_height):
        for j in range(1, matrix_width):
            horizontal_distance = distance_matrix[i][j-1] + \
                                  global_alignment_score[-1][alphabet.index(str_2[j-1])]

            vertical_distance = distance_matrix[i-1][j] + \
                                global_alignment_score[alphabet.index(str_1[i-1])][-1]

            delta = 0 if str_1[i-1] == str_2[j-1] else \
                global_alignment_score[alphabet.index(str_1[i-1])][alphabet.index(str_2[j-1])]

            diagonal_distance = distance_matrix[i-1][j-1] + delta

            distance_matrix[i][j] = min(horizontal_distance,
                                        vertical_distance,
                                        diagonal_distance)

    # at this point, the distance should be the lower right hand corner of the matrix
    return distance_matrix[-1][-1]

