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


def get_kmers(read, k):
    """
    return a set of kmers for the given read

    :param read:
    :return:
    """

    start = 0
    kmers = set()

    for i in range(len(read)-k+1):
        kmers.add(read[start:start+k])
        start += 1

    return kmers


def overlap_all_pairs(reads, k):
    """
    given a set of reads and kmer k, return an overlap map with all reads that contain exact matching

    :param reads:
    :param k:
    :return:
    """

    kmer_map = {}

    for i in range(len(reads)):
        kmers = get_kmers(reads[i], k)
        for kmer in kmers:
            if kmer in kmer_map:
                kmer_map[kmer].add(i)
            else:
                kmer_map[kmer] = set([i])

    matching_pairs = []
    # look up the suffix of this read in the kmer_map
    for i in range(len(reads)):
        read_suffix = reads[i][-k:]

        # if the read's suffix is in the kmer_map, we know there are overlaps
        if read_suffix in kmer_map:

            # go through the reads in the map for this suffix and overlap
            for read_b in kmer_map[read_suffix]:

                # don't match yourself
                if i != read_b:

                    # might help if I called overlap() =P
                    num_overlapping = overlap(reads[i], reads[read_b], k)

                    # only add them in if they are overlapping
                    if num_overlapping > 0:
                        matching_pairs.append((reads[i], reads[read_b]))
    '''
    for matching_pair in matching_pairs:
        print matching_pairs
    '''
    return matching_pairs


def shortest_common_super_string(set_of_strings):
    """
    given a set of strings, brute force the shortest common super string

    SLOW!!!!!!

    :param set_of_strings:
    :return:
    """

    shortest_super_string = None

    for set_of_strings_perms in permutations(set_of_strings):

        # start with first string
        super_string = set_of_strings_perms[0]

        # for every string BUT the first one
        for i in range(len(set_of_strings)-1):

            # get the overlap between the current string and the NEXT string with a minimum of 1
            overlap_length = overlap(set_of_strings_perms[i], set_of_strings_perms[i+1], min_overlap_length=1)

            # concatenate the overlaps
            super_string += set_of_strings_perms[i+1][overlap_length:]

        if shortest_super_string is None or len(super_string) < len(shortest_super_string):
            shortest_super_string = super_string

    return shortest_super_string


def pick_maximal_overlap(reads, k):
    """
    Finds the two reads with the maximal overlap, combines them and returns both reads and the overlap length

    :param reads:
    :param k:
    :return:
    """
    read_a, read_b = None, None

    best_overlap_length = 0

    for a,b in permutations(reads, 2):
        overlap_length = overlap(a, b, min_overlap_length=k)

        # only returns the FIRST overlap length we run into
        if overlap_length > best_overlap_length:
            read_a, read_b = a, b

            best_overlap_length = overlap_length

    return read_a, read_b, best_overlap_length


def greedy_shortest_common_super_string(reads, k):
    """
    Greedy implementation of the shortest common super string, may not return the "shortest" but is much faster than
    "brute force" method.

    :param reads:
    :param k:
    :return:
    """
    read_a, read_b, overlap_length = pick_maximal_overlap(reads, k)

    # remove the two that had the longest overlap
    while overlap_length > 0:
        reads.remove(read_a)
        reads.remove(read_b)

        # now add in the "best overlap" for further consideration
        reads.append(read_a + read_b[overlap_length:])

        # recalculate the reads using the newly modified list
        read_a, read_b, overlap_length = pick_maximal_overlap(reads, k)

    return ''.join(reads)


def deb_bruin_ize(read, k):
    """
    Given a string, create a de-bruijn graph

    :param read:
    :param k:
    :return:
    """
    edges = []
    nodes = set()

    # go through and set up nodes/edges with kmers of size k
    for i in range(len(read) - k + 1):
        edges.append((read[i:i+k-1], read[i+1:i+k]))
        nodes.add(read[i:i+k-1])
        nodes.add(read[i+1:i+k])

    return nodes, edges


def shortest_common_super_string_list(set_of_strings):
    """
    given a set of strings, brute force the shortest common super string and return ALL strings
    who are the length of the shortest common super string

    SLOW!!!!!!

    :param set_of_strings:
    :return:
    """

    shortest_super_string = None
    shortest_super_strings = set()

    for set_of_strings_perms in permutations(set_of_strings):

        # start with first string
        super_string = set_of_strings_perms[0]

        # for every string BUT the first one
        for i in range(len(set_of_strings)-1):

            # get the overlap between the current string and the NEXT string with a minimum of 1
            overlap_length = overlap(set_of_strings_perms[i], set_of_strings_perms[i+1], min_overlap_length=1)

            # concatenate the overlaps
            super_string += set_of_strings_perms[i+1][overlap_length:]

        # add in strings that are less than OR equal to the shortest one
        if shortest_super_string is None or len(super_string) <= len(shortest_super_string):
            shortest_super_strings.add(super_string)

        if shortest_super_string is None or len(super_string) < len(shortest_super_string):
            shortest_super_string = super_string

    shortest_length = len(shortest_super_string)

    remove_vals = []

    # build list of indexes to remove
    for cur_string in shortest_super_strings:
        if len(cur_string) > shortest_length:
            remove_vals.append(cur_string)

    # remove any that are longer than shortest
    for remove_val in remove_vals:
        shortest_super_strings.remove(remove_val)

    return sorted(list(shortest_super_strings))

