"""
Most basic read/alignment algorithms
"""
import sys
sys.path.append("../dna/")

import dna.dna_utils as dnau

def naive_exact(read, sequence):
    """
    Matches exact pattern p in text sequence, returning a list of the occurrences (offests from start of sequence)
    :param read:
    :param sequence:
    :return occurrences:
    :return matches:
    :return mismatches:
    """

    occurrences = []
    matches = 0
    mismatches = 0

    # loop over alignments, make sure to end before length of read
    for i in range(len(sequence) - len(read) + 1):
        match = True

        for j in range(len(read)):
            # if the characters don'sequence match, break out immediately
            if sequence[i+j] != read[j]:
                mismatches += 1
                match = False
                break
            else:
                matches += 1

        if match:
            # all chars matched, add to list
            occurrences.append(i)

    return occurrences, matches, mismatches


def naive_exact_with_rc(read, sequence):
    """
    Matches exact pattern p in text sequence, returning a list of the occurrences (offests from start of sequence) in both 5' and 3' strands
    :param read:
    :param sequence:
    :return occurrences:
    :return matches:
    :return mismatches:
    """

    occurrences = []
    matches_total = 0
    mismatches_total = 0

    occurrences, matches, mismatches = naive_exact(read, sequence)

    matches_total += matches
    mismatches_total +=  mismatches

    occurrences_reverse_compliment, matches, mismatches = naive_exact(dnau.reverse_complement(read), sequence)

    matches_total += matches
    mismatches_total +=  mismatches
    occurrences.extend(occurrences_reverse_compliment)

    # convert occurences to set then back to list to remove any duplicates
    return list(set(occurrences)), matches, mismatches
