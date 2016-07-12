"""
Most basic read/alignment algorithms
"""
import sys
sys.path.append("../dna/")

import dna.dna_utils as dnau


def boyer_moore(p, p_bm, t):
    i = 0
    occurrences = []
    matches = 0
    mismatches = 0

    while i < len(t) - len(p) + 1:
        # ammount we move after compare
        shift = 1
        mismatched = False

        # start at the end
        for j in range(len(p)-1, -1, -1):
            # mismatch!
            if not p[j] == t[i+j]:

                # calculate the bad character rule skip
                skip_bc = p_bm.bad_character_rule(j, t[i+j])

                # calculate the good suffix rule skip
                skip_gs = p_bm.good_suffix_rule(j)

                # skip the largest of skip, BC and GS skips
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                mismatches += 1

                break
        if not mismatched:
            occurrences.append(i)

            # skip the amount in the match skip of bm
            skip_gs = p_bm.match_skip()

            # skip the max of our shift of gs skip
            shift = max(shift, skip_gs)

            matches += 1
        i += shift

    return occurrences, matches, mismatches


def naive_mm_allowed(read, sequence, num_mm_allowed=2):
    """
    Matches exact read in DNA sequence, returning a list of the occurrences (offests from start of sequence), allows up to num_mm_allowed mismatches, with a default of 2
    :param read:
    :param sequence:
    :param num_mm_allowed:
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

        allowed_mm = 0
        for j in range(len(read)):
            # if the characters don'sequence match, break out immediately
            if sequence[i+j] != read[j]:
                mismatches += 1

                allowed_mm += 1

                # allow num_mm_allowed before breaking
                if allowed_mm > num_mm_allowed:
                    match = False
                    break
            else:
                matches += 1

        if match:
            # all chars matched, add to list
            occurrences.append(i)

    return occurrences, matches, mismatches

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

    matches_total = 0
    mismatches_total = 0

    occurrences, matches, mismatches = naive_exact(read, sequence)

    matches_total += matches
    mismatches_total += mismatches

    occurrences_reverse_compliment, matches, mismatches = naive_exact(dnau.reverse_complement(read), sequence)

    matches_total += matches
    mismatches_total += mismatches
    occurrences.extend(occurrences_reverse_compliment)

    # convert occurrences to set then back to list to remove any duplicates
    return list(set(occurrences)), matches, mismatches
