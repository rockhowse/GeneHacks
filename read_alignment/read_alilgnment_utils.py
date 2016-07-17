"""
Most basic read/alignment algorithms
"""
import sys
sys.path.append("../dna/")

import dna.dna_utils as dnau
import boyer_moore as bm


def approximate_match_boyer_moore(read, sequence, num_segments):
    """
    'pigeon hole' matching (approximate matching) using boyer-moore for segments
    :param read:
    :param sequence:
    :param num_segments:
    :return:
    """

    segment_length = int(round(len(read)/(num_segments+1)))
    all_matches = set()

    # go through each segment
    for i in range(num_segments+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(read))

        read_bm = bm.BoyerMoore(read[start:end], alphabet='ACGT')
        matches, _, _ = boyer_moore(read[start:end], read_bm, sequence)

        for match in matches:
            # filter match before start or the match + read len after sequence
            if match < start or match-start+len(read) > len(sequence):
                continue

            mismatches = 0

            # check left hand side matches
            for j in range(0, start):
                if not read[j] == sequence[match-start+j]:
                    mismatches += 1
                    if mismatches > num_segments:
                        break

            # check right hand side matches
            for j in range(end, len(read)):
                if not read[j] == sequence[match-start+j]:
                    mismatches += 1
                    if mismatches > num_segments:
                        break

            if mismatches <= num_segments:
                all_matches.add(match - start)

    return list(all_matches)


def query_k_mer_index(read, sequence, index):
    """
    queries a pre-k-mer indexed sequence with a desired read

    :param read:
    :param sequence:
    :param index:
    :return:
    """

    k = index.k
    offsets = []
    for i in index.query(read):
        # verification
        if read[k:] == sequence[i+k:i+len(read)]:
            offsets.append(i)
    return offsets


def boyer_moore_with_counts(read, p_bm, sequence):
    """
    Uses a pre BoyerMoore() indexed sequence to align a specific read, returns some measurement of work as denoted by
    num_alignments and num_characters_compared

    :param read:
    :param p_bm:
    :param sequence:
    :return:
    """

    i = 0
    occurrences = []
    matches = 0
    mismatches = 0

    num_alignments = 0
    num_characters_compared = 0

    while i < len(sequence) - len(read) + 1:
        # amount we move after compare
        shift = 1
        mismatched = False

        # start at the end
        for j in range(len(read)-1, -1, -1):
            num_characters_compared += 1

            # mismatch!
            if not read[j] == sequence[i+j]:
                # calculate the bad character rule skip
                skip_bc = p_bm.bad_character_rule(j, sequence[i + j])

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

        num_alignments += 1

    return occurrences, matches, mismatches, num_alignments, num_characters_compared


def boyer_moore(read, p_bm, sequence):
    """
    Uses a pre BoyerMoore() indexed sequence to align a specific read

    :param read:
    :param p_bm:
    :param sequence:
    :return:
    """

    i = 0
    occurrences = []
    matches = 0
    mismatches = 0

    while i < len(sequence) - len(read) + 1:
        # amount we move after compare
        shift = 1
        mismatched = False

        # start at the end
        for j in range(len(read)-1, -1, -1):
            # mismatch!
            if not read[j] == sequence[i+j]:

                # calculate the bad character rule skip
                skip_bc = p_bm.bad_character_rule(j, sequence[i + j])

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


def naive_exact_with_counts(read, sequence):
    """
    Matches exact pattern p in text sequence, returning a list of the occurrences (offests from start of sequence)
    additionally, these include counts for the following:

    1. num_alignments
    2. num_characters

    :param read:
    :param sequence:
    :return occurrences:
    :return matches:
    :return mismatches:
    """

    occurrences = []
    matches = 0
    mismatches = 0

    num_alignments = 0
    num_characters_compared = 0

    # loop over alignments, make sure to end before length of read
    for i in range(len(sequence) - len(read) + 1):
        match = True

        for j in range(len(read)):

            num_characters_compared += 1

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

        num_alignments += 1

    return occurrences, matches, mismatches, num_alignments, num_characters_compared


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
    Matches exact pattern p in text sequence, returning a list of the occurrences (offsets from start of sequence) in both 5' and 3' strands
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
