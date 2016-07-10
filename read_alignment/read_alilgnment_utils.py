"""
Most basic read/alignment algorithms
"""

def naive_exact(pattern, sequence):
    """
    Matches exact pattern p in text sequence, returning a list of the occurrences (offests from start of sequence)
    :param pattern:
    :param sequence:
    :return occurrences:
    :return matches:
    :return mismatches:
    """

    occurrences = []
    matches = 0
    mismatches = 0

    # loop over alignments, make sure to end before length of pattern
    for i in range(len(sequence) - len(pattern) + 1):
        match = True

        for j in range(len(pattern)):
            # if the characters don'sequence match, break out immediately
            if sequence[i+j] != pattern[j]:
                mismatches += 1
                match = False
                break
            else:
                matches += 1

        if match:
            # all chars matched, add to list
            occurrences.append(i)

    return occurrences, matches, mismatches
