"""
Most basic read/alignment algorithms
"""

def naive_exact(p, t):
    """
    Matches exact pattern p in text t, returning a list of the occurences (offests from start of t)
    :param p:
    :param t:
    :return occurrences:
    :return matches:
    :return mismatches:
    """

    occurrences = []
    matches = 0
    mismatches = 0

    # loop over alignments, make sure to end before length of pattern
    for i in range(len(t) - len(p) + 1):
        match = True

        for j in range(len(p)):
            # if the characters don't match, break out immediately
            if t[i+j] != p[j]:
                mismatches += 1
                match = False
                break
            else:
                matches += 1

        if match:
            # all chars matched, add to list
            occurrences.append(i)

    return occurrences, matches, mismatches
