"""
Implementation of k-mer index
"""
import bisect


class KMerIndex(object):

    def __init__(self, t, k):
        self.k = k
        self.index = []

        # make sure we only index until the last k-mer
        for i in range(len(t) - k + 1):

            # store as a tuple
            self.index.append((t[i:i+k], i))

        # sort at the end
        self.index.sort()

    def query(self, pattern):
        k_mer = pattern[:self.k]

        i = bisect.bisect_left(self.index, (k_mer, -1))

        hits = []

        while i < len(self.index):

            if self.index[i][0] != k_mer:
                break

            hits.append(self.index[i][1])
            i += 1

        return hits


