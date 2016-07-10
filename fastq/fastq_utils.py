"""
Data structures and functions that help you parse fastq formatted dna sequencing reads
"""


def q_to_phred_33(Q):
    """
    Turn Q into Phred+33 ASCII-encoded quality
    :param Q:
    :return:
    """
    return chr(Q + 33)


def phread_33_to_q(qual):
    """
    Turn Phred+33 ASCII-encoded quality into Q
    :param qual:
    :return:
    """

    return ord(qual) - 33

def read_fastq(file_name):
    """
    Given a fastq filename, will return a list of sequences and their corresponding qualities
    :param file_name:
    :return:
    """
    sequences = []
    qualities = []

    with open(file_name) as in_file:
        while True:
            # skip header for now
            in_file.readline()
            # get the sequence data
            seq = in_file.readline().rstrip()
            # skip the + sign
            in_file.readline()
            # read in quality data
            qual = in_file.readline().rstrip()

            #exit at the end of the file
            if len(seq) == 0:
                break

            sequences.append(seq)
            qualities.append(qual)

    return sequences, qualities


def create_hist(qualities):
    """
    Creates a histagram with frequency of quality scores
    :param qualities:
    :return:
    """

    hist = [0] * 50

    for qual in qualities:
        for phred in qual:
            q = phread_33_to_q(phred)
            hist[q] += 1

    return hist


def find_gc_by_pos(sequences, read_length=100):
    """
    Given a list of reads from a fastq file, get the GC composiiton at a specific position
    :param sequences:
    :param read_length:
    :return:
    """

    gc = [0] * read_length
    totals = [0] * read_length

    for read in sequences:
        for i in range(len(read)):

            # increment GC of read
            if read[i] in ['C', 'G']:
                gc[i] += 1

            # increment totals
            totals[i] += 1

    # calculate the GC ratio for reads
    for i in range(len(gc)):
        # div/0 error
        if totals[i] > 0:
             gc[i] /= float(totals[i])

    return gc