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
