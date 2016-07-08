"""
A set up useful data structures and functions for handling a fasta file containing a single genome
"""


def read_genome(file_name):
    """
    Reads a single genome from a fasta file and returns the gennome as a string
    :param file_name:
    :return:
    """
    genome = ""

    with open(file_name, "r") as in_file:
        for line in in_file:
            if line.startswith(">"):
                continue
            else:
                # take off the newline characater
                genome += line.rstrip()

        return genome