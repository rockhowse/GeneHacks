"""
A set of functiosn useful for processing mult-fasta formatted data
"""

def get_num_records(fasta_file_name):
    """
    counts number of records denoted by lines starting with ">" character
    :return: num_records
    """
    num_records = 0

    with open(fasta_file_name, "r") as in_file:

        for line in in_file:
            if line.startswith(">"):
                num_records += 1

    return num_records
