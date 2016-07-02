"""
A set of functiosn useful for processing mult-fasta formatted data
"""

def get_num_records(fasta_file_name):
    """
    counts number of records denoted by lines starting with ">" character
    :param fasta_file_name:
    :return: num_records, num_errors
    """
    num_records = 0
    num_errors = 0

    with open(fasta_file_name, "r") as in_file:

        for line in in_file:
            if line.startswith(">"):
                # fasta headers shouldn't have a space after the ">" character
                if line[1] == " ":
                    num_errors += 1

                num_records += 1

    return num_records, num_errors


def get_record_headers(fasta_file_name):
    """
    returns a list containing all the headers for each fasta record
    :param fasta_file_name:
    :return: header_records
    """
    header_records = []

    with open(fasta_file_name, "r") as in_file:
        for line in in_file:
            # make sure it doesn't have a space after >
            if line.startswith(">") and line[1] != " ":
                header_records.append(line.strip())

    return header_records
