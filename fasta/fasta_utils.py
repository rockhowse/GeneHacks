"""
A set of functions useful for processing mult-fasta formatted data
"""

class FastaRecord():
    """
    Simple python class for holding fasta records
    """

    def __init__(self, id, header, sequence=""):
        self.id = id
        self.header = header
        self.sequence = sequence

def is_header_line(fasta_line):
    """
    Returns True if line is fasta header line, False otherwise
    :param fasta_line:
    :return: is_fasta_header
    """

    is_fasta_header = False

    if fasta_line.startswith(">") and fasta_line[1] != " ":
        is_fasta_header = True

    return is_fasta_header


def get_num_records(fasta_file_name):
    """
    counts number of records denoted by lines starting with ">" character
    :param fasta_file_name:
    :return num_records, num_errors:
    """
    num_records = 0
    num_errors = 0

    with open(fasta_file_name, "r") as in_file:

        for line in in_file:
            if is_header_line(line):
                # fasta headers shouldn't have a space after the ">" character
                if line[1] == " ":
                    num_errors += 1

                num_records += 1

    return num_records, num_errors


def get_record_headers(fasta_file_name):
    """
    returns a list containing all the headers for each fasta record
    :param fasta_file_name:
    :return header_records:
    """
    header_records = []

    with open(fasta_file_name, "r") as in_file:
        for line in in_file:
            # make sure it doesn't have a space after >
            if is_header_line(line):
                header_records.append(line.strip())

    return header_records


def get_record_id(fasta_line):
    """
    retuns the identifier for this fasta record
    :param fasta_line:
    :return unique_id:
    """

    unique_id = None

    if is_header_line(fasta_line):
        split_line = fasta_line.split("|")

        # trim the > character
        unique_id = split_line[0][1] + "~" + split_line[1] + "~" + split_line[2] + "~" + split_line[3]

    return unique_id


def get_sequences(fasta_file_name):
    """
    Gets a dictionary of FastaRecord objects, each containing the following:
      1. id (key)
      2. header
      3. sequence
    :param fasta_file_name:
    :return sequences:
    """

    sequences = {}

    with open(fasta_file_name, "r") as in_file:
        # as we create the records, we need to append them to the data structure
        cur_fasta_record = None

        for line in in_file:

            # if we have a fasta header, create the key for the dict
            # and use it and the header to init the FastaRecord()
            if is_header_line(line):
                fasta_id = get_record_id(line)
                fasta_header = line.strip()

                cur_fasta_record = FastaRecord(fasta_id, fasta_header)

                if fasta_id in sequences:
                    sequences[fasta_id].append(cur_fasta_record)
                else:
                    sequences[fasta_id] = [cur_fasta_record]

            # otherwise, we assume we are handling the correct FastaRecord()
            # and append the seq data until the next FastaRecord()
            else:
                # append on this sequence line
                cur_fasta_record.sequence += line.strip()

    return sequences


def get_sequence_lengths(fasta_file_name):
    """
    Returns a list of the lengths of the sequences in a fasta file
    :param fasta_file_name:
    :return:
    """

    sequence_lengths = []

    sequences = get_sequences(fasta_file_name)

    for fasta_id in sequences:
        for fasta_record in sequences[fasta_id]:
            sequence_lengths.append(len(fasta_record.sequence))

    return sequence_lengths
