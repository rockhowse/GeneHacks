"""
A set of functions useful for processing mult-fasta formatted data
"""


class CodonInfo:
    """
    Simple python class for holding codon related data parsed from a fasta record
    """
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]

    def __init__(self, codons=None, open_reading_frames=None):
        self.codons = codons
        self.open_reading_frames = open_reading_frames

class FastaRecord():
    """
    Simple python class for holding fasta records
    """

    def __init__(self, id, header, sequence="", frame_1_codon_info=None, frame_2_codon_info=None, frame_3_codon_info=None):
        self.id = id
        self.header = header
        self.sequence = sequence
        self.frame_1_codon_info = frame_1_codon_info
        self.frame_2_codon_info = frame_2_codon_info
        self.frame_3_codon_info = frame_3_codon_info

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
        split_line = fasta_line.split(" ")

        # trim the > character
        unique_id = split_line[0][1:].strip()

    return unique_id


def get_sequences(fasta_file_name, sequence_id=None):
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

                if cur_fasta_record:
                    cur_fasta_record.frame_1_codon_info, \
                    cur_fasta_record.frame_2_codon_info, \
                    cur_fasta_record.frame_3_codon_info = get_codons_from_sequence_three_framed(cur_fasta_record.sequence)

                    # if we have a passed in sequence_id only process this one record then GTFO
                    if sequence_id and sequence_id == fasta_id:
                        break

                fasta_id = get_record_id(line)

                # skip all sequences but this one if it's passed in
                if sequence_id:
                    if fasta_id != sequence_id:
                        continue

                fasta_header = line.strip()

                cur_fasta_record = FastaRecord(fasta_id, fasta_header)

                if fasta_id in sequences:
                    sequences[fasta_id].append(cur_fasta_record)
                else:
                    sequences[fasta_id] = [cur_fasta_record]

            # otherwise, we assume we are handling the correct FastaRecord()
            # and append the seq data until the next FastaRecord()
            else:

                if cur_fasta_record:
                    # append on this sequence line
                    cur_fasta_record.sequence += line.strip()

    # have to do codons for LAST record in the file
    if cur_fasta_record:
        cur_fasta_record.frame_1_codon_info, \
        cur_fasta_record.frame_2_codon_info, \
        cur_fasta_record.frame_3_codon_info = get_codons_from_sequence_three_framed(cur_fasta_record.sequence)

    return sequences


def get_sequence_lengths(fasta_file_name):
    """
    Returns a list of the lengths of the sequences in a fasta file
    :param fasta_file_name:
    :return sequence_lengths:
    """

    sequence_lengths = []

    sequences = get_sequences(fasta_file_name)

    for fasta_id in sequences:
        for fasta_record in sequences[fasta_id]:
            sequence_lengths.append(len(fasta_record.sequence))

    return sequence_lengths


def get_codons_from_sequence(dna_sequence, reading_frame=1):
    """
    Returns a list of all 3 character functions
    :param dna_sequence:
    :return codon_list:
    """

    codons = []

    # have to handle special case of frame 2, 3
    if reading_frame > 1:
        codons.append(dna_sequence[:reading_frame-1])

    # subract 1 because reading frames are 1, 2, 3
    cur_index = reading_frame-1

    while cur_index < len(dna_sequence):
        codons.append(dna_sequence[cur_index:cur_index+3])

        cur_index += 3

    return codons


def get_open_reading_frames_from_codons(codons, frame_num=1):
    """
    Returns a list of open_reading frames, from this list of codons
    :param codons:
    :return:
    """

    open_reading_frames = []
    str_pos = 0

    # get a list of the indexes of ALL start codons
    for start_index in range(len(codons)):
        codon = codons[start_index]

        str_pos += len(codon)

        #codon has to be len 3
        if len(codon) != 3:
            continue

        if codon == CodonInfo.start_codon:
            # increment strpos + 1 to account for DNA positioning (1 based)
            open_reading_frames.append([start_index, start_index, str_pos+1, 3, codon])

    # now go through each potential starting codon and read until we hit an end codon
    for cur_orf in open_reading_frames:

        # start directly after this start_codon
        build_index = cur_orf[0]+1
        found_end_codon = False

        while build_index < len(codons):
            codon = codons[build_index]

            # codon has to be len 3, if not, break from this funciton
            if len(codon) != 3:
                break

            cur_orf[4] += codon

            # if this codon is in the list of stop codons, clear out start codon
            if codon in CodonInfo.stop_codons:
                found_end_codon = True
                cur_orf[3] = len(cur_orf[4])
                cur_orf[1] = build_index
                break

            build_index += 1

        # if we get to here without finding a stop codon... invalid ORF
        if not found_end_codon:
            cur_orf[3] = 0
            cur_orf[4] = None

    derp = 27


    return open_reading_frames

def get_codons_from_sequence_three_framed(dna_sequence):
    """
    Returns frames 1, 2 and 3 given a single sequence
    :param dna_sequence:
    :return codon_list_f1, codon_list_f2, condon_list_f3:
    """

    codon_list_f1 = []
    codon_list_f2 = []
    codon_list_f3 = []

    # get the list of codons for frames 1, 2, and 3
    codon_list_f1 = get_codons_from_sequence(dna_sequence, 1)
    codon_list_f2 = get_codons_from_sequence(dna_sequence, 2)
    codon_list_f3 = get_codons_from_sequence(dna_sequence, 3)

    # use this list of codons to find ORF secions
    codon_orf_f1 = get_open_reading_frames_from_codons(codon_list_f1, 1)
    codon_orf_f2 = get_open_reading_frames_from_codons(codon_list_f2, 2)
    codon_orf_f3 = get_open_reading_frames_from_codons(codon_list_f3, 3)

    # return the codon info for frames 1, 2, and 3 for this single sided DNA sequence
    codon_info_f1 = CodonInfo(codons=codon_list_f1, open_reading_frames=codon_orf_f1)
    codon_info_f2 = CodonInfo(codons=codon_list_f2, open_reading_frames=codon_orf_f2)
    codon_info_f3 = CodonInfo(codons=codon_list_f3, open_reading_frames=codon_orf_f3)

    return codon_info_f1, codon_info_f2, codon_info_f3

def get_overlapping_repeats(dna_sequence, repeat_length):
    """
    Returns a dict containing the repeat_length letter seq as the key with a list containing the start position of each occurrence
    :param dna_sequence:
    :param repeat_length:
    :return overlapping_repeats:
    """

    overlapping_repeats = {}

    # loop through the seq one letter at a time, and generate a dict for all possible repeats
    for start_index in range(len(dna_sequence)):
        potential_repeat = dna_sequence[start_index:start_index+repeat_length]

        if len(potential_repeat) < repeat_length:
            continue

        if potential_repeat in overlapping_repeats:
            overlapping_repeats[potential_repeat].append(start_index+1)
        else:
            overlapping_repeats[potential_repeat] = [start_index+1]

    return overlapping_repeats