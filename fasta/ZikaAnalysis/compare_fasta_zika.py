'''
Compares 4 zika FASTA sequences taken from:

   ~Cuba
   ~Martinique
   ~PuertoRico
   ~USA
   
Data was taken from this study:
http://andersen-lab.com/zika-sequence-local-florida-transmission/
'''

import os

data_dir = "./data/"

# dict of FASTA seq data
sequences = {}

# read in the fasta data with header at 0 and seq at 1
for file_name in os.listdir(data_dir):

    with open(data_dir + file_name, "r") as in_file:

        fasta_data = []
        for line in in_file:
            if line.startswith(">"):
                fasta_data.append(line.strip())
            else:
                fasta_data.append(line.strip())

        sequences[file_name] = fasta_data

# get the lengths of sequences
for fasta_seq in sequences:
    cur_seq = sequences[fasta_seq]
    print (str(len(cur_seq[1])) + "|" + fasta_seq)

# align sequences, two appear to be offset slightly, after this they match up pretty well
try:
    change_seq = 'Martinique_ZF1_001Sa.fasta'
    sequences[change_seq][1] = "_" + sequences[change_seq][1]
except KeyError:
    print("Martinique_ZF1_001Sa data not included")

try:
    change_seq = 'PuertoRico_ZF8_008U.fasta'
    sequences[change_seq][1] = "__" + sequences[change_seq][1]
except KeyError:
    print("PuertoRico_ZF8_008U data not included")

# find the longest seq
longest_seq = 0

# get the longest seq so we can itterate through and compare
for fasta_seq in sequences:
    cur_seq = sequences[fasta_seq]
    seq_len = len(cur_seq[1])

    if seq_len > longest_seq:
        longest_seq = seq_len

    print cur_seq[1]

diff_dict = {}

# iterate through using the longest seq len and collect anywhere there isn't a 4 nucleotide matchup
for cur_seq_index in range(longest_seq):

    char_list = [' ', ' ']

    char_index = 0
    for fasta_seq in sequences:
        cur_seq = sequences[fasta_seq]

        try:
            char_list[char_index] = cur_seq[1][cur_seq_index]
        except IndexError:
            char_list[char_index] = '_'

        char_index += 1

    # see if any chars are different, add this seq_index to the list of differences
    if len(set(char_list)) > 1:
        diff_dict[cur_seq_index] = char_list

    cur_seq_index += 1

# print number of differences
print("num_diff: " + str(len(diff_dict)))

# print out the mis-matched indexes
for cur_index in sorted(diff_dict):
    print (str(cur_index).rjust(5, " ") + "|" + str(diff_dict[cur_index]))







