"""
Testing script for calling fasta_utils functiosn
"""
import fasta_utils as fau

data_dir = "./data/"
test_file_name = "dna.example.fasta"
full_file_name = data_dir + test_file_name

def TestNumRecords():
    num_records = 0

    num_records = fau.get_num_records(full_file_name)

    print(full_file_name + " has: " + str(num_records) + " records.")
    assert(num_records > 0)

TestNumRecords()