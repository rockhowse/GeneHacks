from Bio import SeqIO

# example from http://biopython.org/wiki/SeqIO
handle = open("./data/dna2.fasta", "rU")
for record in SeqIO.parse(handle, "fasta"):
    print(record.id)
handle.close()