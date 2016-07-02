# GeneHacks
Set of tools used for processing Geneomic and Bioinformatic data

# fasta
Package containing utility functions for handling fasta formatted data

It supports the following functions:
* simple FastaRecord containing the following:
  1. identifier (key)
  2. header
  3. sequence (combined)
* get a count of the number of valid and potentially erouneous records
* get a list containing just the record headers
* get a dictionary containing a list FastaRecord() pulled from the file keyed on the identifier in the FastaRecord()
