# Zika Analyis ~ 2016-09-10
Practice analizing .fasta formatted data in Python using zika sequencing data from Anderson Lab

http://andersen-lab.com/zika-sequence-local-florida-transmission/

This was a simple exercise in an attempt to track down the 5 differences between Anderson Lab's more recent 
sequence USA_ZL2_Hu0015 compared to that found in a previous sequence Cuba_ZF10_010U.

It consists of a simple python script that does the following:

* Reads in the .fasta data for both sequences using the header as element 0 and the seq as element 1
* It then goes through each seq and chooses the one with the largest number of base pairs
* Then, it iterates using the largest base pair and keeps track of mismatches and missing nucleotides
* Lastly, it itterates through the list of mismiathes and prints them out sorted by index

The results are shown below:

>10467|USA_ZL2_Hu0015.fasta  
>10463|Cuba_ZF10_010U.fasta  

>num_diff: 9  
>&nbsp;&nbsp;&nbsp;50|['G', 'A']  
>&nbsp;&nbsp;930|['T', 'C']  
>&nbsp;5624|['T', 'C']  
>&nbsp;9117|['C', 'T']  
>&nbsp;9317|['G', 'A']  
>10463|['A', '_']  
>10464|['G', '_']  
>10465|['A', '_']  
>10466|['A', '_']  

You can see that the comparison is very basic, however it does show that the US seq contains another 4 nucleotides and 5 substitutions resulting in a total of 9 differences.