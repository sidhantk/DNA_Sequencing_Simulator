#!/usr/env/bin python

# Name: Sidhant Karamchandani (skaramch)
# Group: Daniel Hwang (djhwang)

'''
Desc: This project will input a fasta file and fragment it using barcodes. It will ask
the user how many big fragments, how many small fragments, and what size fragments they
want. It will then proceed to randomly pick big fragments of the specified size from
the sequences. It will then randomly pick small fragments of the specified size
from the big fragments. It will output the results to two fastq files. One fastq file
will have a list of all the small fragments, and the other fastq file will have a list of
randomly generated barcodes (one barcode per big fragment, so more than one small fragment
will have the same barcode). Each fragment and its corresponding barcode in the other file
will have the same header. The header will say which original sequence the fragment came from,
what number small fragment it is from the big fragment it came from, and the index the small
fragment starts at in its original sequence.

Input: FASTA file with various DNA sequences
Output: two FASTQ files, one with small fragments, one with barcodes
Assumptions: FASTA file will only have DNA sequences, user won't need to use quality values
in the FASTQ files

Pseudocode:
    1) Import modules that are needed, SeqIO and random
    2) Ask user to input the fasta file
    3) Open output fastq files in write mode
    4) Make empty list for sequences.
    5) Set variables for how many big fragments, how many small fragments, and what size fragments the user wants.
    6) Parse the fasta file and add the string of each sequence to the empty list
    7) Set a counter for the total sequence length
    8) Use a for loop to add every sequence length to the total
    9) Make an empty list for the fractions
    10) Use a for loop to get every sequence fraction (sequence length over total length) and add to list
    11) Make empty dictionary
    12) Use for loop to make dictionary. Keys are the sequence fractions and values
    are the actual sequence/sequence number + 1
    13) Sort fraction list numerically
    14) Use for loop to add the list of fractions to each next fraction, so we have a list of fractions from 0 to 1
    15) Set counter, then set while loop that goes while count is less than number of big fragments
    16) Get random float from 0-1 to basically randomly choose which sequence a fragment is chosen from (by percent chance)
    17) Set new counter for another while loop
    18) Make a list with the four base pairs
    19) Make empty lsit for barcodes
    20) Make while loop that occurs 25 times. It will take a random character from base pair list and add it to another list
    21) Turn that barcode list into a string
    22) Make empty list for quality values. For loop for every character in barcode, add 'b' to list. Turn list to string.
    23) For loop to go through every sequence fraction
    24) Make conditional: if random float is less than or equal to sequence fraction, then:
    24) Get the corresponding dictionary value (sequence) for that sequence fraction in list
    25) Generate random integer from 0 to the length of the chosen sequence minus the big fragment length
    26) Get a big fragment from that sequence at the above index to the specified length minus 1
    27) Set 2 more counters
    28) Make while loop that goes while count is less than number of small fragments
    29) Generate random integer from 0 to the big fragment length minus the small fragment length
    30) Get small fragment from big fragment at the above index to the specified length
    31) Make another empty list for quality values
    32) For loop going through small fragment, and appending a 'b' the above list for every base pair in small fragment. Turn list into string
    33) Increase counters
    34) Write results to output files
    35) Increase counters
    36) Make break statment for conditional
    37) Increase counter
    38) Close outputfiles
'''

# import modules needed
from Bio import SeqIO
import random


# tell use to input fasta file
fastafile = raw_input('Enter filename here:')

# make output fastq files
outputfile1 = open('case1read1.fastq','w')
outputfile2 = open('case1read2.fastq','w')

# make list for sequences
seqlist = []

# ask user how many fragments they want
big_frag_number = int(raw_input('How many big fragments do you want? Type number:'))
small_frag_number = int(raw_input('How many small fragments do you want? Type number:'))

# ask user what size fragments they want
big_frag_length = int(raw_input('What size big fragments do you want? Type number:'))
small_frag_length = int(raw_input('What size small fragments do you want? Type number:'))

# parse fasta file
for seq_record in SeqIO.parse(fastafile,'fasta'):
    seqlist.append(str(seq_record.seq))

# set counter for total sequence length
seq_total_len = 0

# for loop to add every sequence length to the total
for seq in seqlist:
    seq_total_len += len(seq)


# empty list
seq_frac_list = []

# for loop to get fraction of each sequence length, and add it to list
for seq in seqlist:
    seq_frac = len(seq)/float(seq_total_len)
    seq_frac_list.append(seq_frac)

# empty dictionary
seq_frac_dict = {}

# make the dictionary with keys of sequence fractions and values of the actual sequences
for seqfrac in range(0, len(seq_frac_list)):
    seq_frac_dict[seq_frac_list[seqfrac]] = (seqlist[seqfrac],seqfrac + 1)


# sort the fraction list numerically
seq_frac_list.sort()

# add the list of fractions to each next fraction, so we have fractions from 0
# to 1, which will be needed when we have a percent chance generator later
runningtotal = 0
cumsum = []
for i in range(0,len(seq_frac_list)):
    runningtotal += seq_frac_list[i]
    cumsum.append(runningtotal)


# set counter
count = 0


# while loop that goes while count is less than number of big fragments
while count < big_frag_number:
    # get random float from 0-1 to basically randomly choose which sequence a fragment is chosen from (by percent chance)
    j = random.random()
    # new counter for while loop under this
    count4 = 1
    # list of base pairs
    base_pair_list = ['A','G','C','T']
    # empty list
    barcode_list = []
    # while loop to generate random barcode
    while count4 < 26:
        # add each random base pair to list
        barcode_list.append(random.choice(base_pair_list))
        # increase count
        count4 += 1
    # turn list into string
    barcode_string = ('').join(barcode_list)
    # make empty list for quality values
    qualityvalues3 = []
    # make quality value for each barcode base pair
    for basepair in range(0,len(barcode_string)):
        qualityvalues3.append('b')
    # turn that list into a string
    qualityvalues4 = ('').join(qualityvalues3)
    # for loop going thru the sequence fraction list
    for seqfrac in range(0, len(seq_frac_list)):
        # conditional if the sequence fraction is less than or equal to j
        if cumsum[seqfrac] >= j:
            # get the corresponding dictionary value (sequence) for that sequence fraction in list
            dict_value = seq_frac_dict[seq_frac_list[seqfrac]][0]
            # get random integer for big fragment
            i = random.randrange(0,(len(seq_frac_dict[seq_frac_list[seqfrac]][0])) - big_frag_length)
            # get a big fragment of that sequence
            big_fragment = dict_value[i:i+(big_frag_length-1)]
            # set second and third counters
            count2 = 0
            count3 = 1
            # while loop that goes while count is less than number of small fragments
            while count2 < small_frag_number:
                # random integer for small fragments
                x = random.randrange(0,(big_frag_length - small_frag_length))
                # get small fragment from big fragment
                small_fragment = big_fragment[x:x+(small_frag_length)]
                # make empty list for quality values
                qualityvalues = []
                # make quality value for each base pair
                for base_pair in range(0,len(small_fragment)):
                    qualityvalues.append('b')
                # turn that list into a string
                qualityvalues2 = ('').join(qualityvalues)
                # increase counter by 1
                count2 += 1
                #print('Frag: Sequence%d_%d_%d/1 \n %s \n + \n %s' %(seq_frac_dict[seq_frac_list[seqfrac]][1], count3, x+i, small_fragment, qualityvalues2))
                outputfile1.write('@Sequence_%d_%d_%d/1\n%s\n+\n%s\n' %(seq_frac_dict[seq_frac_list[seqfrac]][1], count3, x+i, small_fragment, qualityvalues2))
                print('@Sequence_%d_%d_%d/2\n%s\n+\n%s\n' %(seq_frac_dict[seq_frac_list[seqfrac]][1], count3, x+i, barcode_string, qualityvalues4))
                outputfile2.write('@Sequence_%d_%d_%d/2\n%s\n+\n%s\n' %(seq_frac_dict[seq_frac_list[seqfrac]][1], count3, x+i, barcode_string, qualityvalues4))
                count3 += 1
            # break statement
            break
    # increase counter by 1
    count += 1
# close files
outputfile1.close()
outputfile2.close()
