# Name: Daniel Hwang (djhwang)
# Group: Sidhant Karamchandani (skaramch)
# Desc: This program will bin the input file of two fastq file of which one consist
#       of the parsed sequence and another fastq for the barcode.
# Input:
#   read_1.fastq and read 2.fastq
# Processing:
#   1. Import SeqIO and Seq to obtain fastq information. 
#   2. Create dictionary with appropriate key corresponding value such as sequence
#   3. Write n files where n is the number of barcode in read_2 file
# Output:
#   n fastq files (n is number of barcodes) 
# barcode_bin.py


###############################################################################
################################  PSEUDO CODE  ################################
###############################################################################

# 1. import SeqIO, Seq
# 2. Create an empty dictionary and list for both file
# 3. for each file, assign dictionary key with header of the file
# 4. append the sequence to a list and then assign it to appropriate key
# 5. By using the header from the parsed sequence, it is able to sort according
#    to the barcodes

##############################################################################
################################ MAIN PROGRAM ################################
##############################################################################
from Bio import SeqIO
from Bio import Seq

# create empty dictionary and list to append 
bc_dict = {}
bc_ls = [0]
seq_dict = {}
hd_ls = []

#input_one = raw_input("Please enter the fastq file containing the short reads: ")
#input_two = raw_input("Please enter the fastq file containing the barcode reads: ")
# Iteration loop going through barcode sequences (could be user input but for now it is hardcode)
for read_2 in SeqIO.parse("case1read2.fastq", "fastq"):
    #assign the variable with the sequence id
    hd2 =str(read_2.id[:-1])
    #conditional statement where if it will only append the barcode if the sequence are different
    if (str(read_2.seq)) != bc_ls[-1]:
        bc_ls.append(str(read_2.seq))
    #assign the parsed headers as the key of dictionary corresponding to the barcode sequence
    bc_dict[hd2] = (str(read_2.seq))
#remove the first element of the barcode list which is 0
bc_ls = bc_ls[1:]

#Iteration loop giong through the processed sequences(could be user input but for now it is hardcode)
for read_1 in SeqIO.parse("case1read1.fastq", "fastq"):
    #assign the variable with the sequence id
    hd1 =str(read_1.id[:-1])
    #append the header into the list
    hd_ls.append(hd1)
    #assign the parsed headers as the key of dictionary corresponding to the barcode sequence
    seq_dict[hd1] = (str(read_1.seq))

for j in range(0, len(bc_ls)):
    #take the barcode sequence as the filename
    file_name = "Sequence_" + str(j) + ".fastq"
    #binning the different sequence together with the same barcode and then write it to file
    with open(file_name, "w") as fh:
        for head in hd_ls:
            #conditional statement where it will only write to file if the header match
            if bc_ls[j] == bc_dict[head]:
                #Takes the header from parsed sequence and write it on file
                fh.write("@" + head + "1\n")
                #write the corresponding sequence
                fh.write(seq_dict[head] + "\n")
                fh.write("+\n")
                #quality value
                length = len(seq_dict[head])
                i = 0
                while i < length:
                    fh.write("b")
                    i+=1
                fh.write("\n")




