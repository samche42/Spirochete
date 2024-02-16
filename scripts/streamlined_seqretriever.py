import sys
import os
from glob import glob
from Bio import SeqIO

wanted =[]

input_directory = sys.argv[1]
#directory is path to multifastas
nucleotide_file = sys.argv[2]
#nucleotide_file is path the conactenated nucleotide file corresponding to proteins

multifasta_file_list = [os.path.abspath(fn) for fn in glob(input_directory+'/*.aln')]

for file in multifasta_file_list:
        records = []
        wanted = []
        output_file = (os.path.abspath(file))[:-4]+"_nuc.fasta"
#       outputfile named the same as input OG file
        wanted = list(r.id for r in SeqIO.parse(file, "fasta"))
        print("A total of "+str(len(wanted))+" OGs being retrieved for "+str(file)+". Please confirm this is as expected")
#       lists all headers in OG file
        for item in wanted:
                for r in SeqIO.parse(nucleotide_file, "fasta"):
                        if r.id == item:
                                records.append(r)
#       pulls all fastas with headers that correspond to wanted list
        SeqIO.write(records, output_file, "fasta")
#       writes pulled fastas to output files
