#!/usr/bin/env python3

import os
from Bio import SeqIO

input_directory = "./"

first_OG_file = [file for file in os.listdir(input_directory) if file.endswith(".aln")][0]

headers = []
with open(first_OG_file, "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        headers.append(record.description)

genome_list = [s[:s.rfind("_")] + "_" for s in headers]

all_files = [file for file in os.listdir(input_directory) if file.endswith(".aln")]
nucleotide_files = [entry for entry in all_files if "edited_nuc_without_stop" in entry]
protein_files = [entry for entry in all_files if "edited_nuc_without_stop" not in entry]

def add_base_filename(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    modified_lines = []
    base_filename = filename.split('.')[0]
    for line in lines:
        if line.startswith('>'):
            modified_lines.append(f'>{base_filename}_{line[1:]}')
        else:
            modified_lines.append(line)
    with open(filename, 'w') as file:
        file.writelines(modified_lines)

for file in all_files:
	add_base_filename(file)

def sort_seqs(input_file,output_file):
	sequence_dict = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
	sorted_dict = dict(sorted(sequence_dict.items()))
	with open(output_file, "w") as output_handle:
		SeqIO.write(sorted_dict.values(), output_handle, "fasta")

#Create genome files
for genome in genome_list:
	nucleotide_output_file = genome+"multifasta.ffn"
	protein_output_file = genome+"multifasta.faa"
	nucleotide_temp = genome+"_temp.ffn"
	protein_temp = genome+"_temp.faa"
	nucleotide_records = []
	for OG_group in nucleotide_files:
		for r in SeqIO.parse(OG_group, "fasta"):
			if genome in str(r.id):
				nucleotide_records.append(r)
		SeqIO.write(nucleotide_records, nucleotide_temp, "fasta")
	protein_records = []
	for OG_group in protein_files:
		for r in SeqIO.parse(OG_group, "fasta"):
			if genome in str(r.id):
				protein_records.append(r)
		SeqIO.write(protein_records, protein_temp, "fasta")
	#Reorder both files alphabetically so that entries are in the same order
	sort_seqs(nucleotide_temp,nucleotide_output_file)
	sort_seqs(protein_temp,protein_output_file)
	os.remove(nucleotide_temp)
	os.remove(protein_temp)
