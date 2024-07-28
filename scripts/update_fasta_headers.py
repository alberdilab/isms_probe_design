import sys
import pandas as pd
from Bio import SeqIO

def update_fasta_headers(input_fasta, mapping_file, output_fasta):
    header_map = {}
    with open(mapping_file, 'r') as mapfile:
        for line in mapfile.readlines()[1:]:  # Skip the header line
            old_header, new_header = line.strip().split('\t')
            header_map[old_header] = new_header

    sequences = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        record.id = header_map[record.id]
        record.description = header_map[record.id]
        sequences.append(record)
    
    SeqIO.write(sequences, output_fasta, "fasta")

if __name__ == "__main__":
    input_fasta = sys.argv[1]
    mapping_file = sys.argv[2]
    output_fasta = sys.argv[3]
    update_fasta_headers(fasta_file, mapping_file,output_fasta)
