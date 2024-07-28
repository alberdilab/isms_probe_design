import sys
import pandas as pd
from Bio import SeqIO

def update_gtf_file(input_gtf, mapping_file, output_gtf):
    header_map = {}
    with open(mapping_file, 'r') as mapfile:
        for line in mapfile.readlines()[1:]:  # Skip the header line
            old_header, new_header = line.strip().split('\t')
            header_map[old_header] = new_header

    with open(input_gtf, 'r') as infile, open(output_gtf, 'w') as outfile:
        for line in infile:
            for old_header, new_header in header_map.items():
                if old_header in line:
                    line = line.replace(old_header, new_header)
            outfile.write(line)

if __name__ == "__main__":
    input_gtf = sys.argv[1]
    mapping_file = sys.argv[2]
    output_gtf = sys.argv[3]
    update_fasta_headers(input_gtf, mapping_file, output_gtf)
