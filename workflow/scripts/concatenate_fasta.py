import os
import sys

def concatenate_fasta(input_files, output_file):

    # filter for .fa files
    fasta_files = [f for f in input_files if f.endswith('.fna')]

    # raise an error if there are no .fa files
    if not fasta_files:
        raise FileNotFoundError("No FASTA files found in ./genomes/")

    # concatenate each FASTA to output
    with open(output_file, 'w') as outfile:
        for file_path in fasta_files:
            with open(file_path, 'r') as infile:
                for line in infile:
                    if line.strip():  # skip blank lines
                        outfile.write(line)

if __name__ == "__main__":

    input_files = sys.argv[2:]
    output_file = sys.argv[1]

    concatenate_fasta(input_files, output_file)
