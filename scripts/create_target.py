######
# InSituMicrobeSeq probe design pipeline
# Martí Dalmases & Antton Alberdi
# 2024/08/01
# Description: target finder from gtf and fasta files.
######

import argparse
from Bio import SeqIO
import pandas as pd
import os
import glob

def read_gtf(gtf_file):
    # reads the .gtf into a pandas df
    return pd.read_csv(gtf_file, sep='\t', comment='#', header=None,
                       names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

def filter_gtf(df, annotation, exon_only):
    # filters the df based on the annotation and exon flags
    filtered_df = df[df['attribute'].str.contains(annotation)]
    if exon_only:
        filtered_df = filtered_df[filtered_df['feature'] == 'exon']
    return filtered_df

def process_region_mode(gtf_files, annotation, exon_only, output_file):
    # aggregate all filtered dfs
    filtered_dfs = []

    for gtf_file in gtf_files:
        df = read_gtf(gtf_file)
        filtered_df = filter_gtf(df, annotation, exon_only)
        if not filtered_df.empty:
            filtered_dfs.append(filtered_df)

    # concatenate all filtered dfs
    if not filtered_dfs:

        raise LookupError("No matches found for annotation:", annotation)

    else:

        concatenated_df = pd.concat(filtered_dfs)

        # sort by 'seqname' and 'start'
        concatenated_df.sort_values(by=['seqname', 'start'], inplace=True)

        # output to file
        concatenated_df.to_csv(output_file, sep='\t', index=False, header=False)


def process_genome_mode(fasta_file, output_file):
    # if in genome mode, output all provided contigs as targets
    with open(output_file, 'w') as gtf:
        for record in SeqIO.parse(fasta_file, "fasta"):
            contig_id = record.id
            contig_length = len(record.seq)
            gtf_entry = f"{contig_id}\tisms_probe_design\ttarget\t1\t{contig_length}\t.\t+\t.\t.\n"
            gtf.write(gtf_entry)

def main():
    parser = argparse.ArgumentParser(description="Filter GTF files or process genome based on mode",
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-m', '--mode', choices=['region', 'genome'], required=True,
                        help='Mode of operation: \n'
                             'region - Filter GTF files based on annotations\n'
                             'genome - Process a single FASTA file')

    # region mode arguments
    parser.add_argument('-g', '--gtf', nargs='+',
                        help='List of input GTF files or folders containing GTF files (required in region mode)')
    parser.add_argument('-a', '--annotation',
                        help='Annotation to filter for (required in region mode)')
    parser.add_argument('-e', '--exon', action='store_true',
                        help='Filter to include only exons (optional in region mode)')

    # genome mode arguments
    parser.add_argument('-f', '--fasta',
                        help='Input FASTA file (required in genome mode)')

    parser.add_argument('-o', '--output', required=True,
                        help='Output file to write the results')

    args = parser.parse_args()

    if not args.mode:
        parser.print_help()
        parser.exit("Error: --mode is required. Please specify the mode of operation (region or genome).")

    if args.mode == 'region':
        if not args.gtf or not args.annotation:
            parser.error("--gtf and --annotation are required in region mode")

        # Collect all GTF files from the provided paths
        gtf_files = []
        for path in args.gtf:
            if os.path.isdir(path):
                gtf_files.extend(glob.glob(os.path.join(path, '*.gtf')))
            elif os.path.isfile(path) and path.endswith('.gtf'):
                gtf_files.append(path)

        if not gtf_files:
            raise FileNotFoundError("No GTF files found in the provided paths")

        process_region_mode(gtf_files, args.annotation, args.exon, args.output)

    elif args.mode == 'genome':
        if not args.fasta:
            parser.error("--fasta is required in genome mode")

        if not os.path.isfile(args.fasta):
            raise FileNotFoundError(f"FASTA file not found: {args.fasta}")

        process_genome_mode(args.fasta, args.output)

if __name__ == "__main__":
    main()
