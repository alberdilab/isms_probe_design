#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 16:47:23 2019 and Wed Jun 26 10:01:32 2019

@author: hershe

A combined script to calculate k-mer frequencies and append the max k-mer count
for a given probe to the BED row of the probe in a Snakemake pipeline.

Based on kmerFilter.py from OligoMiner and the max k-mer appending script.

"""

import pandas as pd
import numpy as np
import subprocess
import tempfile
import sys
import os

def run_frequency(input_bed_file, output_file, jf_file):

    # Determine the stem of the input filename.
    file_name = os.path.basename(input_bed_file)

    # define mer length as a variable
    mer_length = 18

    # read in BED file
    with open(input_bed_file, 'r') as f:
      file_read = [line.strip() for line in f]

    # Make list to hold the sequences of the imported probes.
    seq_list = []

    # Parse out sequence info from each probe in the input .bed
    # and add to the newly created list.
    for i in range(0, len(file_read), 1):
      seq_list.append('>' + '\n' + file_read[i].split('\t')[3])

    # Generate a randomized number to add to the temporary input/output files
    # required to run Jellyfish.
    random_int = np.random.randint(0, 1000000 + 1)

    # write jf files into temporary directory on cluster node
    with tempfile.TemporaryDirectory() as tmpdir:

        # Create a temporary input .fa file for jellyfish.
        jf_fasta = open('%s/%s_%d_%d_temp.fa' \
                       % (tmpdir, file_name, mer_length, random_int), 'w')

        # Write probe sequences to the .fa file.
        jf_fasta.write('\n'.join(seq_list))
        jf_fasta.close()

        # Create index of kmer positions within each probe to reference later.
        index_list = []
        for i in range(0, len(file_read), 1):
          probe_length = len(seq_list[i]) - 2
          mer_windows = probe_length - mer_length + 1
          for j in range(0, mer_windows, 1):
              index_list.append(file_read[i])

        # Call jellyfish and count kmers in the probe sequences.
        subprocess.call(['jellyfish', 'query', '%s' % jf_file, '-s',
                         '%s/%s_%d_%d_temp.fa' \
                         % (tmpdir, file_name, mer_length, random_int),
                         '-o', '%s/%s_%d_%d_temp.txt' \
                         % (tmpdir, file_name, mer_length, random_int)],
                        stderr=None, shell=False)

        # Open jellyfish output.
        with open('%s/%s_%d_%d_temp.txt' \
            % (tmpdir, file_name, mer_length, random_int), 'r') as jf_data:

            # read in jf kmer counts as string per line
            jf_read = [line.strip() for line in jf_data]

        # create a data frame to store jf results
        jf_result = pd.DataFrame()

        # create column with the index of the probe for each kmer position
        jf_result['probe'] = index_list

        # parse out counts from the jf file read in
        jf_counts = [int(line.split(' ')[1]) for line in jf_read]

        # create column with result (kmer and its count)
        jf_result['count'] = jf_counts

        # calculate max kmer count for each probe
        jf_result = jf_result.groupby('probe')['count'].max().reset_index()

    # Load probes with metadata
    df = pd.read_csv(input_bed_file, sep='\t', header=None)
    df.columns = [
        'chrom',
        'start',
        'stop',
        'parent',
        'Tm',
        'on_target_score',
        'off_target_score',
        'repeat',
    ]

    # Map the max kmer counts to the probes
    max_kmers_dict = dict(zip(jf_result['probe'], jf_result['count']))
    df['max_kmers'] = df.apply(lambda row: max_kmers_dict.get('\t'.join(row.astype(str)), 0), axis=1)

    # Add hard-coded + strand (blockparse output)
    df['strand'] = '+'

    # Sort probes on chromosomal start coordinate
    df = df.sort_values(['chrom', 'start'])

    # Save to disk
    df.to_csv(output_file, sep='\t', index=False, header=False, float_format='%.3f')

if __name__ == "__main__":

    input_file = sys.argv[1]
    jf_file = sys.argv[2]
    output_file = sys.argv[3]

    run_frequency(input_file, output_file, jf_file)
