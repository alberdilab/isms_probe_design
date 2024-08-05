######
# InSituMicrobeSeq probe design pipeline
# Martí Dalmases & Antton Alberdi
# 2024/08/02
# Description: Rank probes based on the on/off-scores, Tm and kmer count
######

import pandas as pd
import sys

# Define weights for each column
weight_on_target_score = 0.5
weight_off_target_score = -1
weight_Tm = 0.2
weight_kmer_count = -0.2

probes_file = sys.argv[1]
output_file = sys.argv[2]

# Define the column names manually since the file has no header
column_names = [
    'chrom', 'start', 'stop', 'parent', 'Tm',
    'on_target_score', 'off_target_score', 'repeat',
    'kmer_count', 'strand'
]

# Define a function to calculate the rank
def calculate_rank(row):
    return (weight_on_target_score * float(row['on_target_score']) +
            weight_off_target_score * float(row['off_target_score']) +
            weight_Tm * float(row['Tm']) +
            weight_kmer_count * float(row['kmer_count']))

def main():
    # Read the TSV file into a DataFrame without a header
    df = pd.read_csv(probes_file, sep='\t', header=None, names=column_names)

    # Add a new column for rank calculation
    df['rank_calculation'] = None

    # Apply the function to each row
    df['rank_calculation'] = df.apply(calculate_rank, axis=1)

    # Sort the DataFrame by the descending rank_calculation column
    df = df.sort_values(by='rank_calculation', ascending=False)

    # Remove the rank_calculation column
    #df = df.drop(columns=['rank_calculation'])

    # Write the DataFrame back to a TSV file
    df.to_csv(output_file, sep='\t', index=False, header=False)

if __name__ == "__main__":
    main()
