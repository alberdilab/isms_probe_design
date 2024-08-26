######
# InSituMicrobeSeq probe design pipeline
# Martí Dalmases & Antton Alberdi
# 2024/08/16
# Description: Select probes for experimentation based on a density and distribution ranking.
######

import os
import sys
import argparse
import glob
import numpy as np
from intervaltree import Interval, IntervalTree
import pandas as pd

## WEIGHTS ##

W_ON_TARGET = 0.5
W_OFF_TARGET = -1
W_TM = 0.2
W_KMER = -0.2
W_DENSITY = 0.3
W_SPARCE = 0.3



## Import files ##

def concat_tsvs(directory):
    all_files = glob.glob(os.path.join(directory, "*.tsv"))
    df_list = []

    for file in all_files:
        # read TSV
        df = pd.read_csv(file, sep="\t", header=None, names=['genome','start','stop','parent','Tm','on_target_score','off_target_score','repeat','max_kmer','strand'])
        # extract filename as tag
        file_name = os.path.basename(file).replace(".tsv", "")
        df['tag'] = file_name
        print("LOADED " + file_name + " FILE")
        # append to the list of dataframes
        df_list.append(df)

    # concatenate
    combined_df = pd.concat(df_list, ignore_index=True)
    return combined_df

def load_priority_file(priority_file):
    priority_df = pd.read_csv(priority_file, sep='\t', header=None, names=['tag', 'fraction'])
    return priority_df

def load_gtf_targets(directory):
    gtf_files = glob.glob(os.path.join(directory, "*.gtf"))
    gtf_data = {}
    for file in gtf_files:
        file_name = os.path.basename(file).replace(".gtf", "")
        with open(file, 'r') as f:
            lines = [line.strip().split('\t') for line in f.readlines() if not line.startswith("#")]
            gtf_data[file_name] = lines  # store lines for each GTF file
    return gtf_data


## run density ##

def build_interval_tree(genome_df):
    tree = IntervalTree()
    for idx, probe in genome_df.iterrows():
        tree[probe['start']:probe['stop'] + 1] = idx  # store the index of the probe in the interval
    print("BUILT TREE")
    return tree

def calculate_density(genome_df, priority_list):

    genome_df['density'] = 0.0
    genome_df['processed'] = False
    priority_map = {tag: i for i, tag in enumerate(priority_list)}

    # Add a priority rank based on the priority list
    genome_df['priority_rank'] = genome_df['tag'].apply(lambda x: priority_map.get(x, len(priority_list)))
    genome_df = genome_df.sort_values(by=['priority_rank', 'start'])

    # build the interval tree
    tree = build_interval_tree(genome_df)

    # calculate density for each probe
    for idx, probe in genome_df.iterrows():
        if not probe['processed']:
            # query the interval tree for overlaps
            overlapping_intervals = tree[probe['start']:probe['stop'] + 1]
            overlapping_indices = [iv.data for iv in overlapping_intervals]

            # extract overlapping probes
            overlapping_probes = genome_df.loc[overlapping_indices]

            # calculate local density
            local_density = len(overlapping_probes)

            # calculate local highest priority
            highest_priority_rank = overlapping_probes['priority_rank'].min()

            # assign density values with priority weight
            for i, (overlap_idx, overlap_probe) in enumerate(overlapping_probes.iterrows()):
                priority_rank = overlap_probe['priority_rank']

                if priority_rank == highest_priority_rank:
                    priority_weight = 1  # highest priority probe gets a value 1
                else:
                    priority_weight = (1 / (priority_rank + 1))  # lower priority probes get a reduced weight

                density_value = (1 / (local_density + 1)) * priority_weight
                genome_df.at[overlap_idx, 'density'] = max(genome_df.at[overlap_idx, 'density'], density_value)
                genome_df.at[overlap_idx, 'processed'] = True

    genome_df = genome_df.drop(columns=['processed', 'priority_rank'])

    return genome_df



## run distribution ##

def calculate_sparcity(target_df, probe_count, gtf_data):
    # open related GTF file based on tag or target
    tag = target_df['tag'].unique()[0]
    if tag not in gtf_data:
        raise ValueError(f"GTF data for {tag} not found. Please make sure all input filenames and target names are correct (case-sensitive)")

    # extract relevant GTF lines
    gtf_lines = gtf_data[tag]

    # calculate the total size of the target
    total_target_size = sum([int(line[4]) - int(line[3]) + 1 for line in gtf_lines])

    # sort by contig
    target_df = target_df.sort_values(by=['genome', 'start', 'stop'])

    # store probe assignments for each contig
    all_contig_dfs = []

    # calculate contig size relative to the total
    for contig, contig_df in target_df.groupby('genome'): # note the 'genome' format is genome1_contig4 for which this actually groups by contigs
        # calculate its lenght from the correspondent GTF line
        contig_gtf = next((line for line in gtf_lines if line[0] == contig), None)

        if not contig_gtf:
            raise ValueError(f"GTF data for {contig} not found in {gtf_lines}. Please make sure all input filenames and target names are correct (case-sensitive)")

        start_pos = int(contig_gtf[3])  # start position from GTF
        stop_pos = int(contig_gtf[4])   # stop position from GTF

        contig_length = stop_pos - start_pos + 1

        # calculate assigned contig probes
        print(type(contig_length), type(total_target_size), type(probe_count))

        num_probes = int(round((contig_length / total_target_size) * probe_count))

        # divide the contig into equal parts based on the number of probes
        split_points = np.linspace(start_pos, stop_pos, num_probes + 1)


        for idx, probe in contig_df.iterrows():
            # find the closest split point for the current probe
            closest_point = min(split_points, key=lambda x: abs(x - probe['start']))

            # apply a logarithmic function based on the distance to the closest split point
            probe_value = -np.log10(abs(probe['start'] - closest_point) + 1)

            # assign distribution value for the probe
            contig_df.at[idx, 'distribution'] = probe_value

        # append the processed contig dataframe to the list
        all_contig_dfs.append(contig_df)

    # combine all contig dataframes back into the target_df
    target_df = pd.concat(all_contig_dfs)

    return target_df

## FINAL WEIGHTED CALC ##

def assign_final_value(target_df):
    target_df['rank'] = (target_df['on_target_score'] * W_ON_TARGET +
                         target_df['off_target_score'] * W_OFF_TARGET +
                         target_df['Tm'] * W_TM +
                         target_df['max_kmer'] * W_KMER +
                         target_df['density'] * W_DENSITY +
                         target_df['distribution'] * W_SPARCE)
    return target_df

def main(input_directory, priority_file, gtf_directory, output_directory, total_probe_count):
    # load priority list
    priority_df = load_priority_file(priority_file)
    priority_map = {row['tag']: row['fraction'] for _, row in priority_df.iterrows()}

    # load target gtfs
    gtf_files = load_gtf_targets(gtf_directory)

    # tag and combine all TSV into one DataFrame
    combined_df = concat_tsvs(input_directory)

    # sort by 'genome'
    combined_df = combined_df.sort_values(by=['genome', 'start', 'stop'])

    # separate by 'genome'
    genome_groups = combined_df.groupby('genome')

    # process each genome for density
    genome_dict = {}
    for genome_name, group in genome_groups:
        genome_df = calculate_density(group, priority_df['tag'])
        genome_dict[genome_name] = genome_df
        print("Calculating density for " + genome_name)

    density_df = pd.concat([df for df in genome_dict.values()])
    print("DENSITY CALCULATED")

    # sort by target
    density_df = density_df.sort_values(by=['tag'])

    # separate by target
    target_groups = density_df.groupby('tag')

    # process each target distribution
    target_dict = {}
    for target_name, group in target_groups:
        target_probe_count = int(priority_map[target_name] * total_probe_count)
        target_df = calculate_sparcity(group, target_probe_count, gtf_files)
        target_df = assign_final_value(target_df)
        target_df = target_df.sort_values(by=['rank'], ascending=False)
        selected_probes = target_df.head(target_probe_count)
        target_dict[target_name] = selected_probes
        print("Calculated sparcity for" + target_name)
    #
    final_df = pd.concat([df for df in target_dict.values()])
    print("CALCULATED SPARCITY")

    # save final result
    output_file = os.path.join(output_directory, f"final_probes.tsv")
    final_df.to_csv(output_file, sep="\t", index=False)

    print(f"Processing complete. {total_probe_count} probes were considered in total.")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Select probes on a ranking based system.",
                                     formatter_class=argparse.RawTextHelpFormatter)

    # define commandline arguments
    parser.add_argument("-i", "--input", required=True, help="Directory containing TSV files to process")
    parser.add_argument("-p", "--priority", required=True, help="Priority file .tsv (contains target name in priority order in the first column and fraction of total probes dedicated in the second column)")
    parser.add_argument("-o", "--output", required=True, help="Output directory for processed files")
    parser.add_argument("-g", "--gtf_directory", required=True, help="Directory containing GTF target files")
    parser.add_argument("-c", "--count", type=int, required=True, help="Number of final probes to consider")

    args = parser.parse_args()

    input_directory = args.input
    priority_file = args.priority
    output_directory = args.output
    total_probe_count = args.count
    gtf_directory = args.gtf_directory


    if not args.count:
        parser.print_help()
        parser.exit("Error: --count is required. Please specify the number of final probes to consider.")

    if not os.path.isdir(input_directory):
        print(f"Error: {input_directory} is not a valid directory.")
        sys.exit(1)

    if not os.path.isfile(priority_file):
        print(f"Error: {priority_file} is not a valid file.")
        sys.exit(1)

    if not os.path.isdir(output_directory):
        print(f"Error: {output_directory} is not a valid directory.")
        sys.exit(1)

    main(input_directory, priority_file, gtf_directory, output_directory, total_probe_count)

