######
# InSituMicrobeSeq probe design pipeline
# MartÃ­ Dalmases & Antton Alberdi
# 2024/08/16
# Description: Select probes for experimentation based on a density and distribution ranking.
######

import os
import sys
import argparse
import glob
import numpy as np
from scipy.stats import gaussian_kde
import multiprocessing as mp
import pandas as pd


num_selected_probes = 0


## WEIGHTS ##

W_ON_TARGET = 0.5
W_OFF_TARGET = -1
W_TM = 0.2
W_KMER = -0.2
W_DENSITY = 0.3
W_SPARCE = 0.3

BANDWIDTH = 40

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


## Multiprocessing Functions ##

def process_genome(genome_args):
    # extract arguments
    group_tuple, priority_df = genome_args
    genome_name, group = group_tuple

    # calculate density
    genome_df = calculate_density(group, priority_df['tag'])
    print("Calculated density for " + genome_name)

    return genome_df

def process_target(target_args):
    # extract arguments
    target_tuple, priority_map, gtf_files, num_selected_probes = target_args
    target_name, group = target_tuple

    # calculate the number of probes proportional to the priority map
    if target_name not in priority_map:
        print(f"ERROR: The target {target_name} does not exist in the priority map.")
        return None

    target_probe_count = int(priority_map[target_name] * total_probe_count)

    if group.shape[0] < target_probe_count:
        print(f"ERROR: The target {target_name} has too few probes. Only {group.shape[0]} will be selected, regardless of quality.")
        target_probe_count = group.shape[0]  # adjust to select available probes

    else:
        print(f"INFO: {target_probe_count} probes will be selected for the target {target_name}")

    # calculate sparsity
    target_df = calculate_sparsity(group, target_probe_count, gtf_files)

    if target_df is None:
        print(f"ERROR: Sparsity calculation failed for target {target_name}.")
        return None

    # final weighted calculation
    target_df = assign_final_value(target_df)

    #sSort by 'rank' and select the top probes
    target_df = target_df.sort_values(by=['rank'], ascending=False)
    selected_probes = target_df.head(target_probe_count)

    return selected_probes



## run density ##

def calculate_density(genome_df, priority_list):

    # assign priority ranks to the probes
    priority_map = {tag: i for i, tag in enumerate(priority_list)}
    genome_df['priority_rank'] = genome_df['tag'].apply(lambda x: priority_map.get(x, len(priority_list)))

    # sort by priority rank and start position
    genome_df = genome_df.sort_values(by=['priority_rank', 'start']).reset_index(drop=True)

    # extract the start point for KDE calculations
    positions = genome_df['start'].values

    # KDE on the probe positions
    kde = gaussian_kde(positions, bw_method=BANDWIDTH / np.std(positions))

    # estimate the density for each probe's start positionProcessing complete. 1000 probes were considered in total.
    densities = kde(positions)

    # assign density values adjusted by priority weight
    genome_df['density'] = 0.0
    for idx, probe in genome_df.iterrows():
        priority_rank = probe['priority_rank']

        if priority_rank == genome_df['priority_rank'].min():
            priority_weight = 1.0  # highest priority gets full weight
        else:
            priority_weight = 1 / (priority_rank + 1)  # lower priority gets reduced weight

        # scale the density by the priority weight
        density_value = densities[idx] * priority_weight
        genome_df.at[idx, 'density'] = density_value

    genome_df = genome_df.drop(columns=['priority_rank'])

    return genome_df



## run distribution ##

def calculate_sparsity(target_df, probe_count, gtf_data):
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

def main(input_directory, priority_file, gtf_directory, output_directory, total_probe_count, num_selected_probes):
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
    genome_groups = list(combined_df.groupby('genome'))

    genome_args = [(group_tuple, priority_df) for group_tuple in genome_groups]

    # process each genome for density
    with mp.Pool(mp.cpu_count()) as pool:
        genome_results = pool.map(process_genome, genome_args)

    pool.close()
    pool.join()

    density_df = pd.concat(genome_results, ignore_index=True)
    print("CALCULATED DENSITY")

    # sort by target
    density_df = density_df.sort_values(by=['tag'])

    # separate by 'target'
    target_groups = density_df.groupby('tag')

    target_args = [(target_tuple, priority_map, gtf_files, num_selected_probes) for target_tuple in target_groups]

    # process each target distribution
    with mp.Pool(mp.cpu_count()) as pool:
        target_results = pool.map(process_target, target_args)

    pool.close()
    pool.join()

    for target in target_results:
        num_selected_probes += target.shape[0]

    final_df = pd.concat(target_results, ignore_index=True)
    print("CALCULATED SPARSITY")

    # save final result
    output_file = os.path.join(output_directory, f"final_probes.tsv")
    final_df.to_csv(output_file, sep="\t", index=False)

    print(f"Processing complete. {num_selected_probes} probes were considered in total.")

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

    main(input_directory, priority_file, gtf_directory, output_directory, total_probe_count, num_selected_probes)
