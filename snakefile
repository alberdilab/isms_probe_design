######
# InSituMicrobeSeq probe design pipeline
# Martí Dalmases & Antton Alberdi
# 2024/07/27
# Description: draft backbone pipeline for probe design. Optimised for slurm-based HPC computation.
# Modified from PaintSHOP pipeline: https://github.com/beliveau-lab/PaintSHOP_pipeline
######

# 1) Clone this repository 'git clone https://github.com/alberdilab/isms_probe_design.git'
# 2) Store the genome sequences in the folder 'genomes' in the working directory with extension .fa, and without any "." or "/" in the file name besides the extension .fa.
# 3) Store the target-defining GTF files in the 'targets' directory with extension .gtf, and without any "." or "/" in the file name besides the extension .gtf.
# 4) Create a screen session and enter into the repository directory
# screen -S isms_probe_design
# cd isms_probe_design
# 5) Launch the snakemake using the following code:
# module purge && module load snakemake/7.20.0 mamba/1.3.1
# snakemake -j 20 --cluster 'sbatch -o logs/{params.jobname}-slurm-%j.out --mem {resources.mem_gb}G --time {resources.time} -c {threads} --job-name={params.jobname} -v'   --use-conda --conda-frontend mamba --conda-prefix conda --latency-wait 600
# 6) Let snakemake to manage the slurm jobs until the 

# List genome and target wildcards
genomes, = glob_wildcards("genomes/{genome}.fa")
targets, = glob_wildcards("targets/{target}.gtf")

# Expand target files
rule all:
    input:
        expand("probes/{target}.tsv", target=targets)

# Merge all input fasta files (bacterial genomes) into a single file.
rule concatenate_fasta:
    input:
        expand("genomes/{genome}.fa", genome=genomes)
    output:
        "pipeline/01_input/allgenomes.fa"
    params:
        jobname="allgenomes.pf"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        python scripts/concatenate_fasta.py {output} {input}
        """

# Create a mapping file of contig headers, with modified new headers in the case of duplications
rule unique_headers:
    input:
        "pipeline/01_input/allgenomes.fa"
    output:
        "pipeline/02_renamed/headers.tsv"
    params:
        jobname="allgenomes.uh"
    threads:
        1
    resources:
        mem_gb=8,
        time=5
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python scripts/unique_headers.py {input} {output}
        """

# If necessary, modify headers to avoid duplicated identities that would affect downstream operations
rule unique_headers_fasta:
    input:
        fasta="pipeline/01_input/allgenomes.fa",
        headers="pipeline/02_renamed/headers.tsv"
    output:
        "pipeline/02_renamed/allgenomes.fa"
    params:
        jobname="allgenomes.uf"
    threads:
        1
    resources:
        mem_gb=8,
        time=5
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python scripts/update_fasta_headers.py {input.fasta} {input.headers} {output}
        """

# Bowtie index renamed fasta file for downstream probe mapping
rule index_fasta:
    input:
        "pipeline/02_renamed/allgenomes.fa"
    output:
        "pipeline/02_renamed/allgenomes.rev.1.bt2"
    params:
        base="pipeline/02_renamed/allgenomes",
        jobname="allgenomes.in"
    threads:
        1
    resources:
        mem_gb=8,
        time=60
    conda:
        "envs/biopython.yaml"
    shell:
        """
        bowtie2-build {input} {params.base}
        """

# Jellyfish count renamed fasta file for downstream kmer count
rule build_jellyfish:
    input:
        "pipeline/02_renamed/allgenomes.fa"
    output:
        "pipeline/02_renamed/allgenomes.jf"
    conda:
        "envs/biopython.yaml"
    params:
        mfree='30G',
        h_rt='3:0:0'
    shell:
        """
        jellyfish count -m 18 -s 3300M -o {output} --out-counter-len 1 -L 2 {input}
        """

# Update the contig id-s in the gtf files using the header mapping file to avoid duplications.
rule unique_ids_gtf:
    input:
        gtf="targets/{target}.gtf",
        headers="pipeline/02_renamed/headers.tsv"
    output:
        "pipeline/02_renamed/{target}.gtf"
    params:
        jobname="{target}.pg"
    threads:
        1
    resources:
        mem_gb=8,
        time=5
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python scripts/update_gtf_ids.py {input.gtf} {input.headers} {output}
        """

# Based on the GTF information, extract target regions from the original fasta files to design probes.
rule extract_targets:
    input:
        fasta="pipeline/02_renamed/allgenomes.fa",
        gtf="pipeline/02_renamed/{target}.gtf"
    output:
        "pipeline/03_extract/{target}.fa"
    params:
        jobname="{target}.ex"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        module load bedtools/2.30.0
        bedtools getfasta -fi {input.fasta} -bed {input.gtf} > {output}
        """

# Generate the initial set of probes to be tested for suitability.
rule generate_probes:
    input:
        "pipeline/03_extract/{target}.fa"
    output:
        "pipeline/04_probes/{target}.fastq"
    params:
        jobname="{target}.gp",
        base="pipeline/04_probes/{target}",
        min_length=36,
        max_length=41,
        min_tm=42,
        max_tm=47
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    conda:
        "envs/oligominer.yaml"
    shell:
        """
        python scripts/blockParse.py  -l {params.min_length} -L {params.max_length} -t {params.min_tm} -T {params.max_tm} -f {input} -o {params.base}
        """

# Align probes against the entire set of genomic sequences (on-targets and off-targets)
rule align_probes:
    input:
        fq="pipeline/04_probes/{target}.fastq",
        ref="pipeline/02_renamed/allgenomes.rev.1.bt2"
    output:
        "pipeline/05_map/{target}.bam"
    params:
        jobname="{target}.mp",
        ref="pipeline/02_renamed/allgenomes"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        module load bowtie2/2.5.2 samtools/1.20
        bowtie2 -x {params.ref} -q {input.fq} --threads {threads} --very-sensitive-local -k 100 | samtools view -bS - > {output}
        """

# Convert alignment output into pairwise file and extract alignment scores
rule alignment_pairwise:
    input:
        "pipeline/05_map/{target}.bam"
    output:
        pair="pipeline/05_map/{target}.pair",
        scores="pipeline/05_map/{target}.txt"
    params:
        jobname="{target}.pa"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    conda:
        "envs/oligominer.yaml"
    shell:
        """
        module load samtools/1.20
        samtools view {input} | sam2pairwise > {output.pair}
        samtools view {input} | awk '{{print $12}}' > {output.scores}
        """

# Use alignment information to identify on-target and off-target matches
rule target_matches:
    input:
        pair="pipeline/05_map/{target}.pair",
        scores="pipeline/05_map/{target}.txt",
        targets="pipeline/02_renamed/{target}.gtf"
    output:
        "pipeline/06_alignments/{target}.csv"
    params:
        jobname="{target}.tm"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python scripts/parse_pairwise.py {input.pair} {input.scores} {input.targets} {output}
        """

# Predict hybridisation probabilities based on machine learning models
rule predict_duplex:
    input:
        probes="pipeline/06_alignments/{target}.csv",
        model="models/42_all_fixed_xgb.pickle.dat",
    output:
        "pipeline/07_predictions/{target}.csv"
    params:
        jobname="{target}.pr"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python scripts/XGB_predict.py {input.probes} {input.model} {output}
        """

# Format probe scores 
rule score_probes:
    input:
        "pipeline/07_predictions/{target}.csv"
    output:
        "pipeline/08_scores/{target}.bed"
    params:
        jobname="{target}.sc"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    conda:
        "envs/biopython.yaml"
    shell:
        """
        python scripts/output_bed.py {input} {output}
        """
# Filter out probes
rule filter_probes:
    input:
        "pipeline/08_scores/{target}.bed"
    output:
        "pipeline/09_kmer/{target}.tsv"
    params:
        jobname="{target}.fi"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        awk '$4 !~ /N/ && $4 !~ /-/' {input} > {output}
        """

# Count kmer frequency and hardcode strands
rule max_kmer:
    input:
        probes="pipeline/09_kmer/{target}.tsv",
        jellyfish="pipeline/02_renamed/allgenomes.jf"
    output:
        "probes/{target}.tsv"
    conda:
        'envs/biopython.yaml'
    params:
        mfree='20G',
        h_rt='3:0:0'
    shell:
        """
        python scripts/kmer_frequency {input.probes} {input.jellyfish} {output}
        """
