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
        "pipeline/input/allgenomes.fa"
    params:
        jobname="allgenomes.pf"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        #pending to add the renaming script
        cat {input} > {output}
        """

# Create a mapping file of contig headers, with modified new headers in the case of duplications
rule unique_headers:
    input:
        "pipeline/input/allgenomes.fa"
    output:
        "pipeline/renamed/headers.tsv"
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
        fasta="pipeline/input/allgenomes.fa",
        headers="pipeline/renamed/headers.tsv"
    output:
        "pipeline/renamed/allgenomes.fa"
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

# Index renamed fasta file for downstream probe mapping
rule index_fasta:
    input:
        "pipeline/renamed/allgenomes.fa"
    output:
        "pipeline/renamed/allgenomes.rev.1.bt2"
    params:
        base="pipeline/renamed/allgenomes",
        jobname="allgenomes.in"
    threads:
        1
    resources:
        mem_gb=8,
        time=60
    shell:
        """
        module load bowtie2/2.5.2
        bowtie2-build {input} {params.base}
        """

# Update the contig id-s in the gtf files using the header mapping file to avoid duplications.
rule unique_ids_gtf:
    input:
        gtf="targets/{target}.gtf",
        headers="pipeline/renamed/headers.tsv"
    output:
        "pipeline/renamed/{target}.gtf"
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
        fasta="pipeline/renamed/allgenomes.fa",
        gtf="pipeline/renamed/{target}.gtf"
    output:
        "pipeline/extract/{target}.fa"
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
        "pipeline/extract/{target}.fa"
    output:
        "pipeline/probes/{target}.fastq"
    params:
        jobname="{target}.gp",
        base="pipeline/probes/{target}",
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
        fq="pipeline/probes/{target}.fastq",
        ref="pipeline/renamed/allgenomes.rev.1.bt2"
    output:
        "pipeline/map/{target}.bam"
    params:
        jobname="{target}.mp",
        ref="pipeline/renamed/allgenomes"
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
        "pipeline/map/{target}.bam"
    output:
        pair="pipeline/map/{target}.pair",
        scores="pipeline/map/{target}.txt"
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
        pair="pipeline/map/{target}.pair",
        scores="pipeline/map/{target}.txt",
        targets="pipeline/renamed/{target}.gtf"
    output:
        "pipeline/alignments/{target}.csv"
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
        probes="pipeline/alignments/{target}.csv",
        model="models/42_all_fixed_xgb.pickle.dat",
    output:
        "pipeline/predictions/{target}.csv"
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
        "pipeline/predictions/{target}.csv"
    output:
        "pipeline/scores/{target}.bed"
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

rule filter_probes:
    input:
        "pipeline/scores/{target}.bed"
    output:
        "probes/{target}.tsv"
    params:
        jobname="{target}.fi"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        #filters probes
        """
