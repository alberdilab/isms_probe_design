######
# InSituMicrobeSeq probe design pipeline
# Martí Dalmases & Antton alberdi
# 2023/07/27
# Description: draft backbone pipeline for probe design
######

# 1) Copy this snakefile to the working directory
# 2) Store the genome sequences in the folder 'genomes' in the working directory with extension .fa, and without any "." or "/" in the file name besides the extension .fa.
# 3) Store the target-defining GTF files in the 'targets' directory with extension .gtf, and without any "." or "/" in the file name besides the extension .gtf.
# 4) Launch the snakemake using the following code:
# snakemake -j 20 --cluster 'sbatch -o logs/{params.jobname}-slurm-%j.out --mem {resources.mem_gb}G --time {resources.time} -c {threads} --job-name={params.jobname} -v'   --use-conda --conda-frontend mamba --conda-prefix conda --latency-wait 600

#List genome and target wildcards
genomes, = glob_wildcards("genomes/{genome}.fa")
targets, = glob_wildcards("targets/{target}.gtf")

#Expand target files
rule all:
    input:
        expand("probes/{target}.tsv", target=targets)

rule prepare_fasta:
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
        #concatenate and rename fasta
        """

rule index_fasta:
    input:
        "pipeline/input/allgenomes.fa"
    output:
        "pipeline/input/allgenomes.bt2"
    params:
        jobname="allgenomes.in"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        #index fasta
        """

rule prepare_gtf:
    input:
        "targets/{target}.gtf"
    output:
        "pipeline/input/{target}.gtf"
    params:
        jobname="{target}.pg"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        #rename gtfs
        """

rule extract_targets:
    input:
        fasta="pipeline/input/allgenomes.fa",
        gtf="pipeline/input/{target}.gtf"
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
        #extract target fasta
        """

rule generate_kmers:
    input:
        "pipeline/extract/{target}.fa"
    output:
        "pipeline/kmers/{target}.jf"
    params:
        jobname="{target}.jf"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        #generate kmer profiles
        """

rule generate_probes:
    input:
        fasta="pipeline/extract/{target}.fa",
        jf="pipeline/kmers/{target}.jf"
    output:
        "pipeline/probes/{target}.fq"
    params:
        jobname="{target}.gp"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        #generate probes with blockparse
        """

rule align_probes:
    input:
        fq="pipeline/probes/{target}.fq",
        ref="pipeline/input/allgenomes.bt2"
    output:
        "pipeline/map/{target}.sam"
    params:
        jobname="{target}.mp"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        #aligns probes to all genomes
        """

rule score_probes:
    input:
        "pipeline/map/{target}.sam"
    output:
        "pipeline/score/{target}.tsv"
    params:
        jobname="{target}.sc"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        #scores probes
        """

rule filter_probes:
    input:
        "pipeline/score/{target}.tsv"
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
