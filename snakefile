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
        expand("results/probes/{target}.tsv", target=targets)

rule prepare_fasta:
    input:
        expand("genomes/{genome}.fa", genome=genomes)
    output:
        "pipeline/input/allgenomes.fa"
    params:
        jobname="{sample}.pf"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        #concatenate and rename fasta
        """
        
rule prepare_gtf:
    input:
        "targets/{target}.gtf"
    output:
        "pipeline/input/{target}.gtf"
    params:
        jobname="{sample}.pg"
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
        fasta="pipeline/input/allgenomes.fa"
        gtf="pipeline/input/{target}.gtf"
    output:
        "pipeline/extract/{target}.fa"
    params:
        jobname="{sample}.ex"
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
        jobname="{sample}.jf"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        #generate kmer profiles
        """

rule generate_kmers:
    input:
        fasta="pipeline/extract/{target}.fa"
        jf="pipeline/kmers/{target}.jf"
    output:
        "pipeline/kmers/{target}.jf"
    params:
        jobname="{sample}.jf"
    threads:
        1
    resources:
        mem_gb=8,
        time=30
    shell:
        """
        #generate kmer profiles
        """
