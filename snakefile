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

rule align_probes:
    input:
        fq="pipeline/probes/{target}.fastq",
        ref="pipeline/renamed/allgenomes.rev.1.bt2"
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
        module load bowtie2/2.5.2 samtools/1.20
        bowtie2 -x {input.ref} -q {input.fq} --threads {threads} --very-sensitive-local -k 100 | samtools sort -@ {threads} -o {output}
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
