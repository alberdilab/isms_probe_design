# InSituMicrobeSeq probe design pipeline

### Installation

1. Make sure you have [conda](https://docs.anaconda.com/miniconda) installed.
2. Make sure you have [mamba](https://github.com/mamba-org/mamba) installed in the base environment of conda.
```
conda install -n base -c conda-forge mamba
```
3. Clone this repository to the desired path in your local environment.
```
cd [the/desired/base/path]
git clone https://github.com/alberdilab/isms_probe_design.git
```
4. Create and activate the isms_probe_design environment.
```
cd isms_probe_design
mamba env create -f environment.yml
conda activate isms_probe_design_env
```

#### Environments

- The environment file `environment.yml` contains the basal tools required to launch the snakemake pipeline and target creator.
- The environment file `envs/python2_env.yml`, which is automatically installed by snakemake, contains python 2 dependencies required by the pipeline.
- The environment file `envs/python3_env.yml`, which is automatically installed by snakemake, contains python 3 dependencies required by the pipeline.

### Define targets

The probe design pipeline requires target regions to be defined using GTF annotation files that reference genomic regions in the considered FASTA files. These target regions can represent entire genomes for the detection of specific taxa or span specific regions across multiple taxa for the detection of specific functions. To streamline the creation of GTF target files, we provide the script create_target.py, which enables the extraction of target regions from FASTA and GTF files.

```
conda activate isms_probe_design_env
python scripts/create_target.py [-h] -m {region,genome} [-g GTF [GTF ...]] [-a ANNOTATION] [-e] [-f FASTA] -o OUTPUT
```
1. Choose mode:
     - -m / --mode: `region` for function detection (GTF) or `genome` for taxonomic detection (FASTA). Mandatory.
  + Region:
     - -g / --gtf: List the input GTF files or the folder containing them.
     - -a / --annotation: Single attribute for filtering.
     - -e / --exon: Filter for exon only regions. Optional.
  + Genome:
     - -f / --fasta: Input FASTA file with the desired genom target.
2. Output:
     - -o/--output: File path to the output file in `./targets/{target}.gtf`. Mandatory
   
#### Genome target example
Create a GTF file targetting the entire genome of ***Escherichia coli***.
```
conda activate isms_probe_design_env
python scripts/create_target.py -m genome \
     -f GCF_000005845.2_ASM584v2_genomic.fna \
     -o targets/escherichia_coli.gtf
```

#### Region target example
Create a GTF file targetting the gene speE within the genomes of ***Escherichia coli***, ***Pseudomonas aeruginosa*** and ***Bacillus subtilis***.
```
conda activate isms_probe_design_env
python scripts/create_target.py -m region \
     -g GCF_000005845.2_ASM584v2_genomic.gtf,GCF_000006765.1_ASM676v1_genomic.gtf,GCF_000009045.1_ASM904v1_genomic.gtf \
     -a speE \
     -o targets/speE.gtf
```

### Prepare input files
Once the reference FASTA files and GTF target files are ready, these need to be stored in specific directories.

- FASTA files (.fa) of all considered genomes must be stored in the `genomes` folder. 
- GTF files (.gtf) containing regions of target sequences must be stored in the `targets` folder.

The mock data files contain 3 fasta and 2 gtf files. The file `target1.gtf` contains two regions from genome2 and genome3, while the file `target2.gtf` contains the entire genome1 as target sequence.

### Run pipeline locally
Snakemake automatically runs jobs locally. On a screen session, launch the snakefile to design the probes using the following commands:
```
screen -S isms_probe_design
conda activate isms_probe_design_env
snakemake \
  --use-conda --conda-frontend mamba --conda-prefix conda \
  --cores 4 \
  --latency-wait 600
```

### Run pipeline in a high computation cluster (using slurm)
Snakemake automatically launches jobs to the slurm queue. On a screen session, launch the snakefile to design the probes using the following commands:
```
screen -S isms_probe_design
cd isms_probe_design
module purge && module load snakemake/7.20.0 mamba/1.3.1
snakemake \
  -j 20 \
  --cluster 'sbatch -o log/{params.jobname}-slurm-%j.out --mem {resources.mem_gb}G --time {resources.time} -c {threads} --job-name={params.jobname} -v' \
  --use-conda --conda-frontend mamba --conda-prefix conda \
  --latency-wait 600
```
