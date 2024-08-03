# InSituMicrobeSeq probe design pipeline
ISMS Probe Design is a bioinformatic pipeline to design probes for in-situ sequencing. The software is based on the [PaintSHOP](https://github.com/beliveau-lab/PaintSHOP_pipeline) pipeline to design oligonucleotides for FISH experiments, and relies on [OligoMiner](https://github.com/beliveau-lab/OligoMiner) for candidate probe identification. It is implemented in [Snakemake](https://snakemake.readthedocs.io/en/stable/) and provided as a standalone repository with all dependencies handled through [conda](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html) environments. The pipeline can be run both locally (suitable for small experiments) or in a remote server (suitable for large experiments) with slurm integration.

- [Installation](#installation): instructions to get isms_probe_design ready.
- [Workflow](#workflow): overview of the procedures.
- [Results](#results): explanation of the final results.
- [Tutorial](#tutorial): example with real data.

## Installation

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
mamba env create -f workflow/envs/environment.yml
conda activate isms_probe_design_env
```

#### Environments

- The environment file `workflow/envs/environment.yml` contains the basal tools required to launch the snakemake pipeline and target creator.
- The environment file `workflow/envs/python2_env.yml`, which is automatically installed by snakemake, contains python 2 dependencies required by the pipeline.
- The environment file `workflow/envs/python3_env.yml`, which is automatically installed by snakemake, contains python 3 dependencies required by the pipeline.

## Workflow

### Define targets

The probe design pipeline requires target regions to be defined using GTF annotation files that reference genomic regions in the considered FASTA files. These target regions can represent entire genomes for the detection of specific taxa or span specific regions across multiple taxa for the detection of specific functions. To streamline the creation of GTF target files, we provide the script create_target.py, which enables the extraction of target regions from FASTA and GTF files.

```
conda activate isms_probe_design_env
python workflow/scripts/create_target.py [-h] -m {region,genome} [-g GTF [GTF ...]] [-a ANNOTATION] [-e] [-f FASTA] -o OUTPUT
```
1. Choose mode:
     - -m / --mode: `region` for function detection (GTF) or `genome` for taxonomic detection (FASTA). Mandatory.
  + Region:
     - -g / --gtf: List the input GTF files or the folder containing them.
     - -a / --annotation: Single attribute for filtering.
  + Genome:
     - -f / --fasta: Input FASTA file with the desired genom target.
2. Output:
     - -o/--output: File path to the output file in `./targets/{target}.gtf`. Mandatory
   
#### Genome target example
Create a GTF file targetting the entire genome of ***Escherichia coli***.
```
conda activate isms_probe_design_env
python workflow/scripts/create_target.py -m genome \
     -f GCF_000005845.2_ASM584v2_genomic.fna \
     -o targets/escherichia_coli.gtf
```

#### Region target example
Create a GTF file targetting the gene speE within the genomes of ***Escherichia coli***, ***Pseudomonas aeruginosa*** and ***Bacillus subtilis***.
```
conda activate isms_probe_design_env
python scripts/create_target.py -m region \
     -g GTF_directory \
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
conda activate isms_probe_design_env
snakemake \
  -j 20 \
  --cluster 'sbatch -o results/log/{params.jobname}-slurm-%j.out --mem {resources.mem_gb}G --time {resources.time} -c {threads} --job-name={params.jobname} -v' \
  --use-conda --conda-frontend mamba --conda-prefix conda \
  --latency-wait 600
```

## Results

Data files generated by the pipeline are stored in the `results` directory. Intermediate files can be found in the `pipeline` folder, while the final set of probes is stored in the `probes` folder. A tsv file containing all the relevant information of the probes is generated from each targets GTF file.

```
- results/probes
 - target1.tsv
 - target2.tsv
```

The final probe set file has the following structure:
| contig | start | end | sequence | Tm | on_target_score | off_target_score | repeat | kmer_count | strand |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| genome1_contig1 | 1 | 2198 | GTTGTGTTGAGAGTGGTGTGAGTGGAACTATCATTCCA | 42.320 | 100.000 | 0.000 | 0 | 0 | + |
| genome1_contig4 | 6770 | 13645 | TGGCCCTTCACCTTTTCAGATGAACCGTAAGCGCTG | 45.900 | 99.077 | 0.000 | 0 | 0 | + |
(...)

- The first three columns define the target region.
- The fourth column is the sequence of the probe.
- The final six columns contain attributes of the probe.

## Tutorial

In this small tutorial four bacterial genomes available at NCBI are used to showcase the design of probes targeting entire genomes and specific genomic regions.

### 1. Get genome data from NCBI

Create a new folder named `data` in the resources directory to store the data downloaded from the NCBI.

```
mkdir resources/data
cd resources/data
```

Use wget to fetch genomic data of ***[Stenotrophomonas rhizophila](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000661955.1/)***, ***[Microbacterium oxydans](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_006540085.1/)***, ***[Xanthomonas retroflexus](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_900143175.1/)*** and ***[Paenibacillus amylolyticus](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_029542105.1/)*** directly from the NCBI FTP.

```
# Get genome data of Stenotrophomonas rhizophila
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/661/955/GCF_000661955.1_ASM66195v1/GCF_000661955.1_ASM66195v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/661/955/GCF_000661955.1_ASM66195v1/GCF_000661955.1_ASM66195v1_genomic.gtf.gz

# Get genome data of Microbacterium oxydans
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/006/540/085/GCF_006540085.1_ASM654008v1/GCF_006540085.1_ASM654008v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/006/540/085/GCF_006540085.1_ASM654008v1/GCF_006540085.1_ASM654008v1_genomic.gtf.gz

# Get genome data of Xanthomonas retroflexus
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/143/175/GCF_900143175.1_ASM90014317v1/GCF_900143175.1_ASM90014317v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/143/175/GCF_900143175.1_ASM90014317v1/GCF_900143175.1_ASM90014317v1_genomic.gtf.gz

# Get genome data of Paenibacillus amylolyticus
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/542/105/GCF_029542105.1_ASM2954210v1/GCF_029542105.1_ASM2954210v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/029/542/105/GCF_029542105.1_ASM2954210v1/GCF_029542105.1_ASM2954210v1_genomic.gtf.gz
```

Decompress all files to get them ready for target GTF generation.

```
gunzip *
cd ../..
```

### 2. Prepare target GTFs

First of all, remove the example files from the `genomes` and `targets` directories.

```
rm resources/genomes/*.fna
rm resources/targets/*.gtf
```

#### 2.1 Genome-level targets

Using the `create_target.py` script, create target files to design probes to detect each genome. The script generates a targets GTF file from the FASTA sequence file.

```
conda activate isms_probe_design_env

# Create targets for Stenotrophomonas rhizophila
python workflow/scripts/create_target.py -m genome \
     -f resources/data/GCF_000661955.1_ASM66195v1_genomic.fna \
     -o resources/targets/stenotrophomonas_rhizophila.gtf

# Create targets for Microbacterium oxydans
python workflow/scripts/create_target.py -m genome \
     -f resources/data/GCF_006540085.1_ASM654008v1_genomic.fna \
     -o resources/targets/microbacterium_oxydans.gtf

# Create targets for Xanthomonas retroflexus
python workflow/scripts/create_target.py -m genome \
     -f resources/data/GCF_900143175.1_ASM90014317v1_genomic.fna \
     -o resources/targets/xanthomonas_retroflexus.gtf

# Create targets for Paenibacillus amylolyticus
python workflow/scripts/create_target.py -m genome \
     -f resources/data/GCF_029542105.1_ASM2954210v1_genomic.fna \
     -o resources/targets/paenibacillus_amylolyticus.gtf
```

Note the script generated single-row GTF files for ***Stenotrophomonas rhizophila*** and ***Paenibacillus amylolyticus***, while multi-row GTF files were generated for ***Microbacterium oxydans*** and ***Xanthomonas retroflexus***. This is because the first two are circularised genomes with a single contig, while the last two are draft genomes with multiple contigs. 

#### 2.2 Region-level targets

Using the `create_target.py` script, create target files to design probes to detect specific regions across genomes. The script parses the annotation identifier and generates the targets GTF file with genomic region information from all candidate genome files.

For instance, biotin synthase BioB is present in ***Stenotrophomonas rhizophila*** and ***Xanthomonas retroflexus***, but not in ***Paenibacillus amylolyticus*** and ***Microbacterium oxydans***.

For instance, methylisocitrate lyase prpB is present in ***Stenotrophomonas rhizophila***,  ***Xanthomonas retroflexus*** and ***Microbacterium oxydans***, but not in ***Paenibacillus amylolyticus***.

```
conda activate isms_probe_design_env

# Create targets for biotin synthase BioB
python workflow/scripts/create_target.py -m region \
     -g resources/data \
     -a 'biotin synthase' \
     -o resources/targets/BioB.gtf

python workflow/scripts/create_target.py -m region \
     -g resources/data \
     -a 'methylisocitrate lyase' \
     -o resources/targets/prpB.gtf
```

### 3. Prepare sequence FASTA files

Sequence FASTA files need to be stored in the `resources/genomes` folder.

```
cp resources/data/*.fna resources/genomes/
```

Now all relevant files should be located in their corresponding paths. Note that different target types (genome and region) can be processed simultaneously by the pipeline.

```
- resources/genomes/
  - GCF_000661955.1_ASM66195v1_genomic.fna
  - GCF_006540085.1_ASM654008v1_genomic.fna
  - GCF_900143175.1_ASM90014317v1_genomic.fna
  - GCF_029542105.1_ASM2954210v1_genomic.fna

- resources/targets/
  - stenotrophomonas_rhizophila.gtf
  - microbacterium_oxydans.gtf
  - xanthomonas_retroflexus.gtf
  - paenibacillus_amylolyticus.gtf
  - BioB.gtf
  - prpB.gtf
```

### 4. Run probe design

```
screen -S isms_probe_design
cd isms_probe_design
conda activate isms_probe_design_env
snakemake \
  -j 20 \
  --cluster 'sbatch -o results/log/{params.jobname}-slurm-%j.out --mem {resources.mem_gb}G --time {resources.time} -c {threads} --job-name={params.jobname} -v' \
  --use-conda --conda-frontend mamba --conda-prefix conda \
  --latency-wait 600
```
