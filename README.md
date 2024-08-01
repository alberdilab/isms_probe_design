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

### Create targets

The provided script `create_target.py` allows for target extraction from FASTA and GTF files. 
```
usage: scripts/create_target.py [-h] -m {region,genome} [-g GTF [GTF ...]] [-a ANNOTATION] [-e] [-f FASTA] -o OUTPUT
```
1. Choose mode:
     - -m / --mode: `region` for function detection (GTF) or `genome` for taxonomic detection (FASTA). Mandatory.
  + Region:
     - -g / --gtf: List the input GTF files or the folder containing them. Mandatory.
     - -a / --annotation: Single attribute for filtering. Mandatory.
     - -e / --exon: Filter for exon only regions. Optional.
  + Genome:
     - -f / --fasta: Input FASTA file with the desired genom target. Mandatory.
2. Output:
     - -o/--output: File path to the output file in `./targets/{target}.gtf`. Mandatory
### Prepare input files
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
