# InSituMicrobeSeq probe design pipeline

### Clone this repository
```
git clone https://github.com/alberdilab/isms_probe_design.git
```

### Prepare input files
- FASTA files (.fa) of all considered genomes must be stored in the `genomes` folder. 
- GTF files (.gtf) containing regions of target sequences must be stored in the `targets` folder.

The mock data files contain 3 fasta and 2 gtf files. The file `target1.gtf` contains two regions from genome2 and genome3, while the file `target2.gtf` contains the entire genome1 as target sequence.

### Run pipeline
On a screen session, launch the snakefile to design the probes
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
