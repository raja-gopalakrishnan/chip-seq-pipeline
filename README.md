# ChIP-Seq analysis pipeline
This pipeline contains the code that was used for analyzing the ChIP-Seq data published in our paper [here](https://www.biorxiv.org/content/10.1101/364521v1). The code has been customized to analyze our specific ChIP-Seq dataset. However, it can be used to analyze any other ChIP-Seq data that uses in-line barcodes for multiplexing libraries and has the same input for all IP samples.

## Getting started

### Required software
- Unix operating system
- git
- conda

### Input files
- FASTQ file (samples multiplexed with in-line barcodes). It is necessary to have input samples for carrying out all the normalization steps. All IP samples from a particular genotype must have the same input. All samples must have two replicates.


## Instructions

### 1. Clone the repository
```
git clone https://github.com/raja-gopalakrishnan/chip-seq-pipeline.git
```

### 2. Edit the config file
Copy the template_config.yaml file to make the config.yaml file
```
cp template_config.yaml config.yaml
```
Open the config.yaml file and edit the following items:
- ```samples``` : Edit sample names, barcode, group, ip, genotype and replicate number. Sample names should follow the ip-genotype-replicate naming system. The group attribute will by 'ip-genotype'.
- ```fastq```: Provide location of fastq file.
- ```normalizations```: Enter samples that each IP library should be normalized by when plotting heatmaps and metagenes. Can enter multiple values for one IP.
- ```comparisons```: Enter samples that each IP library should be normalized by when comparing occupancy within genes across genotypes. Can enter only one value for each IP.
- ```callpeaks```: Enter sample to be normalized by when calling for peaks using macs2. This is usually the input for all samples.
Change other parameters as required.

### 3. Create the environment to run the pipeline
```
# Create an environment where all packages required for running the pipeline are installed
conda env create -f envs/snakemake_chip-seq.yaml

# Activate the environment
source activate snakemake_chip-seq.yaml
```

### 4. Run the pipeline
```
# Do a dry run of the pipeline which will show you the number of jobs that will be submitted
snakemake -np
```

If running the pipeline locally on your computer use the following command to run it:
```
snakemake -p
```

If running the pipeline on an HPC cluster (recommended) that uses the SLURM job scheduler, use the following command:
```
sbatch slurm.sh
```
If using a HPC cluster with a different job scheduler, modify the ```slurm.sh``` and ```cluster.yaml``` files accordingly

## Acknowledgements
- Thank you James Chuang, for providing all the genome annotation files and the R script for making correlation plots!
