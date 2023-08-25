# invert

### Pipeline for quantifying mRNA, cRNA and vRNA transcripts during influenza virus replication.

## Overview
This is a copy of the pipeline used to generate the data in our NHBE manuscript. The original repo with additional details is found at (https://github.com/langloislab/invert). 

## Usage


To run the full pipeline, run from your terminal:
`snakemake --cores 24` (<-- or however many cores you want Snakemake to utilize)

Alternatively, you can use the `workflow/invert.srun` to run on some type of slurm cluster. 

### Data and reference genomes
FASTQ files for the RNAseq are deposited on NCBI under Geo accession number XYZ. They will need to be downloaded, and the `config.yaml` file `in_dir` parameter will need to be updated to match the input file path.

Reference sequences and the annotation gtf file for the Cal09 strain are included in the `ref_files` directory. You will also need a hybrid human-Cal09 reference genome. To create this, see the script at `scripts/compile_genome.srun` and modify for your needs. The reference genome fasta and gtf file paths will need to be updated in the config file as well.


### Running the pipeline

The pipeline requires the following dependencies:
* Snakemake and snakemake-wrapper-utils
* Trimmomatic
* STAR
* fastqc
* cufflinks
* samtools
* pandas

This pipeline runs using snakemake. Each rule contains a conda environment with the software necessary to run it. To run the pipeline, you only need to install snakemake, i.e:
'''
conda env create -n snakemake -c bioconda snakemake
conda activate snakemake
'''

Each rule will create an isolated conda environment automatically as it runs.

The pipeline will output a final csv file with m, c, and vRNA quantifications. These files are used as input to the scripts in `../scripts/NHBE-paper-code.Rmd`. Copies of the final csv files are available for download (see information in `../data/`).