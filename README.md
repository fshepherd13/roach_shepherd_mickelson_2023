# Influenza infection of primary human airway cells
## Roach, Shepherd, Mickelson et al., 2023 (in review)
### 8/25/2023

#### Contents
This is a directory that contains all code and raw data (or directions to storage location for raw data) for published figures in Roach, Shepherd, Mickelson et al. Code to reproduce figures is at `scripts/NHBE paper code.Rmd` and `scripts/NHBE-figure-5-code.Rmd`. A knitted version of the code and subsequent figures is at `scripts/NHBE-paper-code.html`. Some of the raw data needs to be downloaded from FigShare (see below, figures 1B, 3C, and 3D).

Additionally, the directory `invert/` contains the code needed to analyze Influenza replication kinetics. The RNAseq reads for this experiment are deposited on NCBI under PRJNA1010235. It outputs the raw data needed for making figure 3D. 

Finally, the directory `scrnaseq/` contains the code needed to reproduce the mapping step of NHBE single cell RNAseq data. It outputs the raw data needed for making figure 5. Alternatively, the cell by gene matrices output by Cellranger are deposited in GEO under accession numbers GSM7904036 and GSM7904037 and can be used as input for the `scripts/NHBE-figure-5-code.Rmd`. 

#### Figure 1 data:
- 1A bar plot, 1C-1H: Data is located under `data/pre_and_post_infection_phenotyping_data.csv`. 
- 1B: Flow cytometry channel value files need to be downloaded from FigShare at [https://figshare.com/s/b19b3165426c795fddd3](https://figshare.com/s/b19b3165426c795fddd3) into the `data/round_7_channel_value_files` folder. Metadata is located at `figure_1_metadata.csv`.

#### Figure 2 data:
- 2A UMAP: The data for this is the same as figure 2B. 
- 2B,F: Data for alluvial plots are located at `data/nhbe_alluvial_plot_data.xlsx`. 
- 2C: Data is same as figure 1A,C-H. 
- 2D: Data is located at `data/sialic_acid_log_gmfi_values.csv`.
- 2E: Data available by request.

#### Figure 3 data:
- 3A,B: Data available by request.
- 3C: FI values need to be downloaded from FigShare at [https://figshare.com/s/7f48cb8873f0259eab12](https://figshare.com/s/7f48cb8873f0259eab12) into the `data/round_7_influenza_fluorescence_intensity` folder. Metadata is located at `data/figure_3c_metadata`.
- 3D: The main results to generate the graphs need to be downloaded from FigShare at [https://figshare.com/s/17e09862b1c23c0a6631](https://figshare.com/s/17e09862b1c23c0a6631) into the `data/invert_data` folder. Metadata is located at `data/figure_3d_metadata.csv`. To re-run InVERT from the raw data, fastq files can be downloaded via BioProject PRJNA1010235, SRA accession numbers SRR25779437-SRR25779478. The analysis can be run using the code in `invert`. There are two experiments, each should be run separately through the invert pipeline using the sample files under `invert/config`. The csv files in that directory list which SRA accession numbers belong to which experiment. 

#### Figure 4 data:
- 4A: Data available by request.
- 4B,C: Data is located at `data/nhbe_burst_size_data.csv`

#### Figure 5 data:
- Mapped files from Cellranger are also included in GEO under accession numbers GSM7904036 and GSM7904037.They should be downloaded and their paths provided to `Figure-5-code.Rmd` to re-run Seurat.
- Raw scRNAseq fastq files can be downloaded from SRA under accession numbers SRR26842481 and SRR26842482. The mapping steps can be reproduced with the code in `scrnaseq/nhbe_cellranger.sh`.