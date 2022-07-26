# Preparation of Course Materials

This master script describes how to prepare all data objects required for 
the course.

The original data are from the paper:

Caron, M., St-Onge, P., Sontag, T. et al. Single-cell analysis of childhood
leukemia reveals a link between developmental states and ribosomal protein
expression as a source of intra-individual heterogeneity. Sci Rep 10, 8079
(2020). https://doi.org/10.1038/s41598-020-64929-x

Single cell RNA-seq data are available in GEO under GSE132509 or SRA:

https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA548203&o=acc_s%3Aa

# 1. Download raw data from SRA

## a) Get metadata

```
mkdir data
```

The SRA a metadata table was downloaded from the SRA URL above SRA Run
Selector page (URL above).
--> *SraRunTable.txt*

## b) Install SRA tools in a conda environment

We will set up a conda environment to run the SRA tool kit.

You will need to have conda(/miniconda/anaconda) installed and have the `bin``
directory on your PATH.

```
source activate 
conda create -n sra 
conda install -n sra -c bioconda sra-tools 
conda deactivate
```

## c) Download data

### a) Retrieve meta data

The script `00-download_SRA.sh` has been written to run on a cluster using the
SLURM scheduler. It uses a batch array to determine with SRA file to download.
As well as downloading the sra files and extracting the fastq files, it will 
also rename the fastq to the convention expected by CellRanger.

This script should be run in the directory containing the *SraRunTable.txt*
file.

```
sbatch 00-download_SRA.sh
```
--> data/sra/SRR92643*/SRR92643*.sra (12 files)
--> data/fastq/SRR92643*_S1_L001_I1_001.fast.gz (12 files)
--> data/fastq/SRR92643*_S1_L001_R1_001.fast.gz (12 files)
--> data/fastq/SRR92643*_S1_L001_R2_001.fast.gz (12 files)

