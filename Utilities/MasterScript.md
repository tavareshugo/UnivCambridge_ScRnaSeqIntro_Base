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

# 2) CellRanger reference

## a) Full genome reference

We will use the most current reference provided by 10X.

```
mkdir data/references
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -C data/references -xzvf refdata-gex-GRCh38-2020-A.tar.gz
rm -f refdata-gex-GRCh38-2020-A.tar.gz
```
--> data/references/refdata-gex-GRCh38-2020-A

## b) Chr21 reference

The Cell Ranger practical requires a reference that only includes chromosome 21.

```
sbatch 01-prepare_chr21_reference.sh
```

# 3) Install Cell Ranger

NOTE: the URL needs to be updated from the 10X website:

https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

```
wget -O cellranger-7.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.0.0.tar.gz?Expires=1658809505&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NTg4MDk1MDV9fX1dfQ__&Signature=YXaQiedAbBsXkjBHozopnWRWOnMeOTQTLhQ6-Ct~LIE7iKn37s~1egogP4mRvtF0QxH2VXSnGZUTL-1ncY7MzTOlO25ogyOkqIjsCMNs1k-N7hB998UnmsX3X2-km8zGQQENqbOiwvsOJ--dbsdYZn9x0EIrNbg6GUYCFeoI0EhF1XIcpioZ5VKkcX9Q8hqOxEkTspNWT5t7JS3fF1T2c4T8xS2fOTlf-kS9fK3x67~aLY5DhPsd98SZnsiD5VxTN3CdZFKWY~nCUN8Enx7f0Pm9MD8m3T08OyJjr0vG~0X26cqpt9xldWV5zymfczMxgJZ8nVVS-giFJJFgzZlZsw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

mkdir data/software
tar -C data/software -xzvf cellranger-7.0.0.tar.gz
rm -f cellranger-7.0.0.tar.gz
```
--> **data/software/cellranger-7.0.0**
