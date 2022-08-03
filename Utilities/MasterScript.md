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
--> *data/SraRunTable.txt*

NOTE: copy of SraRunTable.txt kept in Git repo for reference, but the scripts
require it in **data**

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

The script `01-download_SRA.sh` has been written to run on a cluster using the
SLURM scheduler. It uses a batch array to determine with SRA file to download.
As well as downloading the sra files and extracting the fastq files, it will 
also rename the fastq to the convention expected by CellRanger.

```
sbatch 01-download_SRA.sh
```
--> data/sra/SRR92643*/SRR92643*.sra (12 files)
--> data/fastq/SRR92643*_S1_L001_I1_001.fast.gz (12 files)
--> data/fastq/SRR92643*_S1_L001_R1_001.fast.gz (12 files)
--> data/fastq/SRR92643*_S1_L001_R2_001.fast.gz (12 files)

# 2) Install Cell Ranger

NOTE: the URL needs to be updated from the 10X website:

https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

```
wget -O cellranger-7.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-7.0.0.tar.gz?Expires=1658809505&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci03LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2NTg4MDk1MDV9fX1dfQ__&Signature=YXaQiedAbBsXkjBHozopnWRWOnMeOTQTLhQ6-Ct~LIE7iKn37s~1egogP4mRvtF0QxH2VXSnGZUTL-1ncY7MzTOlO25ogyOkqIjsCMNs1k-N7hB998UnmsX3X2-km8zGQQENqbOiwvsOJ--dbsdYZn9x0EIrNbg6GUYCFeoI0EhF1XIcpioZ5VKkcX9Q8hqOxEkTspNWT5t7JS3fF1T2c4T8xS2fOTlf-kS9fK3x67~aLY5DhPsd98SZnsiD5VxTN3CdZFKWY~nCUN8Enx7f0Pm9MD8m3T08OyJjr0vG~0X26cqpt9xldWV5zymfczMxgJZ8nVVS-giFJJFgzZlZsw__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

mkdir data/software
tar -C data/software -xzvf cellranger-7.0.0.tar.gz
rm -f cellranger-7.0.0.tar.gz
```
--> **data/software/cellranger-7.0.0**


# 3. Make CellRanger reference - full genome

We will update to use to use the latest Mouse genome/annotation - GRCh38.p13 release
107 and Gencode Human v41.

## a) Retrieve and prepare references

```
sbatch 03_a-download_and_prepare_references.sh
```
--> data/references/gencode.v41.primary_assembly.annotation.filtered.gtf
--> data/references/Homo_sapiens.GRCh38.dna.primary_assembly.gencoded.fa

## b) Make Cell Ranger reference

```
sbatch 03_b-make_Cell_Ranger_reference.sh
```
--> **data/references/refdata-gex-GRCh38.p13-Gencode.v41**

# 4) Run Cell Ranger count

```
sbatch 04-cellranger_count.sh
```

# 5. Make Chr21 reference

For the Cell Ranger exercise we need to run SRR9264343 against the chr21 reference.
We also need to subsample the fastq to 1 million reads

```
sbatch 05-prepare_chr21_reference.sh
```
--> **cellranger_index**
--> gencode.v41.primary_assembly.annotation.chr21.gtf
--> Homo_sapiens.GRCh38.dna.chromosome.21.fa

# 6. Subsample fastq

For the Cell Ranger exercise we need a fastq with just 1 million reads

```
sbatch 06-subsample_fastq.sh
```
--> data/fastq_subsample/SRR9264343_S1_L001_I1_001.fastq.gz
--> data/fastq_subsample/SRR9264343_S1_L001_R1_001.fastq.gz
--> data/fastq_subsample/SRR9264343_S1_L001_R2_001.fastq.gz

# 7. Run Cell Ranger count - exercise

```
sbatch 07-cellranger_exercise.sh
```
--> **data/cellranger_exercise/ETV6_RUNX1_rep1/**

# 8. Make sample sheet

```
08-make_sample_sheet.R
```

