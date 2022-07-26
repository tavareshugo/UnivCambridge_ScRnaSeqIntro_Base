#!/bin/bash
#SBATCH -J mkref
#SBATCH -o 01-prepare_reference.%j.log
#SBATCH -e 01-prepare_reference.%j.err
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH --mem 16G

set -exu

# Cell Ranger software
crDir=`readlink -f data/software/cellranger-7.0.0`
export PATH=${crDir}:${PATH}

mkdir -p data/references
cd data/references

#### Download reference genome ####

# download only chromosome 21 (to reduce size of the data)
ensURL=https://ftp.ensembl.org/pub/release-104
URL=${ensURL}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
wget --no-check-certificate ${URL}
gunzip Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz

#### Download and prepare gene annotation ####

# download from ENSEMBL
URL=${ensURL}/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz
wget --no-check-certificate ${URL}
gunzip Homo_sapiens.GRCh38.104.chr.gtf.gz

grep -E "^21|^#" Homo_sapiens.GRCh38.104.chr.gtf \
  > Homo_sapiens.GRCh38.104.chr21.gtf

#### cellranger indexing ####

cellranger mkref \
  --genome=cellranger_index \
  --fasta=Homo_sapiens.GRCh38.dna.chromosome.21.fa \
  --genes=Homo_sapiens.GRCh38.104.chr21.gtf \
  --nthreads=8 \
  --memgb=16
