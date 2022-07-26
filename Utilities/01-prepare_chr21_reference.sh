#!/bin/bash
#SBATCH -J mkref
#SBATCH -o 01-prepare_reference.%j.log
#SBATCH -e 01-prepare_reference.%j.err
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH -p cclake
#SBATCH --mem-per-cpu=3420MB

set -exu

refDir=data/reference
mkdir -p ${refDir}

#### Download reference genome ####

# download only chromosome 21 (to reduce size of the data)
ensURL=https://ftp.ensembl.org/pub/release-104
URL=${ensURL}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz
wget --no-check-certificate -P ${refDir} ${URL}
gunzip ${refDir}/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz

#### Download and prepare gene annotation ####

# download from ENSEMBL
URL=${ensURL}/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz
wget --no-check-certificate -P ${refDir} ${URL}
gunzip ${refDir}/Homo_sapiens.GRCh38.104.chr.gtf.gz

grep -E "^21|^#" ${refDir}/Homo_sapiens.GRCh38.104.chr.gtf \
  > ${refDir}/Homo_sapiens.GRCh38.104.chr21.gtf

#### cellranger indexing ####

cellranger mkref \
  --nthreads="${SLURM_CPUS_PER_TASK}" \
  --genome=cellranger_index \
  --fasta=${refDir}/Homo_sapiens.GRCh38.dna.chromosome.21.fa \
  --genes=Homo_sapiens.GRCh38.104.chr21.gtf
