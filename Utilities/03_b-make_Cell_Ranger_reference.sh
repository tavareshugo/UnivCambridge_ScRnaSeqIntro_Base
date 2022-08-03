#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -t 2:00:00
#SBATCH --mem=64G
#SBATCH -o cellranger_mkref.%j.out
#SBATCH -e cellranger_mkref.%j.err
#SBATCH -J cellranger_mkref

set -ex

echo "Start: $(date)"

cd data

# Cell Ranger software
crDir=`readlink -f software/cellranger-7.0.0`
export PATH=${crDir}:${PATH}

# references
cd references
genome=refdata-gex-GRCh38.p13-Gencode.v41
fasta="Homo_sapiens.GRCh38.dna.primary_assembly.gencoded.fa"
gtf="gencode.v41.primary_assembly.annotation.filtered.gtf"

# mkref
cellranger mkref --genome=${genome} \
                 --fasta=${fasta} \
                 --genes=${gtf} \
                 --memgb=64 \
                 --nthreads=8
