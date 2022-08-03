#!/bin/bash
#SBATCH -J prepare_chr21_reference
#SBATCH -o prepare_chr21_reference.%j.log
#SBATCH -e prepare_chr21_reference.%j.err
#SBATCH -c 8
#SBATCH -t 12:00:00
#SBATCH --mem 16G

set -exu

# Cell Ranger software
crDir=`readlink -f data/software/cellranger-7.0.0`
export PATH=${crDir}:${PATH}

cd data/references

#### Download reference genome ####

# download only chromosome 21 (to reduce size of the data)
fasta_url="http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz"
fasta_lcl="Homo_sapiens.GRCh38.dna.chromosome.21.fa.gz"
wget --no-check-certificate -O ${fasta_lcl} ${fasta_url} 

# modify contig name

fasta_gencoded="Homo_sapiens.GRCh38.dna.chromosome.21.fa"
zcat ${fasta_lcl} | 
    sed -E -e 's/^>(\S+).*/>\1 \1/' \
           -e 's/^>([0-9]+|[XY]) />chr\1 /' \
           -e 's/^>MT />chrM /' \
    > ${fasta_gencoded}
rm -f ${fasta_lcl}

# filter gtf for chr 21

gtf_filtered="gencode.v41.primary_assembly.annotation.filtered.gtf"
gtf_chr21="gencode.v41.primary_assembly.annotation.chr21.gtf"
grep -E "^chr21|^#" ${gtf_filtered} > ${gtf_chr21}

#### cellranger indexing ####

cellranger mkref \
  --genome=cellranger_index \
  --fasta=${fasta_gencoded} \
  --genes=${gtf_chr21} \
  --nthreads=8 \
  --memgb=16
