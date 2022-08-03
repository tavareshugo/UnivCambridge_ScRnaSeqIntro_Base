#!/bin/bash
#SBATCH -J CellRanger_count 
#SBATCH -o CellRanger_count.%j_%a.out
#SBATCH -e CellRanger_count.%j_%a.err
#SBATCH -a 2-13 
#SBATCH --mincpus 16 
#SBATCH --mem=64G
#SBATCH --time=4-03:22:42

set -ex

echo "Start: $(date)"


# Cell Ranger software
crDir=`readlink -f data/software/cellranger-7.0.0`
export PATH=${crDir}:${PATH}

cd data

# get sample information from the table
ID=`head -n $SLURM_ARRAY_TASK_ID SraRunTable.txt | tail -n 1 | cut -d "," -f 1`

# reference genome
REF="references/refdata-gex-GRCh38.p13-Gencode.v41"

# run cellranger count pipeline
cellranger count \
  --id=${ID} \
  --sample=${ID} \
  --transcriptome=${REF} \
  --fastqs=fastq \
	--localcores=16 \
	--localmem=64

# move to output directory
mkdir -p cellranger/
mv ${ID} cellranger/

rm -f __*.mro

echo "End: $(date)"
