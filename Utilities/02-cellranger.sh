#!/bin/bash
#SBATCH -J CellRanger_count 
#SBATCH -o cellranger_count.%j_%a.log
#SBATCH -a 2-13 
#SBATCH --nodes=1
#SBATCH --mincpus 16 
#SBATCH --mem=64G
#SBATCH --time=2-03:22:42

set -ex

echo "Start: $(date)"

cd data

# Cell Ranger software
crDir=`readlink -f software/cellranger-7.0.0`
export PATH=${crDir}:${PATH}

# get sample information from the table
ID=`head -n $SLURM_ARRAY_TASK_ID SraRunTable.txt | tail -n 1 | cut -d "," -f 1`

# reference genome
REF="references/refdata-gex-GRCh38-2020-A"

# run cellranger count pipeline
cellranger count \
  --id=${ID} \
  --sample=${ID} \
  --transcriptome=${REF} \
  --fastqs=fastq \
	--localcores=16 \
	--localmem=64

# move to output directory
mkdir -p data/cellranger/
mv ${ID} data/cellranger/

rm -f __*.mro

echo "End: $(date)"
