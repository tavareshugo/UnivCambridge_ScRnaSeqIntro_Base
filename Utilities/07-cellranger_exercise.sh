#!/bin/bash
#SBATCH -J CellRanger_count 
#SBATCH -o cellranger_count.%j_%a.log
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

# reference genome
REF="references/cellranger_index"

# run cellranger count pipeline
cellranger count \
  --id=ETV6_RUNX1_rep1 \
  --sample=SRR9264343 \
  --transcriptome=${REF} \
  --fastqs=fastq_subsample \
	--localcores=16 \
	--localmem=64

# move to output directory
mkdir -p data/cellranger_exercise/
mv ETV6_RUNX1_rep1 data/cellranger_exercise/.

rm -f __ETV6_RUNX1_rep1.mro

echo "End: $(date)"
