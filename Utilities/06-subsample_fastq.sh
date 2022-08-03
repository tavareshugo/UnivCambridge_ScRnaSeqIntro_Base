#!/bin/bash
#SBATCH -J SubSampleFQ
#SBATCH -o SubSampleFQ.%j.log
#SBATCH --nodes=1
#SBATCH --mincpus 1 
#SBATCH --mem=4G
#SBATCH --time=2-03:22:42

set -ex

cd data

mkdir -p fastq_subsample
for fqFile in fastq/SRR9264343_S1_L001_*; do
    filename=`basename ${fqFile}`
    zcat ${fqFile} | head -n 4000000 | gzip -c > fastq_subsample/${filename}
done