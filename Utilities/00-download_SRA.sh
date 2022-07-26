#!/bin/bash
#SBATCH -J public_data
#SBATCH -o download_SRA.%j_%a.log
#SBATCH -e download_SRA.%j_%a.err
#SBATCH -c 1
#SBATCH -t 12:00:00
#SBATCH --mem 2G
#SBATCH -a 2-13 # rows in the CSV file with SRA informatioj

set -ex

# activate conda environment
source activate sra

# get sample information from the table
ID=`head -n $SLURM_ARRAY_TASK_ID SraRunTable.txt | tail -n 1 | cut -d "," -f 1`

# create data folders
sraDir=data/sra
mkdir -p ${sraDir}
fqDir=data/fastq
mkdir -p ${fqDir}

# download file
prefetch --max-size 30G -O ${sraDir} ${ID}

# convert to fastq
fastq-dump -O ${fqDir} --gzip --split-files ${sraDir}/${ID}/${ID}.sra

# rename files for cellranger
# 1 = I1
# 2 = R1
# 3 = R2
mv ${fqDir}/${ID}_1.fastq.gz ${fqDir}/${ID}_S1_L001_I1_001.fast.gz
mv ${fqDir}/${ID}_2.fastq.gz ${fqDir}/${ID}_S1_L001_R1_001.fast.gz
mv ${fqDir}/${ID}_3.fastq.gz ${fqDir}/${ID}_S1_L001_R2_001.fast.gz