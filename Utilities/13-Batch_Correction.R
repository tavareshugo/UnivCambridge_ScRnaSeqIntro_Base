#!/usr/bin/env Rscript
#SBATCH -J BatchCorrection
#SBATCH -o BatchCorrection.%j.out
#SBATCH -e BatchCorrection.%j.err
#SBATCH --mincpus 8 
#SBATCH --mem=16G
#SBATCH --time=02:57:42

# This script performs batch correction on the 500 cell per sample data set. 
# For the clustering demonstration we will only use the ETV6-RUNX1 and PBMMC
# samples, so the data set is first subset to these samples. After batch
# correction tSNE and UMAP dimension reduction are run using the corrected data.

library(scater) 
library(scran)
library(batchelor)
library(BiocParallel)

bpp <- MulticoreParam(8)

# Load data
sce <- readRDS("data/R_objects/Caron_dimRed.500.rds")
sce <- sce[,sce$SampleGroup%in%c("ETV6-RUNX1", "PBMMC")]

sce_corrected <- quickCorrect(sce, batch = sce$SampleName)$corrected

reducedDim(sce, "corrected") <- reducedDim(sce_corrected, "corrected")
assay(sce, "reconstructed") <- assay(sce_corrected, "reconstructed")

# Add reduced dimensions plots

set.seed(75390)
sce <- runTSNE(sce, dimred = "corrected", name = "TSNE_corrected")
sce <- runUMAP(sce, dimred = "corrected", name = "UMAP_corrected")

# save object
saveRDS(sce, "data/R_objects/Caron_batch_corrected.500.rds")
