#!/usr/bin/env Rscript
#SBATCH -J ObjectForDI
#SBATCH -o ObjectForDI.%j.out
#SBATCH -e ObjectForDI.%j.err
#SBATCH --mincpus 8 
#SBATCH --mem=16G
#SBATCH --time=02:57:42

# This script starts with the all cells filtered data set. It extracts the two
# technical replicate PBMMC samples - SRR9264351 & SRR9264352. It then
# normalizes the counts, runs dimension reduction and clustering (default of
# clusterCells) on each of these samples separately.

library(scater) 
library(scran)
library(BiocParallel)

bpp <- MulticoreParam(8)

## Full data set 
# Load data
sce <- readRDS("data/R_objects/Caron_filtered.full.rds")

# Subset the object

prepSample <- function(sampleID, sceObj = sce){
    subsce <- sceObj[, sceObj$Sample==sampleID]
    # Remove the per cell QC metrics
    colData(subsce) <- colData(subsce)[,1:5]
    # Remove undetected genes
    keep <- rowSums(counts(subsce)) > 0
    subsce <- subsce[keep,]
    # Normalize the counts
    set.seed(100) 
    clust <- quickCluster(subsce) 
    subsce <- computePooledFactors(subsce, cluster=clust, min.mean=0.1)
    subsce <- logNormCounts(subsce)
    # Identify highly variable genes
    gene_var <- modelGeneVar(subsce)
    hvgs <- getTopHVGs(gene_var, prop=0.1)
    # PCA 
    subsce <- runPCA(subsce, subset_row = hvgs)
    # tSNE
    subsce <- runTSNE(subsce, dimred="PCA", n_dimred=10, BPPARAM = bpp)
    # UMAP
    subsce <- runUMAP(subsce, dimred="PCA", n_dimred=10, BPPARAM = bpp)
    # Clustering
    colLabels(subsce) <- clusterCells(subsce, use.dimred = "PCA")
    subsce
}

pbmmc1a <- prepSample("SRR9264351")
pbmmc1b <- prepSample("SRR9264352")

saveRDS(pbmmc1a, "data/R_objects/PBMMC_1a_dimRed.rds")
saveRDS(pbmmc1b, "data/R_objects/PBMMC_1b_dimRed.rds")




