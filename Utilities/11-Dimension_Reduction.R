#!/usr/bin/env Rscript
#SBATCH -J DimensionReduction
#SBATCH -o DimensionReduction.%j.out
#SBATCH -e DimensionReduction.%j.err
#SBATCH --mincpus 8 
#SBATCH --mem=16G
#SBATCH --time=02:57:42

# This script performs feature selection and then runs PCA,
# tSNE and UMAP.

library(scater) 
library(scran)
library(BiocParallel)

bpp <- MulticoreParam(8)

## Full data set 
# Load data
sce <- readRDS("data/R_objects/Caron_normalized.full.rds")

# use common gene names instead of Ensembl gene IDs
rownames(sce) <- uniquifyFeatureNames(rownames(sce), rowData(sce)$Symbol)

# Identify highly variable genes
gene_var <- modelGeneVar(sce)
hvgs <- getTopHVGs(gene_var, prop=0.1)

# PCA 
sce <- runPCA(sce, subset_row = hvgs)

# tSNE
sce <- runTSNE(sce, dimred="PCA", n_dimred=10, BPPARAM = bpp)

# UMAP

sce <- runUMAP(sce, dimred="PCA", n_dimred=10, BPPARAM = bpp)

# save object

saveRDS(sce, "data/R_objects/Caron_dimRed.full.rds")

## 500 cell/sample data set
# Load data
sce <- readRDS("data/R_objects/Caron_normalized.500.rds")

# use common gene names instead of Ensembl gene IDs
rownames(sce) <- uniquifyFeatureNames(rownames(sce), rowData(sce)$Symbol)

# Identify highly variable genes
gene_var <- modelGeneVar(sce)
hvgs <- getTopHVGs(gene_var, prop=0.1)

# PCA
sce <- runPCA(sce, subset_row = hvgs)

# tSNE
sce <- runTSNE(sce, dimred="PCA", n_dimred=10)

# UMAP
sce <- runUMAP(sce, dimred="PCA", n_dimred=10)

# save object
saveRDS(sce, "data/R_objects/Caron_dimRed.500.rds")
