# Create sce object with 7 samples and just 500 cells per sample.

library(scran)
library(scater)
library(batchelor)
library(bluster)
library(tidyverse)

sce <- readRDS("Robjects/DataIntegration_all_sce_dimred.Rds")
# 17700 genes and 3500 cells
# already has logcounts

# remove genes with 0 counts
keepGenes <- rowSums(counts(sce))>0
sce <- sce[keepGenes]
# 17486 genes

# run batch integration

sce_corrected <- quickCorrect(sce, batch = sce$SampleName)$corrected

reducedDim(sce, "corrected") <- reducedDim(sce_corrected, "corrected")
assay(sce, "reconstructed") <- assay(sce_corrected, "reconstructed")

# Add reduced dimensions plots

set.seed(75390)
sce <- runTSNE(sce, dimred = "corrected", name = "TSNE_corrected")
sce <- runUMAP(sce, dimred = "corrected", name = "UMAP_corrected")

saveRDS(sce, "Robjects/DataIntegration_mnn.Rds")
