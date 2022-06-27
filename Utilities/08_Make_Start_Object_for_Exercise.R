# Create sce object with all samples but just 500 cells be sample.

library(scran)
library(batchelor)
library(tidyverse)

sce <- readRDS("Robjects/Caron_filtered.rds")
sce2 <- readRDS("Robjects/DataIntegration_mnn.out.Rds")

# only keep 500 cells from each sample and discard one of the PBMMC_1 samples.
cells_to_keep <- colData(sce) %>% 
  as.data.frame() %>%
  select(Sample) %>%
  rownames_to_column("Barcode") %>% 
  group_by(Sample) %>%
  slice(sample(n(), 500)) %>%
  filter(Sample!="PBMMC_1a") %>%
  pull(Barcode)

# only keep the same genes as for the materials object
genes <- rownames(sce2)

# subset the object
sce <- sce[genes, cells_to_keep]

# Fix the PBMMC_1b sample name
sce$Sample <- str_remove(sce$Sample, "b$")

# normalise
clust <- quickCluster(sce)
sce <- computePooledFactors(sce, cluster=clust, min.mean=0.1)
sce <- logNormCounts(sce)

# run batch integration

quick.corrected <- quickCorrect(sce, batch = sce$Sample)

mnn_corrected <- quick.corrected$corrected

reducedDim(sce, "corrected") <- reducedDim(mnn_corrected, "corrected")

# Add reduced dimensions plots

set.seed(75390)

