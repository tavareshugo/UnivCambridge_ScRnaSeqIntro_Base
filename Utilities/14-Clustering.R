#!/usr/bin/env Rscript
#SBATCH -J Clustering
#SBATCH -o Clustering.%j.out
#SBATCH -e Clustering.%j.err
#SBATCH --mincpus 8 
#SBATCH --mem=16G
#SBATCH --time=02:57:42

# This script generates the clustered sce objects for the Cluster Marker Genes
# and Differential Expression/Abundance sessions

library(scater) 
library(scran)
library(bluster)

# 500 cells
sce <- readRDS("data/R_objects/Caron_batch_corrected.500.rds")


out <- clusterSweep(reducedDim(sce, "corrected"),
                    BLUSPARAM = NNGraphParam(),
                    k = as.integer(c(10, 15, 20, 25, 30, 35)),
                    cluster.fun = c("walktrap", "louvain", "leiden"),
                    BPPARAM=BiocParallel::MulticoreParam(7))


colData(sce) <- cbind(colData(sce), DataFrame(out$clusters))

colLabels(sce) <- sce$k.25_cluster.fun.leiden

saveRDS(sce, "data/R_objects/Caron_clustered.500.rds")

# For the full data set we want just the PBMMC and ETV6-RUNX1 samples

sce <- readRDS("data/R_objects/Caron_filtered.full.rds")

sce <- sce[, sce$SampleGroup%in%c("PBMMC", "ETV6-RUNX1")]

# Normalise
set.seed(100) 
clust <- quickCluster(sce)
sce <- computePooledFactors(sce, cluster=clust, min.mean=0.1)
sce <- logNormCounts(sce)

# Batch correct
merge_order <- list(list(c("PBMMC_1a", "PBMMC_1b"), "PBMMC_2", "PBMMC_3"),
                    list("ETV6-RUNX1_1","ETV6-RUNX1_2",
                         "ETV6-RUNX1_3", "ETV6-RUNX1_4"))
set.seed(123)
sce_corrected <- quickCorrect(sce, 
                              PARAM = FastMnnParam(merge.order = merge_order),
                              batch = sce$SampleName)$corrected

reducedDim(sce, "corrected") <- reducedDim(sce_corrected, "corrected")
assay(sce, "reconstructed") <- assay(sce_corrected, "reconstructed")

# Cluster
out <- clusterSweep(reducedDim(sce, "corrected"),
                    BLUSPARAM = NNGraphParam(),
                    k = as.integer(c(10, 15, 20, 25, 30, 35)),
                    cluster.fun = c("walktrap", "louvain", "leiden"),
                    BPPARAM=BiocParallel::MulticoreParam(7))


colData(sce) <- cbind(colData(sce), DataFrame(out$clusters))

colLabels(sce) <- sce$k.25_cluster.fun.leiden

# Reduced Dims on Corrected

sce <- runTSNE(sce, dimred="corrected", 
               n_dimred=10, 
               BPPARAM = bpp, 
               name = "TSNE_corrected")
sce <- runUMAP(sce, dimred="corrected", 
               n_dimred=10, 
               BPPARAM = bpp, 
               name = "UMAP_corrected")

saveRDS(sce, "data/R_objects/Caron_clustered.PBMMCandETV6RUNX1.rds")

