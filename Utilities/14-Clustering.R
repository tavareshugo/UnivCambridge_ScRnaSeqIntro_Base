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
