#!/usr/bin/env Rscript
#SBATCH -J Clustering
#SBATCH -o Clustering.%j.out
#SBATCH -e Clustering.%j.err
#SBATCH --mincpus 8 
#SBATCH --mem=16G
#SBATCH --time=02:57:42

# This script generates the clustered sce objects for the
# Differential Expression/Abundance sessions

library(scater)
library(scran)
library(bluster)
library(batchelor)

#sce <- readRDS("data/R_objects/Caron_filtered.full.rds")
sce <- readRDS("R_objects/Caron_filtered.full.rds")
rownames(sce) <- uniquifyFeatureNames(rownames(sce), rowData(sce)$Symbol)


# First dataset ----

# For the first data set we want just the PBMMC and ETV6-RUNX1 samples
sce.1 <- sce[, sce$SampleGroup %in%c ("PBMMC", "ETV6-RUNX1")]

# Normalise
# Note - we should do this per batch and then multiBatchNorm - to fix in future
set.seed(100) 
clust.1 <- quickCluster(sce.1)
sce.1 <- computePooledFactors(sce.1, cluster=clust.1, min.mean=0.1)
sce.1 <- logNormCounts(sce.1)

# Batch correct
merge_order.1 <- list(list("PBMMC_1", "PBMMC_2", "PBMMC_3"),
                    list("ETV6-RUNX1_1","ETV6-RUNX1_2",
                         "ETV6-RUNX1_3", "ETV6-RUNX1_4"))
set.seed(123)
sce_corrected.1 <- quickCorrect(sce.1, 
                              PARAM = FastMnnParam(merge.order = merge_order.1),
                              batch = sce.1$SampleName)$corrected

reducedDim(sce.1, "corrected") <- reducedDim(sce_corrected.1, "corrected")
assay(sce.1, "reconstructed") <- assay(sce_corrected.1, "reconstructed")

# Cluster
out.1 <- clusterSweep(reducedDim(sce.1, "corrected"),
                      BLUSPARAM = NNGraphParam(),
                      k = as.integer(c(60)),
                      cluster.fun = c("leiden"),
                      BPPARAM=BiocParallel::MulticoreParam(7))

colData(sce.1) <- cbind(colData(sce.1), DataFrame(out.1$clusters))
colLabels(sce.1) <- sce.1$k.60_cluster.fun.leiden

# Reduced Dims on Corrected
sce.1 <- runTSNE(sce.1, dimred="corrected", 
                 n_dimred=10, 
                 BPPARAM = BiocParallel::MulticoreParam(7), 
                 name = "TSNE_corrected")
sce.1 <- runUMAP(sce.1, dimred="corrected", 
                 n_dimred=10, 
                 BPPARAM = BiocParallel::MulticoreParam(7), 
                 name = "UMAP_corrected")

# # manual cell type annotation based on these genes
# known_genes <- c(
#   "HBA1", # erythrocytes
#   "CST3", # monocytes
#   "CD3E", # T cells
#   "NKG7", # NK T cells
#   "CD79A",  # B cells
#   "MS4A1" # CD20 B cells
#   )
# plotDots(sce.1,
#          features = known_genes,
#          group = "k.60_cluster.fun.leiden", 
#          block = "SampleGroup",
#          scale = TRUE, center = TRUE, zlim = c(-3, 3))
# plotGroupedHeatmap(sce.1,
#          features = known_genes,
#          group = "label", 
#          block = "SampleGroup",
#          scale = TRUE, center = TRUE, zlim = c(-3, 3))

levels(colLabels(sce.1)) <- c("B (c1)", "B (c2)", "B (c3)", 
                              "T (c4)", "Erythrocytes (c5)", "CD20+ B (c6)", 
                              "B (c7)", "NK T (c8)", "Erythrocytes (c9)", 
                              "Mono (c10)", "B (c11)", "CD20+ B (c12)",
                              "CD20+ B (c13)", "T (c14)", "Erythrocytes (c15)", 
                              "Erythrocytes (c16)", "Mono (c17)")

saveRDS(sce.1, "R_objects/Caron_clustered.PBMMCandETV6RUNX1.rds")


# Second dataset ----

sce.2 <- sce[, sce$SampleGroup %in% c("PRE-T", "HHD")]

# Normalise
set.seed(100) 
clust.2 <- quickCluster(sce.2)
sce.2 <- computePooledFactors(sce.2, cluster=clust.2, min.mean=0.1)
sce.2 <- logNormCounts(sce.2)

# Batch correct
merge_order.2 <- list(list("PRE-T_1", "PRE-T_2"),
                      list("HHD_1", "HHD_2"))
set.seed(123)
sce_corrected.2 <- quickCorrect(sce.2, 
                                PARAM = FastMnnParam(merge.order = merge_order.2),
                                batch = sce.2$SampleName)$corrected

reducedDim(sce.2, "corrected") <- reducedDim(sce_corrected.2, "corrected")
assay(sce.2, "reconstructed") <- assay(sce_corrected.2, "reconstructed")

# Cluster
out.2 <- clusterSweep(reducedDim(sce.2, "corrected"),
                      BLUSPARAM = NNGraphParam(),
                      k = as.integer(c(60)),
                      cluster.fun = c("leiden"),
                      BPPARAM=BiocParallel::MulticoreParam(7))

colData(sce.2) <- cbind(colData(sce.2), DataFrame(out.2$clusters))
colLabels(sce.2) <- sce.2$k.60_cluster.fun.leiden


# Reduced Dims on Corrected
sce.2 <- runTSNE(sce.2, dimred="corrected", 
                 n_dimred=10, 
                 BPPARAM = BiocParallel::MulticoreParam(7), 
                 name = "TSNE_corrected")
sce.2 <- runUMAP(sce.2, dimred="corrected", 
                 n_dimred=10, 
                 BPPARAM = BiocParallel::MulticoreParam(7), 
                 name = "UMAP_corrected")

# # manual cell type annotation based on these genes
# # very rough annotation based on exploratory analysis with plotExpression()
# known_genes <- c(
#   "HBA1", # erythrocytes
#   "CST3", # monocytes
#   "CD3E", # T cells
#   "NKG7", # NK T cells
#   "CD79A",  # B cells
#   "MS4A1" # CD20 B cells
#   )
# plotGroupedHeatmap(sce.2, 
#                    features = known_genes,
#                    group = "k.60_cluster.fun.leiden",
#                    block = "SampleGroup", 
#                    scale = TRUE, center = TRUE, zlim = c(-3, 3))
# plotDots(sce.2, 
#                    features = known_genes,
#                    group = "k.60_cluster.fun.leiden",
#                    block = "SampleGroup", 
#                    scale = TRUE, center = TRUE, zlim = c(-3, 3))

levels(colLabels(sce.2)) <- c("Mono (c1)", "CD20+ B (c2)", 
                              "B (c3)", "B (c4)", "B (c5)", 
                              "T (c6)", "CD20 + B (c7)", 
                              "B (c8)", "B (c9)", "Erythrocytes (c10)")

saveRDS(sce.2, "R_objects/Caron_clustered.PRETandHHD.rds")
