# load libraries ----

library(scater)
library(scran)
library(batchelor)
library(bluster)
library(tidyverse)
library(pheatmap)
library(BiocSingular)


# read the data ----

# SCE objects for the two technical replicates
sce_rep1 <- readRDS("Robjects/BC_sample1_dimred.rds")
sce_rep2 <- readRDS("Robjects/BC_sample2_dimred.rds")

# add information about which replicate each sample is from
# this is added as a new column in the colData DataFrame of the object
colData(sce_rep1)$batch <- "1"
colData(sce_rep2)$batch <- "2"

# DataFrame objects with mean-variance results from modelGeneVar()
gene_var_rep1 <- readRDS("Robjects/BC_dec1_dimred.rds")
gene_var_rep2 <- readRDS("Robjects/BC_dec2_dimred.rds")


# Data preparation - subset common genes ----

# identify genes common to both samples
common_genes <- intersect(rownames(sce_rep1), rownames(sce_rep2))

# Subset the SCE object
sce_rep1 <- sce_rep1[common_genes, ]
sce_rep2 <- sce_rep2[common_genes, ]

# Subset the mean-variance results
gene_var_rep1 <- gene_var_rep1[common_genes, ]
gene_var_rep2 <- gene_var_rep2[common_genes, ]


# Data preparation - rescale size factors ----

# rescale the size factors in each batch to account for sequencing depth differences
# this returns a list with two SCE objects
rescaled_size_factors <- multiBatchNorm(sce_rep1, sce_rep2)

# combine both objects in the list
sce <- cbind(rescaled_size_factors[[1]],
             rescaled_size_factors[[2]])


# Data preparation - select variable genes ----

# summarise the variance estimated across batches
gene_var_combined <- combineVar(gene_var_rep1, gene_var_rep2)

# choose HVGs
hvgs <- gene_var_combined$bio > 0
sum(hvgs) # number of HVGs selected


# Visualise uncorrected data ----

# run PCA - this adds a "PCA" slot to the SCE object
sce <- runPCA(sce, subset_row = hvgs)

# Define cell clusters - this will be covered in detail later
sce$cluster_uncorrected <- clusterCells(sce,
                                        use.dimred = "PCA")

# run and visualise TSNE
sce <- runTSNE(sce, dimred = "PCA", name = "TSNE_uncorrected")
plotReducedDim(sce, dimred = "TSNE_uncorrected",
               colour_by = "batch", text_by = "cluster_uncorrected")

# tabulate cells per cluster
table(Cluster = sce$cluster_uncorrected, Batch = sce$batch)


# Perform MNN correction ----

# run the MNN algorithm with
# d = 50 principal components
# k = 20 nearest neighbours
mnn_corrected <- fastMNN(sce,
                         batch = sce$batch,
                         d = 50, k = 20,
                         subset.row = hvgs)
mnn_corrected

# store the corrected values in a new reducedDim in the original sce object
reducedDim(sce, "corrected") <- reducedDim(mnn_corrected, "corrected")


# Visualise corrected data ----

# Define cell clusters based on the corrected matrix
sce$cluster_corrected <- clusterCells(sce,
                                      use.dimred = "corrected")

# run and visualise TSNE
sce <- runTSNE(sce, dimred = "corrected", name = "TSNE_corrected")
plotReducedDim(sce, dimred = "TSNE_corrected",
               colour_by = "batch", text_by = "cluster_corrected")

# visualise cells per cluster
as.data.frame(table(Cluster = sce$cluster_corrected, Batch = sce$batch)) %>%
  ggplot(aes(Cluster, Freq)) +
  geom_col(aes(fill = Batch)) +
  labs(title = "MNN-corrected data")

as.data.frame(table(Cluster = sce$cluster_uncorrected, Batch = sce$batch)) %>%
  ggplot(aes(Cluster, Freq)) +
  geom_col(aes(fill = Batch)) +
  labs(title = "Uncorrected data")


# Exercise 1 ----

# read ETV6_RUNX1 and PBMMC datasets
sce_all <- readRDS("Robjects/DataIntegration_all_sce_dimred.Rds")

# tabulate the number of cells per sample
table(sce_all$SampleName)

# obtain a batch-corrected SCE object
sce_all_corrected <- quickCorrect(FIXME)$corrected

# add the corrected matrix to the original object - to keep it all together
reducedDim(sce_all, "corrected") <- reducedDim(sce_all_corrected, "corrected")

# update SCE object with a UMAP on corrected data
# name this dimred as "UMAP_corrected" to distinguish it from the original
sce_all <- runUMAP(sce_all,
                   dimred = FIXME,
                   name = FIXME)

# visualise uncorrected UMAP
plotReducedDim(sce_all, dimred = "UMAP", colour_by = "SampleName")

# visualise corrected UMAP
plotReducedDim(sce_all, FIXME)


# Correction Diagnostics ----

# cluster cells
sce_all$cluster_corrected <- clusterCells(sce_all,
                                          use.dimred = "corrected")

# tabulate the number of cells from each batch per cluster
batch_per_cluster <- table(Cluster = sce_all$cluster_corrected, Batch = sce_all$SampleName)

# visualise distribution of samples per cluster
as.data.frame(batch_per_cluster) %>%
  ggplot(aes(Cluster, Freq)) +
  geom_col(aes(fill = Batch)) +
  labs(title = "Corrected data")
