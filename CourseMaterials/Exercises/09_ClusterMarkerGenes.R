# Load packages ----
library(scater)
library(scran)
library(pheatmap)
library(tidyverse) # always load tidyverse after other packages


# Read data ----

# read object (continued from the previous section)
sce <- readRDS("Robjects/Caron_clustering_material.rds")

# visualise cluster assignments on the corrected data
plotReducedDim(sce,
               dimred = "UMAP_corrected",
               colour_by = "louvain15",
               text_by = "louvain15")

# visualise a previously known marker gene (for monocytes)
plotReducedDim(sce,
               dimred = "UMAP_corrected",
               colour_by = "CST3",
               text_by = "louvain15",
               by_exprs_values = "logcounts")


# Marker gene identification ----

# identify marker genes
# by default the function uses "logcounts" as the assay (see help)
markers <- scoreMarkers(sce,
                        groups = sce$louvain15,
                        block = sce$SampleName)

# returns a list of length equal to the number of clusters
markers

# check the result of a particular cluster
markers[[8]]

# extract results for one of the clusters
c8_markers <- markers[["8"]] %>%
  as.data.frame()

# look at top-ranked genes
c8_markers %>%
  select(contains("cohen")) %>%
  filter(rank.logFC.cohen <= 5) %>%
  arrange(rank.logFC.cohen)

# visualise one of the top genes on the MNN-corrected UMAP
plotReducedDim(sce,
               dimred = "UMAP_corrected",
               colour_by = "LYZ",
               text_by = "louvain15")

# visualise the logcounts distribution for this gene
plotExpression(sce,
               features = "LYZ",
               x = "louvain15")


# Exercise ----

# visualise CD3D (T cell marker)
plotReducedDim(sce,
               dimred = "UMAP_corrected",
               colour_by = "CD3D",
               text_by = "louvain15")

plotExpression(sce,
               features = "CD3D",
               x = "louvain15")

# extract results from cluster 4 and convert to data.frame
c4_markers <- FIXME

# filter the data.frame using your choice of ranking statistic
# `rank.logFC.detected` or `rank.logFC.cohen`
# or a combination of both!
c4_markers %>%
  filter(FIXME)

# visualise the expression of genes that seem interesting from your filters
plotExpression(sce,
               features = FIXME,
               x = "louvain15")


# Heatmaps ----

# select some top genes for cluster 8
c8_top10 <- c8_markers %>%
  filter(rank.logFC.cohen <= 10)

# heatmap of expression for each cell
plotHeatmap(sce,
            features = rownames(c8_top10),
            order_columns_by = c("louvain15", "SampleGroup"))

# heatmap of expression with average per cluster
plotGroupedHeatmap(sce,
                   features = rownames(c8_top10),
                   group = "louvain15",
                   block = "SampleGroup")

# heatmap of Z-scores for each cell
plotHeatmap(sce,
            features = rownames(c8_top10),
            order_columns_by = c("louvain15", "SampleGroup"),
            scale = TRUE, zlim = c(-3, 3))

# heatmap of Z-scores averaged per cluster
plotGroupedHeatmap(sce,
                   features = rownames(c8_top10),
                   group = "louvain15",
                   block = "SampleGroup",
                   scale = TRUE, zlim = c(-3, 3))


# LFC threshold ----

# genes with rank 1 in cluster 8
c8_markers %>%
  filter(rank.logFC.cohen == 1) %>%
  select(contains("cohen"))

# plot expression of FCGR3A
plotExpression(sce,
               features = "FCGR3A",
               x = "louvain15")

# run gene marker analysis using a stricter LFC threshold
markers_lfc <- scoreMarkers(sce,
                            groups = sce$louvain15,
                            block = sce$SampleName,
                            lfc = 2)

# extract the results for cluster 8
c8_markers_lfc <- markers_lfc[["8"]] %>% as.data.frame()

# check top 5 ranked genes
c8_markers_lfc %>%
  select(contains("cohen")) %>%
  filter(rank.logFC.cohen <= 5)

# check new rank for FCGR3A
c8_markers_lfc["FCGR3A", c("rank.logFC.cohen")]

# we could have also used other ranks to eliminate this gene in the original analysis
c8_markers %>%
  filter(rank.logFC.cohen == 1) %>%
  select(contains("rank"))
