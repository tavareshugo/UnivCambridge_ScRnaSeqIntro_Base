# load packages 
library(scater)
library(scran)
library(pheatmap)
library(tidyverse)


# Load data 
sce <- readRDS("R_objects/Caron_clustered.500.rds")

# Check rownames - these have been replaced with the gene symbol 
rownames(sce)[11:20]

# Check label column contains correct clustering
all(sce$k.25_cluster.fun.leiden==sce$label)

# Plot UMAP of clusters 
plotReducedDim(sce, 
               dimred = "UMAP_corrected",
               colour_by = "label", 
               text_by = "label")


# Plot CST3 expression on UMAP  
plotReducedDim(sce, 
               dimred = "UMAP_corrected",
               colour_by = "CST3", 
               text_by = "label", 
               by_exprs_values = "reconstructed",
               add_legend = FALSE)


## Score marker genes

markers <- scoreMarkers(sce, 
                        groups = sce$label, 
                        block = sce$SampleName)


# Score markers results for cluster 10 
c10_markers <- as.data.frame(markers[["10"]])
head(c10_markers)

## Selecting top marker genes

# Top marker genes for cluster 10 
c10_markers %>% 
  select(contains("cohen")) %>%
  filter(rank.logFC.cohen <= 5) %>%
  arrange(rank.logFC.cohen)


# Plot LYZ expression
p1 <- plotReducedDim(sce, 
                     dimred = "UMAP_corrected",
                     colour_by = "LYZ", 
                     text_by = "label")

p2 <- plotExpression(sce, features = "LYZ", x = "label")

p1 + p2



## -- Exercise 1 -- ############################################################



################################################################################

## Heatmaps of marker genes

# Cluster 10 top genes 
c10_top_genes <- c10_markers %>% 
  filter(rank.logFC.cohen <= 5)

# cluster 10 heatmap by cells 
plotHeatmap(sce, 
            features = rownames(c10_top_genes),
            order_columns_by = c("label", "SampleGroup"))


# cluster 10 heatmap by block 
plotGroupedHeatmap(sce, 
                   features = rownames(c10_top_genes),
                   group = "label",
                   block = "SampleGroup")


# heatmap with z score 
plotHeatmap(sce, 
            features = rownames(c10_top_genes),
            order_columns_by = c("label", "SampleGroup"),
                   scale = TRUE, 
            zlim = c(-3, 3))

plotGroupedHeatmap(sce, 
                   features = rownames(c10_top_genes),
                   group = "label",
                   block = "SampleGroup",
                   scale = TRUE, 
                   zlim = c(-3, 3))

## Adjusting the log fold change threshold

# c12 top markers 
c12_top_markers <- markers[["12"]] %>% 
  as.data.frame() %>% 
  filter(rank.logFC.cohen <= 2)
c12_top_markers


# c12 flt3 expression 
c12_top_markers["FLT3", ] %>%
  select(min.logFC.cohen, max.logFC.cohen)


# plot flt3 expression 
plotExpression(sce,
               features = "FLT3",
               x = "label")


# score markers with lfc threshold 
markers_lfc <- scoreMarkers(sce,
                           groups = sce$label,
                           block = sce$SampleName,
                           lfc = 2)


# c12 thresholded markers 
c12_markers_lfc <- as.data.frame(markers_lfc[["12"]])

c12_markers_lfc %>%
  select(contains("cohen")) %>%
  filter(rank.logFC.cohen <= 2)


# c12 FLT3 threholded rank 
c12_markers_lfc["FLT3",  c("rank.logFC.cohen")]

