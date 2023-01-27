# Setup & Data ----

# load packages
library(scater)
library(scran)
library(tidyverse)
library(patchwork)

# read single cell object
sce <- readRDS("R_objects/Caron_clustered.500.rds")

# check labels are set to our clusters
all(sce$k.25_cluster.fun.leiden == sce$label)

# visualise UMAP of our clusters
plotReducedDim(sce, 
               dimred = "UMAP_corrected",
               colour_by = "label", 
               text_by = "label")

# visualise specific marker
plotReducedDim(sce, 
               dimred = "UMAP_corrected",
               colour_by = "CST3", 
               text_by = "label", 
               by_exprs_values = "reconstructed",
               add_legend = FALSE)


# Score markers ----

# calculate pairwise marker gene statistics
markers <- scoreMarkers(sce, 
                        groups = sce$label, 
                        block = sce$SampleName)

# extract results for cluster 11
c11_markers <- as.data.frame(markers[["11"]])
head(c11_markers)

# filter markers based on rank statistics
c11_markers %>% 
  select(contains("cohen")) %>%
  filter(rank.logFC.cohen <= 5) %>%
  arrange(rank.logFC.cohen)

# visualise one of the markers
p1 <- plotReducedDim(sce, 
                     dimred = "UMAP_corrected",
                     colour_by = "LYZ", 
                     text_by = "label")

p2 <- plotExpression(sce, features = "LYZ", x = "label")

p1 + p2


# Exercise 1 ----

# CD3D suggests cluster 6 and 7 are T cells
plotReducedDim(sce, 
               dimred = "UMAP_corrected",
               colour_by = "CD3D", 
               text_by = "label")

plotExpression(sce, 
               features = "CD3D", 
               x = "label")

# Confirm this by identifying other genes that differentiate
# these two clusters from the rest of the cells.

# 1. Extract results for cluster 6 and convert it to data.frame
# 2. Filter the data.frame using your choice of ranking statistic -
#  `rank.logFC.detected` or `rank.logFC.cohen` or a combination of both.
# 3. Visualise the expression of genes that seem interesting from your filters.










# Heatmaps ----

# get top-ranked markers for cluster 11
c11_top_genes <- c11_markers %>% 
  filter(rank.logFC.cohen <= 5)

# visualise their expression as a heatmap
plotHeatmap(sce, 
            features = rownames(c11_top_genes),
            order_columns_by = c("label", "SampleGroup"))

# heatmap average per group (cluster)
plotGroupedHeatmap(sce, 
                   features = rownames(c11_top_genes),
                   group = "label",
                   block = "SampleGroup")

# scaled heatmap (z-scores)
plotHeatmap(sce, 
            features = rownames(c11_top_genes),
            order_columns_by = c("label", "SampleGroup"),
            scale = TRUE, 
            center = TRUE,
            zlim = c(-3, 3))

plotGroupedHeatmap(sce, 
                   features = rownames(c11_top_genes),
                   group = "label",
                   block = "SampleGroup",
                   scale = TRUE, 
                   center = TRUE,
                   zlim = c(-3, 3))


# Adjusting log-fold change ----

c12_top_markers <- markers[["12"]] %>% 
  as.data.frame() %>% 
  filter(rank.logFC.cohen <= 2)
c12_top_markers

c12_top_markers["FLT3", ] %>%
  select(min.logFC.cohen, max.logFC.cohen)

plotExpression(sce,
               features = "FLT3",
               x = "label")

markers_lfc <- scoreMarkers(sce,
                           groups = sce$label,
                           block = sce$SampleName,
                           lfc = 2)

c12_markers_lfc <- as.data.frame(markers_lfc[["12"]])

c12_markers_lfc %>%
  select(contains("cohen")) %>%
  filter(rank.logFC.cohen <= 2)

c12_markers_lfc["FLT3",  c("rank.logFC.cohen")]


# Annotation labels ----

# loop through list of marker genes and extract top-ranked gene names
top_markers_all <- lapply(markers, function(x){
  x %>% 
    as.data.frame() %>% 
    filter(rank.logFC.cohen < 10) %>% 
    rownames()
})

# examining this list reveals several known markers of immune cells
top_markers_all

# cell type specific genes
known_genes <- c(
  "HBA1", # erythrocytes
  "CST3", # monocytes
  "CD3E", # T cells
  "NKG7", # NK T cells
  "CD79A",  # B cells
  "MS4A1" # CD20 B cells
  )

# violin plot
plotExpression(sce, x = "label", features = known_genes)

# scaled heatmap of expression
plotGroupedHeatmap(sce, 
                   features = known_genes,
                   group = "label",
                   block = "SampleGroup", 
                   scale = TRUE, center = TRUE, 
                   zlim = c(-3, 3))

# re-label the cells - original cluster in parenthesis
levels(colLabels(sce)) <- c("B (c1)", "B (c2)", 
                            "B (c3)", "B (c4)",
                            "CD20+ B (c5)", 
                            "T (c6)", "NK T (c7)", 
                            "Erythrocytes (c8)", "Erythrocytes (c9)", 
                            "Erythrocytes c(10)",
                            "Monocytes (c11)", "B (c12)")

# visualise UMAP with new labels
plotReducedDim(sce, dimred = "UMAP_corrected", 
               colour_by = "label", text_by = "label")
