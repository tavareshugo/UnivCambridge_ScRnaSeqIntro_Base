# load packages 
library(scater)
library(scran)
library(tidyverse)

# Load data 
sce <- readRDS("R_objects/Caron_clustered.500.rds")

## Score marker genes
markers <- scoreMarkers(sce, 
                        groups = sce$label, 
                        block = sce$SampleName)


# Exercise 1 -----

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
