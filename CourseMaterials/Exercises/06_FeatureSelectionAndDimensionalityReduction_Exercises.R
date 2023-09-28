# run this code to catch up to the demonstration
library(scater)
library(scran)
library(PCAtools)
library(tidyverse)

sce <- readRDS("R_objects/Caron_normalized.500.rds")

rownames(sce) <- uniquifyFeatureNames(rownames(sce), rowData(sce)$Symbol)

gene_var <- modelGeneVar(sce)

hvgs <- getTopHVGs(gene_var, prop=0.1)

sce <- runPCA(sce, subset_row = hvgs)

# Exercise 1

## # Run t-SNE 
## 
## # add the t-SNE result to the reducedDim slot of the SCE object
## # we name this reducedDim "TSNE_perplex50"
## # we set perplexity = 50 (which is the default if we don't specify it)
## # we run t-SNE based on the PCA we ran previously
## # we will use the first 10 principle components

set.seed(123) # set a random seed to ensure reproducibility
sce <- runTSNE(sce,
               name = "TSNE_perplex50",
               perplexity = 50,
               dimred = "PCA",
               n_dimred = 10)
## 
## # Make a custom visualisation using ggcells
ggcells(sce, aes(x = TSNE_perplex50.1, y = TSNE_perplex50.2,
                  colour = SampleName)) +
        geom_point()

## 
## # Part A
## # Re-run the algorithm but change the random seed number.
## # Do the results change dramatically between runs?

YOUR CODE HERE

 ## 
## # Part B
## # Instead of colouring by SampleName, colour by expression of known cell markers
## # CD79A (B cells)
## # CST3 (monocytes)
## # CD3D (T cells)
## # HBA1 (erythrocytes)

YOUR CODE HERE

## 
## # Part C
## # Facet these plots by SampleName to better understand where each marker is mostly expressed

YOUR CODE HERE

## 
## # Part D
## # Explore different perplexity values (for example 5 and 500)
## # Do you get tighter or looser clusters?

YOUR CODE HERE

# Exercise 2

## # Run UMAP 
## 
## # Part A
## # run the UMAP with 50 neighbours
set.seed(123) # set seed for reproducibility
sce <- runUMAP(sce,
               name = "UMAP_neighbors50",
               dimred = "PCA",
               ********YOUR CODE HERE********)

## 
## # Part B
## # visualise the resulting UMAP projection (colour cells by sample)

YOUR CODE HERE

## 
## # Part C
## # run the UMAP with 5 and 500 neighbours and compare the results

YOUR CODE HERE

## 
## # Part D
## # compare the UMAP projection with the t-SNE projections
## # would you prefer one over the other?

YOUR CODE HERE

