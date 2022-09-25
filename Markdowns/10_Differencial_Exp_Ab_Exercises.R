# Catch up by running this

library(scater)
library(scran)
library(batchelor)
library(edgeR)
library(tidyverse)
library(patchwork)
library(DT)
library(bluster)
library(BiocParallel)
library(miloR)

## ----exercise 1---------------------------------------------------------------
## # First load in the other two sample groups
## sce_PRET_HHD <- readRDS("R_objects/Caron_clustered.PRETandHHD.rds")
## 
## # replace the ensembl IDs with gene symbols where possible
## rownames(sce_PRET_HHD) <- uniquifyFeatureNames(rownames(sce_PRET_HHD), rowData(sce_PRET_HHD)$Symbol)
## 
## # check your sce object
## sce_PRET_HHD
## 
## # Part A
## # Re run the analysis using these new samples.
## 
## FIXME
## 
## # Part B
## # Looking at your results, which cluster has the most DEGs?
## 
## FIXME
## 
## # Part C
## # Which genes are sig. DE in that cluster and not the others?
## 
## FIXME
## 

# Code to catch up for exercise 2

milo <- readRDS("R_objects/10_milo_calcNhoodDistance.rds")

da_results <- testNhoods(milo, design = ~ SampleGroup, design.df = milo_design, reduced.dim = "corrected")

milo <- buildNhoodGraph(milo)

da_results <- annotateNhoods(milo, da_results, coldata_col = "label")

da_results$label <- ifelse(da_results$label_fraction < 0.7, "Mixed", da_results$label)

## ----exercise 2---------------------------------------------------------------
## # First load in the other two sample groups
## set.seed(42) # set your random seed for reproducibility
## 
## # Part A
## # rerun the grouping with values for max.lfc.delta of 1, 5, 25 and plot beeswarm for each
## 
## plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = FIXME) , group.by = "NhoodGroup") + ggtitle("max LFC delta=1")
## 
## 
## # Part B
## # Using the max.lfc.delta you preferred from Part A now alter the overlap value trying 1, 3, 10.
## 
## FIXME
## 