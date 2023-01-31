# Setup ----

# load packages
library(scater)
library(scran)
library(batchelor)
library(bluster)
library(edgeR)
library(tidyverse)

# load the SCE object
sce <- readRDS("R_objects/Caron_clustered.PBMMCandETV6RUNX1.rds")

# plot UMAP done on the batch-corrected data
plotReducedDim(sce, dimred = "UMAP_corrected", 
               colour_by = "label", 
               text_by = "label")


# Create pseudo-bulk samples -----

# tabulate number of cells per label + sample combination
table(sce$label, sce$SampleName)

# sum counts across cells - by label and sample
summed <- aggregateAcrossCells(sce, 
                               ids = colData(sce)[, c("label", "SampleName")])

# the output is a new SCE object with aggregated counts matrix
summed


# Run DE analysis -----

# "ncells" is added to our colData
colData(summed)[, c("SampleName", "ncells")]

# use this to filter our pseudo-bulk object
summed <- summed[, summed$ncells >= 20]

# perform differential analysis for each label
# this uses edgeR package behind the scenes.
de_results <- pseudoBulkDGE(summed, 
                            label = summed$label,
                            design = ~ SampleGroup, 
                            coef = "SampleGroupPBMMC",
                            condition = summed$SampleName)

# the model matrix being used by edgeR
head(model.matrix(~ SampleGroup, data = colData(summed)))

# the results come as a list for each cell label
de_results

# extract one of the tables from the list
b_c1_res <- de_results[[1]]

# extract the edgeR object used for differential expression
b_c1_dgelist <- metadata(b_c1_res)$y

# plot mean-CV relationship across genes
plotBCV(b_c1_dgelist)

# plot mean-difference plot to assess library size normalisation
# we expect simmetric distribution centered around zero
par(mfrow=c(2,4))
for (i in seq_len(ncol(b_c1_dgelist))) {
  plotMD(b_c1_dgelist, column=i)
  abline(h = 0, col = "salmon", lwd = 2)
}

# plot MDS 
plotMDS(cpm(b_c1_dgelist, log = TRUE), 
        col = ifelse(b_c1_dgelist$samples$SampleGroup == "PBMMC", "tomato", "steelblue"))


# Explore DEGs ----

# identify DEGs based on FDR threshold
is_de <- decideTestsPerLabel(de_results, threshold = 0.05)

# this outputs a large matrix
is_de[1:10, 1:5]

# summarise the results
summarizeTestsPerLabel(is_de)

# filter the results table for DEGs
b_c1_res %>% 
  as.data.frame() %>% 
  arrange(FDR) %>%
  head()

# remove size factors from the object (otherwise plotting function complains)
# and add normalised log counts to the object
sizeFactors(summed) <- NULL
summed <- logNormCounts(summed)

# visualise from our summed single cell object
plotExpression(summed, 
               features = "HTR1F",
               x = "SampleName", colour_by="SampleGroup", 
               other_fields="label") + 
  facet_wrap(~ label) +
  scale_x_discrete(guide = guide_axis(angle = 45))


# Exercise -----

# First load in the other two sample groups
sce_PRET_HHD <- readRDS("R_objects/Caron_clustered.PRETandHHD.rds")

# check your sce object
sce_PRET_HHD

# Part A
# Run the analysis using these new samples

FIXME

# Part B
# Determine which cluster has the most DEGs

FIXME

# Part C
# Visualise the expression of one of the top genes for that cluster

FIXME
