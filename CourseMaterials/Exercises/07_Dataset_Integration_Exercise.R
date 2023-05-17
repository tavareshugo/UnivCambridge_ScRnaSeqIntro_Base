# Dataset Integration Exercise worksheet
library(scater)
library(scran)
library(batchelor)
library(bluster)
library(pheatmap)
library(magrittr)
set.seed(1704)

setwd("~/Course_Materials")

# load the data
sce_all <- readRDS("R_objects/Caron_dimRed.500.rds")

# tabulate the number of cells per sample
table(sce_all$SampleName)

### Add your own code to complete the following steps

# obtain a batch-corrected SCE object

# add the corrected matrix to the original object - to keep it all together

# add a tSNE using the corrected data

# visualise uncorrected tSNE

# visualise corrected tSNE