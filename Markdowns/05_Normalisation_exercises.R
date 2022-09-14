# Normalisation Practical

# In the demonstration we ran the normalisation using just 500 cells per sample.
# For this practical you will carry out normalisation, but this time using the 
# entire data set. Some of the commands will take a little longer to run this 
# time.

# load_packages 
library(scater)
library(scran)
library(tidyverse)
library(BiocParallel)
 
bpp <- MulticoreParam(7)

## Prepare the data object

# load_data 
sce <- readRDS("R_objects/Caron_filtered.full.rds")


# subset_data 
etvr1 <- which(sce$SampleName=="ETV6-RUNX1_1")
sce <- sce[, etvr1]

## -- Exercise 1 -- ############################################################

# Now that we have reduced the number of cells in the data set, it may be that
# there are genes in the object that have not been detected in any of the 
# remaining cells. Filter the object to remove any genes that have not been 
# detected.

# PUT YOUR CODE HERE

# Q1. How many cells and genes are there in the ETV6_RUNX1_1 data set?

# PUT YOUR CODE HERE

# Q2. In how many cells has the gene ORC1 been detected with at least 1 UMI?

# PUT YOUR CODE HERE

################################################################################

## Normalisation by deconvolution

## -- Exercise 2 -- ############################################################

# Now normalise the data set using the deconvolution method. You will need to

# 1. cluster the cells - remember to set a seed so that your results are
#    reproducible
# 2. compute size factors using the deconvolution method  
# 3. log normalize the counts by applying the size factors  
# 4. check that the single cell experiment object contains a new "logcounts" 
#    assay

# PUT YOUR CODE HERE

################################################################################

## Normalisation with variance stabilising transformation

# extract_counts_matrix 
counts <- counts(sce)

## -- Exercise 3 -- ############################################################

# Now use the `vst` function to estimate model parameters and perform the
# variance stabilizing transformation.


# PUT YOUR CODE HERE

################################################################################

# check_vst_model 
print(vst_out$model_str)


## -- Exercise 4 -- ############################################################

# Use the `plot_model` function to inspect the genes 'RPL10' and 'FTL' to see if
# the fitted relationship between cell total UMI and gene expression, and to 
# check the residuals.

# PUT YOUR CODE HERE

################################################################################

# residual_mean_histogram
ggplot(vst_out$gene_attr, aes(residual_mean)) +
  geom_histogram(binwidth=0.01)


# residual_variance_histogram 
ggplot(vst_out$gene_attr, aes(residual_variance)) +
  geom_histogram(binwidth=0.1) +
  geom_vline(xintercept=1, color='red') +
  xlim(0, 10)


# residual_variance_v_mean_plot 
ggplot(vst_out$gene_attr,
       aes(log10(gmean), residual_variance)) +
       geom_point(alpha=0.3, shape=16) +
       geom_density_2d(size = 0.3)

