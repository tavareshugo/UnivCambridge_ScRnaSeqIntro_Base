###############################################################################
# load packages
library(knitr)
library(scater)
library(scran)
library(sctransform)
library(tidyverse)
library(BiocParallel)
library(patchwork)
#################################################################################


#################################################################################

# set number of cores for processing
bpp <- MulticoreParam(7) 

# load data
sce <- readRDS("R_objects/Caron_filtered.500.rds")

# check data object
sce

# Check number of cells per sample
table(sce$SampleName)

#################################################################################

#################################################################################
# PBMMC_1 cell UMI counts distribution before normalization

# Subset the data we need from the sce object
oneSamTab <- colData(sce) %>% 
  as.data.frame() %>% 
  filter(SampleName == "PBMMC_1") %>% 
  dplyr::select(SampleName,Barcode, sum) %>% 
  mutate(cell_num = 1:n())

# Make a bar plot of the UMI counts per cell
p_before_nom <- ggplot(data=oneSamTab, aes(x=cell_num, y=sum)) +
  geom_bar(stat = 'identity') +
  labs( x= 'Cell Index',
        y='Cell UMI counts',
        title = "PBMMC_1: Before Normalization" ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size=20, color = 'red')
  )

# Display the plot
p_before_nom
#################################################################################



#################################################################################
# Deconvolution Normalisation
# Deconvolution : Compute size factors
set.seed(100)

# get clusters - deconvolution pools are made only within clusters
clust <- quickCluster(sce,BPPARAM=bpp)

# how many clusters do we have and how big are they?
table(clust)

# We have just the raw counts assay at the moment
assayNames(sce)

# We have no size factors for the cells at this stage
sizeFactors(sce)

# compute pooled size factors
sce <- computePooledFactors(sce,
			 clusters = clust,
			 min.mean = 0.1,
			 BPPARAM = bpp)

# We still have just the raw counts assay
assayNames(sce)

# We now have size factors which we can extract and summarise
deconv.sf <- sizeFactors(sce)
summary(deconv.sf)

#################################################################################


#################################################################################
# deconvolution vs library_size size factors

# Calculate library size factors
lib.sf <- librarySizeFactors(sce)

# Combine the size factor results with the sample group names and make a 
# scatter plot of the results
# This shows how the library size and deconvolution size factors differ
data.frame(LibrarySizeFactors = lib.sf, 
           DeconvolutionSizeFactors = deconv.sf,
			     SampleGroup = sce$SampleGroup) %>%
  ggplot(aes(x=LibrarySizeFactors, y=DeconvolutionSizeFactors)) +
      geom_point(aes(col=SampleGroup)) +
      geom_abline(slope = 1, intercept = 0)


#################################################################################

#################################################################################
# Deconvolution: scale and transform raw counts

# We still only have a raw counts assay
assayNames(sce)

# Use the deconvolution scale factors to generate an assay of normalised counts
# called "logcounts"
sce <- logNormCounts(sce)

assayNames(sce)

# PBMMC_1 cell UMI counts distribution after normalization
# Get the normalised counts, without doing the log transformation and
# sum them for each cell
norm_counts <- logNormCounts(sce,transform='none' ) %>% 
  assay('normcounts') %>% 
  colSums()

# Make a tibble dataframe of the total normalised counts per cell
norm_counts <- tibble(Barcode=names(norm_counts),
                      normCounts = log2(norm_counts)
                      )
head(norm_counts)

# Add the raw counts sum, cell numbers and SampleName to the tibble
norm_counts <- inner_join(norm_counts, oneSamTab, by='Barcode')

# Plot the raw and normalised cell UMI counts
p_after_norm <- ggplot(data=norm_counts, aes(x=cell_num, y=normCounts)) +
  geom_bar(stat = 'identity') +
  labs( x= 'Cell Index',
        y='Normalized Cell UMI counts',
        title = "PBMMC_1:After Normalization" ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size=20, color = 'red')
  )

p_before_nom + p_after_norm

# Let's separate out the scaling normalisation from the log transformation
# What do the un-normalised data look if we log them?
p_before_nom_nolog <- ggplot(data=oneSamTab, aes(x=cell_num, y=log2(sum))) +
  geom_bar(stat = 'identity') +
  labs( x= 'Cell Index',
        y='Cell UMI counts',
        title = "Logged raw counts" ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size=20, color = 'red')
  )

p_before_nom_nolog + p_after_norm

#################################################################################

#################################################################################
# Deconvolution: the effect of log transformation on mean and variance correlation for genes

# The log transformation is meant to reduce the correlation between mean and variance 
# for genes - has this worked?

# We can look at the relationship between the mean gene expression and variance
# for raw UMI counts, scaled counts and scaled, logged counts

# mean and variance for raw counts
mean <- rowMeans(assay(sce, "counts"))
var <- rowVars(assay(sce, "counts"))

# There is a strong linear relationship between mean and variance
plot(log(mean), log(var))
abline(a=1, b=1, col="red")

# Mean and variance for scaled counts
mean_scaled <- logNormCounts(sce,transform='none' ) %>% 
  assay('normcounts') %>% 
  rowMeans()
var_scaled <- logNormCounts(sce,transform='none' ) %>% 
  assay('normcounts') %>% 
  rowVars()

plot(log(mean_scaled), log(var_scaled))
abline(a=1, b=1, col="red")

# Mean and variance for scaled, log transformed counts
mean_norm <- rowMeans(assay(sce, "logcounts"))
var_norm <- rowVars(assay(sce, "logcounts"))

plot(mean_norm, var_norm)
abline(a=1, b=1, col="red")

# We see that the log transformation removes a large part of the relationship between
# mean and variance for gene expression values

# save sce object after normalisation
saveRDS(sce, "results/caron_normalised.rds")

#################################################################################

#################################################################################
#################################################################################
# Exercise 1
# Exercise: apply the deconvolution normalization on a single sample: ETV6-RUNX1_1 (aka GSM3872434).
#################################################################################
#################################################################################

#################################################################################
# sctransform: Variant Stabilising Transformation
# extract counts 
counts <- counts(sce)
class(counts)
dim(sce)
#################################################################################

#################################################################################
# add gene attributes 
gene_attr <- data.frame(mean = rowMeans(counts), 
                        detection_rate = rowMeans(counts > 0),
                        var = rowVars(counts)) %>% 
  mutate(log_mean = log10(mean)) %>% 
  mutate(log_var = log10(var))

dim(gene_attr)
head(gene_attr)
#################################################################################

#################################################################################
# add cell attributes
cell_attr <- data.frame(n_umi = colSums(counts),
                        n_gene = colSums(counts > 0))

dim(cell_attr)
head(cell_attr)
#################################################################################

#################################################################################
# over dispersion plot
ggplot(gene_attr, aes(log_mean, log_var)) + 
  geom_point(alpha=0.3, shape=16) +
  geom_abline(intercept = 0, slope = 1, color='red')
#################################################################################

#################################################################################
# detection rate plot
x = seq(from = -3, to = 2, length.out = 1000)
poisson_model <- data.frame(log_mean = x,
			    detection_rate = 1 - dpois(0, lambda = 10^x))

ggplot(gene_attr, aes(log_mean, detection_rate)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_line(data=poisson_model, color='red') +
  theme_gray(base_size = 8)
#################################################################################

#################################################################################
# n_gene vs n_umi (Total UMI counts in a cell)
ggplot(cell_attr, aes(n_umi, n_gene)) + 
  geom_point(alpha=0.3, shape=16) + 
  geom_density_2d(size = 0.3)
#################################################################################


#################################################################################
# sctransform (VST run)
set.seed(44)
vst_out <- vst(umi = counts,
               latent_var = c('log_umi'),
               return_gene_attr = TRUE,
               return_cell_attr = TRUE,
               verbosity = 2
  )
#################################################################################


#################################################################################
# sctransform plot model parameters
plot_model_pars(vst_out, verbosity = 1)
#################################################################################

#################################################################################
# sctransform model
vst_out$model_str
#################################################################################


#################################################################################
# Inspect model using two genes
rowData(sce) %>%
	as.data.frame %>%
	filter(Symbol %in% c('RPL10', 'HBB')) %>%
  select(ID, Symbol, Type, Chromosome)

ensId <- rowData(sce) %>%
	as.data.frame %>%
	filter(Symbol %in% c('RPL10', 'HBB')) %>%
  pull("ID")

plot_model(x = vst_out,
           umi = counts,
           goi = ensId,
           plot_residual = TRUE)

#################################################################################

#################################################################################
# properties of transformed data
# residual mean is centered around 0
ggplot(vst_out$gene_attr, aes(x = residual_mean)) +
	geom_histogram(binwidth=0.01)
#################################################################################

#################################################################################
# properties of transformed data
# residual variance is centered around 1
ggplot(vst_out$gene_attr, aes(residual_variance)) +
	geom_histogram(binwidth=0.1) +
	geom_vline(xintercept=1, color='red') +
	xlim(0, 10)
#################################################################################

#################################################################################
# properties of transformed data
# mean UMI vs variance plot
ggplot(vst_out$gene_attr, aes(x = log10(gmean), y = residual_variance)) +
       geom_point(alpha=0.3, shape=16)

#################################################################################

#################################################################################
# top genes with high residual variance
vst_out$gene_attr %>%
  arrange(desc(residual_variance)) %>% 
	top_n(n = 10) %>%
	mutate(across(where(is.numeric), round, 2)) %>% 
  rownames_to_column("ID") %>%
  left_join(as.data.frame(rowData(sce))[,c("ID", "Symbol")], "ID")
#################################################################################

#################################################################################
# add_vst to sce
keepGenes <- rownames(sce)%in%rownames(vst_out$y)
sce <- sce[keepGenes,]
vstMat <- as(vst_out$y[rownames(sce),], "dgCMatrix")

assay(sce, "sctrans_norm", withDimnames=FALSE) <- vstMat
#################################################################################


#################################################################################
#################################################################################
# Exercise 2
# Exercise: apply the sctransform VST normalisation on a single sample: ETV6-RUNX1_1 (aka SRR9264343).
#################################################################################
#################################################################################
