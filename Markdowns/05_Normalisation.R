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
# load data
bpp <- MulticoreParam(7) # set number of cores

sce <- readRDS("R_objects/Caron_filtered.500.rds")
sce
table(sce$SampleName)
#################################################################################

#################################################################################
# PBMMC_1 sample cell UMI counts distribution before normalization
oneSamTab <- colData(sce) %>% 
  as.data.frame() %>% 
  filter(SampleName == "PBMMC_1") %>% 
  dplyr::select(SampleName,Barcode, sum) %>% 
  mutate(cell_num = 1:n())

p_before_nom <- ggplot(data=oneSamTab, aes(x=cell_num, y=sum)) +
  geom_bar(stat = 'identity') +
  labs( x= 'Cell Index',
        y='Cell UMI counts',
        title = "PBMMC_1: Before Normalization" ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size=20, color = 'red')
  )

p_before_nom
#################################################################################



#################################################################################
# Deconvolution Normalisation
# Deconvolution : Compute size factors
set.seed(100)
# get clusters 
clust <- quickCluster(sce,BPPARAM=bpp)
table(clust)

assayNames(sce)

sizeFactors(sce)
# compute pooled size factors
sce <- computePooledFactors(sce,
			 clusters = clust,
			 min.mean = 0.1,
			 BPPARAM = bpp)

assayNames(sce)
deconv.sf <- sizeFactors(sce)
summary(deconv.sf)
#################################################################################


#################################################################################
# deconvolution vs library_size size factors
lib.sf <- librarySizeFactors(sce)

data.frame(LibrarySizeFactors = lib.sf, 
           DeconvolutionSizeFactors = deconv.sf,
			     SampleGroup = sce$SampleGroup) %>%
  ggplot(aes(x=LibrarySizeFactors, y=DeconvolutionSizeFactors)) +
      geom_point(aes(col=SampleGroup)) +
      geom_abline(slope = 1, intercept = 0)

#################################################################################

#################################################################################
# Deconvolution: scale and transform raw counts and save object
assayNames(sce)

sce <- logNormCounts(sce)

assayNames(sce)

# PBMMC_1 sample cell UMI counts distribution after normalization
norm_counts <- logNormCounts(sce,transform='none' ) %>% 
  assay('normcounts') %>% 
  colSums()

norm_counts <- tibble(Barcode=names(norm_counts),
                      normCounts = log2(norm_counts)
                      )
head(norm_counts)

norm_counts <- inner_join(norm_counts, oneSamTab, by='Barcode')


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

p_after_norm

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
