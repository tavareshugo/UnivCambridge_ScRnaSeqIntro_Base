## ----load_packages-------------------------------------------------------------
library(scater)
library(scran)
library(PCAtools)
library(tidyverse)

## ----read_data----------------------------------------------------------------
sce <- readRDS("R_objects/Caron_normalized.500.rds")
sce


## ----set_rownames_to_symbols---------------------------------------------------------------
rownames(sce) <- uniquifyFeatureNames(rownames(sce), rowData(sce)$Symbol)


## ----mean_variance_model-------------------------------------------------------------------
gene_var <- modelGeneVar(sce)

gene_var


## ----mean_variance_plot--------------------------------------------------------------------
gene_var %>% 
  as.data.frame() %>% 
  ggplot(aes(mean, total)) +
  geom_point() +
  geom_line(aes(y = tech), colour = "dodgerblue", size = 1) +
  labs(x = "Mean of log-expression", y = "Variance of log-expression")

plot(gene_var$mean,gene_var$total)
lines(gene_var$mean,gene_var$tech, type = "p", col = "green")
dev.off()


## ----get_hvgs------------------------------------------------------------------------------
hvgs <- getTopHVGs(gene_var, prop=0.1)
length(hvgs)
hvgs[1:10]


## ----plot_count_HVGtop20-------------------------------------------------------------------
plotExpression(sce, features = hvgs[1:20], point_alpha = 0.05)


## ----pca_computation-----------------------------------------------------------------------
sce <- runPCA(sce, subset_row = hvgs)
sce


## ----view_reducedDim-----------------------------------------------------------------------
reducedDim(sce, "PCA")[1:10, 1:5]


## ----scree_plot----------------------------------------------------------------------------
percent.var <- attr(reducedDim(sce), "percentVar")
plot(percent.var, log="y", xlab="PC", ylab="Variance explained (%)")


## ----pca_plot_1----------------------------------------------------------------------------
plotReducedDim(sce, dimred = "PCA", colour_by = "SampleName")


## ----sce_pca_plotReducedDim----------------------------------------------------------------
plotReducedDim(sce, dimred = "PCA", ncomponents = 3, colour_by = "SampleName")


## ----ggcells-------------------------------------------------------------------------------
ggcells(sce, aes(x = PCA.1, y = PCA.2, colour = SampleName)) +
  geom_point(size = 0.5) +
  facet_wrap(~ SampleName) +
  labs(x = "PC1", y = "PC2", colour = "Sample")


## ----explanatory_variables-----------------------------------------------------------------
explain_pcs <- getExplanatoryPCs(sce,
                                variables = c("sum",
                                              "detected",
                                              "SampleGroup",
                                              "SampleName",
                                              "subsets_Mito_percent")
                                )

plotExplanatoryPCs(explain_pcs/100)


## ----explanatory_variables_correlation_plot------------------------------------------------
plotExplanatoryVariables(sce,
                         variables = c(
                           "sum",
                           "detected",
                           "SampleGroup",
                           "SampleName",
                           "subsets_Mito_percent"
                         ))




## ----choose_elbow--------------------------------------------------------------------------
chosen_elbow <- findElbowPoint(percent.var)
chosen_elbow


## ----plot_eblow_choice---------------------------------------------------------------------
plot(percent.var)
abline(v=chosen_elbow, col="dodgerblue")


## ----denoisePCA-------------------------------------------------------------
sce.denoised <- denoisePCA(sce, technical = gene_var)


## ----denoisePCA_number_of_PCs--------------------------------------------------------------
ncol(reducedDim(sce.denoised, "PCA"))


## ----basic_tsne----------------------------------------------------------------------------
sce <- runTSNE(sce)

plotTSNE(sce, colour_by = "SampleName")

## ----Exercise 1-----------------------------------------------------------------
## # Run t-SNE ----
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
## FIXME
## 
## # Part B
## # Instead of colouring by SampleName, colour by expression of known cell markers
## # CD79A (B cells)
## # CST3 (monocytes)
## # CD3D (T cells)
## # HBA1 (erythrocytes)
## FIXME
## 
## # Part C
## # Facet these plots by SampleName to better understand where each marker is mostly expressed
## FIXME
## 
## # Part D
## # Explore different perplexity values (for example 5 and 500)
## # Do you get tighter or looser clusters?
## FIXME






 
 















## ----Exercise 1 solution--------------------------------------------------------
## # Run t-SNE ----
## 
## # add the t-SNE result to the reducedDim slot of the SCE object
## # we name this reducedDim "TSNE_perplex50"
## # we set perplexity = 50 (which is the default if we don't specify it)
## # we run t-SNE based on the PCA we ran previously
## set.seed(123) # set a random seed to ensure reproducibility
## sce <- runTSNE(sce,
##                name = "TSNE_perplex50",
##                perplexity = 50,
##                dimred = "PCA",
##                n_dimred = 10)
## 
## # Make a custom visualisation using ggcells
## ggcells(sce, aes(x = TSNE_perplex50.1, y = TSNE_perplex50.2,
##                  colour = SampleName)) +
##   geom_point()
## 
## # Re-run the algorithm but change the random seed number.
## # Do the results change dramatically between runs?
## set.seed(321)
## sce <- runTSNE(sce,
##                name = "TSNE_perplex50_seed321",
##                perplexity = 50,
##                dimred = "PCA",
##                n_dimred = 10)
## 
## ggcells(sce, aes(x = TSNE_perplex50_seed321.1, y = TSNE_perplex50_seed321.2,
##                  colour = SampleName)) +
##   geom_point()
## 
## 
## # Modify the visualisation to colour the points based on logcounts of known cell markers
## # CD79A (B cells)
## # CST3 (monocytes)
## # CD3D (T cells)
## # HBA1 (erythrocytes)
## ggcells(sce, aes(x = TSNE_perplex50_seed321.1, y = TSNE_perplex50_seed321.2,
##                  colour = CD79A)) +
##   geom_point()
## 
## # Facet these plots by SampleName to better understand where each marker is mostly expressed
## ggcells(sce, aes(x = TSNE_perplex50_seed321.1, y = TSNE_perplex50_seed321.2,
##                  colour = CD79A)) +
##   geom_point() +
##   facet_wrap(~ SampleName)
## 
## # Explore different perplexity values (for example 5 and 500)
## # Do you get tighter or looser clusters?
## set.seed(321)
## sce <- runTSNE(sce,
##                name = "TSNE_perplex5",
##                perplexity = 5,
##                dimred = "PCA",
##                n_dimred = 10)
## sce <- runTSNE(sce,
##                name = "TSNE_perplex500",
##                perplexity = 500,
##                dimred = "PCA",
##                n_dimred = 10)
## 
## # visualise
## ggcells(sce, aes(x = TSNE_perplex5.1, y = TSNE_perplex5.2,
##                  colour = SampleName)) +
##   geom_point() +
##   labs(title = "Perplexity 5")
## 
## ggcells(sce, aes(x = TSNE_perplex50.1, y = TSNE_perplex50.2,
##                  colour = SampleName)) +
##   geom_point() +
##   labs(title = "Perplexity 50")
## 
## ggcells(sce, aes(x = TSNE_perplex500.1, y = TSNE_perplex500.2,
##                  colour = SampleName)) +
##   geom_point() +
##   labs(title = "Perplexity 500")


## ----runUMAP-------------------------------------------------------------------------------
set.seed(123)
sce <- runUMAP(sce)

plotUMAP(sce, colour_by = "SampleName")


## ----Exercise 2-----------------------------------------------------------------
## # Run UMAP ----
## 
## # Part A
## # run the UMAP with 50 neighbours
 set.seed(123) # set seed for reproducibility
 sce <- runUMAP(sce,
                name = "UMAP_neighbors50",
                dimred = "PCA",
                FIXME)
## 
## # Part B
## # visualise the resulting UMAP projection (colour cells by sample)
## FIXME
## 
## # Part C
## # run the UMAP with 5 and 500 neighbours and compare the results
## FIXME
## 
## # Part D
## # compare the UMAP projection with the t-SNE projections
## # would you prefer one over the other?
## FIXME


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
## ----Exercise 2 solution--------------------------------------------------------
## # Run UMAP ----
## 
## # run the UMAP with 50 neighbours
## set.seed(123) # set seed for reproducibility
## sce <- runUMAP(sce,
##                name = "UMAP_neighbors50",
##                dimred = "PCA",
##                n_neighbors = 50)
## 
## # visualise the resulting UMAP projection
## # colour cells by sample
## ggcells(sce, aes(x = UMAP_neighbors50.1, y = UMAP_neighbors50.2,
##                  colour = SampleName)) +
##   geom_point()
## 
## # run the UMAP with 5 and 500 neighbours and compare the results
## set.seed(123) # set seed for reproducibility
## sce <- runUMAP(sce,
##                name = "UMAP_neighbors5",
##                dimred = "PCA",
##                n_neighbors = 5)
## sce <- runUMAP(sce,
##                name = "UMAP_neighbors500",
##                dimred = "PCA",
##                n_neighbors = 500)
## 
## ggcells(sce, aes(x = UMAP_neighbors5.1, y = UMAP_neighbors5.2,
##                  colour = SampleName)) +
##   geom_point() +
##   labs(title = "Neighbours = 5")
## ggcells(sce, aes(x = UMAP_neighbors500.1, y = UMAP_neighbors500.2,
##                  colour = SampleName)) +
##   geom_point() +
##   labs(title = "Neighbours = 500")
## 
## # compare the UMAP projection with the t-SNE projections
## # would you prefer one over the other?
## ggcells(sce, aes(x = TSNE_perplex50.1, y = TSNE_perplex50.2,
##                  colour = SampleName)) +
##   geom_point()
## ggcells(sce, aes(x = UMAP_neighbors50.1, y = UMAP_neighbors50.2,
##                  colour = SampleName)) +
##   geom_point()

