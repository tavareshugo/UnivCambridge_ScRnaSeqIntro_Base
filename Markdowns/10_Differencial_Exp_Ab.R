## ----load_packages-----------------------------------------------------------
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


## ----read_data--------------------------------------------------------------

sce <- readRDS("R_objects/Caron_clustered.PBMMCandETV6RUNX1.rds")

rownames(sce) <- uniquifyFeatureNames(rownames(sce), rowData(sce)$Symbol)

sce


## ----plot_TSNE_1-------------------------------------------------------------------------

plotReducedDim(sce, dimred = "TSNE_corrected", colour_by = "label")



## ----aggregate_cells---------------------------------------------------------------------

summed <- aggregateAcrossCells(sce, 
                               id=colData(sce)[,c("label", "SampleName")])
summed


## ----choose_cluster_and_subset-----------------------------------------------------------

current <- summed[,summed$label == "1"]

colData(current)



## ----make_edgeR_object-------------------------------------------------------------------

countsToUse <- counts(current)
colnames(countsToUse) <- colData(current)$SampleName

y <- DGEList(countsToUse, samples=colData(current))
y



## ----discard_low_cells-------------------------------------------------------------------
discarded <- current$ncells < 20
y <- y[,!discarded]
summary(discarded)


## ----remove_low_genes--------------------------------------------------------------------

keep <- filterByExpr(y, group=current$SampleGroup)
y <- y[keep,]
summary(keep)



## ----calc_norm_factors-------------------------------------------------------------------

y <- calcNormFactors(y)
y$samples



## ----plot_MD-----------------------------------------------------------------------------
par(mfrow=c(2,4))
for (i in seq_len(ncol(y))) {
    plotMD(y, column=i)
}


## ----plot_MDS----------------------------------------------------------------------------
y$samples$SampleGroup <- factor(y$samples$SampleGroup)
limma::plotMDS(cpm(y, log=TRUE), col = as.numeric(y$samples$SampleGroup))


## ----design_model------------------------------------------------------------------------
design <- model.matrix(~factor(SampleGroup), y$samples)
design


## ----estimate_dispersions----------------------------------------------------------------
y <- estimateDisp(y, design)


## ----plot_BCV----------------------------------------------------------------------------
plotBCV(y)


## ----fit_QL------------------------------------------------------------------------------
fit <- glmQLFit(y, design, robust=TRUE)


## ----plot_QL_dispersion------------------------------------------------------------------
plotQLDisp(fit)


## ----perform_DEA-------------------------------------------------------------------------
results <- glmQLFTest(fit, coef = 2)
summary(decideTests(results))


## ----get_results-------------------------------------------------------------------------
topTags(results)


## ----filter_summed_counts----------------------------------------------------------------
summed.filt <- summed[,summed$ncells >= 20]


## ----get_common_metadata-----------------------------------------------------------------

targets <- colData(sce)[!duplicated(sce$SampleName),] %>%
  data.frame()


## ----design_matrix_for_looping-----------------------------------------------------------

design <- model.matrix(~factor(SampleGroup), data=targets)
rownames(design) <- targets$SampleName


## ----run_DEA_loop------------------------------------------------------------------------
summed.filt$SampleGroup <- factor(summed.filt$SampleGroup)

de.results <- pseudoBulkDGE(summed.filt, 
    label = summed.filt$label,
    design = ~SampleGroup,
    coef = "SampleGroupPBMMC",
    condition = summed.filt$SampleName
)

de.results

de.results[[1]]


## ----per_label_DEGs----------------------------------------------------------------------
is.de <- decideTestsPerLabel(de.results, threshold=0.05)
summarizeTestsPerLabel(is.de)


## ----upregulated_genes-------------------------------------------------------------------

up.de <- is.de > 0 & !is.na(is.de)
head(sort(rowMeans(up.de), decreasing=TRUE), 10)


## ----downregulated_genes-----------------------------------------------------------------

down.de <- is.de < 0 & !is.na(is.de)
head(sort(rowMeans(down.de), decreasing=TRUE), 10)


## ----de_in_1_but_not_in_rest-------------------------------------------------------------
remotely.de <- decideTestsPerLabel(de.results, threshold=0.5)
not.de <- remotely.de==0 | is.na(remotely.de)

cx <- "1"

other.labels <- setdiff(colnames(not.de), cx)

unique.degs <- is.de[,cx]!=0 & rowMeans(not.de[,other.labels])==1
unique.degs <- names(which(unique.degs))
head(unique.degs)


## ----------------------------------------------------------------------------------------
top_gene <- rownames(de.results[[1]][order(de.results[[1]]$FDR),])[1]

sizeFactors(summed.filt) <- NULL

plotExpression(logNormCounts(summed.filt), 
    features=top_gene,
    x="SampleName", colour_by="SampleName", 
    other_fields="label") + 
    facet_wrap(~label) +
  ggtitle(top_gene)


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


















## ----exercise1_solution, eval=FALSE------------------------------------------------------
## # First load in the other two sample groups
 sce_PRET_HHD <- readRDS("R_objects/Caron_clustered.PRETandHHD.rds")
## 
## # replace the ensembl IDs with gene symbols where possible
 rownames(sce_PRET_HHD) <- uniquifyFeatureNames(rownames(sce_PRET_HHD), rowData(sce_PRET_HHD)$Symbol)
## 
## # check your sce object
 sce_PRET_HHD
## 
## # Part A
## # Re run the analysis using these new samples.
 summed_PRET_HHD <- aggregateAcrossCells(sce_PRET_HHD,
                                id=colData(sce_PRET_HHD)[,c("label", "SampleName")])
 summed_PRET_HHD
## 
 summed_PRET_HHD.filt <- summed_PRET_HHD[,summed_PRET_HHD$ncells >= 20]
## 
 targets_PRET_HHD <- colData(sce_PRET_HHD)[!duplicated(sce_PRET_HHD$SampleName),] %>%
   data.frame()
## 
 design_PRET_HHD <- model.matrix(~factor(SampleGroup), data=targets_PRET_HHD)
 rownames(design_PRET_HHD) <- targets_PRET_HHD$SampleName
## 
 summed_PRET_HHD.filt$SampleGroup <- factor(summed_PRET_HHD.filt$SampleGroup)
## 
 de.results_PRET_HHD <- pseudoBulkDGE(summed_PRET_HHD.filt,
                             label = summed_PRET_HHD.filt$label,
                             design = ~SampleGroup,
                             coef = "SampleGroupPRE-T",
                             condition = summed_PRET_HHD.filt$SampleName)
 
## # Part B
## # Looking at your results, which cluster has the most DEGs?
## 
 is.de_PRET_HHD <- decideTestsPerLabel(de.results_PRET_HHD, threshold=0.05)
 summarizeTestsPerLabel(is.de_PRET_HHD)
## 
## # Part C
## # Which genes are sig. DE in that cluster and not the others?
## 
 remotely.de_PRET_HHD <- decideTestsPerLabel(de.results_PRET_HHD, threshold=0.5)
 not.de_PRET_HHD <- remotely.de_PRET_HHD==0 | is.na(remotely.de_PRET_HHD)
## 
 cz <- "2"
## 
other.labels_PRET_HHD <- setdiff(colnames(not.de_PRET_HHD), cz)
## 
 unique.degs_PRET_HHD <- is.de_PRET_HHD[,cz]!=0 & rowMeans(not.de_PRET_HHD[,other.labels_PRET_HHD])==1
 unique.degs_PRET_HHD <- names(which(unique.degs_PRET_HHD))
 head(unique.degs_PRET_HHD)


## ----multi_cluster_DE--------------------------------------------------------------------
summed.sub <- summed[,summed$label %in% c("3", "4")]

between.res <- pseudoBulkDGE(summed.sub,
    label=rep("dummy", ncol(summed.sub)),
    design=~factor(SampleName) + factor(label),
    coef="factor(label)4")[[1]]

table(Sig=between.res$FDR <= 0.05, Sign=sign(between.res$logFC))



## ----multi_cluster_DE_results------------------------------------------------------------
between.res[order(between.res$FDR),]


## ----multi_cluster_DE_plots--------------------------------------------------------------
summed.sub <- logNormCounts(summed.sub, size.factors=NULL)
plotExpression(summed.sub, 
    features=head(rownames(between.res)[order(between.res$FDR)]),
    x="label", 
    colour_by=I(factor(summed.sub$SampleName)))


## Differential Abundance


## ----make_milo_object--------------------------------------------------------------------

milo <- Milo(sce)

milo


## ----build_Graph, message=FALSE----------------------------------------------------------

milo <- buildGraph(milo, k = 60, d = 30, reduced.dim = "corrected", BPPARAM = MulticoreParam(7))


## ----make_neighbourhoods, message=FALSE--------------------------------------------------

milo <- makeNhoods(milo, prop = 0.1, k = 60, d=30, refined = TRUE, reduced_dims = "corrected")


## ----plot_neighbourhood_sizes, message=FALSE---------------------------------------------

plotNhoodSizeHist(milo)


## ----count_cells, message=FALSE----------------------------------------------------------

milo <- countCells(milo, meta.data = data.frame(colData(milo)), sample="SampleName")

head(nhoodCounts(milo))


## ----milo_design, message=FALSE----------------------------------------------------------

milo_design <- data.frame(colData(milo))[,c("SampleName", "SampleGroup")]
milo_design <- distinct(milo_design)
rownames(milo_design) <- milo_design$SampleName

milo_design


## ----calcultate_neighbourhood_distance, message=FALSE------------------------------------

# milo <- calcNhoodDistance(milo, d=30, reduced.dim = "corrected")
milo <- readRDS("R_objects/10_milo_calcNhoodDistance.rds")


## ----test_neighbourhoods, message=FALSE--------------------------------------------------

da_results <- testNhoods(milo, design = ~ SampleGroup, design.df = milo_design, reduced.dim = "corrected")

da_results %>%
  arrange(SpatialFDR) %>%
  head()


## ----plot_DA_pvalue_histogram, message=FALSE---------------------------------------------

ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)


## ----plot_DA_volcano, message=FALSE------------------------------------------------------
ggplot(da_results, aes(logFC, -log10(SpatialFDR))) + 
  geom_point() +
  geom_hline(yintercept = 1) 


## ----build_neighbourhood_graph, message=FALSE--------------------------------------------
milo <- buildNhoodGraph(milo)



## ----plot_umap_neighbourhood_combined, message=FALSE-------------------------------------

umap_plot <- plotReducedDim(milo, dimred = "UMAP_corrected", colour_by="label", text_by = "label")

nh_graph_plot <- plotNhoodGraphDA(milo, da_results, layout="UMAP_corrected",alpha=0.05)

umap_plot + nh_graph_plot +
  plot_layout(guides="collect")


## ----annotate_neighbourhoods, message=FALSE----------------------------------------------

da_results <- annotateNhoods(milo, da_results, coldata_col = "label")
head(da_results)


## ----plot_label_fraction, message=FALSE--------------------------------------------------
ggplot(da_results, aes(label_fraction)) + geom_histogram(bins=50)


## ----excluded_mixed_neighbourhoods, message=FALSE----------------------------------------

da_results$label <- ifelse(da_results$label_fraction < 0.7, "Mixed", da_results$label)

head(da_results)


## ----beeswarm_plot_by_label, message=FALSE-----------------------------------------------

plotDAbeeswarm(da_results, group.by = "label")






## ----autogroup, message=FALSE------------------------------------------------------------

da_results <- groupNhoods(milo, da_results, max.lfc.delta = 10, overlap = 1)
head(da_results)

length(unique(da_results$NhoodGroup))

## ----grouped_umap, message=FALSE---------------------------------------------------------
plotNhoodGroups(milo, da_results, layout="UMAP_corrected") 


## ----grouped_beeswarm, message=FALSE-----------------------------------------------------
plotDAbeeswarm(da_results, "NhoodGroup")



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














## ----exercise2_solution------------------------------------------------------
## # First load in the other two sample groups
 set.seed(42) # set your random seed for reproducibility
## 
## # Part A
## # rerun the grouping with values for max.lfc.delta of 1, 5, 25 and plot beeswarm for each
## 
 plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 1) , 
                group.by = "NhoodGroup") + ggtitle("max LFC delta=1")
## 
 plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 5) , 
                group.by = "NhoodGroup") + ggtitle("max LFC delta=5")
## 
 plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 25) , 
                group.by = "NhoodGroup") + ggtitle("max LFC delta=25")
## 
## # Part B
## # Using the max.lfc.delta you preferred from Part A now alter the overlap value trying 1, 3, 10.
## 
 plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 5, overlap = 1) , 
                group.by = "NhoodGroup") + ggtitle("overlap=1")
## 
 plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 5, overlap = 3) , 
                group.by = "NhoodGroup") + ggtitle("overlap=3")
## 
 plotDAbeeswarm(groupNhoods(milo, da_results, max.lfc.delta = 5, overlap = 10) , 
                group.by = "NhoodGroup") + ggtitle("overlap=10")
## 
## # The final code
## 
 set.seed(42)
 da_results <- groupNhoods(milo, da_results, max.lfc.delta = 5, overlap=1)
 plotNhoodGroups(milo, da_results, layout="UMAP_corrected")
 plotDAbeeswarm(da_results, group.by = "NhoodGroup")
## 


## ----get_hvgs, message=FALSE-------------------------------------------------------------
set.seed(101)
dec <- modelGeneVar(milo)
hvgs <- getTopHVGs(dec, n=2000)
head(hvgs)


## ----get_neighbourhood_group_markers, message=FALSE--------------------------------------
set.seed(42)
nhood_markers <- findNhoodGroupMarkers(milo, da_results, subset.row = hvgs, 
                                       aggregate.samples = TRUE, sample_col = "SampleName")
head(nhood_markers)


## ----group_9_markers, message=FALSE------------------------------------------------------
gr9_markers <- nhood_markers[c("logFC_9", "adj.P.Val_9","GeneID")] 
colnames(gr9_markers) <- c("logFC", "adj.P.Val", "GeneID")
head(gr9_markers[order(gr9_markers$adj.P.Val), ])


## ----plot_group_markers, message=FALSE---------------------------------------------------
markers <- nhood_markers$GeneID[nhood_markers$adj.P.Val_9 < 0.01 
                                & nhood_markers$logFC_9 > 0]

plotNhoodExpressionGroups(milo, da_results, features=intersect(rownames(milo), markers[1:10]),
                          subset.nhoods = da_results$NhoodGroup %in% c('6','9'), 
                          scale=TRUE,
                          grid.space = "fixed")



