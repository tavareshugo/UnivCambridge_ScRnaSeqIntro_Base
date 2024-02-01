# load packages 
library(scater)
library(scran)
library(bluster)
library(cluster)
library(igraph)
library(pheatmap)
library(patchwork)
library(tidyverse)

setwd("~/Course_Materials")

# load data 
sce <- readRDS("R_objects/Caron_batch_corrected.500.rds")

# check samples 
table(sce$SampleName)

###############################
# Clustering with the walktrap method
#####################################

# cluster cells using default clustering parameters and the batch-corrected reduced dimensions
# The default algorithm for clusterCells constructs a SNN graph with k = 10 and 
# edge weighting by rank. It then uses walktrap to identify communities. 
clustering1 <- clusterCells(sce, use.dimred="corrected", full=TRUE)

# table of clusters 
table(clustering1$clusters)

## NOTE: we can plot the kNN graph, but we wont go through that here, please see
## the main materials if you are interested to know how. Generally, instead, we
## use a tSNE or UMAP to visualise the clusters.

# plot tSNE with default clusters 
sce$Clusters1 <- clustering1$clusters
plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               colour_by="Clusters1",
               text_by = "Clusters1")


# Generate an alternative clustering using Shared Nearest Neighbours with k = 15 
# to generate the graph and the walktrap algorithm to identify communities
sce$walktrap15 <- clusterCells(sce, 
                           use.dimred = "corrected", 
                           BLUSPARAM = SNNGraphParam(k = 15, 
                                                     cluster.fun = "walktrap"))


# A heatmap showing the (logged) numbers of cells in each cluster for each dataset with 
# walktrap k = 15 
# This gives us an overview of how well each cluster is represented across the samples and the replicates
w15_table <- log(table(sce$walktrap15, sce$SampleName)+1)
pheatmap(w15_table, cluster_rows = TRUE, cluster_cols = FALSE)

# plot tSNE walktrap k15 
plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               colour_by="walktrap15", 
               text_by = "walktrap15")


# plot tSNE walktrap k15 by samplegroup 
plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               colour_by="walktrap15", 
               text_by = "walktrap15",
               other_fields = list("SampleGroup")) +
  facet_wrap("SampleGroup")

###############################
# Clustering with the Louvain method
#####################################

# cluster cells - louvain k15 
sce$louvain15 <- clusterCells(sce, 
                           use.dimred = "corrected", 
                           BLUSPARAM = SNNGraphParam(k = 15, 
                                                     cluster.fun = "louvain"))

# plot tSNE louvain k15 
plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               colour_by = "louvain15", 
               text_by = "louvain15")


# plot tSNE louvain k15 by samplegroup 
plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               colour_by="louvain15", 
               text_by = "louvain15",
               other_fields = list("SampleGroup")) +
  facet_wrap("SampleGroup")


## -- Exercise 1 -- ############################################################


################################################################################

####################################
# Assessing cluster behaviour
###################################

# calculate silhouette widths for the Leiden, k=20 clustering
sil.approx <- approxSilhouette(reducedDim(sce, "corrected"),
                               clusters=sce$leiden20)

# Examine the silouette data
sil.approx

# silhouette width plot with beeswarm 
plotSilBeeswarm <- function(silDat){
  silTab <- silDat %>% 
    as.data.frame() %>% 
    mutate(closestCluster = ifelse(width > 0, cluster, other) %>% factor())
  
  plt <- silTab %>% 
      ggplot(aes(x=cluster, y=width, colour=closestCluster)) +
        ggbeeswarm::geom_quasirandom(method="smiley", alpha=0.6) +
        theme_bw()
  
  plt <- plt + scale_color_manual(
    values = scater:::.get_palette("tableau20"),
    name = "closestCluster")
  
  plt
}

p1 <- plotSilBeeswarm(sil.approx)
p2 <- plotReducedDim(sce, 
                     dimred = "TSNE_corrected", 
                     colour_by="leiden20", 
                     text_by = "leiden20")
p1 + p2


# silhouette width grid 
plotSilGrid <- function(silDat){
  silDat %>% 
    as.data.frame() %>% 
    mutate(closestCluster = ifelse(width > 0, cluster, other) %>% factor()) %>% 
    count(cluster, closestCluster,  name="olap") %>% 
    group_by(cluster) %>% 
    mutate(total  = sum(olap)) %>% 
    mutate(proportion = olap / total) %>% 
    mutate(proportion = ifelse(cluster == closestCluster, proportion, -proportion)) %>% 
    ggplot(aes(x = cluster, y = closestCluster)) +
      geom_tile(aes(fill = proportion)) +
      geom_text(aes(label = olap), size=5) +
      scale_fill_gradientn(colors = c("#fc8d59", "#ffffbf", "#91cf60"),
                            limits = c(-1, 1)) +
      geom_vline(xintercept=seq(0.5, 30.5, by=1)) +
      geom_hline(yintercept=seq(0.5, 30.5, by=1), colour="lightgrey", linetype=2) +
      guides(fill = "none") +
      theme(
          aspect.ratio = 1,
          panel.background = element_blank())
}

# Save Leiden 20 sil grid plot
le20_3 <- plotSilGrid(sil.approx)

# plot silhouette width for walktrap k=15 clustering 
sil.approx <- approxSilhouette(reducedDim(sce, "corrected"),
                               clusters=sce$walktrap15)

wp1 <- plotSilBeeswarm(sil.approx)

wp2 <- plotReducedDim(sce, 
                     dimred = "TSNE_corrected", 
                     colour_by="walktrap15", 
                     text_by = "walktrap15")

wp3 <- plotSilGrid(sil.approx)

wp1 + wp2 + wp3


# plot silhouette width for louvain k=15 clustering 
sil.approx <- approxSilhouette(reducedDim(sce, "corrected"),
                               clusters=sce$louvain15)

lp1 <- plotSilBeeswarm(sil.approx)

lp2 <- plotReducedDim(sce, 
                     dimred = "TSNE_corrected", 
                     colour_by="louvain15", 
                     text_by = "louvain15")

lp3 <- plotSilGrid(sil.approx)

lp1 + lp2 + lp3

le20_3 + wp3+ lp3

# Look at the concordance of different clusterings using the jaccard index 
jacc.mat <- linkClustersMatrix(sce$louvain15, sce$walktrap15)
rownames(jacc.mat) <- paste("Louvain", rownames(jacc.mat))
colnames(jacc.mat) <- paste("Walktrap", colnames(jacc.mat))
pheatmap(jacc.mat, color=viridis::viridis(100), cluster_cols=FALSE, cluster_rows=FALSE)


# Use clustersweep to examine a range clustering parameters
# clusterSweepwith k = 5, 10, 15, 20, 25 and walktrap method
out <- clusterSweep(reducedDim(sce, "corrected"),
                    BLUSPARAM = NNGraphParam(),
                    k = as.integer(c(5, 10, 15, 20, 25)),
                    cluster.fun = "walktrap",
                    BPPARAM=BiocParallel::MulticoreParam(7))

# make a data frame of clustersweep results
df <- as.data.frame(out$parameters)

# get the number of clusters
df$num.clusters <- apply(out$clusters, 2, max)

# get the mean silhouette width
getMeanSil <- function(cluster) {
    sil <- approxSilhouette(reducedDim(sce, "corrected"), cluster)
    mean(sil$width)
}
df$silhouette <- map_dbl(as.list(out$clusters), getMeanSil)

# Generate line plots of the cluster number and mean silouette width
# for different values of k
nclPlot <- ggplot(df, aes(x = k, y = num.clusters)) + 
                  geom_line(lwd=2)
silPlot <- ggplot(df, aes(x = k, y = silhouette)) + 
                  geom_line(lwd=2)
nclPlot + silPlot


## -- Exercise 2 -- ############################################################


################################################################################

# add_clusterSweep output to sce 
colData(sce) <- cbind(colData(sce), DataFrame(out$clusters))

# set labels to our favourite clustering
colLabels(sce) <- sce$k.25_cluster.fun.leiden

# plot tSNE leiden k25 
plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               colour_by = "label", 
               text_by = "label") +
  ggtitle("Leiden k=25 clusters")

# switch the rownames in the SCE object to be gene symbols
# (and make sure they are unique!) 
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)

# plot B cell marker expression on TSNE 
plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               by_exprs_values = "logcounts",
               colour_by = "CD79A",
               text_by = "label", other_fields = list("SampleGroup"))


# plot expression B cell marker as a violin plot 
plotExpression(sce, 
               exprs_values = "logcounts",
               x = "label", 
               colour_by = "label",
               features=c("CD79A"))

# plot tSNE monocyte markers 
plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               by_exprs_values = "logcounts",
               colour_by = "LYZ",
               text_by = "label")

# plot expression monocyte markers 
plotExpression(sce, 
               exprs_values = "logcounts",
               x = "label", 
               colour_by = "label",
               features=c("LYZ"))

