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

# cluster cells default 
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


# cluster cells - walktrap k15 
sce$walktrap15 <- clusterCells(sce, 
                           use.dimred = "corrected", 
                           BLUSPARAM = SNNGraphParam(k = 15, 
                                                     cluster.fun = "walktrap"))


# heatmap walktrap k15 
table(sce$walktrap15, sce$SampleName) %>% 
  pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)


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


# calculate silhouette widths 
sil.approx <- approxSilhouette(reducedDim(sce, "corrected"),
                               clusters=sce$leiden20)


# silhouette width beeswarm 
plotSilBeeswarm <- function(silDat){
  silTab <- silDat %>% 
    as.data.frame() %>% 
    mutate(closestCluster = ifelse(width > 0, cluster, other) %>% factor())
  
  plt <- silTab %>% 
      ggplot(aes(x=cluster, y=width, colour=closestCluster)) +
        ggbeeswarm::geom_quasirandom(method="smiley", alpha=0.6) +
        theme_bw()
  
  plt <- scater:::.resolve_plot_colours(plt, silTab$closestCluster, "closestCluster")
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
plotSilGrid(sil.approx)


# plot silhouette width walktrap k15 
sil.approx <- approxSilhouette(reducedDim(sce, "corrected"),
                               clusters=sce$walktrap15)

wp1 <- plotSilBeeswarm(sil.approx)

wp2 <- plotReducedDim(sce, 
                     dimred = "TSNE_corrected", 
                     colour_by="walktrap15", 
                     text_by = "walktrap15")

wp3 <- plotSilGrid(sil.approx)

wp1 + wp2 + wp3


# plot silhouette width louvain k15 
sil.approx <- approxSilhouette(reducedDim(sce, "corrected"),
                               clusters=sce$louvain15)

lp1 <- plotSilBeeswarm(sil.approx)

lp2 <- plotReducedDim(sce, 
                     dimred = "TSNE_corrected", 
                     colour_by="louvain15", 
                     text_by = "louvain15")

lp3 <- plotSilGrid(sil.approx)

lp1 + lp2 + lp3


# pairwise modularity 
walktrap15 <- clusterCells(sce, 
                           use.dimred = "corrected", 
                           BLUSPARAM = SNNGraphParam(k = 15, 
                                                     cluster.fun = "walktrap"),
                           full = TRUE)
g <- walktrap15$objects$graph
ratio <- pairwiseModularity(g, walktrap15$clusters, as.ratio=TRUE)

hm1 <- pheatmap(log2(ratio+1),
         cluster_rows=FALSE, 
         cluster_cols=FALSE,
         color=colorRampPalette(c("white", "blue"))(100))


# plot modularity v silhouette widths 
wp4 <- ggplotify::as.ggplot(hm1)
wp2 + wp3 + wp4


# plot modularity network 
cluster.gr <- igraph::graph_from_adjacency_matrix(log2(ratio+1),
                                                  mode="upper", 
                                                  weighted=TRUE, diag=FALSE)

set.seed(11001010)
plot(cluster.gr, 
     edge.width=igraph::E(cluster.gr)$weight*5,
     layout=igraph::layout_with_lgl)


# jaccard index 
jacc.mat <- linkClustersMatrix(sce$louvain15, sce$walktrap15)
rownames(jacc.mat) <- paste("Louvain", rownames(jacc.mat))
colnames(jacc.mat) <- paste("Walktrap", colnames(jacc.mat))
pheatmap(jacc.mat, color=viridis::viridis(100), cluster_cols=FALSE, cluster_rows=FALSE)


# clusterSweep 
out <- clusterSweep(reducedDim(sce, "corrected"),
                    BLUSPARAM = NNGraphParam(),
                    k = as.integer(c(5, 10, 15, 20, 25)),
                    cluster.fun = "walktrap",
                    BPPARAM=BiocParallel::MulticoreParam(7))

# assess clustersweep 
df <- as.data.frame(out$parameters)

# get the number of clusters
df$num.clusters <- apply(out$clusters, 2, max)

# get the mean silhouette width
getMeanSil <- function(cluster) {
    sil <- approxSilhouette(reducedDim(sce, "corrected"), cluster)
    mean(sil$width)
}
df$silhouette <- map_dbl(as.list(out$clusters), getMeanSil)

nclPlot <- ggplot(df, aes(x = k, y = num.clusters)) + 
                  geom_line(lwd=2)
silPlot <- ggplot(df, aes(x = k, y = silhouette)) + 
                  geom_line(lwd=2)
nclPlot + silPlot

# clusterSweep Jaccard index 
jacc.mat <- linkClustersMatrix(out$clusters$k.15_cluster.fun.walktrap, 
                               out$clusters$k.25_cluster.fun.walktrap)
rownames(jacc.mat) <- paste("Walktrap_15", rownames(jacc.mat))
colnames(jacc.mat) <- paste("Walktrap_25", colnames(jacc.mat))
pheatmap(jacc.mat, 
         color = viridis::viridis(100), 
         cluster_cols = FALSE, 
         cluster_rows = FALSE)

## -- Exercise 2 -- ############################################################


################################################################################

# add_clusterSweep output to sce 
colData(sce) <- cbind(colData(sce), DataFrame(out$clusters))

# set labels 
colLabels(sce) <- sce$k.25_cluster.fun.leiden

# plot tSNE leiden k25 
plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               colour_by = "label", 
               text_by = "label") +
  ggtitle("Leiden k=25 clusters")


# symbols to rownames 
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)


# plot tSNE B cell markers 
p1 <- plotReducedDim(sce, 
               dimred = "TSNE_corrected", 
               by_exprs_values = "logcounts",
               colour_by = "MS4A1",
               text_by = "label")
p2 <- plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               by_exprs_values = "logcounts",
               colour_by = "CD79A",
               text_by = "label")
p1 + p2


# plot expression B cell markers 
plotExpression(sce, 
               exprs_values = "logcounts",
               x = "label", 
               colour_by = "label",
               features=c("MS4A1", "CD79A"))


# plot tSNE monocyte markers 
p1 <- plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               by_exprs_values = "logcounts",
               colour_by = "FCGR3A",
               text_by = "label")
p2 <- plotReducedDim(sce, 
               dimred = "TSNE_corrected", 
               by_exprs_values = "logcounts",
               colour_by = "MS4A7",
               text_by = "label")
p1 + p2


# plot expression monocyte markers 
plotExpression(sce, 
               exprs_values = "logcounts",
               x = "label", 
               colour_by = "label",
               features=c("FCGR3A", "MS4A7"))


# plot tSNE dendritic cell markers 
p1 <- plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               by_exprs_values = "logcounts",
               colour_by = "FCER1A",
               text_by = "label")
p2 <- plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               by_exprs_values = "logcounts",
               colour_by = "CST3",
               text_by = "label")
p1 + p2


# plot expression dendritic cell markers 
plotExpression(sce, 
               exprs_values = "logcounts",
               x = "label", 
               colour_by = "label",
               features=c("FCER1A", "CST3"))

