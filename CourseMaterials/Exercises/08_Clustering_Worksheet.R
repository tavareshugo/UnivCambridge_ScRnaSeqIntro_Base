# Load packages
library(scater) # scRnaSeq QC
library(scran) # scRnaSeq normalisation
library(bluster) # scRnaSeq clustering
library(cluster) # for silhouette
library(igraph) # for graph-based clustering and plotting networks
library(pheatmap) # for heatmap plotting
library(patchwork) # to combine plots
library(tidyverse) # data wrangling and plotting (ggplot2)

# Load data

sce <- readRDS("Robjects/DataIntegration_mnn.Rds")

#-------------------------------------------------------------------------------

# Run graph based clustering - default parameters

clustering1 <- clusterCells(sce, use.dimred="corrected", full=TRUE)

# How many clusters and cells per cluster?

table(clustering1$clusters)

#-------------------------------------------------------------------------------

# extract the graph
snn.gr <- clustering1$objects$graph

# Add Sample group to vertices (nodes, ie cells)
V(snn.gr)$SampleGroup <- as.character(colData(sce)$SampleGroup)

# pick 1000 nodes randomly
set.seed(1423)
selectedNodes <- sample(3500, 1000)

# subset graph for these 1000 randomly chosen nodes
snn.gr.subset <- subgraph(snn.gr, selectedNodes)

# set colors for clusters
grps <-  V(snn.gr.subset)$SampleGroup
cols <- c("dodgerblue", "lightyellow")[as.numeric(factor(grps))]
names(cols) <- grps

# plot graph
plot.igraph(snn.gr.subset,
            layout = layout_with_fr(snn.gr.subset),
            vertex.size = 3, 
            vertex.label = NA,
            vertex.color = cols,
            frame.color = cols,
            main = "default parameters"
)

# add legend
legend('bottomright',
       legend=unique(names(cols)),
       pch=21,
       pt.bg=unique(cols),
       pt.cex=1, cex=.6, bty="n", ncol=1)
#-------------------------------------------------------------------------------

# Visualise the clusters on a tSNE plot
sce$Clusters1 <- clustering1$clusters
plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               colour_by="Clusters1",
               text_by = "Clusters1")

#-------------------------------------------------------------------------------

# The Walktrap method with k = 15

sce$walktrap15 <- clusterCells(sce, 
                               use.dimred = "corrected", 
                               BLUSPARAM = SNNGraphParam(k = 15, 
                                                      cluster.fun = "walktrap"))

table(sce$walktrap15)

## Visualise the number of cells from each sample in each cluster 
table(sce$walktrap15, sce$SampleName) %>% 
  pheatmap(cluster_rows = FALSE, cluster_cols = FALSE)

plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               colour_by="walktrap15", 
               text_by = "walktrap15",
               other_fields = list("SampleGroup")) +
  facet_wrap(vars(SampleGroup))

#-------------------------------------------------------------------------------

# The Louvain method

sce$louvain15 <- clusterCells(sce, 
                           use.dimred = "corrected", 
                           BLUSPARAM = SNNGraphParam(k = 15, 
                                                     cluster.fun = "louvain"))

table(sce$louvain15)

plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               colour_by="louvain15", 
               text_by = "louvain15",
               other_fields = list("SampleGroup")) +
  facet_wrap("SampleGroup")

#-------------------------------------------------------------------------------

# The Leiden method

## EXERCISE 1 ##

#-------------------------------------------------------------------------------

# Assessing cluster behaviour

## Silhouette width

sil.approx <- approxSilhouette(reducedDim(sce, "corrected"), 
                               clusters=sce$leiden20)
sil.approx

### Visualise the results in as a beeswarm plot

plotSilBeeswarm <- function(silDat){
  silTab <- silDat %>% 
    as.data.frame() %>% 
    mutate(closestCluster = ifelse(width > 0, cluster, other) %>% factor())
  
  plt <- silTab %>% 
      ggplot(aes(x=cluster, y=width, colour=closestCluster)) +
        ggbeeswarm::geom_quasirandom(method="smiley", alpha=0.6) +
        theme_bw()
  
  plt <- scater:::.resolve_plot_colours(plt, silTab$closestCluster, "closestCluster")
}

p1 <- plotSilBeeswarm(sil.approx)
p2 <- plotReducedDim(sce, 
                     dimred = "TSNE_corrected", 
                     colour_by="leiden20", 
                     text_by = "leiden20")
p1 + p2


### Visualise results on a grid 

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

### Again with Walktrap + k = 15

sil.approx <- approxSilhouette(reducedDim(sce, "corrected"),
                               clusters=sce$walktrap15)

wp1 <- plotSilBeeswarm(sil.approx)

wp2 <- plotReducedDim(sce, 
                     dimred = "TSNE_corrected", 
                     colour_by="walktrap15", 
                     text_by = "walktrap15")

wp3 <- plotSilGrid(sil.approx)

wp1 + wp2 + wp3

## Again with Louvain + k = 15

sil.approx <- approxSilhouette(reducedDim(sce, "corrected"),
                               clusters=sce$louvain15)

lp1 <- plotSilBeeswarm(sil.approx)

lp2 <- plotReducedDim(sce, 
                     dimred = "TSNE_corrected", 
                     colour_by="louvain15", 
                     text_by = "louvain15")

lp3 <- plotSilGrid(sil.approx)

lp1 + lp2 + lp3

#-------------------------------------------------------------------------------

##  Modularity to assess clusters quality

walktrap15 <- clusterCells(sce, 
                           use.dimred = "corrected", 
                           BLUSPARAM = SNNGraphParam(k = 15, 
                                                     cluster.fun = "walktrap"),
                           full = TRUE)
g <- walktrap15$objects$graph
ratio <- pairwiseModularity(g, walktrap15$clusters, as.ratio=TRUE)

### Visualise modularity on heatmap

hm1 <- pheatmap(log2(ratio+1), 
                cluster_rows=FALSE, 
                cluster_cols=FALSE,
                color=colorRampPalette(c("white", "blue"))(100))

### Compare to the equivalent plot based on silhouette widths

wp4 <- ggplotify::as.ggplot(hm1)

wp2 + wp3 + wp4

### Visualise as a network graph

cluster.gr <- igraph::graph_from_adjacency_matrix(log2(ratio+1),
                                                  mode="upper", 
                                                  weighted=TRUE, diag=FALSE)

set.seed(11001010)
plot(cluster.gr, 
     edge.width=igraph::E(cluster.gr)$weight*5,
     layout=igraph::layout_with_lgl)

## Comparing two sets of clusters

jacc.mat <- linkClustersMatrix(sce$louvain15, sce$walktrap15)
rownames(jacc.mat) <- paste("Louvain", rownames(jacc.mat))
colnames(jacc.mat) <- paste("Walktrap", colnames(jacc.mat))
pheatmap(jacc.mat, 
         color=viridis::viridis(100), 
         cluster_cols=FALSE, 
         cluster_rows=FALSE)

#-------------------------------------------------------------------------------

# Cluster sweep

out <- clusterSweep(reducedDim(sce, "corrected"),
                    BLUSPARAM = NNGraphParam(),
                    k = as.integer(c(5, 10, 15, 20, 25)),
                    cluster.fun = "walktrap",
                    BPPARAM=BiocParallel::MulticoreParam(7))


out$clusters[,1:4]
out$parameters

## Assess all clusters using number of clusters and mean silhouette width

df <- as.data.frame(out$parameters)

### get the number of clusters

df$num.clusters <- apply(out$clusters, 2, max)

### get the mean silhouette width

getMeanSil <- function(cluster) {
    sil <- approxSilhouette(reducedDim(sce, "corrected"), cluster)
    mean(sil$width)
}
df$silhouette <- map_dbl(as.list(out$clusters), getMeanSil)

### plot scores

nclPlot <- ggplot(df, aes(x = k, y = num.clusters)) + 
                  geom_line(lwd=2)
silPlot <- ggplot(df, aes(x = k, y = silhouette)) + 
                  geom_line(lwd=2)
nclPlot + silPlot

## Compare clusterings using Jaccard index

jacc.mat <- linkClustersMatrix(out$clusters$k.10_cluster.fun.walktrap, 
                               out$clusters$k.15_cluster.fun.walktrap)
rownames(jacc.mat) <- paste("Walktrap_10", rownames(jacc.mat))
colnames(jacc.mat) <- paste("Walktrap_15", colnames(jacc.mat))
pheatmap(jacc.mat, color=viridis::viridis(100), cluster_cols=FALSE, cluster_rows=FALSE)


jacc.mat <- linkClustersMatrix(out$clusters$k.20_cluster.fun.walktrap, 
                               out$clusters$k.25_cluster.fun.walktrap)
rownames(jacc.mat) <- paste("Walktrap_20", rownames(jacc.mat))
colnames(jacc.mat) <- paste("Walktrap_25", colnames(jacc.mat))
pheatmap(jacc.mat, color=viridis::viridis(100), cluster_cols=FALSE, cluster_rows=FALSE)

#-------------------------------------------------------------------------------

## EXERCISE 2 ##


#-------------------------------------------------------------------------------

# Add `clusterSweep` results to SCE object.

colData(sce) <- cbind(colData(sce), DataFrame(out$clusters))

#-------------------------------------------------------------------------------

# Finalise clustering selection

colLabels(sce) <- sce$k.20_cluster.fun.walktrap

#-------------------------------------------------------------------------------

# Expression of known marker genes

plotReducedDim(sce, 
               dimred = "TSNE_corrected",
               colour_by = "label", 
               text_by = "label") +
  ggtitle("Walktrap clusters")

## We will switch the rownames to be gene symbols

rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)


## B-cells markers


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



plotExpression(sce, 
               exprs_values = "logcounts",
               x = "label", 
               colour_by = "label",
               features=c("MS4A1", "CD79A"))


## Monocyte markers


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



plotExpression(sce, 
               exprs_values = "logcounts",
               x = "label", 
               colour_by = "label",
               features=c("FCGR3A", "MS4A7"))


## Dendritic cell markers


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



plotExpression(sce, 
               exprs_values = "logcounts",
               x = "label", 
               colour_by = "label",
               features=c("FCER1A", "CST3"))

