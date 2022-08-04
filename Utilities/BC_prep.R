library(scater)
library(scran)
library(batchelor)
library(bluster)
library(tidyverse)
library(pheatmap)
library(clustree)
library(Cairo)
library(BiocSingular)
library(cowplot)

##### single samples

sce.sample.1 <- readRDS("../Robjects/BC_sample1.rds")
sce.sample.2 <- readRDS("../Robjects/BC_sample2.rds")

sce.sample.1 <- logNormCounts(sce.sample.1)
sce.sample.2 <- logNormCounts(sce.sample.2)

dec.sample.1 <- modelGeneVar(sce.sample.1)
dec.sample.2 <- modelGeneVar(sce.sample.2)

hvgs.sample.1 <- getTopHVGs(dec.sample.1, prop=0.1)
hvgs.sample.2 <- getTopHVGs(dec.sample.2, prop=0.1)

sce.sample.1 <- runPCA(sce.sample.1, subset_row = hvgs.sample.1)
sce.sample.2 <- runPCA(sce.sample.2, subset_row = hvgs.sample.2)

sce.sample.1 <- runTSNE(sce.sample.1, dimred = "PCA")
sce.sample.2 <- runTSNE(sce.sample.2, dimred = "PCA")

s1.snn.gr <- buildSNNGraph(sce.sample.1, use.dimred="PCA")
s1.clusters <- igraph::cluster_walktrap(s1.snn.gr)$membership
colLabels(sce.sample.1) <- factor(s1.clusters)
s2.snn.gr <- buildSNNGraph(sce.sample.2, use.dimred="PCA")
s2.clusters <- igraph::cluster_walktrap(s2.snn.gr)$membership
colLabels(sce.sample.2) <- factor(s2.clusters)

saveRDS(sce.sample.1, "../Robjects/BC_sample1_dimred.rds")
saveRDS(sce.sample.2, "../Robjects/BC_sample2_dimred.rds")

saveRDS(dec.sample.1, "../Robjects/BC_dec1_dimred.rds")
saveRDS(dec.sample.2, "../Robjects/BC_dec2_dimred.rds")

###### multiple samples

all.sce <- readRDS("../Robjects/DataIntegration_all_sce.Rds")

list2env(all.sce, globalenv())
all.sce <- cbind(`ETV6-RUNX1_1`,`ETV6-RUNX1_2`,`ETV6-RUNX1_3`,`ETV6-RUNX1_4`,
                 PBMMC_1b,PBMMC_2,PBMMC_3)

all.sce <- logNormCounts(all.sce)
all.dec <- modelGeneVar(all.sce)
all.hvgs <- getTopHVGs(all.dec, prop = 0.1)
all.sce <- runPCA(all.sce, subset_row = all.hvgs)
all.sce <- runTSNE(all.sce, dimred = "PCA")

saveRDS(all.sce, "../Robjects/DataIntegration_all_sce_dimred.Rds")

plotTSNE(all.sce, colour_by="SampleName")

#####

all.sce <- lapply(all.sce, logNormCounts)

all.dec <- lapply(all.sce, modelGeneVar)
all.hvgs <- lapply(all.dec, getTopHVGs, prop=0.1)

set.seed(10000)
all.sce <- mapply(FUN=runPCA,
                  x=all.sce,
                  subset_row=all.hvgs,
                  MoreArgs=list(ncomponents=25,
                                BSPARAM=RandomParam()),
                  SIMPLIFY=FALSE)

set.seed(100000)
all.sce <- lapply(all.sce, runTSNE, dimred="PCA")

saveRDS(all.sce,"../Robjects/DataIntegration_all_sce_dimred.Rds")

all.dec <- setNames(all.dec, c("ETV6-RUNX1_1_dec","ETV6-RUNX1_2_dec","ETV6-RUNX1_3_dec","ETV6-RUNX1_4_dec","PBMMC_1b_dec","PBMMC_2_dec","PBMMC_3_dec"))

saveRDS(all.dec,"../Robjects/DataIntegration_all_dec_dimred.Rds")

##### other idea

all.sce <- readRDS("../Robjects/postQC_caron_allcells.rds")

a <- colData(all.sce) %>%
  as.data.frame() %>%
  separate(Sample, sep = "_", into = c("Dataset","no"), remove = FALSE) %>%
  select(-no) %>%
  DataFrame()
colData(all.sce) <- a

Caron <- all.sce[,all.sce$Dataset == "ETV6-RUNX1" | all.sce$Dataset == "PBMMC"]

Caron <- logNormCounts(Caron)
Caron.dec <- modelGeneVar(Caron)
Caron.hvgs <- getTopHVGs(Caron.dec, prop = 0.1)
Caron <- runPCA(Caron, subset_row = Caron.hvgs)
Caron <- runTSNE(Caron, dimred = "PCA")

#saveRDS(all.sce, "../Robjects/DataIntegration_all_sce_dimred.Rds")
saveRDS(Caron, "../Robjects/07_semiprocessedCaronSamples.rds")

plotTSNE(Caron, colour_by="Sample")
