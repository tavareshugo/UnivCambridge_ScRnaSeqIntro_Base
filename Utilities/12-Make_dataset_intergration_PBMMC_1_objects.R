#!/usr/bin/env Rscript
#SBATCH -J ObjectForDI
#SBATCH -o ObjectForDI.%j.out
#SBATCH -e ObjectForDI.%j.err
#SBATCH --mincpus 8 
#SBATCH --mem=16G
#SBATCH --time=02:57:42

# This script loads the two technical replicate PBMMC samples - SRR9264351 &
# SRR9264352. It then filteres the cells, normalizes the counts, runs dimension
# reduction and clustering (default of clusterCells) on each of these samples
# separately.

library(DropletUtils)
library(scater) 
library(scran)
library(BiocParallel)
library(rtracklayer)
library(tidyverse)

bp.params <- MulticoreParam(workers = 8)

prepSample <- function(sampleID){
    # load data
    cr_dir <- str_c("data/cellranger/", 
                    sampleID, 
                    "/outs/filtered_feature_bc_matrix")
    names(cr_dir) <- sampleID
    sce <- read10xCounts(cr_dir, col.names=TRUE, BPPARAM = bp.params)

    # add meta data 
    sce$Barcode <- rownames(colData(sce))
    sce$SampleName <- ifelse(sampleID == "SRR9264351", "PBMMC_1a", "PBMMC_1b") 
    sce$SampleGroup <- "PBMMC"

    # Keep only detected genes
    detected_genes <- rowSums(counts(sce)) > 0
    sce <- sce[detected_genes,]

    # Annotate genes

    # or just get the info from the GTF
    gtf <- "data/references/refdata-gex-GRCh38.p13-Gencode.v41/genes/genes.gtf.gz"
    annot <- readGFF(gtf, 
                    columns="seqid", 
                    tags="gene_id", 
                    filter=list(type="gene")) %>%
        rename(ID = gene_id, Chromosome = seqid)

    rowData(sce) <- merge(rowData(sce), annot, by = "ID", sort=FALSE)
    rownames(rowData(sce)) <- rowData(sce)$ID


    # qc addPerCellQC
    is.mito <- rowData(sce)$Chromosome%in%c("chrM", "MT")
    sce <- addPerCellQC(sce, subsets=list(Mito=is.mito), BPPARAM = bp.params)


    # adpative thresholds with Sample as batch
    cell_qc_results <- quickPerCellQC(colData(sce),
                                    percent_subsets=c("subsets_Mito_percent"),
                                    batch=sce$Sample)

    sce <- sce[, !cell_qc_results$discard]

    # Remove the per cell QC metrics
    colData(sce) <- colData(sce)[,1:4]

    # Remove undetected genes
    keep <- rowSums(counts(sce)) > 0
    sce <- sce[keep,]

    # Normalize the counts
    set.seed(100) 
    clust <- quickCluster(sce) 
    sce <- computePooledFactors(sce, cluster=clust, min.mean=0.1)
    sce <- logNormCounts(sce)

    # Identify highly variable genes
    gene_var <- modelGeneVar(sce)
    hvgs <- getTopHVGs(gene_var, prop=0.1)

    # Dimension reduction
    sce <- runPCA(sce, subset_row = hvgs)
    sce <- runTSNE(sce, dimred="PCA", n_dimred=10, BPPARAM = bp.params)
    sce <- runUMAP(sce, dimred="PCA", n_dimred=10, BPPARAM = bp.params)

    # Clustering
    colLabels(sce) <- clusterCells(sce, use.dimred = "PCA")
    sce
}

pbmmc1a <- prepSample("SRR9264351")
pbmmc1b <- prepSample("SRR9264352")

saveRDS(pbmmc1a, "data/R_objects/PBMMC_1a_dimRed.rds")
saveRDS(pbmmc1b, "data/R_objects/PBMMC_1b_dimRed.rds")




