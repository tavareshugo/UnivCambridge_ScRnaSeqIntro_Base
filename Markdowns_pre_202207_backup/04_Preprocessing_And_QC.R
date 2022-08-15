#############################################################################
## Load libraries
library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)
library(patchwork)
library(ggvenn)
#############################################################################


#############################################################################
## Read sample sheet
samplesheet <- read_tsv("Data/sample_sheet.tsv")
#############################################################################

samplesheet %>%
	as.data.frame() %>%
	datatable(rownames = FALSE, options = list(dom="tpl", nrows=20))


## ----------------------------------------------------------------------------------
bp.params <- MulticoreParam(workers = 7)


## ----example_load------------------------------------------------------------------
sample.path <- "CellRanger_Outputs/SRR9264343/outs/filtered_feature_bc_matrix/"
sce.sing <- read10xCounts(sample.path, col.names=TRUE, BPPARAM = bp.params)
sce.sing


## ----example_dim-------------------------------------------------------------------
dim(counts(sce.sing))


## ----countsMatrix------------------------------------------------------------------
counts(sce.sing)[1:10, 1:10]


## ----example_row-------------------------------------------------------------------
rowData(sce.sing)


## ----example_col-------------------------------------------------------------------
colData(sce.sing)


## ----------------------------------------------------------------------------------
colnames(counts(sce.sing))[1:6]


## ----echo=FALSE--------------------------------------------------------------------
genesPerCell <- colSums(counts(sce.sing) > 0)


## ----properties_distribNbGenesDetected---------------------------------------------
genesPerCell <- colSums(counts(sce.sing) > 0)
summary(genesPerCell)

plot(density(genesPerCell), main="", xlab="Genes per cell")


## ----------------------------------------------------------------------------------
tmpCounts <- counts(sce.sing)[,1:1000]

plot(rowSums(tmpCounts),
     rowMeans(tmpCounts > 0),
     log = "x",
     xlab="total number of UMIs",
     ylab="proportion of cells expressing the gene"
)
rm(tmpCounts)


## ----fig.width = 8, fig.height = 12------------------------------------------------
rel_expression <- t( t(counts(sce.sing)) / colSums(counts(sce.sing))) * 100
rownames(rel_expression) <- rowData(sce.sing)$Symbol
most_expressed <- sort(rowSums( rel_expression ),T)[20:1]
plot_data <- as.matrix(t(rel_expression[names(most_expressed),]))

boxplot(plot_data, cex=0.1, las=1, xlab="% total count per cell", horizontal=TRUE)


## ----echo=FALSE--------------------------------------------------------------------
rm(rel_expression, plot_data)
ncolRaw <- ncol(sce.sing)
rm(sce.sing)


## ----make_file_list----------------------------------------------------------------
samples_list <- samplesheet %>% 
    group_by(SampleGroup) %>%  
    slice(1) %>%  
    pull(SampleId)
list_of_files <- str_c("CellRanger_Outputs/", 
                       samples_list, 
                       "/outs/filtered_feature_bc_matrix")
names(list_of_files) <- samples_list
list_of_files


## ----load_data_sets----------------------------------------------------------------
sce <- read10xCounts(list_of_files, col.names=TRUE, BPPARAM = bp.params)
sce


## ----------------------------------------------------------------------------------
colData(sce)


## ----dataSets_addSampleSheet-------------------------------------------------------
colData(sce) <- colData(sce) %>% 
    as.data.frame() %>%
    rownames_to_column("RowName") %>% 
    mutate(SampleNum = str_extract(RowName, "^[0-9]+")) %>%
    mutate(Barcode = str_replace(Barcode, "1$", SampleNum)) %>%
    left_join(samplesheet, by=c(Sample="SampleId")) %>%
    rename(SampleId=Sample) %>% 
    rename(Sample=SampleName) %>%    
#     mutate(Sample = case_when(
#         SampleId == "SRR9264351" ~ str_c(Sample, "a"),
#         SampleId == "SRR9264352" ~ str_c(Sample, "b"),
#         TRUE ~ Sample)) %>% 
    column_to_rownames("RowName") %>% 
    select(Sample, Barcode, SampleId, SampleGroup, DatasetName) %>%
    DataFrame()


## ----------------------------------------------------------------------------------
colData(sce)


## ----detected_genes----------------------------------------------------------------
detected_genes <- rowSums(counts(sce)) > 0
table(detected_genes)


## ----remove_undetected_genes-------------------------------------------------------
sce <- sce[detected_genes,]


## ----annotate_genes----------------------------------------------------------------
ah <- AnnotationHub()
ens.mm.98 <- query(ah, c("Homo sapiens", "EnsDb", 98))[[1]] 

genes <- rowData(sce)$ID
gene_annot <- AnnotationDbi::select(ens.mm.98, 
                                    keys = genes,
                                    keytype = "GENEID",
                                    columns = c("GENEID", "SEQNAME")) %>%
    set_names(c("ID", "Chromosome"))
rowData(sce) <- merge(rowData(sce), gene_annot, by = "ID", sort=FALSE)
rownames(rowData(sce)) <- rowData(sce)$ID

rowData(sce)


## ----qc_addPerCellQC---------------------------------------------------------------
is.mito <- which(rowData(sce)$Chromosome=="MT")

sce <- addPerCellQC(sce, subsets=list(Mito=is.mito), BPPARAM = bp.params)


## ----qc_addPerCellQCTab, eval=TRUE-------------------------------------------------
colData(sce)


## ---- fig.width=12, fig.height=4---------------------------------------------------
plotColData(sce, x="Sample", y="sum",other_fields="SampleGroup") + 
    facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
    scale_y_log10() + 
    ggtitle("Total count")


## ---- fig.width=12, fig.height=4---------------------------------------------------
plotColData(sce, x="Sample", y="detected", other_fields="SampleGroup") + 
    facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
    scale_y_log10() + 
    ggtitle("Detected features")


## ---- fig.width=12, fig.height=4---------------------------------------------------
plotColData(sce, x="Sample", y="subsets_Mito_percent", other_fields="SampleGroup") + 
    facet_wrap(~SampleGroup, nrow=1, scales = "free_x") +
    ggtitle("Mito percent")


## ---- fig.width=10, fig.height=6---------------------------------------------------
colData(sce) %>% 
    as.data.frame() %>% 
    arrange(subsets_Mito_percent) %>% 
    ggplot(aes(x = sum, y = detected)) +
      geom_point(aes(colour = subsets_Mito_percent > 10)) + 
      facet_wrap(vars(SampleGroup))


## ----adapThresTab_libSize----------------------------------------------------------
low_lib_size <- isOutlier(sce$sum, log=TRUE, type="lower")
table(low_lib_size)


## ----adapThresVal_libSize----------------------------------------------------------
attr(low_lib_size, "thresholds")


## ---- fig.width=12, fig.height=5---------------------------------------------------
colData(sce)$low_lib_size <- low_lib_size
plotColData(sce, 
            x="Sample", 
            y="sum",
            other_fields="SampleGroup", 
            colour_by = "low_lib_size") + 
    facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
    scale_y_log10() + 
    labs(y = "Total count", title = "Total count") +
    guides(colour=guide_legend(title="Discarded"))


## ----adapThresTab_detected---------------------------------------------------------
low_n_features <- isOutlier(sce$detected, log=TRUE, type="lower")
table(low_n_features)


## ----adapThresVal_detected---------------------------------------------------------
attr(low_n_features, "thresholds")[1]


## ---- fig.width=12, fig.height=5---------------------------------------------------
colData(sce)$low_n_features <- low_n_features
plotColData(sce, 
            x="Sample", 
            y="detected",
            other_fields="SampleGroup", 
            colour_by = "low_n_features") + 
    facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
    scale_y_log10() + 
    labs(y = "Genes detected", title = "Genes detected") +
    guides(colour=guide_legend(title="Discarded"))


## ----adapThresTab_mito-------------------------------------------------------------
high_Mito_percent <- isOutlier(sce$subsets_Mito_percent, type="higher")
table(high_Mito_percent)


## ----adapThresVal_mito-------------------------------------------------------------
attr(high_Mito_percent, "thresholds")[2]


## ---- fig.width=12, fig.height=5---------------------------------------------------
colData(sce)$high_Mito_percent <- high_Mito_percent
plotColData(sce,  
            x="Sample",
            y="subsets_Mito_percent",
            other_fields="SampleGroup",
            colour_by = "high_Mito_percent") + 
    facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
    labs(y = "Percentage mitochondrial UMIs",
         title = "Mitochondrial UMIs") +
    guides(colour=guide_legend(title="Discarded"))


## ----adapThres_summary-------------------------------------------------------------
data.frame(`Library Size` = sum(low_lib_size),
           `Genes detected` = sum(low_n_features),
           `Mitochondrial UMIs` = sum(high_Mito_percent),
           Total = sum(low_lib_size | low_n_features | high_Mito_percent))


## ----adapThres_quickPerCellQC------------------------------------------------------
cell_qc_results <- quickPerCellQC(colData(sce),
			  percent_subsets=c("subsets_Mito_percent"))
colSums(as.data.frame(cell_qc_results))


## ----quickPerCellQC_batch_compute--------------------------------------------------
batch.cell_qc_results <- quickPerCellQC(colData(sce),
                                percent_subsets=c("subsets_Mito_percent"),
                                batch=sce$DatasetName)
colSums(as.data.frame(batch.cell_qc_results))


## ----------------------------------------------------------------------------------
all.thresholds <- tibble(`Batch`="All",
       `Library Size`=attr(cell_qc_results$low_lib_size, "thresholds")[1],
       `Genes detected`=attr(cell_qc_results$low_n_features, "thresholds")[1],
       `Mitochondrial UMIs`=attr(cell_qc_results$high_subsets_Mito_percent, "thresholds")[2])


tibble(`Batch`=names(attr(batch.cell_qc_results$low_lib_size, "thresholds")[1,]),
       `Library Size`=attr(batch.cell_qc_results$low_lib_size, "thresholds")[1,],
       `Genes detected`=attr(batch.cell_qc_results$low_n_features, "thresholds")[1,],
       `Mitochondrial UMIs`=attr(batch.cell_qc_results$high_subsets_Mito_percent, "thresholds")[2,]) %>% 
    bind_rows(all.thresholds) %>% 
    mutate(across(where(is.numeric), round, digits=2)) %>% 
    datatable(rownames = FALSE, options = list(dom="t"))


## ----quickPerCellQC_batch_replace--------------------------------------------------
sce$low_lib_size <- batch.cell_qc_results$low_lib_size
sce$low_n_features <- batch.cell_qc_results$low_n_features
sce$high_Mito_percent <- batch.cell_qc_results$high_subsets_Mito_percent
sce$discard <- batch.cell_qc_results$discard


## ---- fig.width=12, fig.height=4---------------------------------------------------
plotColData(sce, 
            x="Sample", 
            y="sum",
            other_fields="SampleGroup", 
            colour_by = "low_lib_size") + 
    facet_wrap(vars(SampleGroup), nrow=1, scales = "free_x") + 
    scale_y_log10() + 
    labs(y = "Total count", title = "Total count") +
    guides(colour=guide_legend(title="Discarded"))


## ---- fig.width=12, fig.height=4---------------------------------------------------
plotColData(sce, 
            x="Sample", 
            y="detected",
            other_fields="SampleGroup", 
            colour_by = "low_n_features") + 
    facet_wrap(vars(SampleGroup), nrow=1, scales = "free_x") + 
    scale_y_log10() + 
    labs(y = "Genes detected", title = "Genes detected") +
    guides(colour=guide_legend(title="Discarded"))


## ---- fig.width=12, fig.height=4---------------------------------------------------
plotColData(sce, 
        x="Sample", 
        y="subsets_Mito_percent",
        other_fields="SampleGroup", 
        colour_by = "high_Mito_percent") + 
    facet_wrap(vars(SampleGroup), nrow=1, scales = "free_x") + 
    labs(y = "Percentage mitochondrial UMIs",
         title = "Mitochondrial UMIs") +
    guides(colour=guide_legend(title="Discarded"))


## ---- fig.width=12, fig.height=8---------------------------------------------------
libDat <- tibble(`All together`=cell_qc_results$low_lib_size, 
                 `By batch`=batch.cell_qc_results$low_lib_size,
                 Batch=sce$Sample)
    
ph1 <- libDat %>% 
    dplyr::filter(Batch=="HCA") %>% 
    ggvenn(show_percentage = FALSE) +
        labs(title="Library Size - HCA")
pc1 <- libDat %>% 
    dplyr::filter(Batch=="Caron") %>% 
    ggvenn(show_percentage = FALSE) +
        labs(title="Library Size - Caron")

nGenDat <- tibble(`All together`=cell_qc_results$low_n_features, 
                  `By batch`=batch.cell_qc_results$low_n_features,
                 Batch=sce$DatasetName)
ph2 <- nGenDat %>% 
    dplyr::filter(Batch=="HCA") %>% 
        ggvenn(show_percentage = FALSE) +
            labs(title="Genes detected - HCA")
pc2 <- nGenDat %>% 
    dplyr::filter(Batch=="Caron") %>% 
           ggvenn(show_percentage = FALSE) +
            labs(title="Genes detected - Caron")


mitDat <- tibble(`All together`=cell_qc_results$high_subsets_Mito_percent, 
       `By batch`=batch.cell_qc_results$high_subsets_Mito_percent,
                 Batch=sce$DatasetName)
ph3 <- mitDat %>% 
    dplyr::filter(Batch=="HCA") %>% 
        ggvenn(show_percentage = FALSE) +
            labs(title="Mitochondrial UMIs - HCA")
pc3 <- mitDat %>% 
    dplyr::filter(Batch=="Caron") %>% 
           ggvenn(show_percentage = FALSE) +
            labs(title="Mitochondrial UMIs - Caron")

(pc1 + pc2 + pc3) / (ph1 + ph2 + ph3)


## ----------------------------------------------------------------------------------
sce.filtered <- sce[, !sce$discard]


## ----------------------------------------------------------------------------------
colData(sce.filtered) <- colData(sce.filtered)[,1:3]
sce.filtered <- addPerCellQC(sce.filtered, BPPARAM = bp.params)


## ----------------------------------------------------------------------------------
sce.BM1 <- sce[ , sce$Sample == "PBMMC_1"]


## ----------------------------------------------------------------------------------
library(robustbase)
stats <- cbind(log10(sce.BM1$sum),
               log10(sce.BM1$detected),
               sce.BM1$subsets_Mito_percent)

outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
multi.outlier <- isOutlier(outlying, type = "higher")
summary(multi.outlier)


## ----------------------------------------------------------------------------------
sce.BM1$log10sum <- log10(sce.BM1$sum)
sce.BM1$log10detected <- log10(sce.BM1$detected)
sce.BM1 <- runColDataPCA(sce.BM1, 
                     variables=list("log10sum", "log10detected", "subsets_Mito_percent"),
                     outliers=TRUE,
			         BPPARAM = bp.params)


## ----------------------------------------------------------------------------------
head(reducedDim(sce.BM1))


## ----------------------------------------------------------------------------------
summary(sce.BM1$outlier)


## ----fig.width = 12, fig.height = 4------------------------------------------------
plotColData(sce, 
            x="sum", 
            y="subsets_Mito_percent", 
            other_fields="Sample",
            colour_by="discard") +
    facet_wrap(~Sample, ncol=5, scale="free_x")


## ----qc_addPerFeatureQC, eval=TRUE-------------------------------------------------
sce <- addPerFeatureQC(sce, BPPARAM = bp.params)
rowData(sce)


## ----sparsity_compute--------------------------------------------------------------
colData(sce)$cell_sparsity <- 1 - (colData(sce)$detected / nrow(sce))
rowData(sce)$gene_sparsity <- (100 - rowData(sce)$detected) / 100


## ----------------------------------------------------------------------------------
hist(sce$cell_sparsity, breaks=50, col="grey80", xlab="Cell sparsity", main="")


## ----------------------------------------------------------------------------------
hist(rowData(sce)$gene_sparsity, breaks=50, col="grey80", xlab="Gene sparsity", main="")


## ----sparsity_filter---------------------------------------------------------------
sparse.cells <- sce$cell_sparsity > 0.99
mito.cells <- sce$subsets_Mito_percent > 10

min.cells <- 1 - (20 / ncol(sce))
sparse.genes <- rowData(sce)$gene_sparsity > min.cells


## ----------------------------------------------------------------------------------
table(sparse.genes)


## ----------------------------------------------------------------------------------
table(sparse.cells, mito.cells)


## ----------------------------------------------------------------------------------
sessionInfo()

