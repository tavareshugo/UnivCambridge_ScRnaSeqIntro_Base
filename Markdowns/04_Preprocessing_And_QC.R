# load packages 
library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)
library(patchwork)
library(ggvenn)

# load samplesheet 
samplesheet <- read_tsv("Data/sample_sheet.tsv")


# view metadata 
samplesheet %>%
	as.data.frame() %>%
	datatable(rownames = FALSE, options = list(dom="tpl", nrows=20))


# set up parallelisation 
bp.params <- MulticoreParam(workers = 7)


# load a single sample 
sample.path <- "Data/CellRanger_Outputs/SRR9264343/outs/filtered_feature_bc_matrix/"
sce.sing <- read10xCounts(sample.path, col.names=TRUE, BPPARAM = bp.params)
sce.sing


# check dimensions 
dim(counts(sce.sing))


# the counts matrix 
counts(sce.sing)[1:10, 1:10]


# feature metadata 
rowData(sce.sing)

rownames(counts(sce.sing))[1:6]


# droplet metadata 
colData(sce.sing)

colnames(counts(sce.sing))[1:6]


# genes_per_cell 
genesPerCell <- colSums(counts(sce.sing) > 0)


# plot_genes_per_cell 
genesPerCell <- colSums(counts(sce.sing) > 0)
plot(density(genesPerCell), main="", xlab="Genes per cell")


# expression_v_detected 
plot(rowSums(counts(sce.sing)) / rowSums(counts(sce.sing) > 0),
     rowMeans(counts(sce.sing) > 0),
     log = "x",
     xlab="Mean UMIs per cell",
     ylab="proportion of cells expressing the gene"
)


# top_20_genes 
rel_expression <- t( t(counts(sce.sing)) / colSums(counts(sce.sing))) * 100
rownames(rel_expression) <- rowData(sce.sing)$Symbol
most_expressed <- sort(rowSums( rel_expression ), decreasing = T)[20:1]
plot_data <- as.matrix(t(rel_expression[names(most_expressed),]))

boxplot(plot_data, cex=0.1, las=1, xlab="% total count per cell", horizontal=TRUE)

# make_file_list 
list_of_files <- samplesheet %>% 
    group_by(SampleGroup) %>%  
    slice(1) %>%  
    mutate(File=str_c("Data/CellRanger_Outputs/", 
                       Sample, 
                       "/outs/filtered_feature_bc_matrix")) %>%
    pull(File, Sample)
list_of_files


# load_data_sets 
sce <- read10xCounts(list_of_files, col.names=TRUE, BPPARAM = bp.params)
sce


# four_samples_colData 
colData(sce)


# add_metadata 
sce$Barcode <- rownames(colData(sce))
colData(sce) <- merge(colData(sce), samplesheet, by="Sample", sort=FALSE)
rownames(colData(sce)) <- sce$Barcode


# detected_genes 
detected_genes <- rowSums(counts(sce)) > 0
table(detected_genes)


# remove_undetected_genes 
sce <- sce[detected_genes,]


# annotate_genes 
ah <- AnnotationHub()
ens.hs.107<- query(ah, c("Homo sapiens", "EnsDb", 107))[[1]] 

genes <- rowData(sce)$ID
gene_annot <- AnnotationDbi::select(ens.hs.107, 
                                    keys = genes,
                                    keytype = "GENEID",
                                    columns = c("GENEID", "SEQNAME")) %>%
    set_names(c("ID", "Chromosome"))
rowData(sce) <- merge(rowData(sce), gene_annot, by = "ID", sort=FALSE)
rownames(rowData(sce)) <- rowData(sce)$ID

rowData(sce)


# add_per_cell_QC 
is.mito <- which(rowData(sce)$Chromosome=="MT")

sce <- addPerCellQC(sce, subsets=list(Mito=is.mito), BPPARAM = bp.params)


# cell_qc_table 
colData(sce)


# plot_total_counts 
plotColData(sce, x="SampleName", y="sum") + 
    scale_y_log10() + 
    ggtitle("Total count")


# plot_detected_genes 
plotColData(sce, x="SampleName", y="detected") + 
    scale_y_log10() + 
    ggtitle("Detected features")


# plot_MT_content 
plotColData(sce, x="SampleName", y="subsets_Mito_percent") + 
    ggtitle("Mito percent")


# plot_genes_v_library_size 
colData(sce) %>% 
    as.data.frame() %>% 
    arrange(subsets_Mito_percent) %>% 
    ggplot(aes(x = sum, y = detected)) +
      geom_point(aes(colour = subsets_Mito_percent > 10)) + 
      facet_wrap(vars(SampleGroup))


# outliers_library_size 
low_lib_size <- isOutlier(sce$sum, log=TRUE, type="lower")
table(low_lib_size)


# libSize_thresholds 
attr(low_lib_size, "thresholds")


# plot_library_size_filtering 
colData(sce)$low_lib_size <- low_lib_size
plotColData(sce, x="SampleName", y="sum", colour_by = "low_lib_size") + 
    scale_y_log10() + 
    labs(y = "Total count", title = "Total count") +
    guides(colour=guide_legend(title="Discarded"))


# outliers_detected_genes 
low_n_features <- isOutlier(sce$detected, log=TRUE, type="lower")
table(low_n_features)


# detected_genes_thresholds 
attr(low_n_features, "thresholds")[1]


# plot_detected_genes_filtering 
colData(sce)$low_n_features <- low_n_features
plotColData(sce, x="SampleName", y="detected", colour_by = "low_n_features") + 
    scale_y_log10() + 
    labs(y = "Genes detected", title = "Genes detected") +
    guides(colour=guide_legend(title="Discarded"))


# outlier_MT_content 
high_Mito_percent <- isOutlier(sce$subsets_Mito_percent, type="higher")
table(high_Mito_percent)


# MT_content_thresholds 
attr(high_Mito_percent, "thresholds")[2]


# plot_MT_content_filtering 
colData(sce)$high_Mito_percent <- high_Mito_percent
plotColData(sce,
            x="SampleName",
            y="subsets_Mito_percent",
            colour_by = "high_Mito_percent") + 
    labs(y = "Percentage mitochondrial UMIs",
         title = "Mitochondrial UMIs") +
    guides(colour=guide_legend(title="Discarded"))


# outliers_summary_table 
data.frame(`Library Size` = sum(low_lib_size),
           `Genes detected` = sum(low_n_features),
           `Mitochondrial UMIs` = sum(high_Mito_percent),
           Total = sum(low_lib_size | low_n_features | high_Mito_percent))


# quickPerCellQC 
cell_qc_results <- quickPerCellQC(colData(sce),
			  percent_subsets=c("subsets_Mito_percent"))
cell_qc_results %>%
  as.data.frame() %>% 
  mutate(SampleName=colData(sce)$SampleName) %>% 
  group_by(SampleName)  %>%
  summarise(across(contains("_"), sum))


# quickPerCellQC_batch 
batch.cell_qc_results <- quickPerCellQC(colData(sce),
                                percent_subsets=c("subsets_Mito_percent"),
                                batch=sce$Sample)

batch.cell_qc_results %>%
  as.data.frame() %>% 
  mutate(SampleName=colData(sce)$SampleName) %>% 
  group_by(SampleName)  %>%
  summarise(across(contains("_"), sum))


# compare_thresholds 
all.thresholds <- tibble(`SampleName`="All",
       `Library Size`=attr(cell_qc_results$low_lib_size, "thresholds")[1],
       `Genes detected`=attr(cell_qc_results$low_n_features, "thresholds")[1],
       `Mitochondrial UMIs`=attr(cell_qc_results$high_subsets_Mito_percent, "thresholds")[2])


tibble(`Sample`=names(attr(batch.cell_qc_results$low_lib_size, "thresholds")[1,]),
       `Library Size`=attr(batch.cell_qc_results$low_lib_size, "thresholds")[1,],
       `Genes detected`=attr(batch.cell_qc_results$low_n_features, "thresholds")[1,],
       `Mitochondrial UMIs`=attr(batch.cell_qc_results$high_subsets_Mito_percent, "thresholds")[2,]) %>% 
    left_join(samplesheet) %>% 
    select(SampleName, `Library Size`, `Genes detected`, `Mitochondrial UMIs`) %>% 
    bind_rows(all.thresholds) %>% 
    mutate(across(where(is.numeric), round, digits=2)) %>% 
    datatable(rownames = FALSE, options = list(dom="t"))


# replace_filters_in_sce 
sce$low_lib_size <- batch.cell_qc_results$low_lib_size
sce$low_n_features <- batch.cell_qc_results$low_n_features
sce$high_Mito_percent <- batch.cell_qc_results$high_subsets_Mito_percent
sce$discard <- batch.cell_qc_results$discard


# plot_library_size_batch_filters 
plotColData(sce, x="SampleName", y="sum", colour_by = "low_lib_size") + 
    scale_y_log10() + 
    labs(y = "Total count", title = "Total count") +
    guides(colour=guide_legend(title="Discarded"))


# plot_detected_genes_batch_filters 
plotColData(sce, x="SampleName", y="detected", colour_by = "low_n_features") + 
    scale_y_log10() + 
    labs(y = "Genes detected", title = "Genes detected") +
    guides(colour=guide_legend(title="Discarded"))


# plot_MT_content_batch_filters 
plotColData(sce, 
        x="Sample", 
        y="subsets_Mito_percent",
        colour_by = "high_Mito_percent") + 
    labs(y = "Percentage mitochondrial UMIs",
         title = "Mitochondrial UMIs") +
    guides(colour=guide_legend(title="Discarded"))


# filtering_venns 
pc1 <- tibble(`All together`=cell_qc_results$low_lib_size, 
              `By batch`=batch.cell_qc_results$low_lib_size) %>% 
           ggvenn(show_percentage = FALSE) +
               labs(title="Library Size")

pc2 <- tibble(`All together`=cell_qc_results$low_n_features, 
              `By batch`=batch.cell_qc_results$low_n_features) %>% 
           ggvenn(show_percentage = FALSE) +
               labs(title="Genes detected")

pc3 <- tibble(`All together`=cell_qc_results$high_subsets_Mito_percent, 
                 `By batch`=batch.cell_qc_results$high_subsets_Mito_percent) %>% 
           ggvenn(show_percentage = FALSE) +
               labs(title="Mitochondrial UMIs")

pc1 + pc2 + pc3


# remove_cells 
sce.filtered <- sce[, !sce$discard]


# rerun_per_cell_qc 
colData(sce.filtered) <- colData(sce.filtered)[,1:3]
sce.filtered <- addPerCellQC(sce.filtered, BPPARAM = bp.params)


# subset_PBMMC 
sce.BM1 <- sce[ , sce$SampleName == "PBMMC_1a"]


# outlyingness_filter 
library(robustbase)
stats <- cbind(log10(sce.BM1$sum),
               log10(sce.BM1$detected),
               sce.BM1$subsets_Mito_percent)

outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
multi.outlier <- isOutlier(outlying, type = "higher")
summary(multi.outlier)


# PCA_filter 
sce.BM1$log10sum <- log10(sce.BM1$sum)
sce.BM1$log10detected <- log10(sce.BM1$detected)
sce.BM1 <- runColDataPCA(sce.BM1, 
                     variables=list("log10sum", 
                                    "log10detected", 
                                    "subsets_Mito_percent"),
                     outliers=TRUE,
			         BPPARAM = bp.params)


# pca_results 
head(reducedDim(sce.BM1))


# outlyingness_results 
summary(sce.BM1$outlier)


# plot_MT_content_v_library_size 
plotColData(sce, 
            x="sum", 
            y="subsets_Mito_percent", 
            other_fields="SampleName",
            colour_by="discard") +
    facet_wrap(~SampleName, ncol=5, scale="free_x")


# add_per_feature_QC 
sce <- addPerFeatureQC(sce, BPPARAM = bp.params)
rowData(sce)


# compute_sparsity 
colData(sce)$cell_sparsity <- 1 - (colData(sce)$detected / nrow(sce))
rowData(sce)$gene_sparsity <- (100 - rowData(sce)$detected) / 100


# plot_cell_sparsity_histogram 
hist(sce$cell_sparsity, breaks=50, col="grey80", xlab="Cell sparsity")


# plot_gene_sparsity_histogram 
hist(rowData(sce)$gene_sparsity, breaks=50, col="grey80", xlab="Gene sparsity")


# filter_by_sparsity 
sparse.cells <- sce$cell_sparsity > 0.99
mito.cells <- sce$subsets_Mito_percent > 10

min.cells <- 1 - (20 / ncol(sce))
sparse.genes <- rowData(sce)$gene_sparsity > min.cells


# sparsity_gene_filter_results 
table(sparse.genes)


# sparsity_cell_filter_results 
table(sparse.cells, mito.cells)


#  
sessionInfo()

