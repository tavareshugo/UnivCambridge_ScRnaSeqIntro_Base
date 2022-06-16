#############################################################################
## 2 Load packages
library(DropletUtils)
library(scater)
library(ensembldb)
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)
library(patchwork)
library(ggvenn)
library(DT)
#############################################################################


#############################################################################
## 3 Reading CellRanger output into R
# Read sample sheet
samplesheet <- read_tsv("Data/sample_sheet.tsv")

# display samplesheet
samplesheet %>%
  datatable(rownames = FALSE, options = list(dom="tpl", nrows=20))
#############################################################################

#############################################################################
## 3.1.1 Parallelisation
# Parallelisation package BiocParallel
bp.params <- MulticoreParam(workers = 7)
#############################################################################

#############################################################################
## 3.1.2 Loading a single sample
sample.path <- "CellRanger_Outputs/SRR9264343/outs/filtered_feature_bc_matrix/"
sce.sing <- read10xCounts(sample.path, col.names=TRUE, BPPARAM = bp.params)
sce.sing
#############################################################################

#############################################################################
## 3.1.4 The counts matrix
dim(counts(sce.sing))
# View first 10 roes and columns of a matrix
counts(sce.sing)[1:10, 1:10]
#############################################################################

#############################################################################
## 3.1.5 Features
rowData(sce.sing)
#############################################################################

#############################################################################
## 3.1.6 Droplet annotation
colData(sce.sing)
# View first 6 column names of a count matrix
colnames(counts(sce.sing))[1:6]
#############################################################################



#############################################################################
## 4 Properties of scRNA-seq data
# Number of genes detected per cell
genesPerCell <- colSums(counts(sce.sing) > 0)
summary(genesPerCell)

# density plot
plot(density(genesPerCell), main="", xlab="Genes per cell")

# Total UMI for a gene versus the number of times detected
tmpCounts <- counts(sce.sing)[,1:1000]

plot(rowSums(tmpCounts),
     rowMeans(tmpCounts > 0),
     log = "x",
     xlab="total number of UMIs / gene",
     ylab="proportion of cells expressing the gene"
)
rm(tmpCounts)

# Distribution of counts for a gene across cells
rel_expression <- t( t(counts(sce.sing)) / colSums(counts(sce.sing))) * 100
rownames(rel_expression) <- rowData(sce.sing)$Symbol
most_expressed <- sort(rowSums( rel_expression ),T)[20:1]
plot_data <- as.matrix(t(rel_expression[names(most_expressed),]))

boxplot(plot_data, cex=0.1, las=1, xlab="% total count per cell", horizontal=TRUE)


# clear objects
rm(rel_expression, plot_data)
ncolRaw <- ncol(sce.sing)
rm(sce.sing)
#############################################################################
#--------- End of single sample exploration ------------------------------- 



#############################################################################
## 5.1 Load multiple samples
# create files list
samples_list <- samplesheet %>% 
  group_by(SampleGroup) %>%  
  slice(1) %>%  
  pull(SampleId)
list_of_files <- str_c("CellRanger_Outputs/", 
                       samples_list, 
                       "/outs/filtered_feature_bc_matrix")
names(list_of_files) <- samples_list
list_of_files

# read 10x data into R and create a sce object
sce <- read10xCounts(list_of_files, col.names=TRUE, BPPARAM = bp.params)
sce

# 5.2 Modify the droplet annotation
# explore col. data before modification
colData(sce)

# modify barcode names
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

# explore col. data after modification
colData(sce)



# 5.3 Undetected genes
# remove undetected genes
detected_genes <- rowSums(counts(sce)) > 0
table(detected_genes)
sce <- sce[detected_genes,]


## 5.4 Annotate genes
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


# 5.5 Add per cell QC metrics
# Is gene a MT gene?
is.mito <- which(rowData(sce)$Chromosome=="MT")

# add per cell QC metrics col. data
sce <- addPerCellQC(sce, subsets=list(Mito=is.mito), BPPARAM = bp.params)

# explore col. data
colData(sce)

# QC metric distribution
# Per cell UMT counts distribution / sample
plotColData(sce, x="Sample", y="sum",other_fields="SampleGroup") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  ggtitle("Total count")

# Per cell gene counts distribution / sample
plotColData(sce, x="Sample", y="detected", other_fields="SampleGroup") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  ggtitle("Detected features")

# Mito. gene percentage distrubution / sample
plotColData(sce, x="Sample", y="subsets_Mito_percent", other_fields="SampleGroup") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") +
  ggtitle("Mito percent")

#  library size vs numbers of genes detected
colData(sce) %>% 
  as.data.frame() %>% 
  arrange(subsets_Mito_percent) %>% 
  ggplot(aes(x = sum, y = detected)) +
  geom_point(aes(colour = subsets_Mito_percent > 10)) + 
  facet_wrap(vars(SampleGroup))


# 5.7 Identification of low-quality cells with adaptive thresholds
# adaptive thresholds for total UMI counts
low_lib_size <- isOutlier(sce$sum, log=TRUE, type="lower")
table(low_lib_size)

# get adaptive thresholds value for UMT counts
attr(low_lib_size, "thresholds")


# add low_lib_size filter column to col. data
colData(sce)$low_lib_size <- low_lib_size

# view the effect of the filtering using plotColData
plotColData(sce, 
            x="Sample", 
            y="sum",
            other_fields="SampleGroup", 
            colour_by = "low_lib_size") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  labs(y = "Total count", title = "Total count") +
  guides(colour=guide_legend(title="Discarded"))


# 5.7.2 Number of genes
# detect outliers based on number of genes detected
low_n_features <- isOutlier(sce$detected, log=TRUE, type="lower")
table(low_n_features)

# get adaptive thresholds value for detected genes
attr(low_n_features, "thresholds")[1]

# add low_n_features filter column to col. data
colData(sce)$low_n_features <- low_n_features

# view the effect of the filtering using plotColData
plotColData(sce, 
            x="Sample", 
            y="detected",
            other_fields="SampleGroup", 
            colour_by = "low_n_features") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  labs(y = "Genes detected", title = "Genes detected") +
  guides(colour=guide_legend(title="Discarded"))



# 5.7.3 Mitochondrial content
# detect outliers based on Mitochondrial content
high_Mito_percent <- isOutlier(sce$subsets_Mito_percent, type="higher")
table(high_Mito_percent)


# get adaptive thresholds value for Mitochondrial content
attr(high_Mito_percent, "thresholds")[2]


# add high_Mito_percent filter column to col. data
colData(sce)$high_Mito_percent <- high_Mito_percent

# view the effect of the filtering using plotColData.
plotColData(sce,  
            x="Sample",
            y="subsets_Mito_percent",
            other_fields="SampleGroup",
            colour_by = "high_Mito_percent") + 
  facet_wrap(~SampleGroup, nrow=1, scales = "free_x") + 
  labs(y = "Percentage mitochondrial UMIs",
       title = "Mitochondrial UMIs") +
  guides(colour=guide_legend(title="Discarded"))


# 5.7.4 Summary of discarded cells
data.frame(`Library Size` = sum(low_lib_size),
           `Genes detected` = sum(low_n_features),
           `Mitochondrial UMIs` = sum(high_Mito_percent),
           Total = sum(low_lib_size | low_n_features | high_Mito_percent))


# 5.7.5 All three filter steps at once
# The three steps above can be run in one go
cell_qc_results <- quickPerCellQC(colData(sce),
                                  percent_subsets=c("subsets_Mito_percent"))
colSums(as.data.frame(cell_qc_results))


# 5.7.6 Assumptions

# 5.7.7 Considering experimental factors when filtering
batch.cell_qc_results <- quickPerCellQC(colData(sce),
                                        percent_subsets=c("subsets_Mito_percent"),
                                        batch=sce$DatasetName)

# get summary of  filter data
colSums(as.data.frame(batch.cell_qc_results))


# thresholds for each metric differ between 
# the batch-wise analysis and the analysis using all samples.
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


# replace the columns in the droplet annotation with these new filters 
sce$low_lib_size <- batch.cell_qc_results$low_lib_size
sce$low_n_features <- batch.cell_qc_results$low_n_features
sce$high_Mito_percent <- batch.cell_qc_results$high_subsets_Mito_percent
sce$discard <- batch.cell_qc_results$discard


# visualise how the new filters look using violin plots
# cell UMI counts / sample
plotColData(sce, 
            x="Sample", 
            y="sum",
            other_fields="SampleGroup", 
            colour_by = "low_lib_size") + 
  facet_wrap(vars(SampleGroup), nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  labs(y = "Total count", title = "Total count") +
  guides(colour=guide_legend(title="Discarded"))


# Detected genes per cell / sample
plotColData(sce, 
            x="Sample", 
            y="detected",
            other_fields="SampleGroup", 
            colour_by = "low_n_features") + 
  facet_wrap(vars(SampleGroup), nrow=1, scales = "free_x") + 
  scale_y_log10() + 
  labs(y = "Genes detected", title = "Genes detected") +
  guides(colour=guide_legend(title="Discarded"))


# Mito percent per cell / sample
plotColData(sce, 
            x="Sample", 
            y="subsets_Mito_percent",
            other_fields="SampleGroup", 
            colour_by = "high_Mito_percent") + 
  facet_wrap(vars(SampleGroup), nrow=1, scales = "free_x") + 
  labs(y = "Percentage mitochondrial UMIs",
       title = "Mitochondrial UMIs") +
  guides(colour=guide_legend(title="Discarded"))


#  venn diagrams below show how the number of discarded droplets in HCA and Caron 
# have changed for each filter in comparison to when the MAD 
# filtering was applied across all samples.
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



# Filtering out poor quality cells
sce.filtered <- sce[, !sce$discard]


# recalculate them QC metrics using filtered data
colData(sce.filtered) <- colData(sce.filtered)[,1:3]
sce.filtered <- addPerCellQC(sce.filtered, BPPARAM = bp.params)


# 6 QC and filtering by combining the metrics.
sce.BM1 <- sce[ , sce$Sample == "PBMMC_1"]


# Using “outlyingness”
library(robustbase)
stats <- cbind(log10(sce.BM1$sum),
               log10(sce.BM1$detected),
               sce.BM1$subsets_Mito_percent)

outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
multi.outlier <- isOutlier(outlying, type = "higher")
summary(multi.outlier)

#################################################################################
# The below sections (PCA and sparsity) is not covered, however, there are materials available
# When you have free time, explore
#################################################################################

