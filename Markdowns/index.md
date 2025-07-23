---
title: "Single-cell RNA-seq Analysis"
output:
  html_document:
    css: ../css/boxes.css
    includes:
      in_header: ../css/navbar.html
---

## Outline
This workshop is designed for biologists seeking to learn standard single-cell RNA-seq analysis techniques. 
The course emphasises droplet-based assays using 10X Genomics technology and covers the use of the Cell Ranger pipeline for read alignment and gene-level quantification. 
Participants will learn how to import count data into R, perform quality control, normalisation, data integration, clustering, and identify cluster marker genes. 
The workshop also addresses differential expression and abundance analyses. 
Additionally, attendees will gain experience in generating common visualisations, including t-SNE, UMAP, and violin plots, to interpret gene expression data.

## Topics

(Note: session times are approximate)

| Time|Topic|Links|
|----:|-----|:----|
|15 min|**Welcome!**||
|45 min|**Introduction to single-cell technologies** |[Slides]()|
|15 min|**Preamble: data set and workflow** |[Slides](02_PreambleSlides.html)|
|1 h 30 min|**Library structure, cellranger for alignment and cell calling** |[Slides](03_CellRangerSlides.html)<br>[Demonstration](03_CellRanger.html)|
|2 h|**QC and exploratory analysis** |[Slides](04_QualityControlSlides.html)<br>[Demonstration](04_Preprocessing_And_QC.html)<br>[Exercises](04_Preprocessing_And_QC.Exercise.html)|
|2 h|**Normalisation** |[Slides](05_NormalisationSlides.html)<br>[Demonstration](05_Normalisation.html)<br>[Exercises](05_Normalisation_exercises.html)|
|2 h|**Feature selection and dimensionality reduction**|[Slides](06_FeatureSelectionAndDimensionalityReduction_slides.html)<br>[Practical](06_FeatureSelectionAndDimensionalityReduction.html)|
|2 h|**Batch correction and data set integration**|[Slides](07_DataIntegrationAndBatchCorrectionSlides)<br>[Demonstration](07_Dataset_Integration.html)|
|1 h 30 min|**Cell clustering**|[Slides](08_ClusteringSlides.html)<br>[Demonstration](08_Clustering.html)|
|1 h|**Identification of cluster marker genes**|[Slides](09_ClusterMarkerGenes.html)<br>[Demonstration](09_Cluster_Marker_Genes.html)|
|1 h 30 min|**Differential expression analysis**|[Demonstration](10_Differential_Expression.html)|
|1 h 30 min|**Differential abundance analysis**|[Demonstration](11_Differential_Abundance.html)|

<!-- https://stackoverflow.com/a/58338258 -->
<style>
table th:first-of-type {
    width: 10%;
}
table th:nth-of-type(2) {
    width: 70%;
}
table th:nth-of-type(3) {
    width: 20%;
}
</style>

## Data

* The course data is based on '[CaronBourque2020](https://www.nature.com/articles/s41598-020-64929-x)'
  relating to pediatric leukemia, with four sample types, including:
  * pediatric Bone Marrow Mononuclear Cells (PBMMCs)
  * three tumour types: ETV6-RUNX1, HHD, PRE-T  
* The data used in the course can be [downloaded from Dropbox](https://www.dropbox.com/scl/fo/9x1rg6qxqw5crq2vtb1ho/AMawguf1kqRYQQs-qZPhFZA?rlkey=6y4w1skyjzpq36zfyocis24t6&st=yyv0kusl&dl=0). Please note that:
  * these data have been processed for teaching purposes and are therefore not suitable for research use;
  * all the data is provided on our training machines, you don't need to download it to attend the course.


## Additional Resources

* [Modern Statistics for Modern Biology](https://www.huber.embl.de/msmb/08-chap.html) - a fantastic resource to go deeper into some of the topics covered in this course. In particular chapter 8 of this book covers the statistical models for count data used in RNA-seq.
* [Orchestrating Single-Cell Analysis with Bioconductor](https://bioconductor.org/books/release/OSCA/)
* [Single-cell best practices](https://www.sc-best-practices.org/preamble.html)


## Acknowledgements

Much of the material in this course has been derived from the demonstrations found in
[OSCA book](https://bioconductor.org/books/release/OSCA/)
and the [Hemberg Group course materials](https://www.singlecellcourse.org/). Additional material concerning `miloR` has been based on the [demonstration from the Marioni Lab](https://marionilab.github.io/miloR/articles/milo_demo.html).

The materials have been contributed to by many individuals over the last few years, including:

Abigail Edwards, Ashley D Sawle, Chandra Chilamakuri, Kamal Kishore, Stephane Ballereau, Zeynep Kalendar Atak, Hugo Tavares, Jon Price, Katarzyna Kania, Roderik Kortlever, Adam Reid, Tom Smith

Apologies if we have missed anyone!