#!/bin/bash

set -eou pipefail

# Render the index page
Rscript -e "rmarkdown::render('Markdowns/index.md', output_format = 'html_document', output_dir = 'docs')"

# Render all Markdown files as HTML documents
for i in Markdowns/03_CellRanger.Rmd \
         Markdowns/04_Preprocessing_And_QC.Rmd \
         Markdowns/04_Preprocessing_And_QC.Exercise.Rmd \
         Markdowns/05_Normalisation.Rmd \
         Markdowns/05_Normalisation_exercises.Rmd \
         Markdowns/06_FeatureSelectionAndDimensionalityReduction.Rmd \
         Markdowns/07_Dataset_Integration.Rmd \
         Markdowns/08_Clustering.Rmd \
         Markdowns/09_Cluster_Marker_Genes.Rmd \
         Markdowns/10_Differential_Expression.Rmd \
         Markdowns/11_Differential_Abundance.Rmd
do
  echo "Rendering $i as HTML document..."
  Rscript -e "rmarkdown::render('$i', output_format = 'html_document', output_dir = 'docs')"
done

# Render all Slides files as ioslides presentations
for i in Slides/00_Day1_Recap.Rmd \
         Slides/02_PreambleSlides.Rmd \
         Slides/03_CellRangerSlides.Rmd \
         Slides/04_QualityControlSlides.Rmd \
         Slides/05_NormalisationSlides.Rmd \
         Slides/06_FeatureSelectionAndDimensionalityReduction_slides.Rmd \
         Slides/07_DataIntegrationAndBatchCorrectionSlides.Rmd \
         Slides/08_ClusteringSlides.Rmd \
         Slides/09_ClusterMarkerGenes.Rmd \
         Slides/10_DifferentialAnalysis.Rmd
do
  echo "Rendering $i as ioslides presentation..."
  Rscript -e "rmarkdown::render('$i', output_format = 'ioslides_presentation', output_dir = 'docs')"
done

echo "All files rendered successfully!"
