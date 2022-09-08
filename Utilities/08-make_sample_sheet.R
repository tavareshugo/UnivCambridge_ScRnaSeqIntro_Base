#!/usr/bin/env Rscript

library(tidyverse)

# We are going to omit the technical replicate for PBMMC_1 - SRR9264352

read_csv("data/SraRunTable.txt") %>%
    select(Sample = Run, SampleGroup = source_name) %>%
    filter(Sample != "SRR9264352") %>%
    arrange(Sample) %>%
    group_by(SampleGroup) %>%
    mutate(Rep = as.numeric(as.factor(Sample))) %>%
    ungroup() %>%
    mutate(SampleName = str_c(SampleGroup, "_", Rep)) %>%
    select(Sample, SampleName, SampleGroup) %>%
    write_tsv("data/sample_sheet.tsv")
