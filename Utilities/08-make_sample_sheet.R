#!/usr/bin/env Rscript

library(tidyverse)

read_csv("data/SraRunTable.txt") %>%
    select(Sample = Run, SampleGroup = source_name) %>%
    arrange(Sample) %>%
    mutate(BioSample = Sample) %>%
    mutate(BioSample = ifelse(BioSample=="SRR9264352", 
                              "SRR9264351", 
                              BioSample)) %>%
    group_by(SampleGroup) %>%
    mutate(Rep = as.numeric(as.factor(BioSample))) %>%
    ungroup() %>%
    mutate(BioSample = str_c(SampleGroup, "_", Rep)) %>%
    mutate(SampleName = case_when(
           Sample == "SRR9264351" ~ str_c(BioSample, "a"),
           Sample == "SRR9264352" ~ str_c(BioSample, "b"),
           TRUE ~ BioSample)) %>%
    select(Sample, BioSample, SampleName, SampleGroup) %>%
    write_tsv("data/sample_sheet.tsv")
