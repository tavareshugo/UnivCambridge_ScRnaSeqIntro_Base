#!/usr/bin/env Rscript

library(tidyverse)

read_csv("data/SraRunTable.txt") %>%
    select(SampleID = Run, SampleGroup = source_name) %>%
    arrange(SampleID) %>%
    mutate(BioSample = SampleID) %>%
    mutate(BioSample = ifelse(BioSample=="SRR9264352", 
                              "SRR9264351", 
                              BioSample)) %>%
    group_by(SampleGroup) %>%
    mutate(Rep = as.numeric(as.factor(BioSample))) %>%
    ungroup() %>%
    mutate(BioSample = str_c(SampleGroup, "_", Rep)) %>%
    mutate(SampleName = case_when(
           SampleID == "SRR9264351" ~ str_c(BioSample, "a"),
           SampleID == "SRR9264352" ~ str_c(BioSample, "b"),
           TRUE ~ BioSample)) %>%
    select(SampleID, BioSample, SampleName, SampleGroup) %>%
    write_tsv("data/sample_sheet.tsv")
