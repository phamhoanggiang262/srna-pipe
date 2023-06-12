#!/usr/bin/env Rscript

library(optparse)
library(readr)
library(tidyverse)

option_list <- list(
  make_option(
    c("-s", "--samplesheet"),
    type = "character",
    default = NULL,
    help = " CSV-format sample sheet file."),
  make_option(
    c("-p", "--path"),
    type = "character",
    default = NULL,
    help = " CSV-format sample sheet file.")
)

opt_parser <- OptionParser(option_list = option_list);
opt <- parse_args(opt_parser);

# read sample sheet file
samplesheet <- read.csv(opt$samplesheet)

#rename column name
samplesheet <- rename(samplesheet, FileName = fastq)

#split file name from path
samplesheet$FileName <- basename(samplesheet$FileName)

#set to absolute path
samplesheet$FileName <- file.path(opt$path, paste0(samplesheet$SampleName, ".bam"))


#move FileName in the front
samplesheet %>%  relocate(FileName, .before = SampleName) -> samplesheet

write.table(samplesheet, file="sampleInfo.csv", sep=",", row.names = F, col.names = T)


