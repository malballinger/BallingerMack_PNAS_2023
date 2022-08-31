#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script examines sample-to-sample distances with hierarchical clustering,
# and generates Figure S14 in BallingerMack_2022.

# Main Result: unlike liver tissue, BAT samples do not cluster by sex

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clears R's "memory"
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

set.seed(19910118)

source("./code/postprocess_RNAseq/clean_parentalReadCounts.R") # where matrices are generated

##############################################################
# Construct and plot distance heatmaps
##############################################################

# set color scheme
colors <- colorRampPalette(rev(brewer.pal(9, "Greys")))(200) # this (number) extends the number of colors from 9 to (number)

# rename sampleIDs (from "XXXX.counts" to "new-names")
tissue_diffs_rows <- c("Brazil-Male-Warm", "Brazil-Male-Warm", "Brazil-Male-Warm", "Brazil-Male-Warm", "Brazil-Male-Warm", "Brazil-Male-Warm",
                       "Brazil-Female-Warm", "Brazil-Female-Warm", "Brazil-Female-Warm", "Brazil-Female-Warm", "Brazil-Female-Warm", "Brazil-Female-Warm",
                       "Brazil-Male-Cold", "Brazil-Male-Cold", "Brazil-Male-Cold", "Brazil-Male-Cold", "Brazil-Male-Cold", "Brazil-Male-Cold",
                       "Brazil-Female-Cold", "Brazil-Female-Cold", "Brazil-Female-Cold", "Brazil-Female-Cold", "Brazil-Female-Cold", "Brazil-Female-Cold",
                       "NewYork-Male-Warm", "NewYork-Male-Warm", "NewYork-Male-Warm", "NewYork-Male-Warm", "NewYork-Male-Warm", "NewYork-Male-Warm",
                       "NewYork-Female-Warm", "NewYork-Female-Warm", "NewYork-Female-Warm", "NewYork-Female-Warm", "NewYork-Female-Warm", "NewYork-Female-Warm",
                       "NewYork-Male-Cold", "NewYork-Male-Cold", "NewYork-Male-Cold", "NewYork-Male-Cold", "NewYork-Male-Cold", "NewYork-Male-Cold",
                       "NewYork-Female-Cold", "NewYork-Female-Cold", "NewYork-Female-Cold", "NewYork-Female-Cold", "NewYork-Female-Cold", "NewYork-Female-Cold")

# construct DESeq datasets
dds_liver <- DESeqDataSetFromMatrix(countData = all_liver_cts_mtrx,
                                    colData = all_liver_sampleinfo,
                                    design = ~ population + sex + temperature)

dds_BAT <- DESeqDataSetFromMatrix(countData = all_BAT_cts_mtrx,
                                  colData = all_BAT_sampleinfo,
                                  design = ~ population + sex + temperature)
# create pheatmap function
pheatmapping <- function(dds)
{
  vsd <- vst(dds, blind = TRUE) # calculate the across-all-samples variability
  
  sampleDists <- dist(t(assay(vsd))) # calculate distances
  
  sampleDistMatrix <- as.matrix(sampleDists) # create distance matrix
  colnames(sampleDistMatrix) <- NULL
  
  pheatmap(sampleDistMatrix, # construct heatmap
           clustering_distance_rows=sampleDists,
           labels_row=tissue_diffs_rows,
           clustering_distance_cols=sampleDists,
           col=colors, scale = "none", border_color = NA, fontsize = 14)
}

# plot
all_sex_temp_liver <- pheatmapping(dds = dds_liver)
ggsave("results/figures/liver_pheatmap.pdf", plot = all_sex_temp_liver,
       height = 9, width = 7.5)

all_sex_temp_BAT <- pheatmapping(dds = dds_BAT)
ggsave("results/figures/BAT_pheatmap.pdf", plot = all_sex_temp_BAT,
       height = 9, width = 7.5)


