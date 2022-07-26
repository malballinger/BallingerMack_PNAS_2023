#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots expression level variation across genotype and environment (via PCA), for each
# tissue and sex separately. The male-specific plots are used in the main manuscript (Fig1D),
# while female-specific plots are included in the supplement (Figure S2) of BallingerMack_2022.

# Major Result(s): genotype explains the most variation

##############################################################
# Required packages
##############################################################

#rm(list = ls())
library(tidyverse)
library(cowplot)
library(DESeq2)

set.seed(19910118)

source("./code/postprocess_RNAseq/clean_parentalReadCounts.R") # where matrices are generated

##############################################################
# Construct DESeq datasets
##############################################################

## males only
# liver
dds_liver_males <- DESeqDataSetFromMatrix(countData = males_liver_cts_mtrx,
                                          colData = sampleinfo_males_liver,
                                          design = ~ population + temperature)

vsd_liver_males <- vst(dds_liver_males, blind = TRUE) # calculate the across-all-samples variability

# BAT
dds_BAT_males <- DESeqDataSetFromMatrix(countData = males_BAT_cts_mtrx,
                                          colData = sampleinfo_males_BAT,
                                          design = ~ population + temperature)

vsd_BAT_males <- vst(dds_BAT_males, blind = TRUE) # calculate the across-all-samples variability

## females only
# liver
dds_liver_females <- DESeqDataSetFromMatrix(countData = females_liver_cts_mtrx,
                                          colData = sampleinfo_females_liver,
                                          design = ~ population + temperature)

vsd_liver_females <- vst(dds_liver_females, blind = TRUE) # calculate the across-all-samples variability

# BAT
dds_BAT_females <- DESeqDataSetFromMatrix(countData = females_BAT_cts_mtrx,
                                        colData = sampleinfo_females_BAT,
                                        design = ~ population + temperature)

vsd_BAT_females <- vst(dds_BAT_females, blind = TRUE) # calculate the across-all-samples variability


##############################################################
# PCA functions
##############################################################

plotPCA_12 <- function (object, intgroup = c("population", "temperature"),
                        ntop = 500, returnData = FALSE, shape) 
{
  rv <- rowVars(assay(object))
  pick <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[pick, ]))
  percentVar <- signif(as.double(100*(pca$sdev^2/sum(pca$sdev^2))), digits = 3)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                  intgroup.df, name = colData(object)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1,2]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", fill = "interaction(population, temperature)")) +
    geom_point(aes_string(color = "interaction(population, temperature)"), stroke = 0.5, size = 2.5, alpha = 0.9, shape = shape) +
    scale_fill_manual(breaks = c("Brazil.Cold", "Brazil.Warm", "NewYork.Cold", "NewYork.Warm"),
                      values = c("#FFFFFF", "#E09832", "#FFFFFF", "#1683B5")) +
    scale_color_manual(breaks = c("Brazil.Cold", "Brazil.Warm", "NewYork.Cold", "NewYork.Warm"),
                       values = c("#E09832", "#E09832", "#1683B5", "#1683B5"),
                       guides) +
    xlab(paste0("PC1 (",percentVar[1],"%)")) +
    ylab(paste0("PC2 (",percentVar[2],"%)")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = 'none',
          axis.text = element_text (size = 10, family = "sans"),
          axis.title.y = element_text(vjust = 1, size = 10.5, family = "sans"),
          axis.title.x = element_text(vjust = -1, size = 10.5, family = "sans"))
}

plotPCA_13 <- function (object, intgroup = c("population", "temperature"),
                              ntop = 500, returnData = FALSE, shape) 
{
  rv <- rowVars(assay(object))
  pick <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[pick, ]))
  percentVar <- signif(as.double(100*(pca$sdev^2/sum(pca$sdev^2))), digits = 3)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, name = colData(object)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1,3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC3", fill = "interaction(population, temperature)")) +
    geom_point(aes_string(color = "interaction(population, temperature)"), stroke = 1.75, size = 4.75, alpha = 0.9, shape = shape) +
    scale_fill_manual(breaks = c("Brazil.Cold", "Brazil.Warm", "NewYork.Cold", "NewYork.Warm"),
                      values = c("#FFFFFF", "#E09832", "#FFFFFF", "#1683B5")) +
    scale_color_manual(breaks = c("Brazil.Cold", "Brazil.Warm", "NewYork.Cold", "NewYork.Warm"),
                       values = c("#E09832", "#E09832", "#1683B5", "#1683B5"),
                       guides) +
    xlab(paste0("PC1 (",percentVar[1],"%)")) +
    ylab(paste0("PC3 (",percentVar[3],"%)")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = 'none',
          axis.text = element_text (size = 12.5, family = "sans"),
          axis.title.y = element_text(vjust = 1, size = 14, family = "sans"),
          axis.title.x = element_text(vjust = -1, size = 14, family = "sans"))
}

plotPCA_13_supp <- function (object, intgroup = c("population", "temperature"),
                        ntop = 500, returnData = FALSE, shape) 
{
  rv <- rowVars(assay(object))
  pick <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[pick, ]))
  percentVar <- signif(as.double(100*(pca$sdev^2/sum(pca$sdev^2))), digits = 3)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, name = colData(object)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1,3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC3", fill = "interaction(population, temperature)")) +
    geom_point(aes_string(color = "interaction(population, temperature)"), stroke = 0.5, size = 2.5, alpha = 0.9, shape = shape) +
    scale_fill_manual(breaks = c("Brazil.Cold", "Brazil.Warm", "NewYork.Cold", "NewYork.Warm"),
                      values = c("#FFFFFF", "#E09832", "#FFFFFF", "#1683B5")) +
    scale_color_manual(breaks = c("Brazil.Cold", "Brazil.Warm", "NewYork.Cold", "NewYork.Warm"),
                       values = c("#E09832", "#E09832", "#1683B5", "#1683B5"),
                       guides) +
    xlab(paste0("PC1 (",percentVar[1],"%)")) +
    ylab(paste0("PC3 (",percentVar[3],"%)")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = 'none',
          axis.text = element_text (size = 10, family = "sans"),
          axis.title.y = element_text(vjust = 1, size = 10.5, family = "sans"),
          axis.title.x = element_text(vjust = -1, size = 10.5, family = "sans"))
}

### execute function and save resulting plots

## males
# liver
males_liver_PCA_12 <- plotPCA_12(object = vsd_liver_males, shape = 24)
#ggsave("results/figures/males_liver_PCA12.pdf", plot = males_liver_PCA_12, height = 1.75, width = 2, units = "in")

males_liver_PCA_13 <- plotPCA_13(object = vsd_liver_males, shape = 21)
ggsave("results/figures/males_liver_PCA13.pdf", plot = males_liver_PCA_13, height = 3.25, width = 3.55, units = "in") #1.75 and 2

# BAT
males_BAT_PCA_12 <- plotPCA_12(object = vsd_BAT_males, shape = 21)
#ggsave("results/figures/males_BAT_PCA12.pdf", plot = males_BAT_PCA_12, height = 1.75, width = 2) # 1.75 by 2

males_BAT_PCA_13 <- plotPCA_13(object = vsd_BAT_males, shape = 21)
ggsave("results/figures/males_BAT_PCA13.pdf", plot = males_BAT_PCA_13, height = 3.25, width = 3.55) # 1.75 by 2


## females
# liver
females_liver_PCA_12 <- plotPCA_12(object = vsd_liver_females, shape = 24)
#ggsave("results/figures/females_liver_PCA12.pdf", plot = females_liver_PCA_12, height = 1.75, width = 2, units = "in")

females_liver_PCA_13 <- plotPCA_13_supp(object = vsd_liver_females, shape = 24)
#ggsave("results/figures/females_liver_PCA13.pdf", plot = females_liver_PCA_13, height = 1.75, width = 2, units = "in") #1.75 and 2

# BAT
females_BAT_PCA_12 <- plotPCA_12(object = vsd_BAT_females, shape = 21)
#ggsave("results/figures/females_BAT_PCA12.pdf", plot = females_BAT_PCA_12, height = 1.75, width = 2) # 1.75 by 2

females_BAT_PCA_13 <- plotPCA_13_supp(object = vsd_BAT_females, shape = 21)
#ggsave("results/figures/females_BAT_PCA13.pdf", plot = females_BAT_PCA_13, height = 1.75, width = 2) # 1.75 by 2

## plot figure S2
figS2 <- cowplot::plot_grid(males_BAT_PCA_12, males_liver_PCA_12, females_BAT_PCA_12, females_liver_PCA_12, females_BAT_PCA_13, females_liver_PCA_13,
                   nrow = 3, ncol = 2)
ggsave("results/figures/FigS2_PCA.pdf", plot = figS2, height = 5.5, width = 4)
