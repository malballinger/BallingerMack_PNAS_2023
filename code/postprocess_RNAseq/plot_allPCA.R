#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots expression level variation across tissue-type, genotype, sex, and environment (via PCA),
# and generates Figure S1 in BallingerMack_2022.

# Major Result(s): tissue-type explains the most variation

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clears R's "memory"
library(tidyverse)
library(DESeq2)
library(cowplot)

set.seed(19910118)

source("./code/postprocess_RNAseq/clean_parentalReadCounts.R") # where matrices are generated

##############################################################
# Construct DESeq datasets
##############################################################

## all individuals

dds_all <- DESeqDataSetFromMatrix(countData = all_cts_mtrx,
                                  colData = SampleInfo,
                                  design = ~ population + sex + temperature + condition)

vsd_all <- vst(dds_all, blind = TRUE) # calculate the across-all-samples variability

## separate out liver and BAT

# liver-only
dds_liver <- DESeqDataSetFromMatrix(countData = all_liver_cts_mtrx,
                                    colData = all_liver_sampleinfo,
                                    design = ~ population + sex + temperature)

vsd_liver <- vst(dds_liver, blind = TRUE) # calculate the across-all-samples variability

# BAT-only
dds_BAT <- DESeqDataSetFromMatrix(countData = all_BAT_cts_mtrx,
                                  colData = all_BAT_sampleinfo,
                                  design = ~ population + sex + temperature)

vsd_BAT <- vst(dds_BAT, blind = TRUE) # calculate the across-all-samples variability

##############################################################
# PCA functions
##############################################################

plotPCA_all_12 <- function (object, intgroup = c("population", "sex", "temperature"),
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
  ggplot(data = d, aes_string(x = "PC1", y = "PC2", size = "sex", fill = "interaction(population, temperature)")) +
    geom_point(aes_string(color = "interaction(population, temperature)"), stroke = 0.8, shape = shape, alpha = 0.9) +
    scale_size_manual(breaks = c("M", "F"),
                      values = c(4,1.5)) +
    scale_fill_manual(breaks = c("Brazil.Cold", "Brazil.Warm", "NewYork.Cold", "NewYork.Warm"),
                      values = c("#FFFFFF", "#E09832", "#FFFFFF", "#1683B5")) +
    scale_color_manual(breaks = c("Brazil.Cold", "Brazil.Warm", "NewYork.Cold", "NewYork.Warm"),
                       values = c("#E09832", "#E09832", "#1683B5", "#1683B5")) +
    xlab(paste0("PC1 (",percentVar[1],"%)")) +
    ylab(paste0("PC2 (",percentVar[2],"%)")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = 'none',
          axis.text = element_text(family = "sans"),
          axis.title = element_text(size = 9, family = "sans"))
}

plotPCA_all_13 <- function (object, intgroup = c("population", "sex", "temperature"),
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
  ggplot(data = d, aes_string(x = "PC1", y = "PC3", size = "sex", fill = "interaction(population, temperature)")) +
    geom_point(aes_string(color = "interaction(population, temperature)"), stroke = 0.8, shape = shape, alpha = 0.9) +
    scale_size_manual(breaks = c("M", "F"),
                      values = c(4,1.5)) +
    scale_fill_manual(breaks = c("Brazil.Cold", "Brazil.Warm", "NewYork.Cold", "NewYork.Warm"),
                      values = c("#FFFFFF", "#E09832", "#FFFFFF", "#1683B5")) +
    scale_color_manual(breaks = c("Brazil.Cold", "Brazil.Warm", "NewYork.Cold", "NewYork.Warm"),
                       values = c("#E09832", "#E09832", "#1683B5", "#1683B5")) +
    xlab(paste0("PC1 (",percentVar[1],"%)")) +
    ylab(paste0("PC3 (",percentVar[3],"%)")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = 'none',
          axis.text = element_text(family = "sans"),
          axis.title = element_text(size = 9, family = "sans"))
}

##############################################################
# Make all sub plots
##############################################################

# FigS1A

pcaData_all <- plotPCA(vsd_all, intgroup = c("population", "sex", "temperature", "condition"), returnData = TRUE)
percentVar_all <- signif(as.double(100*attr(pcaData_all, "percentVar")), digits = 3)

all_plot_12 <-
  ggplot(pcaData_all, aes(PC1, PC2, shape = condition, size = sex,
                          fill = interaction(population,temperature), color = interaction(population,temperature))) +
  geom_point(stroke = 0.8, alpha = 0.9) +
  scale_size_manual(breaks = c("M", "F"),
                    values = c(5.5,1.5),
                    labels = c("Male", "Female")) +#,
                    #guide = guide_legend(title = "sex",
                                         #override.aes = list(color = "black",
                                                             #shape = 0,
                                                             #size = c(6,1)))) +
  scale_shape_manual(breaks = c("BAT", "LIVER"),
                     values = c(21, 24)) +#,
                     #guide = guide_legend(title = "tissue",
                                          #override.aes = list(color = "black",
                                                              #stroke = 1,
                                                              #shape = c(16,17)))) +
  scale_fill_manual(breaks = c("Brazil.Cold", "Brazil.Warm", "NewYork.Cold", "NewYork.Warm"),
                    values = c("#FFFFFF", "#E09832", "#FFFFFF", "#1683B5")) +#,
                    #labels = c("Brazil-Cold", "Brazil-Warm", "NewYork-Cold", "NewYork-Warm"),
                    #guide = guide_legend(title = "population - environment",
                                         #override.aes = list(color = c("#E09832", "#875B1E", "#388BFF", "#193D70"),
                                                             #shape = 3))) +
  scale_color_manual(breaks = c("Brazil.Cold", "Brazil.Warm", "NewYork.Cold", "NewYork.Warm"),
                     values = c("#E09832", "#E09832", "#1683B5", "#1683B5")) +
  xlab(paste0("PC1 (",percentVar_all[1],"%)")) +
  ylab(paste0("PC2 (",percentVar_all[2],"%)")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = 'none',
        #legend.position = c(0.5,0.5),
        #legend.background = element_blank(),
        #legend.box.background = element_rect(color="black"),
        #legend.key = element_blank(),
        #legend.key.height = unit(0.4, "cm"),
        #legend.key.width = unit(0.6, "cm"),
        #legend.spacing.x = unit(0.5, "mm"),
        #legend.spacing.y = unit(0.5, "mm"),
        #legend.title = element_text(size = 8, face = "italic"),
        #legend.text = element_text(size = 7),
        axis.text = element_text(family = "sans"),
        axis.title = element_text(family = "sans"))


# FigS1B
all_liver_PCA_12 <- plotPCA_all_12(object = vsd_liver, shape = 24)
all_BAT_PCA_12 <- plotPCA_all_12(object = vsd_BAT, shape = 21)

# FigS1C
all_liver_PCA_13 <- plotPCA_all_13(object = vsd_liver, shape = 24)
all_BAT_PCA_13 <- plotPCA_all_13(object = vsd_BAT, shape = 21)


##############################################################
# Plot Figure S1
##############################################################
# combine liver and BAT plots
all_liver_BAT_123 <- cowplot::plot_grid(all_liver_PCA_12, all_BAT_PCA_12, all_liver_PCA_13, all_BAT_PCA_13,
                                        nrow = 2, ncol = 2)
# combine all plots together
figS1 <- cowplot::plot_grid(all_plot_12, all_liver_BAT_123, nrow = 2)
ggsave("results/figures/FigS1_allPCA.pdf", plot = figS1, height = 7, width = 5)
