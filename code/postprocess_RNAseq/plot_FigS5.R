#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script filters and subsets counts metadata files and
# generates matrices that are used to generates Figure S5.

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clear R's environment
library(tidyverse)
library(here)

##############################################################
# Import data
##############################################################

F1s_Counts_liver1 <- read_table(here("./data/raw/ReadCounts/F1_liver_warm_malefemale.txt"), col_names = TRUE, col_types = NULL)
F1s_Counts_liver2 <- read_table(here("./data/raw/ReadCounts/F1_liver_cold_malefemale.txt"), col_names = TRUE, col_types = NULL)
F1s_Counts_bat <- read_table(here("./data/raw/ReadCounts/F1_BAT_allreads_warmcold_malefemale.txt"), col_names = TRUE, col_types = NULL)
Parents_Counts_Metadata <- read_table(here("./data/raw/ReadCounts/all_parents_counts.txt"), col_names = TRUE, col_types = NULL)

Counts_Metadata <- full_join(F1s_Counts_liver1, F1s_Counts_liver2, by = "gene") %>%
                  full_join(., F1s_Counts_bat, by = "gene") %>%
                  full_join(., Parents_Counts_Metadata, by = "gene")

males_SampleInfo <- read_delim(here("./data/raw/ReadCounts/males_parents_F1s_sample_info.txt"))
colnames(males_SampleInfo)[1] <- "sampleID"

females_SampleInfo <- read_delim(here("./data/raw/ReadCounts/females_parents_F1s_sample_info.txt"))
colnames(females_SampleInfo)[1] <- "sampleID"

##############################################################
# Subset data
##############################################################

## males-liver only
males_liver_counts <- Counts_Metadata %>%
  dplyr::select(gene, starts_with(c("079", "082", "086", "090", "097", "098", "078", "083", "088", "093", "095", "099",
                                    "MBSD_002", "MBSD_003", "MBSD_009", "MBSD_012", "MBSD_013", "MBSD_022",
                                    "MBSD_005", "MBSD_007", "MBSD_011", "MBSD_015", "MBSD_020", "MBSD_024",
                                    "MBSD_055", "MBSD_061", "MBSD_065", "MBSD_072", "MBSD_075", "MBSD_077",
                                    "MBSD_057", "MBSD_059", "MBSD_063", "MBSD_067", "MBSD_069", "MBSD_074")))
                
males_liver_sampleinfo <- males_SampleInfo %>%
  filter(condition == "LIVER")

## males-BAT only
males_bat_counts <- Counts_Metadata %>%
  dplyr::select(gene, starts_with(c("177", "182", "187", "192", "194", "198", "196", "197", "189", "185", "178", "181",
                                    "MBSD_101", "MBSD_102", "MBSD_108", "MBSD_111", "MBSD_112", "MBSD_121",
                                    "MBSD_104", "MBSD_106", "MBSD_110", "MBSD_114", "MBSD_119", "MBSD_123",
                                    "MBSD_154", "MBSD_160", "MBSD_164", "MBSD_171", "MBSD_174", "MBSD_176",
                                    "MBSD_156", "MBSD_158", "MBSD_162", "MBSD_166", "MBSD_168", "MBSD_173")))

males_bat_sampleinfo <- males_SampleInfo %>%
  filter(condition == "BAT")



## females-liver only
females_liver_counts <- Counts_Metadata %>%
  dplyr::select(gene, starts_with(c("080", "089", "091", "096", "081", "084", "085", "087", "092", "094",
                                    "MBSD_004", "MBSD_006", "MBSD_008", "MBSD_010", "MBSD_017", "MBSD_019",
                                    "MBSD_001", "MBSD_014", "MBSD_016", "MBSD_018", "MBSD_021", "MBSD_023",
                                    "MBSD_054", "MBSD_056", "MBSD_062", "MBSD_068", "MBSD_073", "MBSD_076",
                                    "MBSD_058", "MBSD_060", "MBSD_064", "MBSD_066", "MBSD_070", "MBSD_071")))

females_liver_sampleinfo <- females_SampleInfo %>%
  filter(condition == "LIVER")

## females-BAT only
females_bat_counts <- Counts_Metadata %>%
  dplyr::select(gene, starts_with(c("191", "193", "180", "183", "184", "186", "179", "188", "190", "195",
                                    "MBSD_103", "MBSD_105", "MBSD_107", "MBSD_109", "MBSD_116", "MBSD_118",
                                    "MBSD_100", "MBSD_113", "MBSD_115", "MBSD_117", "MBSD_120", "MBSD_122",
                                    "MBSD_153", "MBSD_155", "MBSD_161", "MBSD_167", "MBSD_172", "MBSD_175",
                                    "MBSD_157", "MBSD_159", "MBSD_163", "MBSD_165", "MBSD_169", "MBSD_170")))

females_bat_sampleinfo <- females_SampleInfo %>%
  filter(condition == "BAT")


##############################################################
# Convert tibbles into matrices
##############################################################

tibb_to_mtrx <- function(tibb)
{
  
  genes <- pull(tibb, gene) # extract the column of all the genes
  matrix <- as.matrix(tibb[,-1]) # create matrix of tibble by omitting first column (which has the genes)
  rownames(matrix) <- genes # put genes back into the now matrix (but the genes column does not have a header)
  
  return(matrix)
  
}

males_liver_cts_mtrx <- tibb_to_mtrx(tibb = males_liver_counts) # males_liver_sampleinfo
males_BAT_cts_mtrx <- tibb_to_mtrx(tibb = males_bat_counts) # males_bat_sampleinfo
females_liver_cts_mtrx <- tibb_to_mtrx(tibb = females_liver_counts) # females_liver_sampleinfo
females_BAT_cts_mtrx <- tibb_to_mtrx(tibb = females_bat_counts) # females_bat_sampleinfo

##############################################################
# Construct DESeq datasets
##############################################################

# males-liver-only
dds_liver_males <- DESeqDataSetFromMatrix(countData = males_liver_cts_mtrx,
                                    colData = males_liver_sampleinfo,
                                    design = ~ population + temperature)

vsd_liver_males <- vst(dds_liver_males, blind = TRUE) # calculate the across-all-samples variability

# males-BAT-only
dds_BAT_males <- DESeqDataSetFromMatrix(countData = males_BAT_cts_mtrx,
                                  colData = males_bat_sampleinfo,
                                  design = ~ population + temperature)

vsd_BAT_males <- vst(dds_BAT_males, blind = TRUE) # calculate the across-all-samples variability



# females-liver-only
dds_liver_females <- DESeqDataSetFromMatrix(countData = females_liver_cts_mtrx,
                                          colData = females_liver_sampleinfo,
                                          design = ~ population + temperature)

vsd_liver_females <- vst(dds_liver_females, blind = TRUE) # calculate the across-all-samples variability

# females-BAT-only
dds_BAT_females <- DESeqDataSetFromMatrix(countData = females_BAT_cts_mtrx,
                                        colData = females_bat_sampleinfo,
                                        design = ~ population + temperature)

vsd_BAT_females <- vst(dds_BAT_females, blind = TRUE) # calculate the across-all-samples variability


##############################################################
# PCA functions
##############################################################

plotPCA_all_12 <- function (object, intgroup = c("population", "temperature"),
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
    geom_point(aes_string(color = "interaction(population, temperature)"), stroke = 0.8, shape = shape, alpha = 0.9) +
    scale_size_manual(breaks = c("M", "F"),
                      values = c(4,1.5)) +
    scale_fill_manual(breaks = c("Brazil.Cold", "Brazil.Warm", "NewYork.Cold", "NewYork.Warm", "F1.Cold", "F1.Warm"),
                      values = c("#FFFFFF", "#E09832", "#FFFFFF", "#1683B5", "#FFFFFF", "darkgray")) +
    scale_color_manual(breaks = c("Brazil.Cold", "Brazil.Warm", "NewYork.Cold", "NewYork.Warm", "F1.Cold", "F1.Warm"),
                       values = c("#E09832", "#E09832", "#1683B5", "#1683B5", "darkgray", "darkgray")) +
    xlab(paste0("PC1 (",percentVar[1],"%)")) +
    ylab(paste0("PC2 (",percentVar[2],"%)")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = 'none',
          axis.text = element_text(family = "sans"),
          axis.title = element_text(size = 9, family = "sans"))
}

males_liver_PCA_12 <- plotPCA_all_12(object = vsd_liver_males, shape = 24)
males_bat_PCA_12 <- plotPCA_all_12(object = vsd_BAT_males, shape = 21)

females_liver_PCA_12 <- plotPCA_all_12(object = vsd_liver_females, shape = 24)
females_bat_PCA_12 <- plotPCA_all_12(object = vsd_BAT_females, shape = 21)

# combine all plots together
figS5 <- cowplot::plot_grid(males_bat_PCA_12, males_liver_PCA_12,
                            females_bat_PCA_12, females_liver_PCA_12,
                            nrow = 2, ncol = 2)
ggsave(here("./results/figures/FigS5_PCA.pdf"), plot = figS5, height = 4, width = 4)
