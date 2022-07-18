#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots differential expression in parental samples, with specific focus on patterns of GxS.
# This script generates Figure S4C in BallingerMack_2022.

# Main Result(s): BAT harbors *very little* GxS; NY shows more sex-specific expression than BZ in liver 

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clear R's environment
library(tidyverse)
library(ggtext)
library(DESeq2)

set.seed(19910118)

source("./code/postprocess_RNAseq/get_DE_parents_analysis.R") # where DESeq data frames are generated

##############################################################
# Generate GxS scatter plots
##############################################################

# create function
GxS <- function(resNY, resBZ, resGxS)
{

  # filter to get mean of 10 reads across all samples
  DE_base_resNY <- resNY %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    filter(baseMean >= 10)
  
  DE_base_resBZ <- resBZ %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    filter(baseMean >= 10)
  
  # merge datasets
  merge_resBZNY <- dplyr::full_join(DE_base_resBZ, DE_base_resNY, by = "gene") %>%
    dplyr::select(
      "gene",
      "lfc.BZ" = "log2FoldChange.x",
      "lfcSE.BZ" = "lfcSE.x",
      "padj.BZ" = "padj.x",
      "lfc.NY" = "log2FoldChange.y",
      "lfcSE.NY" = "lfcSE.y",
      "padj.NY" = "padj.y",
      )
  
  n_total_genes <- nrow(merge_resBZNY)
  n_sig_DE_NY_genes <- merge_resBZNY %>% filter(padj.NY < 0.05) %>% nrow()
  n_sig_DE_BZ_genes <- merge_resBZNY %>% filter(padj.BZ < 0.05) %>% nrow()
  
  # add a column indicating whether BZ and NY groups are in the same direction
  merge_resBZNY_direction <-  merge_resBZNY %>%
    mutate(direction = case_when((sign(lfc.BZ)) == (sign(lfc.NY)) ~ 1, # same direction
                                 (sign(lfc.BZ)) != (sign(lfc.NY)) ~ 0)) # opposite direction
  
  # shrink plotting window by assigning large LFC to 5's
  merge_resBZNY_direction$lfc.NY[merge_resBZNY_direction$lfc.NY>5]<-5
  merge_resBZNY_direction$lfc.BZ[merge_resBZNY_direction$lfc.BZ>5]<-5
  merge_resBZNY_direction$lfc.NY[merge_resBZNY_direction$lfc.NY<(-5)]<-(-5)
  merge_resBZNY_direction$lfc.BZ[merge_resBZNY_direction$lfc.BZ<(-5)]<-(-5)
  
  # categorize genes by direction and/or significance
  sigNY = dplyr::filter(merge_resBZNY_direction, padj.NY < 0.05 & padj.BZ >= 0.05)
  sigBZ = dplyr::filter(merge_resBZNY_direction, padj.BZ < 0.05 & padj.NY >= 0.05)
  sigsame = dplyr::filter(merge_resBZNY_direction, direction == 1 & padj.NY < 0.05 & padj.BZ < 0.05)
  sigopps = dplyr::filter(merge_resBZNY_direction, direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05)
  nonsig = dplyr::filter(merge_resBZNY_direction, padj.NY >= 0.05 & padj.BZ >= 0.05)
  
  # collate all genes that show evidence of GxS
  # first, genes identified as GxS by DESEQ (very stringent)
  DESEQ_GxS <- resGxS %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    filter(baseMean >= 10, padj < 0.05)
  # save a list of GxS genes defined by DESEQ
  DESEQ_GxS_list <- DESEQ_GxS %>% pull(gene)
  
  # second, genes I identify as GxS from categorization above + DESEQ genes
  GxS <- merge_resBZNY_direction %>%
    filter((padj.NY < 0.05 & padj.BZ >= 0.05) |
             (padj.BZ < 0.05 & padj.NY >= 0.05) |
             (direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) |
             (gene %in% DESEQ_GxS_list))
  
  # save a list of all GxS genes
  GxS_list <- GxS %>% pull(gene)
  
  # sample sizes of categories:
  n_sigNY <- nrow(sigNY)
  n_sigsame <- nrow(sigsame)
  n_sigBZ <- nrow(sigBZ)
  n_sigopps <- nrow(sigopps)
  n_GxS <- nrow(GxS)
  
  # plot
  plot <- merge_resBZNY_direction %>%
    ggplot(aes(x = lfc.NY, y = lfc.BZ)) +
    geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
    geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
    geom_abline(slope = 1, linetype = "dotted", size = 0.25) +
    geom_point(data = nonsig, aes(x = lfc.NY, y = lfc.BZ), fill = "lightgray", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
    geom_point(data = sigNY, aes(x = lfc.NY, y = lfc.BZ), fill =  "#1683B5", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
    geom_point(data = sigBZ, aes(x = lfc.NY, y = lfc.BZ), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
    geom_point(data = sigopps, aes(x = lfc.NY, y = lfc.BZ), fill = "#000000", color = "000000", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
    geom_point(data = sigsame, aes(x = lfc.NY, y = lfc.BZ), fill = "#B22222", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
    scale_x_continuous(limits = c(-5,5), expand = c(0.01,0.075)) +
    scale_y_continuous(limits = c(-5,5),expand = c(0.01,0.075)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 7.5, family = "sans"),
          axis.title = element_blank())
    #labs(x = "Log<sub>2</sub> Fold Change in NY<br>(male vs female)",
         #y = "Log<sub>2</sub> Fold Change in BZ<br>(male vs female)")
  
  return(list(n_total_genes = n_total_genes, n_sig_DE_NY_genes = n_sig_DE_NY_genes, n_sig_DE_BZ_genes = n_sig_DE_BZ_genes,
              n_sigNY = n_sigNY, n_sigBZ = n_sigBZ, n_sigsame = n_sigsame,
              n_sigopps = n_sigopps, n_GxS = n_GxS, GxS_list = GxS_list, plot = plot))
}

# execute function and save resulting plots
warm_liver <- GxS(resNY = res_NY_liver_warm_MvF, resBZ = res_BZ_liver_warm_MvF, resGxS = res_sex_liver_warm_GxS)
plot_warm_liver <- warm_liver$plot
warm_liver_GxS_list <- warm_liver$GxS_list
ggsave("results/figures/sex_DE_warm_liver.pdf", plot = plot_warm_liver, height = 2, width = 2.1)

warm_BAT <- GxS(resNY = res_NY_BAT_warm_MvF, resBZ = res_BZ_BAT_warm_MvF, resGxS = res_sex_BAT_warm_GxS)
plot_warm_BAT <- warm_BAT$plot
warm_BAT_GxS_list <- warm_BAT$GxS_list
ggsave("results/figures/sex_DE_warm_BAT.pdf", plot = plot_warm_BAT, height = 2, width = 2.1)

cold_liver <- GxS(resNY = res_NY_liver_cold_MvF, resBZ = res_BZ_liver_cold_MvF, resGxS = res_sex_liver_cold_GxS)
plot_cold_liver <- cold_liver$plot
cold_liver_GxS_list <- cold_liver$GxS_list
ggsave("results/figures/sex_DE_cold_liver.pdf", plot = plot_cold_liver, height = 2, width = 2.1)

cold_BAT <- GxS(resNY = res_NY_BAT_cold_MvF, resBZ = res_BZ_BAT_cold_MvF, resGxS = res_sex_BAT_cold_GxS)
plot_cold_BAT <- cold_BAT$plot
cold_BAT_GxS_list <- cold_BAT$GxS_list
ggsave("results/figures/sex_DE_cold_BAT.pdf", plot = plot_cold_BAT, height = 2, width = 2.1)




##############################################################
# Statistical Analysis - Chi-square Test
##############################################################

## Is the number of DEG in BZ significantly more/less than the number of DEG in NY?

# create function
x2_setup <- function(n_NY, all_subNY, n_BZ, all_subBZ) #
{
  data <- matrix(c(n_NY, all_subNY, n_BZ, all_subBZ), ncol = 2, byrow = TRUE)
  colnames(data) <- c("DE", "not_DE")
  rownames(data) <- c("NewYork", "Brazil")
  data <- as.table(data)
}

warm_liver_x2 <- x2_setup(n_NY = warm_liver$n_sig_DE_NY_genes, all_subNY = (warm_liver$n_total_genes)-(warm_liver$n_sig_DE_NY_genes),
                          n_BZ = warm_liver$n_sig_DE_BZ_genes, all_subBZ = (warm_liver$n_total_genes)-(warm_liver$n_sig_DE_BZ_genes))

cold_liver_x2 <- x2_setup(n_NY = cold_liver$n_sig_DE_NY_genes, all_subNY = (cold_liver$n_total_genes)-(cold_liver$n_sig_DE_NY_genes),
                          n_BZ = cold_liver$n_sig_DE_BZ_genesZ, all_subBZ = (cold_liver$n_total_genes)-(cold_liver$n_sig_DE_BZ_genes))

warm_BAT_x2 <- x2_setup(n_NY = warm_BAT$n_sig_DE_NY_genes, all_subNY = (warm_BAT$n_total_genes)-(warm_BAT$n_sig_DE_NY_genes),
                          n_BZ = warm_BAT$n_sig_DE_BZ_genes, all_subBZ = (warm_BAT$n_total_genes)-(warm_BAT$n_sig_DE_BZ_genes))

cold_BAT_x2 <- x2_setup(n_NY = cold_BAT$n_sig_DE_NY_genes, all_subNY = (cold_BAT$n_total_genes)-(cold_BAT$n_sig_DE_NY_genes),
                          n_BZ = cold_BAT$n_sig_DE_BZ_genes, all_subBZ = (cold_BAT$n_total_genes)-(cold_BAT$n_sig_DE_BZ_genes))

# perform chi-square tests
warm_liver_x2_test <- chisq.test(warm_liver_x2, correct = FALSE)
cold_liver_x2_test <- chisq.test(cold_liver_x2, correct = FALSE)
warm_BAT_x2_test <- chisq.test(warm_BAT_x2, correct = FALSE)
cold_BAT_x2_test <- chisq.test(cold_BAT_x2, correct = FALSE)

# (I used correct=FALSE since there are large sample sizes (cell >= 5 observations))

