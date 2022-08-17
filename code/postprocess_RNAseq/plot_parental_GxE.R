#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots differential expression in parental samples, with specific focus on patterns of GxE.
# This script generates Figures 2B (males) and S3C (females) in BallingerMack_2022.

# Main Result(s): BZ harbors more GxE than NY; canalization of NY suggests adaptation to cold environment

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
# Generate GxE scatter plots
##############################################################

# create function
GxE <- function(resNY, resBZ, resGxE)
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
  
  # add a column indicating whether warm and cold groups are in the same direction
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
  
  # collate all genes that show evidence of GxE
  # first, genes identified as GxE by DESEQ (very stringent)
  DESEQ_GxE <- resGxE %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    filter(baseMean >= 10, padj < 0.05)
  # save a list of GxE genes defined by DESEQ
  DESEQ_GxE_list <- DESEQ_GxE %>% pull(gene)
  
  # second, genes I identify as GxE from categorization above + DESEQ genes
  GxE <- merge_resBZNY_direction %>%
    filter((padj.NY < 0.05 & padj.BZ >= 0.05) |
             (padj.BZ < 0.05 & padj.NY >= 0.05) |
             (direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) |
             (gene %in% DESEQ_GxE_list))
  # save a list of all GxE genes
  GxE_list <- GxE %>% pull(gene)
  
  # sample sizes of categories:
  n_sigNY <- nrow(sigNY)
  list_sigNY <- sigNY %>% pull(gene)
  n_sigsame <- nrow(sigsame)
  n_sigBZ <- nrow(sigBZ)
  n_sigopps <- nrow(sigopps)
  n_GxE <- nrow(GxE)
  
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
    #labs(x = "Log<sub>2</sub> Fold Change in NY<br>(warm vs cold)",
         #y = "Log<sub>2</sub> Fold Change in BZ<br>(warm vs cold)")
  
  return(list(n_total_genes = n_total_genes, n_sig_DE_NY_genes = n_sig_DE_NY_genes, n_sig_DE_BZ_genes = n_sig_DE_BZ_genes,
              n_sigNY = n_sigNY, n_sigBZ = n_sigBZ, n_sigsame = n_sigsame, merge_resBZNY = merge_resBZNY,
              n_sigopps = n_sigopps, n_GxE = n_GxE, GxE_list = GxE_list, plot = plot))
}

# execute function and save resulting plots
male_liver_gxe <- GxE(resNY = res_NY_males_liver_WvC, resBZ = res_BZ_males_liver_WvC, resGxE = res_males_liver_GxE)
plot_male_liver_gxe <- male_liver_gxe$plot
male_liver_GxE_list <- male_liver_gxe$GxE_list
n_GxE_genes_male_liver <- male_liver_gxe$n_GxE
n_total_genes_male_liver <- male_liver_gxe$n_total_genes
n_total_genes_male_liver_NY <- male_liver_gxe$n_sig_DE_NY_genes
n_total_genes_male_liver_BZ <- male_liver_gxe$n_sig_DE_BZ_genes
ggsave("results/figures/males_DE_plasticity_liver.pdf", plot = plot_male_liver_gxe, height = 2, width = 2.1)

male_BAT_gxe <- GxE(resNY = res_NY_males_BAT_WvC, resBZ = res_BZ_males_BAT_WvC, resGxE = res_males_BAT_GxE)
plot_male_BAT_gxe <- male_BAT_gxe$plot
male_BAT_GxE_list <- male_BAT_gxe$GxE_list
n_GxE_genes_male_BAT <- male_BAT_gxe$n_GxE
n_total_genes_male_BAT <- male_BAT_gxe$n_total_genes
n_total_genes_male_BAT_NY <- male_BAT_gxe$n_sig_DE_NY_genes
n_total_genes_male_BAT_BZ <- male_BAT_gxe$n_sig_DE_BZ_genes
ggsave("results/figures/males_DE_plasticity_BAT.pdf", plot = plot_male_BAT_gxe, height = 2, width = 2.1)

female_liver_gxe <- GxE(resNY = res_NY_females_liver_WvC, resBZ = res_BZ_females_liver_WvC, resGxE = res_females_liver_GxE)
plot_female_liver_gxe <- female_liver_gxe$plot
female_liver_GxE_list <- female_liver_gxe$GxE_list
n_GxE_genes_female_liver <- female_liver_gxe$n_GxE
n_total_genes_female_liver <- female_liver_gxe$n_total_genes
n_total_genes_female_liver_NY <- female_liver_gxe$n_sig_DE_NY_genes
n_total_genes_female_liver_BZ <- female_liver_gxe$n_sig_DE_BZ_genes
ggsave("results/figures/females_DE_plasticity_liver.pdf", plot = plot_female_liver_gxe, height = 2, width = 2.1)

female_BAT_gxe <- GxE(resNY = res_NY_females_BAT_WvC, resBZ = res_BZ_females_BAT_WvC, resGxE = res_females_BAT_GxE)
plot_female_BAT_gxe <- female_BAT_gxe$plot
female_BAT_GxE_list <- female_BAT_gxe$GxE_list
n_GxE_genes_female_BAT <- female_BAT_gxe$n_GxE
n_total_genes_female_BAT <- female_BAT_gxe$n_total_genes
n_total_genes_female_BAT_NY <- female_BAT_gxe$n_sig_DE_NY_genes
n_total_genes_female_BAT_BZ <- female_BAT_gxe$n_sig_DE_BZ_genes
ggsave("results/figures/females_DE_plasticity_BAT.pdf", plot = plot_female_BAT_gxe, height = 2, width = 2.1)




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

male_liver_gxe_x2 <- x2_setup(n_NY = male_liver_gxe$n_sig_DE_NY_genes, all_subNY = (male_liver_gxe$n_total_genes)-(male_liver_gxe$n_sig_DE_NY_genes),
                              n_BZ = male_liver_gxe$n_sig_DE_BZ_genes, all_subBZ = (male_liver_gxe$n_total_genes)-(male_liver_gxe$n_sig_DE_BZ_genes))

male_BAT_gxe_x2 <- x2_setup(n_NY = male_BAT_gxe$n_sig_DE_NY_genes, all_subNY = (male_BAT_gxe$n_total_genes)-(male_BAT_gxe$n_sig_DE_NY_genes),
                            n_BZ = male_BAT_gxe$n_sig_DE_BZ_genes, all_subBZ = (male_BAT_gxe$n_total_genes)-(male_BAT_gxe$n_sig_DE_BZ_genes))

female_liver_gxe_x2 <- x2_setup(n_NY = female_liver_gxe$n_sig_DE_NY_genes, all_subNY = (female_liver_gxe$n_total_genes)-(female_liver_gxe$n_sig_DE_NY_genes),
                                n_BZ = female_liver_gxe$n_sig_DE_BZ_genes, all_subBZ = (female_liver_gxe$n_total_genes)-(female_liver_gxe$n_sig_DE_BZ_genes))

female_BAT_gxe_x2 <- x2_setup(n_NY = female_BAT_gxe$n_sig_DE_NY_genes, all_subNY = (female_BAT_gxe$n_total_genes)-(female_BAT_gxe$n_sig_DE_NY_genes),
                              n_BZ = female_BAT_gxe$n_sig_DE_BZ_genes, all_subBZ = (female_BAT_gxe$n_total_genes)-(female_BAT_gxe$n_sig_DE_BZ_genes))

# perform chi-square tests
male_liver_gxe_x2_test <- chisq.test(male_liver_gxe_x2, correct = FALSE)
male_BAT_gxe_x2_test <- chisq.test(male_BAT_gxe_x2, correct = FALSE)
female_liver_gxe_x2_test <- chisq.test(female_liver_gxe_x2, correct = FALSE)
female_BAT_gxe_x2_test <- chisq.test(female_BAT_gxe_x2, correct = FALSE)

# (I used correct=FALSE since there are large sample sizes (cell >= 5 observations))

