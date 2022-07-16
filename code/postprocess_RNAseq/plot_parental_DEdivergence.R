#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots differential expression in parental samples, with specific focus on patterns of divergence.
# This script generates Figures 1E (males), S3A (females), and S4A (sex) in BallingerMack_2022.

# Major Result(s): divergence is main gene expression pattern: expression divergence between NY and BZ is
# concordant across both environments and sexes

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
# Expression divergence (NY vs BZ) across environments (warm vs cold)
##############################################################

# create function
env_divergence <- function(resW, resC)
{
  # filter to get mean of 10 reads across all samples (per gene)
  DE_base_resW <- resW %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    filter(baseMean >= 10)
  
  DE_base_resC <- resC %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    filter(baseMean >= 10)
  
  # merge datasets
  merge_resWC <- dplyr::full_join(DE_base_resW, DE_base_resC, by = c("gene")) %>%
    dplyr::select(
      "gene",
      "lfc.w" = "log2FoldChange.x",
      "lfcSE.w" = "lfcSE.x",
      "padj.w" = "padj.x",
      "lfc.c" = "log2FoldChange.y",
      "lfcSE.c" = "lfcSE.y",
      "padj.c" = "padj.y"
      )
  n_total_genes <- nrow(merge_resWC)
  n_sig_genes <- merge_resWC %>% filter(padj.w < 0.05 | padj.c < 0.05) %>% nrow()
  
  # add a column indicating whether warm and cold groups are in the same direction
  merge_resWC_direction <-  merge_resWC %>%
    mutate(direction = case_when((sign(lfc.w)) == (sign(lfc.c)) ~ 1, # same direction
                                 (sign(lfc.w)) != (sign(lfc.c)) ~ 0)) # opposite direction
  
  # shrink plotting window by assigning large LFC to 10's
  merge_resWC_direction$lfc.w[merge_resWC_direction$lfc.w>10]<-10
  merge_resWC_direction$lfc.c[merge_resWC_direction$lfc.c>10]<-10
  merge_resWC_direction$lfc.w[merge_resWC_direction$lfc.w<(-10)]<-(-10)
  merge_resWC_direction$lfc.c[merge_resWC_direction$lfc.c<(-10)]<-(-10)
  
  # categorize genes by direction and/or significance
  sigC = dplyr::filter(merge_resWC_direction, padj.c < 0.05 & padj.w >= 0.05)
  sigconserveWvC = dplyr::filter(merge_resWC_direction, direction == 1 & padj.w < 0.05 & padj.c < 0.05)
  sigW = dplyr::filter(merge_resWC_direction, padj.w < 0.05 & padj.c >= 0.05)
  sigoppsWvC = dplyr::filter(merge_resWC_direction, direction == 0 & padj.w < 0.05 & padj.c < 0.05)
  nonsigWvC = dplyr::filter(merge_resWC_direction, padj.w >= 0.05 & padj.c >= 0.05)
  
  # sample sizes of categories:
  n_sigC <- nrow(sigC)
  n_sigconserveWvC <- nrow(sigconserveWvC)
  n_sigW <- nrow(sigW)
  n_sigoppsWvC <- nrow(sigoppsWvC)
  
  # plot
  plot <- merge_resWC_direction %>%
    ggplot(aes(x = lfc.w, y = lfc.c)) +
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_abline(slope = 1, linetype = "dotted") +
    geom_point(data = sigconserveWvC, aes(x = lfc.w, y = lfc.c), fill = "#44B379", color = "grey100", stroke = 0.4, size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
    geom_point(data = nonsigWvC, aes(x = lfc.w, y = lfc.c), fill = "darkgray", color = "grey100", stroke = 0.4, size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
    geom_point(data = sigW, aes(x = lfc.w, y = lfc.c), fill = "#E8E430", color = "grey100", stroke = 0.4, size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
    geom_point(data = sigC, aes(x = lfc.w, y = lfc.c), fill = "#452D72", color = "grey100", stroke = 0.4, size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
    geom_point(data = sigoppsWvC, aes(x = lfc.w, y = lfc.c), fill = "#000000", color = "000000", size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
    scale_x_continuous(limits = c(-10,10), expand = c(0.01,0.075)) +
    scale_y_continuous(limits = c(-10,10), expand = c(0.01,0.075)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text (size = 12.5, family = "sans"),
          axis.title = element_blank())
    #       axis.title.y = element_markdown(vjust = 1, size = 14, family = "sans"),
    #       axis.title.x = element_markdown(vjust = -1.5, size = 14, family = "sans")) +
    # labs(x = "Log<sub>2</sub> Fold Change in warm<br>(NY vs BZ)",
    #      y = "Log<sub>2</sub> Fold Change in cold<br>(NY vs BZ)")
  
  return(list(n_total_genes = n_total_genes, n_sig_genes = n_sig_genes,
              n_sigW = n_sigW, n_sigC = n_sigC, n_sigconserveWvC = n_sigconserveWvC,
              n_sigoppsWvC = n_sigoppsWvC, plot = plot))
}

# execute function and save resulting plots
male_livers <- env_divergence(resW = res_warm_males_liver_NYvsBZ, resC = res_cold_males_liver_NYvsBZ)
plot_males_liver <- male_livers$plot
ggsave("results/figures/males_DE_divergence_liver.pdf", plot = plot_males_liver, height = 3.25, width = 3.5)

male_BAT <- env_divergence(resW = res_warm_males_BAT_NYvsBZ, resC = res_cold_males_BAT_NYvsBZ)
plot_males_BAT <- male_BAT$plot
ggsave("results/figures/males_DE_divergence_BAT.pdf", plot = plot_males_BAT, height = 3.25, width = 3.5)

female_liver <- env_divergence(resW = res_warm_females_liver_NYvsBZ, resC = res_cold_females_liver_NYvsBZ)
plot_females_liver <- female_liver$plot
ggsave("results/figures/females_DE_divergence_liver.pdf", plot = plot_females_liver, height = 2, width = 2.1)

female_BAT <- env_divergence(resW = res_warm_females_BAT_NYvsBZ, resC = res_cold_females_BAT_NYvsBZ)
plot_females_BAT <- female_BAT$plot
ggsave("results/figures/females_DE_divergence_BAT.pdf", plot = plot_females_BAT, height = 2, width = 2.1)




##############################################################
# Sex expression patterns (males vs females) - divergence
##############################################################

# create function
sex_divergence <- function(resM, resF)
{
  # filter to get mean of 10 reads across all samples (per gene)
  DE_base_resM <- resM %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    filter(baseMean >= 10)
  
  DE_base_resF <- resF %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    filter(baseMean >= 10)
  
  # merge datasets
  merge_resMF <- dplyr::full_join(DE_base_resM, DE_base_resF, by = c("gene")) %>%
    dplyr::select(
      "gene",
      "lfc.m" = "log2FoldChange.x",
      "lfcSE.m" = "lfcSE.x",
      "padj.m" = "padj.x",
      "lfc.f" = "log2FoldChange.y",
      "lfcSE.f" = "lfcSE.y",
      "padj.f" = "padj.y"
    )
  n_total_genes <- nrow(merge_resMF)
  n_sig_genes <- merge_resMF %>% filter(padj.m < 0.05 | padj.f < 0.05) %>% nrow()
  
  # add a column indicating whether males and females are in the same direction
  merge_resMF_direction <-  merge_resMF %>%
    mutate(direction = case_when((sign(lfc.m)) == (sign(lfc.f)) ~ 1, # same direction
                                 (sign(lfc.m)) != (sign(lfc.f)) ~ 0)) # opposite direction
  
  # shrink plotting window by assigning large LFC to 10's
  merge_resMF_direction$lfc.m[merge_resMF_direction$lfc.m>10]<-10
  merge_resMF_direction$lfc.f[merge_resMF_direction$lfc.f>10]<-10
  merge_resMF_direction$lfc.m[merge_resMF_direction$lfc.m<(-10)]<-(-10)
  merge_resMF_direction$lfc.f[merge_resMF_direction$lfc.f<(-10)]<-(-10)
  
  # categorize genes by direction and/or significance
  sigM = dplyr::filter(merge_resMF_direction, padj.m < 0.05 & padj.f >= 0.05)
  sigconserveMvF = dplyr::filter(merge_resMF_direction, direction == 1 & padj.m < 0.05 & padj.f < 0.05)
  sigF = dplyr::filter(merge_resMF_direction, padj.f < 0.05 & padj.m >= 0.05)
  sigoppsMvF = dplyr::filter(merge_resMF_direction, direction == 0 & padj.m < 0.05 & padj.f < 0.05)
  nonsigMvF = dplyr::filter(merge_resMF_direction, padj.m >= 0.05 & padj.f >= 0.05)
  
  # sample sizes of categories:
  n_sigM <- nrow(sigM)
  n_sigconserveMvF <- nrow(sigconserveMvF)
  n_sigF <- nrow(sigF)
  n_sigoppsMvF <- nrow(sigoppsMvF)
  
  # plot
  plot <- merge_resMF_direction %>%
    ggplot(aes(x = lfc.m, y = lfc.f)) +
    geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
    geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
    geom_abline(slope = 1, linetype = "dotted", size = 0.25) +
    geom_point(data = sigconserveMvF, aes(x = lfc.m, y = lfc.f), fill = "#458A00", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
    geom_point(data = nonsigMvF, aes(x = lfc.m, y = lfc.f), fill = "darkgray", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
    geom_point(data = sigM, aes(x = lfc.m, y = lfc.f), fill = "#05BFFF", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
    geom_point(data = sigF, aes(x = lfc.m, y = lfc.f), fill = "#E93E8B", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
    geom_point(data = sigoppsMvF_warm_liver, aes(x = lfc.m, y = lfc.f), fill = "#000000", color = "000000", size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
    scale_x_continuous(limits = c(-10,10), expand = c(0.01,0.075)) +
    scale_y_continuous(limits = c(-10,10), expand = c(0.01,0.075)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 7.5, family = "sans"),
          axis.title = element_blank())
  
  return(list(n_total_genes = n_total_genes, n_sig_genes = n_sig_genes,
              n_sigM = n_sigM, n_sigF = n_sigF, n_sigconserveMvF = n_sigconserveMvF,
              n_sigoppsMvF = n_sigoppsMvF, plot = plot))

}

# execute function and save resulting plots
warm_liver_sex <- sex_divergence(resM = res_males_liver_warm_NYvBZ, resF = res_females_liver_warm_NYvBZ)
plot_warm_liver_sex <- warm_liver_sex$plot
ggsave("results/figures/MvF_warm_liver_DE_divergence.pdf", plot = plot_warm_liver_sex, height = 2, width = 2.1)

cold_liver_sex <- sex_divergence(resM = res_males_liver_cold_NYvBZ, resF = res_females_liver_cold_NYvBZ)
plot_cold_liver_sex <- cold_liver_sex$plot
ggsave("results/figures/MvF_cold_liver_DE_divergence.pdf", plot = plot_cold_liver_sex, height = 2, width = 2.1)

warm_BAT_sex <- sex_divergence(resM = res_males_BAT_warm_NYvBZ, resF = res_females_BAT_warm_NYvBZ)
plot_warm_BAT_sex <- warm_BAT_sex$plot
ggsave("results/figures/MvF_warm_BAT_DE_divergence.pdf", plot = plot_warm_BAT_sex, height = 2, width = 2.1)

cold_BAT_sex <- sex_divergence(resM = res_males_BAT_cold_NYvBZ, resF = res_females_BAT_cold_NYvBZ)
plot_cold_BAT_sex <- cold_BAT_sex$plot
ggsave("results/figures/MvF_cold_BAT_DE_divergence.pdf", plot = plot_cold_BAT_sex, height = 2, width = 2.1)

