#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script explores broad patterns of genotype (G), sex (S),
# and genotype-by-sex interactions (GxS) on differentially expressed genes
# in parental samples. Liver, BAT, warm, and cold groups were explored separately.
# This script generates Figure S4B in BallingerMack_2022.

# Main Result: both tissues show strong genotype effects; liver shows greater sex-effects

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clear R's environment
library(tidyverse)
library(DESeq2)
library(ggtext)
library(ggtern)

source("./code/postprocess_RNAseq/get_DE_parents_analysis.R") # where DESeq data frames are generated
source("./code/postprocess_RNAseq/plot_parental_GxS.R") # where GxS datasets are defined

##############################################################
# Generate ternaries
##############################################################

# create function
ternary_sex <- function(resG, resS, resGS, GxS_list)
{
  # Genotype
  # filter to get mean of 10 reads across all samples (per gene)
  DE_base_resG <- resG %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    filter(baseMean >= 10) %>%
    mutate(abslfc = abs(log2FoldChange))
  
  # Sex
  # filter to get mean of 10 reads across all samples (per gene)
  DE_base_resS <- resS %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    filter(baseMean >= 10) %>%
    mutate(abslfc = abs(log2FoldChange))
  
  # GxS
  # filter to get mean of 10 reads across all samples
  DE_base_resGS <- resGS %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble() %>%
    filter(baseMean >= 10) %>%
    mutate(abslfc = abs(log2FoldChange))
  
  # merge datasets
  merge_resGS <- dplyr::full_join(DE_base_resG, DE_base_resS, by = "gene") %>%
    dplyr::select("gene",
           "lfc.g" = "log2FoldChange.x",
           "padj.g" = "padj.x",
           "abslfc.g" = "abslfc.x",
           "lfc.s" = "log2FoldChange.y",
           "padj.s" = "padj.y",
           "abslfc.s" = "abslfc.y"
           )
  
  mergeall_resGS <- dplyr::full_join(merge_resGS, DE_base_resGS, by = "gene") %>%
    dplyr::select("gene",
           "lfc.g", "padj.g", "abslfc.g",
           "lfc.s", "padj.s", "abslfc.s",
           "lfc.gxs" = "log2FoldChange",
           "padj.gxs" = "padj",
           "abslfc.gxs" = "abslfc"
           )
  
  n_total_genes <- nrow(mergeall_resGS)
  n_sig_genes <- mergeall_resGS %>% filter(padj.g < 0.05 | padj.s < 0.05 | padj.gxs < 0.05) %>% nrow()
  
  # effect sizes
  # How many sig. DEG for genotype only, with no influence of sex?
  g <- mergeall_resGS %>%
    filter(padj.g < 0.05 &
           (!gene %in% GxS_list)) # this list of genes has been generated in ./code/postprocess_RNAseq/plot_parental_GxS.R
  n_g_genes <- nrow(g)
  effectsize_g <- g %>% summarise(effectsize = mean(abslfc.g))
  
  # How many sig. DEG for sex only, with no genotype influence?
  s <- mergeall_resGS %>%
    filter(padj.s < 0.05 &
           (!gene %in% GxS_list))
  n_s_genes <- nrow(s)
  effectsize_s <- s %>% summarise(effectsize = mean(abslfc.s))
  
  # assign each gene to a main category for plotting purposes
  sig_groups <- mergeall_resGS %>%
    filter((padj.g < 0.05) | (padj.s < 0.05) | (gene %in% GxS_list)) %>%
    mutate(main_category = case_when(
                                      padj.g < 0.05 & (!gene %in% GxS_list) ~ 1, # g = 1
                                      padj.s < 0.05 & (!gene %in% GxS_list) ~ 2, # s = 2
                                      gene %in% GxS_list ~ 3)) %>% # gxs = 3
    mutate(totalvar = (abslfc.g + abslfc.s + abslfc.gxs)) %>%
    mutate(g_propvar = (abslfc.g/totalvar),
           s_propvar = (abslfc.s/totalvar),
           gxs_propvar = (abslfc.gxs/totalvar))
  n_sig_genes_2 <- nrow(sig_groups) # should be same as n_sig_genes
    
  # plot parental expression patterns of parents as a ternary plot
  plot <- sig_groups %>%
    ggtern(aes(x=g_propvar, y=gxs_propvar, z=s_propvar)) +
    geom_point(fill = "black", color = "black", size = 2, stroke = 0.5, shape = 21, alpha = 0.1) +
    geom_mask() +
    scale_L_continuous(breaks=c(.25,.50,.75,1)) +
    scale_R_continuous(breaks=c(.25,.50,.75,1)) +
    scale_T_continuous(breaks=c(.25,.50,.75,1)) +
    theme_bw(base_size=6)+
    theme_showarrows() +
    theme(axis.text = element_text(size = 10, family = "sans"),
          axis.title = element_text(size = 12, family = "sans"),
          aspect.ratio = 1/1,
          tern.axis.arrow = element_line(color = "black", size = 1)) +
    labs(x = "G",
         xarrow = NULL,
         y = "GxS",
         yarrow = NULL,
         z = "S",
         zarrow = NULL) #+
    # geom_density_tern(color="#884C48",linetype=1, alpha=1, size = 1, n = 100, bins = 5,
    #                   stat = "DensityTern", base = "ilr", bdl = 0.010, na.rm = FALSE)
  
  return(list(n_total_genes = n_total_genes, n_sig_genes = n_sig_genes, n_sig_genes_2 = n_sig_genes_2, plot = plot,
              n_g_genes = n_g_genes, effectsize_g = effectsize_g, n_s_genes = n_s_genes, effectsize_s = effectsize_s))

}

# execute function and save resulting plots
warm_liver_tern <- ternary_sex(resG = res_sex_liver_warm_G, resS = res_sex_liver_warm_S,
                      resGS = res_sex_liver_warm_GxS, GxS_list = warm_liver_GxS_list)
warm_liver_tern_effectsize_g = warm_liver_tern$effectsize_g
warm_liver_tern_effectsize_s = warm_liver_tern$effectsize_s
warm_liver_tern_plot <- warm_liver_tern$plot
ggsave("results/figures/GxS_ternary_liver_warm.pdf", plot = warm_liver_tern_plot, height = 3, width = 3.5)


warm_BAT_tern <- ternary_sex(resG = res_sex_BAT_warm_G, resS = res_sex_BAT_warm_S,
                         resGS = res_sex_BAT_warm_GxS, GxS_list = warm_BAT_GxS_list)
warm_BAT_tern_effectsize_g = warm_BAT_tern$effectsize_g
warm_BAT_tern_effectsize_s = warm_BAT_tern$effectsize_s
warm_BAT_tern_plot <- warm_BAT_tern$plot
ggsave("results/figures/GxS_ternary_BAT_warm.pdf", plot = warm_BAT_tern_plot, height = 3, width = 3.5)


cold_liver_tern <- ternary_sex(resG = res_sex_liver_cold_G, resS = res_sex_liver_cold_S,
                           resGS = res_sex_liver_cold_GxS, GxS_list = cold_liver_GxS_list)
cold_liver_tern_effectsize_g = cold_liver_tern$effectsize_g
cold_liver_tern_effectsize_s = cold_liver_tern$effectsize_s
cold_liver_tern_plot <- cold_liver_tern$plot
ggsave("results/figures/GxS_ternary_liver_cold.pdf", plot = cold_liver_tern_plot, height = 3, width = 3.5)


cold_BAT_tern <- ternary_sex(resG = res_sex_BAT_cold_G, resS = res_sex_BAT_cold_S,
                         resGS = res_sex_BAT_cold_GxS, GxS_list = cold_BAT_GxS_list)
cold_BAT_tern_effectsize_g = cold_BAT_tern$effectsize_g
cold_BAT_tern_effectsize_s = cold_BAT_tern$effectsize_s
cold_BAT_tern_plot <- cold_BAT_tern$plot
ggsave("results/figures/GxS_ternary_BAT_cold.pdf", plot = cold_BAT_tern_plot, height = 3, width = 3.5)

