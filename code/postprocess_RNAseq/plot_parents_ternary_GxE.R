#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script explores broad patterns of the effects of genotype (g), environment (e),
# and genotype-by-environment interactions (GxE) on differentially expressed genes
# in parental samples. Liver, BAT, males, and females were explored separately.
# This script generates Figures 2A (males) and S3B (females) in BallingerMack_2022.

# Major Result(s): BAT shows more environmental effects (plasticity) than liver; both tissues show strong genotype effects

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clear R's environment
library(tidyverse)
library(DESeq2)
library(ggtext)
library(ggtern)

source("./code/postprocess_RNAseq/get_DE_parents_analysis.R") # where DESeq data frames are generated
source("./code/postprocess_RNAseq/plot_parental_GxE.R") # where GxE datasets are defined

##############################################################
# Generate ternaries
##############################################################

# create function
ternary_env <- function(resG, resE, resGE, GxE_list)
{
# Genotype
# filter to get mean of 10 reads across all samples (per gene)
DE_base_resG <- resG %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

# Environment
# filter to get mean of 10 reads across all samples (per gene)
DE_base_resE <- resE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

# GxE
# filter to get mean of 10 reads across all samples
DE_base_resGxE <- resGE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

# merge datasets
merge_resGE <- dplyr::full_join(DE_base_resG, DE_base_resE, by = "gene") %>%
  dplyr::select("gene",
         "lfc.g" = "log2FoldChange.x",
         "padj.g" = "padj.x",
         "abslfc.g" = "abslfc.x",
         "lfc.e" = "log2FoldChange.y",
         "padj.e" = "padj.y",
         "abslfc.e" = "abslfc.y"
         )

mergeall_resGE <- dplyr::full_join(merge_resGE, DE_base_resGxE, by = "gene") %>%
  dplyr::select("gene",
         "lfc.g", "padj.g", "abslfc.g",
         "lfc.e", "padj.e", "abslfc.e",
         "lfc.gxe" = "log2FoldChange",
         "padj.gxe" = "padj",
         "abslfc.gxe" = "abslfc"
         )

n_total_genes <- nrow(mergeall_resGE)
n_sig_genes <- mergeall_resGE %>% filter(padj.g < 0.05 | padj.e < 0.05 | padj.gxe < 0.05) %>% nrow()

# effect sizes
# How many sig. DEG for genotype only, with no environmental influence?
g <- mergeall_resGE %>%
  filter(padj.g < 0.05 &
         (!gene %in% GxE_list)) # this list of genes has been generated in ./code/postprocess_RNAseq/plot_males_parents_GxE.R
n_g_genes <- nrow(g)
effectsize_g <- g %>% summarise(effectsize = mean(abslfc.g))

# How many sig. DEG for env only, with no genotype influence?
e <- mergeall_resGE %>%
  filter(padj.e < 0.05 &
         (!gene %in% GxE_list))
n_e_genes <- nrow(e)
effectsize_e <- e %>% summarise(effectsize = mean(abslfc.e))

# assign each gene to a main category for plotting purposes
sig_groups <- mergeall_resGE %>%
  filter((padj.g < 0.05) | (padj.e < 0.05) | (gene %in% GxE_list)) %>%
  mutate(main_category = case_when(
                                    padj.g < 0.05 & (!gene %in% GxE_list) ~ 1, # g_liver = 1
                                    padj.e < 0.05 & (!gene %in% GxE_list) ~ 2, # e_liver = 1
                                    gene %in% GxE_list ~ 3)) %>% # gxe_liver = 3
  mutate(totalvar = (abslfc.g + abslfc.e + abslfc.gxe)) %>%
  mutate(g_propvar = (abslfc.g/totalvar),
         e_propvar = (abslfc.e/totalvar),
         gxe_propvar = (abslfc.gxe/totalvar))
n_sig_genes_2 <- nrow(sig_groups) # should be same as n_sig_genes

# plot parental expression patterns of parents as a ternary plot
plot <- sig_groups %>%
  ggtern(aes(x=g_propvar, y=gxe_propvar, z=e_propvar)) +
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
       y = "GxE",
       yarrow = NULL,
       z = "E",
       zarrow = NULL) #+
  # geom_density_tern(color="#884C48",linetype=1, alpha=1, size = 1, n = 100, bins = 5,
  #                   stat = "DensityTern", base = "ilr", bdl = 0.010, na.rm = FALSE)

return(list(n_total_genes = n_total_genes, n_sig_genes = n_sig_genes, n_sig_genes_2 = n_sig_genes_2, plot = plot,
            n_g_genes = n_g_genes, effectsize_g = effectsize_g, n_e_genes = n_e_genes, effectsize_e = effectsize_e))

}

# execute function and save resulting plots
males_liver_tern <- ternary_env(resG = res_males_liver_G, resE = res_males_liver_E,
                                resGE = res_males_liver_GxE, GxE_list = male_liver_GxE_list)
males_liver_tern_plot <- males_liver_tern$plot
males_liver_tern_n_g_genes <- males_liver_tern$n_g_genes
males_liver_tern_n_e_genes <- males_liver_tern$n_e_genes
males_liver_tern_effectsize_g <- males_liver_tern$effectsize_g
males_liver_tern_effectsize_e <- males_liver_tern$effectsize_e
ggsave("results/figures/males_GxE_ternary_liver.pdf", plot = males_liver_tern_plot, height = 3, width = 3.5)


males_BAT_tern <- ternary_env(resG = res_males_BAT_G, resE = res_males_BAT_E,
                              resGE = res_males_BAT_GxE, GxE_list = male_BAT_GxE_list)
males_BAT_tern_plot <- males_BAT_tern$plot
males_BAT_tern_n_g_genes <- males_BAT_tern$n_g_genes
males_BAT_tern_n_e_genes <- males_BAT_tern$n_e_genes
males_BAT_tern_effectsize_g <- males_BAT_tern$effectsize_g
males_BAT_tern_effectsize_e <- males_BAT_tern$effectsize_e
ggsave("results/figures/males_GxE_ternary_BAT.pdf", plot = males_BAT_tern_plot, height = 3, width = 3.5)


females_liver_tern <- ternary_env(resG = res_females_liver_G, resE = res_females_liver_E,
                                  resGE = res_females_liver_GxE, GxE_list = female_liver_GxE_list)
females_liver_tern_plot <- females_liver_tern$plot
females_liver_tern_n_g_genes <- females_liver_tern$n_g_genes
females_liver_tern_n_e_genes <- females_liver_tern$n_e_genes
females_liver_tern_effectsize_g <- females_liver_tern$effectsize_g
females_liver_tern_effectsize_e <- females_liver_tern$effectsize_e
ggsave("results/figures/females_GxE_ternary_liver.pdf", plot = females_liver_tern_plot, height = 3, width = 3.5)


females_BAT_tern <- ternary_env(resG = res_females_BAT_G, resE = res_females_BAT_E,
                                resGE = res_females_BAT_GxE, GxE_list = female_BAT_GxE_list)
females_BAT_tern_plot <- females_BAT_tern$plot
females_BAT_tern_n_g_genes <- females_BAT_tern$n_g_genes
females_BAT_tern_n_e_genes <- females_BAT_tern$n_e_genes
females_BAT_tern_effectsize_g <- females_BAT_tern$effectsize_g
females_BAT_tern_effectsize_e <- females_BAT_tern$effectsize_e
ggsave("results/figures/females_GxE_ternary_BAT.pdf", plot = females_BAT_tern_plot, height = 3, width = 3.5)

