#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script explores broad patterns of the effects of genotype (g), environment (e),
# and genotype-by-environment interactions (GxE) on differentially expressed genes
# in male parental samples. Liver and BAT were explored separately.
# This script generates Figure 2A in BallingerMack_2022.

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clear R's environment
library(tidyverse)
library(DESeq2)
library(ggtext)
library(ggtern)

source("./code/get_DE_males_parents.R") # where DESeq datasets are generated
source("./code/plot_males_parents_GxE.R") # where GxE datasets are defined

##############################################################
# Liver expression - effects of G, E, and GxE
##############################################################

## Genotype
# filter to get mean of 10 reads across all samples (per gene)
DE_base_geno_males_liver <- res_males_liver_G %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## Environment
# filter to get mean of 10 reads across all samples (per gene)
DE_base_env_males_liver <- res_males_liver_E %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## GxE
# filter to get mean of 10 reads across all samples
DE_base_GxE_males_liver <- res_males_liver_GxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## merge
merge_ge_males_liver <- dplyr::full_join(DE_base_geno_males_liver, DE_base_env_males_liver, by = "gene") %>%
  dplyr::select("gene",
         "lfc.g" = "log2FoldChange.x",
         "padj.g" = "padj.x",
         "abslfc.g" = "abslfc.x",
         "lfc.e" = "log2FoldChange.y",
         "padj.e" = "padj.y",
         "abslfc.e" = "abslfc.y"
         )
# 14,028 genes

mergeall_males_liver <- dplyr::full_join(merge_ge_males_liver, DE_base_GxE_males_liver, by = "gene") %>%
  dplyr::select("gene",
         "lfc.g", "padj.g", "abslfc.g",
         "lfc.e", "padj.e", "abslfc.e",
         "lfc.gxe" = "log2FoldChange",
         "padj.gxe" = "padj",
         "abslfc.gxe" = "abslfc"
         )
# 14,028 genes


## effect sizes

# How many sig. DE genes for genotype only, with no environmental influence?
g_males_liver <- mergeall_males_liver %>%
  filter(padj.g < 0.05 &
         (!gene %in% liver_males_plast_GxE_list)) # this list of genes (liver_males_plast_GxE_list) has been generated in ./code/plot_males_parents_GxE.R
# 5,422 at FDR < 0.05
effectsize_g_males_liver <- g_males_liver %>% summarise(effectsize = mean(abslfc.g)) # 0.802

# How many sig. DE genes for env only, with no geno inlfuence?
e_males_liver <- mergeall_males_liver %>%
  filter(padj.e < 0.05 &
         (!gene %in% liver_males_plast_GxE_list))
# 553 at FDR < 0.05
effectsize_e_males_liver <- e_males_liver %>% summarise(effectsize = mean(abslfc.e)) # 0.456

# Double-check that number of sig. DE genes for GxE is same as before
#nrow(liver_males_plast_GxE) # yes
# 817 at FDR < 0.05


## assign each gene to a main category for plotting purposes

sig_males_liver <- mergeall_males_liver %>%
  filter((padj.g < 0.05) | (padj.e < 0.05) | (gene %in% liver_males_plast_GxE_list)) %>%
  mutate(main_category = case_when(padj.g < 0.05 & (!gene %in% liver_males_plast_GxE_list) ~ 1, # g_liver = 1
                             padj.e < 0.05 & (!gene %in% liver_males_plast_GxE_list) ~ 2, # e_liver = 1
                             gene %in% liver_males_plast_GxE_list ~ 3)) %>% # gxe_liver = 3
  mutate(totalvar = (abslfc.g + abslfc.e + abslfc.gxe)) %>%
  mutate(g_propvar = (abslfc.g/totalvar),
         e_propvar = (abslfc.e/totalvar),
         gxe_propvar = (abslfc.gxe/totalvar))
# 6,506 genes
  
## plot parental expression patterns of parents as a terny plot
liver_males_tern <-
  ggtern(data=sig_males_liver, aes(x=g_propvar, y=gxe_propvar, z=e_propvar)) +
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
ggsave("results/figures/males_GxE_ternary_liver.pdf", plot = liver_males_tern,
       height = 3, width = 3.5)

# liver_males_tern_talk <- ggtern(data=sig_males_liver, aes(x=g_propvar, y=gxe_propvar, z=e_propvar)) +
#   geom_point(fill = "black", color = "grey100", size = 2, stroke = 0.5, shape = 21, alpha = 0.8) +
#   theme_bw(base_size=6)+
#   theme(axis.text = element_text(size = 12)) +
#   labs(x = "",
#        y="",
#        z="")



##############################################################
# BAT expression - effects of G, E, and GxE
##############################################################

## Genotype
# filter to get mean of 10 reads across all samples
DE_base_geno_BAT <- res_males_BAT_G %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## Environment
# filter to get mean of 10 reads across all samples
DE_base_env_BAT <- res_males_BAT_E %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## GxE
# filter to get mean of 10 reads across all samples
DE_base_GxE_males_BAT <- res_males_BAT_GxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

# merge
merge_ge_males_BAT <- dplyr::full_join(DE_base_geno_BAT, DE_base_env_BAT, by = "gene") %>%
  dplyr::select("gene",
         "lfc.g" = "log2FoldChange.x",
         "padj.g" = "padj.x",
         "abslfc.g" = "abslfc.x",
         "lfc.e" = "log2FoldChange.y",
         "padj.e" = "padj.y",
         "abslfc.e" = "abslfc.y"
  )
# 14,176 genes

mergeall_males_BAT <- dplyr::left_join(merge_ge_males_BAT, DE_base_GxE_males_BAT, by = "gene") %>%
  dplyr::select("gene",
         "lfc.g", "padj.g", "abslfc.g",
         "lfc.e", "padj.e", "abslfc.e",
         "lfc.gxe" = "log2FoldChange",
         "padj.gxe" = "padj",
         "abslfc.gxe" = "abslfc"
  )
# 14,176 genes


# How many sig. DE genes for geno only, with no environmental influence?
g_males_BAT <- mergeall_males_BAT %>%
  filter(padj.g < 0.05 &
           (!gene %in% BAT_males_plast_GxE_list)) # this list of genes (BAT_males_plast_GxE_list) has been generated in ./code/plot_males_parents_GxE.R
# 4,922 at FDR < 0.05
effectsize_g_males_BAT <- g_males_BAT %>% summarise(effectsize = mean(abslfc.g)) # 0.791

# How many sig. DE genes for env only, with no geno inlfuence?
e_males_BAT <- mergeall_males_BAT %>%
  filter(padj.e < 0.05 &
           (!gene %in% BAT_males_plast_GxE_list))
# 1,180 at FDR < 0.05
effectsize_e_males_BAT <- e_males_BAT %>% summarise(effectsize = mean(abslfc.e)) # 0.479

# How many sig. DE genes for GxE?
nrow(BAT_males_plast_GxE)
# 1,425 at FDR < 0.05

# keep only significant genes (in either g, e, or gxe) and make new columns

sig_males_BAT <- mergeall_males_BAT %>%
  filter((padj.g < 0.05) | (padj.e < 0.05) | (gene %in% BAT_males_plast_GxE_list)) %>%
  #replace(is.na(.), 0) %>%
  mutate(main_category = case_when(padj.g < 0.05 & (!gene %in% BAT_males_plast_GxE_list) & padj.e >= 0.05 ~ 1, # g_liver = 1
                                   padj.e < 0.05 & (!gene %in% BAT_males_plast_GxE_list) & padj.g >= 0.05 ~ 2, # e_liver = 1
                                   gene %in% BAT_males_plast_GxE_list ~ 3)) %>% # gxe_liver = 3
  mutate(totalvar = (abslfc.g + abslfc.e + abslfc.gxe)) %>%
  mutate(g_propvar = (abslfc.g/totalvar),
         e_propvar = (abslfc.e/totalvar),
         gxe_propvar = (abslfc.gxe/totalvar))
# 6,971 genes

### plot parental expression patterns of parents - LIVER
BAT_males_tern <-
  ggtern(data=sig_males_BAT, aes(x=g_propvar, y=gxe_propvar, z=e_propvar)) +
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
  # geom_density_tern(color="#FCECAA",linetype=1, alpha=1, size = 01, n = 100, bins = 5,
  #                   stat = "DensityTern", base = "ilr", bdl = 0.010, na.rm = FALSE)
ggsave("results/figures/males_GxE_ternary_BAT.pdf", plot = BAT_males_tern,
       height = 3, width = 3.5)

