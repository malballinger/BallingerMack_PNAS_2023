#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script explores broad patterns of genotype (G), sex (S),
# and genotype-by-sex interactions (GxS) on differentially expressed genes
# in parental samples. Liver, BAT, warm, and cold groups were explored separately.
# This script generates Figure S4 in BallingerMack_2022.

# Main Result: both tissues show strong genotype effects, and BAT harbors very little GxS

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clear R's environment
library(tidyverse)
library(DESeq2)
library(ggtext)
library(ggtern)

source("./code/get_DE_parents_analysis.R") # where DESeq data frames are generated
source("./code/plot_parental_GxS.R") # where GxE datasets are defined

##############################################################
# Liver expression (warm) - effects of G, S, and GxS
##############################################################

## Genotype
# filter to get mean of 10 reads across all samples (per gene)
DE_base_geno_sex_liver_warm <- res_sex_liver_warm_G %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## Sex
# filter to get mean of 10 reads across all samples (per gene)
DE_base_sex_liver_warm <- res_sex_liver_warm_S %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## GxS
# filter to get mean of 10 reads across all samples
DE_base_GxS_liver_warm <- res_sex_liver_warm_GxS %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## merge
merge_gs_liver_warm <- dplyr::full_join(DE_base_geno_sex_liver_warm, DE_base_sex_liver_warm, by = "gene") %>%
  dplyr::select("gene",
         "lfc.g" = "log2FoldChange.x",
         "padj.g" = "padj.x",
         "abslfc.g" = "abslfc.x",
         "lfc.s" = "log2FoldChange.y",
         "padj.s" = "padj.y",
         "abslfc.s" = "abslfc.y"
         )
#nrow(merge_gs_liver_warm) 14,153 genes

mergeall_sex_liver_warm <- dplyr::full_join(merge_gs_liver_warm, DE_base_GxS_liver_warm, by = "gene") %>%
  dplyr::select("gene",
         "lfc.g", "padj.g", "abslfc.g",
         "lfc.s", "padj.s", "abslfc.s",
         "lfc.gxs" = "log2FoldChange",
         "padj.gxs" = "padj",
         "abslfc.gxs" = "abslfc"
         )
#nrow(mergeall_sex_liver_warm) 14,153 genes


## effect sizes

# How many sig. DEG for genotype only, with no influence of sex?
g_liver_warm <- mergeall_sex_liver_warm %>%
  filter(padj.g < 0.05 &
         (!gene %in% sex_warm_liver_GxS_list)) # this list of genes has been generated in ./code/plot_parental_GxS.R
#nrow(g_liver_warm) 5,307 at FDR < 0.05
#effectsize_g_liver_warm <- g_liver_warm %>% summarise(effectsize = mean(abslfc.g)) # 0.782

# How many sig. DEG for sex only, with no genotype influence?
s_liver_warm <- mergeall_sex_liver_warm %>%
  filter(padj.s < 0.05 &
         (!gene %in% sex_warm_liver_GxS_list))
#nrow(s_liver_warm) 714 at FDR < 0.05
#effectsize_s_liver_warm <- s_liver_warm %>% summarise(effectsize = mean(abslfc.s)) # 0.844

# Double-check that number of sig. DEG for GxS is same as before (in ./code/plot_parental_GxS.R)
#nrow(sex_warm_liver_GxS) # yes
# 908 at FDR < 0.05


## assign each gene to a main category for plotting purposes

sig_sex_liver_warm <- mergeall_sex_liver_warm %>%
  filter((padj.g < 0.05) | (padj.s < 0.05) | (gene %in% sex_warm_liver_GxS_list)) %>%
  mutate(main_category = case_when(padj.g < 0.05 & (!gene %in% sex_warm_liver_GxS_list) ~ 1, # g_liver_warm = 1
                             padj.s < 0.05 & (!gene %in% sex_warm_liver_GxS_list) ~ 2, # s_liver_warm = 2
                             gene %in% sex_warm_liver_GxS_list ~ 3)) %>% # gxs_liver_warm = 3
  mutate(totalvar = (abslfc.g + abslfc.s + abslfc.gxs)) %>%
  mutate(g_propvar = (abslfc.g/totalvar),
         s_propvar = (abslfc.s/totalvar),
         gxs_propvar = (abslfc.gxs/totalvar))
#nrow(sig_sex_liver_warm) 6,583 genes
  
## plot parental expression patterns of parents as a ternary plot
sex_liver_warm_tern <- sig_sex_liver_warm %>%
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
ggsave("results/figures/GxS_ternary_liver_warm.pdf", plot = sex_liver_warm_tern,
       height = 3, width = 3.5)


##############################################################
# BAT expression (warm) - effects of G, S, and GxS
##############################################################

## Genotype
# filter to get mean of 10 reads across all samples (per gene)
DE_base_geno_sex_BAT_warm <- res_sex_BAT_warm_G %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## Sex
# filter to get mean of 10 reads across all samples (per gene)
DE_base_sex_BAT_warm <- res_sex_BAT_warm_S %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## GxS
# filter to get mean of 10 reads across all samples (per gene)
DE_base_GxS_BAT_warm <- res_sex_BAT_warm_GxS %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## merge
merge_gs_BAT_warm <- dplyr::full_join(DE_base_geno_sex_BAT_warm, DE_base_sex_BAT_warm, by = "gene") %>%
  dplyr::select("gene",
                "lfc.g" = "log2FoldChange.x",
                "padj.g" = "padj.x",
                "abslfc.g" = "abslfc.x",
                "lfc.s" = "log2FoldChange.y",
                "padj.s" = "padj.y",
                "abslfc.s" = "abslfc.y"
  )
#nrow(merge_gs_BAT_warm) 14,162 genes

mergeall_sex_BAT_warm <- dplyr::full_join(merge_gs_BAT_warm, DE_base_GxS_BAT_warm, by = "gene") %>%
  dplyr::select("gene",
                "lfc.g", "padj.g", "abslfc.g",
                "lfc.s", "padj.s", "abslfc.s",
                "lfc.gxs" = "log2FoldChange",
                "padj.gxs" = "padj",
                "abslfc.gxs" = "abslfc"
  )
#nrow(mergeall_sex_BAT_warm) 14,162 genes


## effect sizes

# How many sig. DEG for genotype only, with no influence of sex?
g_BAT_warm <- mergeall_sex_BAT_warm %>%
  filter(padj.g < 0.05 &
           (!gene %in% sex_warm_BAT_GxS_list)) # this list of genes has been generated in ./code/plot_parental_GxS.R
#nrow(g_BAT_warm) 5,647 at FDR < 0.05
#effectsize_g_BAT_warm <- g_BAT_warm %>% summarise(effectsize = mean(abslfc.g)) # 0.731

# How many sig. DEG for sex only, with no genotype influence?
s_BAT_warm <- mergeall_sex_BAT_warm %>%
  filter(padj.s < 0.05 &
           (!gene %in% sex_warm_BAT_GxS_list))
#nrow(s_BAT_warm) 19 at FDR < 0.05
#effectsize_s_BAT_warm <- s_BAT_warm %>% summarise(effectsize = mean(abslfc.s)) # 3.98

# Double-check that number of sig. DEG for GxS is same as before (in ./code/plot_parental_GxS.R)
#nrow(sex_warm_BAT_GxS) # yes
# 13 at FDR < 0.05


## assign each gene to a main category for plotting purposes

sig_sex_BAT_warm <- mergeall_sex_BAT_warm %>%
  filter((padj.g < 0.05) | (padj.s < 0.05) | (gene %in% sex_warm_BAT_GxS_list)) %>%
  mutate(main_category = case_when(padj.g < 0.05 & (!gene %in% sex_warm_BAT_GxS_list) ~ 1, # g_liver_warm = 1
                                   padj.s < 0.05 & (!gene %in% sex_warm_BAT_GxS_list) ~ 2, # s_liver_warm = 2
                                   gene %in% sex_warm_BAT_GxS_list ~ 3)) %>% # gxs_liver_warm = 3
  mutate(totalvar = (abslfc.g + abslfc.s + abslfc.gxs)) %>%
  mutate(g_propvar = (abslfc.g/totalvar),
         s_propvar = (abslfc.s/totalvar),
         gxs_propvar = (abslfc.gxs/totalvar))
#nrow(sig_sex_BAT_warm) 5,673 genes

## plot parental expression patterns of parents as a ternary plot
sex_BAT_warm_tern <- sig_sex_BAT_warm %>%
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
ggsave("results/figures/GxS_ternary_BAT_warm.pdf", plot = sex_BAT_warm_tern,
       height = 3, width = 3.5)


##############################################################
# Liver expression (cold) - effects of G, S, and GxS
##############################################################

## Genotype
# filter to get mean of 10 reads across all samples (per gene)
DE_base_geno_sex_liver_cold <- res_sex_liver_cold_G %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## Sex
# filter to get mean of 10 reads across all samples (per gene)
DE_base_sex_liver_cold <- res_sex_liver_cold_S %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## GxS
# filter to get mean of 10 reads across all samples
DE_base_GxS_liver_cold <- res_sex_liver_cold_GxS %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## merge
merge_gs_liver_cold <- dplyr::full_join(DE_base_geno_sex_liver_cold, DE_base_sex_liver_cold, by = "gene") %>%
  dplyr::select("gene",
                "lfc.g" = "log2FoldChange.x",
                "padj.g" = "padj.x",
                "abslfc.g" = "abslfc.x",
                "lfc.s" = "log2FoldChange.y",
                "padj.s" = "padj.y",
                "abslfc.s" = "abslfc.y"
  )
#nrow(merge_gs_liver_cold) 14,124 genes

mergeall_sex_liver_cold <- dplyr::full_join(merge_gs_liver_cold, DE_base_GxS_liver_cold, by = "gene") %>%
  dplyr::select("gene",
                "lfc.g", "padj.g", "abslfc.g",
                "lfc.s", "padj.s", "abslfc.s",
                "lfc.gxs" = "log2FoldChange",
                "padj.gxs" = "padj",
                "abslfc.gxs" = "abslfc"
  )
#nrow(mergeall_sex_liver_cold) 14,124 genes


## effect sizes

# How many sig. DEG for genotype only, with no influence of sex?
g_liver_cold <- mergeall_sex_liver_cold %>%
  filter(padj.g < 0.05 &
           (!gene %in% sex_cold_liver_GxS_list)) # this list of genes has been generated in ./code/plot_parental_GxS.R
#nrow(g_liver_cold)  # 4,930 genes at FDR < 0.05
#effectsize_g_liver_cold <- g_liver_cold %>% summarise(effectsize = mean(abslfc.g)) # 0.852

# How many sig. DEG for sex only, with no genotype influence?
s_liver_cold <- mergeall_sex_liver_cold %>%
  filter(padj.s < 0.05 &
           (!gene %in% sex_cold_liver_GxS_list))
#nrow(s_liver_cold) # 978 genes at FDR < 0.05
#effectsize_s_liver_cold <- s_liver_cold %>% summarise(effectsize = mean(abslfc.s)) # 0.994

# Double-check that number of sig. DEG for GxS is same as before (in ./code/plot_parental_GxS.R)
#nrow(sex_cold_liver_GxS) # yes
# 947 at FDR < 0.05


## assign each gene to a main category for plotting purposes

sig_sex_liver_cold <- mergeall_sex_liver_cold %>%
  filter((padj.g < 0.05) | (padj.s < 0.05) | (gene %in% sex_cold_liver_GxS_list)) %>%
  mutate(main_category = case_when(padj.g < 0.05 & (!gene %in% sex_cold_liver_GxS_list) ~ 1, # g_liver_warm = 1
                                   padj.s < 0.05 & (!gene %in% sex_cold_liver_GxS_list) ~ 2, # s_liver_warm = 2
                                   gene %in% sex_cold_liver_GxS_list ~ 3)) %>% # gxs_liver_warm = 3
  mutate(totalvar = (abslfc.g + abslfc.s + abslfc.gxs)) %>%
  mutate(g_propvar = (abslfc.g/totalvar),
         s_propvar = (abslfc.s/totalvar),
         gxs_propvar = (abslfc.gxs/totalvar))
#nrow(sig_sex_liver_cold) 6,402 genes

## plot parental expression patterns of parents as a ternary plot
sex_liver_cold_tern <- sig_sex_liver_cold %>%
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
ggsave("results/figures/GxS_ternary_liver_cold.pdf", plot = sex_liver_cold_tern,
       height = 3, width = 3.5)


##############################################################
# BAT expression (cold) - effects of G, S, and GxS
##############################################################

## Genotype
# filter to get mean of 10 reads across all samples (per gene)
DE_base_geno_sex_BAT_cold <- res_sex_BAT_cold_G %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## Sex
# filter to get mean of 10 reads across all samples (per gene)
DE_base_sex_BAT_cold <- res_sex_BAT_cold_S %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## GxS
# filter to get mean of 10 reads across all samples (per gene)
DE_base_GxS_BAT_cold <- res_sex_BAT_cold_GxS %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## merge
merge_gs_BAT_cold <- dplyr::full_join(DE_base_geno_sex_BAT_cold, DE_base_sex_BAT_cold, by = "gene") %>%
  dplyr::select("gene",
                "lfc.g" = "log2FoldChange.x",
                "padj.g" = "padj.x",
                "abslfc.g" = "abslfc.x",
                "lfc.s" = "log2FoldChange.y",
                "padj.s" = "padj.y",
                "abslfc.s" = "abslfc.y"
  )
#nrow(merge_gs_BAT_cold) 15,913 genes

mergeall_sex_BAT_cold <- dplyr::full_join(merge_gs_BAT_cold, DE_base_GxS_BAT_cold, by = "gene") %>%
  dplyr::select("gene",
                "lfc.g", "padj.g", "abslfc.g",
                "lfc.s", "padj.s", "abslfc.s",
                "lfc.gxs" = "log2FoldChange",
                "padj.gxs" = "padj",
                "abslfc.gxs" = "abslfc"
  )
#nrow(mergeall_sex_BAT_cold) 15,913 genes


## effect sizes

# How many sig. DEG for genotype only, with no influence of sex?
g_BAT_cold <- mergeall_sex_BAT_cold %>%
  filter(padj.g < 0.05 &
           (!gene %in% sex_cold_BAT_GxS_list)) # this list of genes has been generated in ./code/plot_parental_GxS.R
#nrow(g_BAT_cold) 5,454 at FDR < 0.05
#effectsize_g_BAT_cold <- g_BAT_cold %>% summarise(effectsize = mean(abslfc.g)) # 0.875

# How many sig. DEG for sex only, with no genotype influence?
s_BAT_cold <- mergeall_sex_BAT_cold %>%
  filter(padj.s < 0.05 &
           (!gene %in% sex_cold_BAT_GxS_list))
#nrow(s_BAT_cold) 27 at FDR < 0.05
#effectsize_s_BAT_cold <- s_BAT_cold %>% summarise(effectsize = mean(abslfc.s)) # 2.95

# Double-check that number of sig. DEG for GxS is same as before (in ./code/plot_parental_GxS.R)
#nrow(sex_cold_BAT_GxS) # yes
# 18 at FDR < 0.05


## assign each gene to a main category for plotting purposes

sig_sex_BAT_cold <- mergeall_sex_BAT_cold %>%
  filter((padj.g < 0.05) | (padj.s < 0.05) | (gene %in% sex_cold_BAT_GxS_list)) %>%
  mutate(main_category = case_when(padj.g < 0.05 & (!gene %in% sex_cold_BAT_GxS_list) ~ 1, # g_liver_warm = 1
                                   padj.s < 0.05 & (!gene %in% sex_cold_BAT_GxS_list) ~ 2, # s_liver_warm = 2
                                   gene %in% sex_cold_BAT_GxS_list ~ 3)) %>% # gxs_liver_warm = 3
  mutate(totalvar = (abslfc.g + abslfc.s + abslfc.gxs)) %>%
  mutate(g_propvar = (abslfc.g/totalvar),
         s_propvar = (abslfc.s/totalvar),
         gxs_propvar = (abslfc.gxs/totalvar))
#nrow(sig_sex_BAT_cold) 5,494 genes

## plot parental expression patterns of parents as a ternary plot
sex_BAT_cold_tern <- sig_sex_BAT_cold %>%
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
ggsave("results/figures/GxS_ternary_BAT_cold.pdf", plot = sex_BAT_cold_tern,
       height = 3, width = 3.5)

