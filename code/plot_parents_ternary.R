#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script explores broad patterns of the effects of genotype (g), environment (e),
# and genotype-by-environment interactions (GxE) on differentially expressed genes
# in parental samples. Liver, BAT, males, and females were explored separately.
# This script generates Figures 2A (males) and S3B (females) in BallingerMack_2022.

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clear R's environment
library(tidyverse)
library(DESeq2)
library(ggtext)
library(ggtern)

source("./code/get_DE_parents_analysis.R") # where DESeq data frames are generated
source("./code/plot_parental_GxE.R") # where GxE datasets are defined

##############################################################
# Male liver expression - effects of G, E, and GxE
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
#nrow(merge_ge_males_liver) 14,028 genes

mergeall_males_liver <- dplyr::full_join(merge_ge_males_liver, DE_base_GxE_males_liver, by = "gene") %>%
  dplyr::select("gene",
         "lfc.g", "padj.g", "abslfc.g",
         "lfc.e", "padj.e", "abslfc.e",
         "lfc.gxe" = "log2FoldChange",
         "padj.gxe" = "padj",
         "abslfc.gxe" = "abslfc"
         )
#nrow(mergeall_males_liver) 14,028 genes


## effect sizes

# How many sig. DEG for genotype only, with no environmental influence?
g_males_liver <- mergeall_males_liver %>%
  filter(padj.g < 0.05 &
         (!gene %in% liver_males_plast_GxE_list)) # this list of genes (liver_males_plast_GxE_list) has been generated in ./code/plot_males_parents_GxE.R
#nrow(g_males_liver) 5,422 at FDR < 0.05
#effectsize_g_males_liver <- g_males_liver %>% summarise(effectsize = mean(abslfc.g)) # 0.802

# How many sig. DEG for env only, with no genotype influence?
e_males_liver <- mergeall_males_liver %>%
  filter(padj.e < 0.05 &
         (!gene %in% liver_males_plast_GxE_list))
#nrow(e_males_liver) 553 at FDR < 0.05
#effectsize_e_males_liver <- e_males_liver %>% summarise(effectsize = mean(abslfc.e)) # 0.456

# Double-check that number of sig. DEG for GxE is same as before (in ./code/plot_parental_GxE.R)
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
#nrow(sig_males_liver) 6,506 genes
  
## plot parental expression patterns of parents as a ternary plot
liver_males_tern <- sig_males_liver %>%
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
ggsave("results/figures/males_GxE_ternary_liver.pdf", plot = liver_males_tern,
       height = 3, width = 3.5)


##############################################################
# Male BAT expression - effects of G, E, and GxE
##############################################################

## Genotype
# filter to get mean of 10 reads across all samples
DE_base_geno_males_BAT <- res_males_BAT_G %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## Environment
# filter to get mean of 10 reads across all samples
DE_base_env_males_BAT <- res_males_BAT_E %>%
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
merge_ge_males_BAT <- dplyr::full_join(DE_base_geno_males_BAT, DE_base_env_males_BAT, by = "gene") %>%
  dplyr::select("gene",
         "lfc.g" = "log2FoldChange.x",
         "padj.g" = "padj.x",
         "abslfc.g" = "abslfc.x",
         "lfc.e" = "log2FoldChange.y",
         "padj.e" = "padj.y",
         "abslfc.e" = "abslfc.y"
  )
#nrow(merge_ge_males_BAT) 14,176 genes

mergeall_males_BAT <- dplyr::left_join(merge_ge_males_BAT, DE_base_GxE_males_BAT, by = "gene") %>%
  dplyr::select("gene",
         "lfc.g", "padj.g", "abslfc.g",
         "lfc.e", "padj.e", "abslfc.e",
         "lfc.gxe" = "log2FoldChange",
         "padj.gxe" = "padj",
         "abslfc.gxe" = "abslfc"
  )
#nrow(mergeall_males_BAT) 14,176 genes


# How many sig. DEG for genotype only, with no environmental influence?
g_males_BAT <- mergeall_males_BAT %>%
  filter(padj.g < 0.05 &
           (!gene %in% BAT_males_plast_GxE_list)) # this list of genes (BAT_males_plast_GxE_list) has been generated in ./code/plot_males_parents_GxE.R
#nrow(g_males_BAT) 4,922 at FDR < 0.05
#effectsize_g_males_BAT <- g_males_BAT %>% summarise(effectsize = mean(abslfc.g)) # 0.791

# How many sig. DEG for env only, with no genotype influence?
e_males_BAT <- mergeall_males_BAT %>%
  filter(padj.e < 0.05 &
           (!gene %in% BAT_males_plast_GxE_list))
#nrow(e_males_BAT) 1,180 at FDR < 0.05
#effectsize_e_males_BAT <- e_males_BAT %>% summarise(effectsize = mean(abslfc.e)) # 0.479

# # Double-check that number of sig. DEG for GxE is same as before (in ./code/plot_parental_GxE.R)
#nrow(BAT_males_plast_GxE) # yes
# 1,425 at FDR < 0.05


## assign each gene to a main category for plotting purposes

sig_males_BAT <- mergeall_males_BAT %>%
  filter((padj.g < 0.05) | (padj.e < 0.05) | (gene %in% BAT_males_plast_GxE_list)) %>%
  mutate(main_category = case_when(padj.g < 0.05 & (!gene %in% BAT_males_plast_GxE_list) & padj.e >= 0.05 ~ 1, # g_liver = 1
                                   padj.e < 0.05 & (!gene %in% BAT_males_plast_GxE_list) & padj.g >= 0.05 ~ 2, # e_liver = 1
                                   gene %in% BAT_males_plast_GxE_list ~ 3)) %>% # gxe_liver = 3
  mutate(totalvar = (abslfc.g + abslfc.e + abslfc.gxe)) %>%
  mutate(g_propvar = (abslfc.g/totalvar),
         e_propvar = (abslfc.e/totalvar),
         gxe_propvar = (abslfc.gxe/totalvar))
#nrow(sig_males_BAT) 6,971 genes

### plot parental expression patterns of parents - BAT
BAT_males_tern <- sig_males_BAT %>%
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
  # geom_density_tern(color="#FCECAA",linetype=1, alpha=1, size = 01, n = 100, bins = 5,
  #                   stat = "DensityTern", base = "ilr", bdl = 0.010, na.rm = FALSE)
ggsave("results/figures/males_GxE_ternary_BAT.pdf", plot = BAT_males_tern,
       height = 3, width = 3.5)




##############################################################
# Female liver expression - effects of G, E, and GxE
##############################################################

## Genotype
# filter to get mean of 10 reads across all samples (per gene)
DE_base_geno_females_liver <- res_females_liver_G %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## Environment
# filter to get mean of 10 reads across all samples (per gene)
DE_base_env_females_liver <- res_females_liver_E %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## GxE
# filter to get mean of 10 reads across all samples
DE_base_GxE_females_liver <- res_females_liver_GxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## merge
merge_ge_females_liver <- dplyr::full_join(DE_base_geno_females_liver, DE_base_env_females_liver, by = "gene") %>%
  dplyr::select("gene",
                "lfc.g" = "log2FoldChange.x",
                "padj.g" = "padj.x",
                "abslfc.g" = "abslfc.x",
                "lfc.e" = "log2FoldChange.y",
                "padj.e" = "padj.y",
                "abslfc.e" = "abslfc.y"
  )
#nrow(merge_ge_females_liver) 14,163 genes

mergeall_females_liver <- dplyr::full_join(merge_ge_females_liver, DE_base_GxE_females_liver, by = "gene") %>%
  dplyr::select("gene",
                "lfc.g", "padj.g", "abslfc.g",
                "lfc.e", "padj.e", "abslfc.e",
                "lfc.gxe" = "log2FoldChange",
                "padj.gxe" = "padj",
                "abslfc.gxe" = "abslfc"
  )
#nrow(mergeall_females_liver) 14,163 genes


## effect sizes

# How many sig. DEG for genotype only, with no environmental influence?
g_females_liver <- mergeall_females_liver %>%
  filter(padj.g < 0.05 &
           (!gene %in% liver_females_plast_GxE_list)) # this list of genes (liver_males_plast_GxE_list) has been generated in ./code/plot_males_parents_GxE.R
#nrow(g_females_liver) 4,968 at FDR < 0.05
#effectsize_g_females_liver <- g_females_liver %>% summarise(effectsize = mean(abslfc.g)) # 0.854

# How many sig. DEG for env only, with no genotype influence?
e_females_liver <- mergeall_females_liver %>%
  filter(padj.e < 0.05 &
           (!gene %in% liver_females_plast_GxE_list))
#nrow(e_females_liver) 607 at FDR < 0.05
#effectsize_e_females_liver <- e_females_liver %>% summarise(effectsize = mean(abslfc.e)) # 0.417

# Double-check that number of sig. DEG for GxE is same as before (in ./code/plot_parental_GxE.R)
#nrow(liver_females_plast_GxE) # yes
# 649 at FDR < 0.05


## assign each gene to a main category for plotting purposes

sig_females_liver <- mergeall_females_liver %>%
  filter((padj.g < 0.05) | (padj.e < 0.05) | (gene %in% liver_females_plast_GxE_list)) %>%
  mutate(main_category = case_when(padj.g < 0.05 & (!gene %in% liver_females_plast_GxE_list) ~ 1, # g_liver = 1
                                   padj.e < 0.05 & (!gene %in% liver_females_plast_GxE_list) ~ 2, # e_liver = 1
                                   gene %in% liver_females_plast_GxE_list ~ 3)) %>% # gxe_liver = 3
  mutate(totalvar = (abslfc.g + abslfc.e + abslfc.gxe)) %>%
  mutate(g_propvar = (abslfc.g/totalvar),
         e_propvar = (abslfc.e/totalvar),
         gxe_propvar = (abslfc.gxe/totalvar))
#nrow(sig_females_liver) 5,957 genes

## plot parental expression patterns of parents as a ternary plot
liver_females_tern <- sig_females_liver %>%
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
ggsave("results/figures/females_GxE_ternary_liver.pdf", plot = liver_females_tern,
       height = 3, width = 3.5)


##############################################################
# Female BAT expression - effects of G, E, and GxE
##############################################################

## Genotype
# filter to get mean of 10 reads across all samples
DE_base_geno_females_BAT <- res_females_BAT_G %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## Environment
# filter to get mean of 10 reads across all samples
DE_base_env_females_BAT <- res_females_BAT_E %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

## GxE
# filter to get mean of 10 reads across all samples
DE_base_GxE_females_BAT <- res_females_BAT_GxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  mutate(abslfc = abs(log2FoldChange))

# merge
merge_ge_females_BAT <- dplyr::full_join(DE_base_geno_females_BAT, DE_base_env_females_BAT, by = "gene") %>%
  dplyr::select("gene",
                "lfc.g" = "log2FoldChange.x",
                "padj.g" = "padj.x",
                "abslfc.g" = "abslfc.x",
                "lfc.e" = "log2FoldChange.y",
                "padj.e" = "padj.y",
                "abslfc.e" = "abslfc.y"
  )
#nrow(merge_ge_females_BAT) 14,261 genes

mergeall_females_BAT <- dplyr::left_join(merge_ge_females_BAT, DE_base_GxE_females_BAT, by = "gene") %>%
  dplyr::select("gene",
                "lfc.g", "padj.g", "abslfc.g",
                "lfc.e", "padj.e", "abslfc.e",
                "lfc.gxe" = "log2FoldChange",
                "padj.gxe" = "padj",
                "abslfc.gxe" = "abslfc"
  )
#nrow(mergeall_females_BAT) 14,261 genes


# How many sig. DEG for genotype only, with no environmental influence?
g_females_BAT <- mergeall_females_BAT %>%
  filter(padj.g < 0.05 &
           (!gene %in% BAT_females_plast_GxE_list)) # this list of genes (BAT_males_plast_GxE_list) has been generated in ./code/plot_males_parents_GxE.R
#nrow(g_females_BAT) 4,420 at FDR < 0.05
#effectsize_g_females_BAT <- g_females_BAT %>% summarise(effectsize = mean(abslfc.g)) # 0.811

# How many sig. DEG for environment only, with no genotype influence?
e_females_BAT <- mergeall_females_BAT %>%
  filter(padj.e < 0.05 &
           (!gene %in% BAT_females_plast_GxE_list))
#nrow(e_females_BAT) 867 at FDR < 0.05
#effectsize_e_females_BAT <- e_females_BAT %>% summarise(effectsize = mean(abslfc.e)) # 0.499

# # Double-check that number of sig. DEG for GxE is same as before (in ./code/plot_parental_GxE.R)
#nrow(BAT_females_plast_GxE) # yes
# 1,425 at FDR < 0.05


## assign each gene to a main category for plotting purposes

sig_females_BAT <- mergeall_females_BAT %>%
  filter((padj.g < 0.05) | (padj.e < 0.05) | (gene %in% BAT_females_plast_GxE_list)) %>%
  mutate(main_category = case_when(padj.g < 0.05 & (!gene %in% BAT_females_plast_GxE_list) & padj.e >= 0.05 ~ 1, # g_liver = 1
                                   padj.e < 0.05 & (!gene %in% BAT_females_plast_GxE_list) & padj.g >= 0.05 ~ 2, # e_liver = 1
                                   gene %in% BAT_females_plast_GxE_list ~ 3)) %>% # gxe_liver = 3
  mutate(totalvar = (abslfc.g + abslfc.e + abslfc.gxe)) %>%
  mutate(g_propvar = (abslfc.g/totalvar),
         e_propvar = (abslfc.e/totalvar),
         gxe_propvar = (abslfc.gxe/totalvar))
#nrow(sig_females_BAT) 6,559 genes

### plot parental expression patterns of parents - BAT
BAT_females_tern <- sig_females_BAT %>%
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
# geom_density_tern(color="#FCECAA",linetype=1, alpha=1, size = 01, n = 100, bins = 5,
#                   stat = "DensityTern", base = "ilr", bdl = 0.010, na.rm = FALSE)
ggsave("results/figures/females_GxE_ternary_BAT.pdf", plot = BAT_females_tern,
       height = 3, width = 3.5)

