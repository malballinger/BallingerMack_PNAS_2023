#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots patterns of concordant (same direction) and discordant (opposite direction) plasticity
# in genes showing cis-regulatory divergence (ASE). This re-analysis allows us to ask about independence regarding
# plasticity and genetic divergence.
# This script generates Figures S9 in BallingerMack_2023.

# Main Result: 

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clear R's environment
library(tidyverse)
library(DESeq2)
library(ggtext)
library(car)
library(performance)
library(glue)
library(here)

set.seed(1991118)

source("./code/postprocess_RNAseq/get_DE_parents_analysis.R") # where DESeq data frames are generated

##############################################################
# Male liver expression
##############################################################

## Genes showing divergence
DE_base_evol_liver_males <- res_warm_males_liver_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_evol_liver_males) # 14,028 genes

## Genes showing plasticity in Brazil
DE_base_BZ_plast_liver_males <- res_BZ_males_liver_CvW %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_plast_liver_males)   # 14,028 genes

## merge evol div to Brazil plast (shared genes)
merge_EvP_liver_males <- dplyr::inner_join(DE_base_evol_liver_males, DE_base_BZ_plast_liver_males, by = "gene") %>%
  dplyr::select("gene",
         "lfc.e" = "log2FoldChange.x",
         "padj.e" = "padj.x",
         "lfc.p" = "log2FoldChange.y",
         "padj.p" = "padj.y")
#nrow(merge_EvP_liver_males) 14,028 genes


## add a column indicating whether evolved and plastic expression are in the same or opposite direction
merge_EvP_liver_males_direction <-  merge_EvP_liver_males %>%
  mutate(direction = case_when((sign(lfc.e)) == (sign(lfc.p)) ~ 1, # same direction (e.g., higher in NY warm and BZ cold)
                               (sign(lfc.e)) != (sign(lfc.p)) ~ 0)) # opposite direction

## pull out categories of genes that we are interested in
# FDR < 0.05
sigevol_liver_males_05 = dplyr::filter(merge_EvP_liver_males_direction, padj.e < 0.05 & padj.p >= 0.05) # 3,864
sigplas_liver_males_05 = dplyr::filter(merge_EvP_liver_males_direction, padj.p < 0.05 & padj.e >= 0.05) # 410
sigsame_liver_males_05 = dplyr::filter(merge_EvP_liver_males_direction, direction == 1 & padj.e < 0.05 & padj.p < 0.05) # 296
sigopps_liver_males_05 = dplyr::filter(merge_EvP_liver_males_direction, direction == 0 & padj.e < 0.05 & padj.p < 0.05) # 111
sigboth_liver_males_05 = dplyr::filter(merge_EvP_liver_males_direction, padj.e < 0.05 & padj.p < 0.05) # 407



### genes that are ASE outliers (NHVT), how many show adaptive and non-adaptive plasticity?
# get list of ASEoutliers
ASE_outlier <- read_delim(here("./data/processed/ASEoutlier_NHVT.txt"))
ASE_outlier_list <- ASE_outlier %>% pull(gene)

sigsame_ASE_outlier_liver_males_05 <- sigsame_liver_males_05 %>% filter(gene %in% ASE_outlier_list) # 4
sigopps_ASE_outlier_liver_males_05 <- sigopps_liver_males_05 %>% filter(gene %in% ASE_outlier_list) # 5

merge_EvP_liver_males_direction %>% filter(gene %in% ASE_outlier_list)

# plot candidate genes
#res_warm_males_liver_NYvsBZ %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% View()
level_order <- c("Warm", "Cold")
sigopps_gene <- plotCounts(deseq_males_liver, gene = "ENSMUSG00000071176", intgroup = c("temperature", "population"), returnData = TRUE)
ggplot(data=sigopps_gene, aes(x = factor(temperature, level = level_order), y=count, color = population)) + #, group = population)) +
geom_boxplot(position = "identity", fill = NA) + labs(x = "Temp", y = "Normalized Counts")

## sigsame, where NY is canalized
# ENSMUSG00000026718
# ENSMUSG00000040272
## sigsame, where both BZ and NY are same plast direction
# ENSMUSG00000032279
# ENSMUSG00000074264 (amy1)

## sigopps, where NY is canalized
# ENSMUSG00000001119

## sany sigopps show higher warm BZ and higher cold BZ (but NY shows same increased plasticity)


### genes that are regulated in cis, what patterns of expression plasticity do they exhibit?
## cis-effects in liver (warm)
liver_warm_male_metadata <- read_delim(here("data/processed/Liver.MALE.WARM.categories.forplot.txt")) %>%
  select("gene" = "...1",
         "lfc.parental" = "log2FoldChange...2",
         "padj.parental" = "padj...3",
         "lfc.F1" = "log2FoldChange...4",
         "padj.F1" = "padj...5",
         "padj.PvF1" = "padj...6",
         "category" = "Amb&Conserved")

cis_liver_warm_male = dplyr::filter(liver_warm_male_metadata, category == "CIS_ONLY") # 431 genes

## cis-effects in liver (cold)
liver_cold_male_metadata <- read_delim(here("data/processed/Liver.MALE.COLD.categories.forplot.txt")) %>%
  select("gene" = "...1",
         "lfc.parental" = "log2FoldChange...2",
         "padj.parental" = "padj...3",
         "lfc.F1" = "log2FoldChange...4",
         "padj.F1" = "padj...5",
         "padj.PvF1" = "liverpvals_coldmale_adj",
         "category" = "Amb&Conserved")

cis_liver_cold_male = dplyr::filter(liver_cold_male_metadata, category == "CIS_ONLY") # 450 genes

## merge
merge_cis_warm_cold_liver <- dplyr::full_join(cis_liver_warm_male, cis_liver_cold_male, by = "gene") %>%
                                dplyr::select("gene") # 614 genes
# save a list of cis-genes
merge_cis_warm_cold_liver_list <- merge_cis_warm_cold_liver %>% pull(gene)

### filter male genes to only include cis-genes
cis_merge_EvP_liver_males_direction <- merge_EvP_liver_males_direction %>% filter(gene %in% merge_cis_warm_cold_liver_list) # 613 genes

## pull out categories of genes that we are interested in
# FDR < 0.05
sigevol_liver_males_cis_05 = dplyr::filter(cis_merge_EvP_liver_males_direction, padj.e < 0.05 & padj.p >= 0.05) # 474
sigplas_liver_males_cis_05 = dplyr::filter(cis_merge_EvP_liver_males_direction, padj.p < 0.05 & padj.e >= 0.05) # 3
sigsame_liver_males_cis_05 = dplyr::filter(cis_merge_EvP_liver_males_direction, direction == 1 & padj.e < 0.05 & padj.p < 0.05) # 36
sigopps_liver_males_cis_05 = dplyr::filter(cis_merge_EvP_liver_males_direction, direction == 0 & padj.e < 0.05 & padj.p < 0.05) # 28
sigboth_liver_males_cis_05 = dplyr::filter(cis_merge_EvP_liver_males_direction, padj.e < 0.05 & padj.p < 0.05) # 64


plast_outliers_liver <- cis_merge_EvP_liver_males_direction %>% filter(gene %in% ASE_outlier_list)

sigsame_BAT_males_ASE_05 = dplyr::filter(plast_outliers_liver, direction == 1 & padj.e < 0.05 & padj.p < 0.05) # 3
sigopps_BAT_males_ASE_05 = dplyr::filter(plast_outliers_liver, direction == 0 & padj.e < 0.05 & padj.p < 0.05) # 3



### plot adaptive and non-adaptive plasticity
liver_cis_adap_plasticity <- cis_merge_EvP_liver_males_direction %>%
  ggplot(aes(x = lfc.p, y = lfc.e)) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.25) +
  geom_point(data = sigsame_liver_males_cis_05, aes(x = lfc.p, y = lfc.e), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_point(data = sigopps_liver_males_cis_05, aes(x = lfc.p, y = lfc.e), fill = "#000000", color = "000000", stroke = 0.175, size = 1.5, shape = 21, alpha = 0.8, show.legend = TRUE) +
  scale_x_continuous(limits = c(-4,2.15), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-4,4.25), expand = c(0.01,0.075)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans")) +
  labs(x = "Log<sub>2</sub> Fold Change in BZ (cold vs warm)",
       y = "Log<sub>2</sub> Fold Change in Warm (NY vs BZ)")
#3BAF85
ggsave("results/figures/cis_adap_plast_liver.pdf", plot = liver_cis_adap_plasticity, height = 2, width = 2.1)


### Correlation Test
##### liver
# check normality of liver data
cis_liver.model_males_full <- lm(lfc.e ~ lfc.p, data = sigboth_liver_males_cis_05)
#check_model(cis_liver.model_males) # data are not normally distributed - use Spearman correlation

r.obs_liver_males_cis_05 <- cor(x = sigboth_liver_males_cis_05$lfc.p, y = sigboth_liver_males_cis_05$lfc.e, method = "spearman")
P.cor_liver_males_cis_05 <- cor.test(x = sigboth_liver_males_cis_05$lfc.p, y = sigboth_liver_males_cis_05$lfc.e, method = "spearman")$p.value

# permutation function
cor.perm_liver_males_cis <- function (data = data, nperm = 10000)
{
  random.subsample_liver_males <- dplyr::slice_sample(data, n=64) # randomly subsample the number of sigboth
  r.per_liver_males <- sapply (1:nperm, FUN = function (i) cor (x = random.subsample_liver_males$lfc.p, y = sample (random.subsample_liver_males$lfc.e)))
  r.per_liver_males <- c(r.per_liver_males, r.obs_liver_males_cis_05)
  hist (r.per_liver_males, breaks = 100, xlab = "Correlation coefficient", main=NULL,
        col = "white")
  abline (v = r.obs_liver_males_cis_05, col = 'red', lty = "dotted")
  P.per_liver_males <- sum (abs (r.per_liver_males) >= abs (r.obs_liver_males_cis_05))/(nperm + 1) 
  return (list (r.obs_liver_males_cis_05 = r.obs_liver_males_cis_05, P.cor_liver_males_cis_05 = P.cor_liver_males_cis_05,
                P.per_liver_males = P.per_liver_males, r.per_liver_males = r.per_liver_males))
}

liver_perm_males_cis_05 <- cor.perm_liver_males_cis(data = merge_EvP_liver_males) # test using all cis genes
liver_p_perm_males_cis_05 <- liver_perm_males_cis_05$P.per_liver_males # P = 0.04 (barely)

# make liver_perm a df so you can plot hist in ggplot
df_liver_perm_males_cis_05 <- data.frame(liver_perm_males_cis_05)
# make histogram
liver_perm_males_cis_05.inset <- df_liver_perm_males_cis_05 %>%
  ggplot(aes(x=r.per_liver_males)) +
  geom_histogram(fill = "white", color = "gray", bins = 100, binwidth = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 400)) +
  scale_x_continuous(expand = c(0,0), limits = c(NA, 0.45)) +
  geom_vline(xintercept = r.obs_liver_males_cis_05, color = "#E09832", linetype = "solid", size = 0.5) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 5, family = "sans"),
        axis.text.y = element_text(size = 5, family = "sans"),
        plot.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Correlation Coefficient",
       y = "Count")
ggsave("results/figures/males_cis_05_liver_inset.pdf", plot = liver_perm_males_cis_05.inset,
       height = 1.5, width = 2.25, units = "cm")


### determine if any of these categories are ASE outliers

### genes that are ASE outliers (NHVT), how many show adaptive and non-adaptive plasticity?
# get list of ASEoutliers
ASE_outlier <- read_delim(here("./data/processed/ASEoutlier_NHVT.txt"))
ASE_outlier_list <- ASE_outlier %>% pull(gene)

sigsame_ASE_outlier_liver_males_05 <- sigsame_liver_males_05 %>% filter(gene %in% ASE_outlier_list) # 4
sigopps_ASE_outlier_liver_males_05 <- sigopps_liver_males_05 %>% filter(gene %in% ASE_outlier_list) # 5

plast_outliers <- merge_EvP_liver_males_direction %>% filter(gene %in% ASE_outlier_list)

sigsame_liver_males_ASE_05 = dplyr::filter(plast_outliers, direction == 1 & padj.e < 0.05 & padj.p < 0.05) # 4
sigopps_liver_males_ASE_05 = dplyr::filter(plast_outliers, direction == 0 & padj.e < 0.05 & padj.p < 0.05) # 5
sigboth_liver_males_ASE_05 = dplyr::filter(plast_outliers, padj.e < 0.05 & padj.p < 0.05) # 9





##############################################################
# Male BAT expression
##############################################################

### New York vs Brazil (in the warm -- evolved differences)
# filter to get mean of 10 reads across all samples
DE_base_evol_BAT_males <- res_warm_males_BAT_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_evol_BAT_males) # 14,176 genes

### Brazil plasticity (cold vs warm):
# filter to get mean of 10 reads across all samples
DE_base_BZ_plast_BAT_males <- res_BZ_males_BAT_CvW %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_plast_BAT_males) # 14,176 genes

### merge evol DE to plast Brazil DE in BAT
merge_EvP_BAT_males <- dplyr::inner_join(DE_base_evol_BAT_males, DE_base_BZ_plast_BAT_males, by = "gene") %>%
  dplyr::select("gene",
         "lfc.e" = "log2FoldChange.x",
         "padj.e" = "padj.x",
         "lfc.p" = "log2FoldChange.y",
         "padj.p" = "padj.y")
#nrow(merge_EvP_BAT_males) # 14,176 genes


## add a column indicating whether evolved and plastic expression are in the same or opposite direction
merge_EvP_BAT_males_direction <-  merge_EvP_BAT_males %>%
  mutate(direction = case_when((sign(lfc.e)) == (sign(lfc.p)) ~ 1, # same direction (e.g., higher in NY warm and BZ cold)
                               (sign(lfc.e)) != (sign(lfc.p)) ~ 0)) # opposite direction

## pull out categories of genes that we are interested in
# FDR < 0.05
sigevol_BAT_males_05 = dplyr::filter(merge_EvP_BAT_males_direction, padj.e < 0.05 & padj.p >= 0.05) # 3,215
sigplas_BAT_males_05 = dplyr::filter(merge_EvP_BAT_males_direction, padj.p < 0.05 & padj.e >= 0.05) # 836
sigsame_BAT_males_05 = dplyr::filter(merge_EvP_BAT_males_direction, direction == 1 & padj.e < 0.05 & padj.p < 0.05) # 295
sigopps_BAT_males_05 = dplyr::filter(merge_EvP_BAT_males_direction, direction == 0 & padj.e < 0.05 & padj.p < 0.05) # 2013
sigboth_BAT_males_05 = dplyr::filter(merge_EvP_BAT_males_direction, padj.e < 0.05 & padj.p < 0.05) # 498


## genes with cis-effects in BAT (warm)

BAT_warm_male_metadata <- read_delim(here("data/processed/BAT.MALE.WARM.categories.forplot.txt")) %>%
  select("gene" = "...1",
         "lfc.parental" = "log2FoldChange...2",
         "padj.parental" = "padj...3",
         "lfc.F1" = "log2FoldChange...4",
         "padj.F1" = "padj...5",
         "padj.PvF1" = "batpadj",
         "category" = "Amb&Conserved")

cis_BAT_warm_male = dplyr::filter(BAT_warm_male_metadata, category == "CIS_ONLY") # 452 genes
#cistrans_warm_BAT = dplyr::filter(BAT_warm_male_metadata, category == "Cis&Trans") #378 genes

#merge_cis_cistrans_warm_BAT <- dplyr::full_join(cis_warm_BAT, cistrans_warm_BAT, by = "gene") %>%
#                                dplyr::select("gene") # 830 genes

## genes with cis-effects in BAT (cold)

BAT_cold_male_metadata <- read_delim(here("data/processed/BAT.MALE.COLD.categories.forplot.txt")) %>%
  select("gene" = "...1",
         "lfc.parental" = "log2FoldChange...2",
         "padj.parental" = "padj...3",
         "lfc.F1" = "log2FoldChange...4",
         "padj.F1" = "padj...5",
         "padj.PvF1" = "batpvals_coldmale_adj",
         "category" = "Amb&Conserved")

cis_BAT_cold_male = dplyr::filter(BAT_cold_male_metadata, category == "CIS_ONLY") # 478 genes
#cistrans_cold_BAT = dplyr::filter(BAT_cold_male_metadata, category == "Cis&Trans") #312 genes

#merge_cis_cistrans_cold_BAT <- dplyr::full_join(cis_cold_BAT, cistrans_cold_BAT, by = "gene") %>%
#                                dplyr::select("gene") # 790 genes


## merge
merge_cis_warm_cold_BAT <- dplyr::full_join(cis_BAT_warm_male, cis_BAT_cold_male, by = "gene")

# save a list of cis-genes
merge_cis_warm_cold_BAT_list <- merge_cis_warm_cold_BAT %>% pull(gene)

### filter male genes to only include cis-genes
cis_merge_EvP_BAT_males_direction <- merge_EvP_BAT_males_direction %>% filter(gene %in% merge_cis_warm_cold_BAT_list) # 660 genes

## pull out categories of genes that we are interested in
## (FDR < 0.05)
sigevol_BAT_males_cis_05 = dplyr::filter(cis_merge_EvP_BAT_males_direction, padj.e < 0.05 & padj.p >= 0.05) # 464 genes
sigplas_BAT_males_cis_05 = dplyr::filter(cis_merge_EvP_BAT_males_direction, padj.p < 0.05 & padj.e >= 0.05) # 23 genes
sigsame_BAT_males_cis_05 = dplyr::filter(cis_merge_EvP_BAT_males_direction, direction == 1 & padj.e < 0.05 & padj.p < 0.05) # 51 genes
sigopps_BAT_males_cis_05 = dplyr::filter(cis_merge_EvP_BAT_males_direction, direction == 0 & padj.e < 0.05 & padj.p < 0.05) # 41 genes
sigboth_BAT_males_cis_05 = dplyr::filter(cis_merge_EvP_BAT_males_direction, padj.e < 0.05 & padj.p < 0.05) # 92 genes


### plot adaptive and non-adaptive plasticity
BAT_cis_adap_plasticity <- cis_merge_EvP_BAT_males_direction %>%
  ggplot(aes(x = lfc.p, y = lfc.e)) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.25) +
  geom_point(data = sigsame_BAT_males_cis_05, aes(x = lfc.p, y = lfc.e), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_point(data = sigopps_BAT_males_cis_05, aes(x = lfc.p, y = lfc.e), fill = "#000000", color = "000000", stroke = 0.175, size = 1.5, shape = 21, alpha = 0.8, show.legend = TRUE) +
  scale_x_continuous(limits = c(-4,2.15), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-4,4.25), expand = c(0.01,0.075)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans")) +
  labs(x = "Log<sub>2</sub> Fold Change in BZ (cold vs warm)",
       y = "Log<sub>2</sub> Fold Change in Warm (NY vs BZ)")
#3BAF85
ggsave("results/figures/cis_adap_plast_BAT.pdf", plot = BAT_cis_adap_plasticity, height = 2, width = 2.1)



### Correlation Test
##### liver
# check normality of liver data
cis_BAT.model_males_05 <- lm(lfc.e ~ lfc.p, data = sigboth_BAT_males_cis_05)
#check_model(cis_BAT.model_males_full) # data are not normally distributed - use Spearman correlation

r.obs_BAT_males_cis_05 <- cor(x = sigboth_BAT_males_cis_05$lfc.p, y = sigboth_BAT_males_cis_05$lfc.e, method = "spearman")
P.cor_BAT_males_cis_05 <- cor.test(x = sigboth_BAT_males_cis_05$lfc.p, y = sigboth_BAT_males_cis_05$lfc.e, method = "spearman")$p.value

# permutation function
cor.perm_BAT_males_cis <- function (data = data, nperm = 10000)
{
  random.subsample_BAT_males <- dplyr::slice_sample(data, n=92) # randomly subsample the number of sigboth
  r.per_BAT_males <- sapply (1:nperm, FUN = function (i) cor (x = random.subsample_BAT_males$lfc.p, y = sample (random.subsample_BAT_males$lfc.e)))
  r.per_BAT_males <- c(r.per_BAT_males, r.obs_BAT_males_cis_05)
  hist (r.per_BAT_males, breaks = 100, xlab = "Correlation coefficient", main=NULL,
        col = "white")
  abline (v = r.obs_BAT_males_cis_05, col = 'red', lty = "dotted")
  P.per_BAT_males <- sum (abs (r.per_BAT_males) >= abs (r.obs_BAT_males_cis_05))/(nperm + 1) 
  return (list (r.obs_BAT_males_cis_05 = r.obs_BAT_males_cis_05, P.cor_BAT_males_cis_05 = P.cor_BAT_males_cis_05,
                P.per_BAT_males = P.per_BAT_males, r.per_BAT_males = r.per_BAT_males))
}

BAT_perm_males_cis_05 <- cor.perm_BAT_males_cis(data = merge_EvP_BAT_males) # test using all cis genes #cis_merge_EvP_BAT_males_direction
BAT_p_perm_males_cis_05 <- BAT_perm_males_cis_05$P.per_BAT_males # P = 0.7

# make BAT_perm a df so you can plot hist in ggplot
df_BAT_perm_males_cis_05 <- data.frame(BAT_perm_males_cis_05)
# make histogram
BAT_perm_males_cis_05.inset <- df_BAT_perm_males_cis_05 %>%
  ggplot(aes(x=r.per_BAT_males)) +
  geom_histogram(fill = "white", color = "gray", bins = 100, binwidth = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 400)) +
  scale_x_continuous(expand = c(0,0), limits = c(NA, 0.45)) +
  geom_vline(xintercept = r.obs_BAT_males_cis_05, color = "#E09832", linetype = "solid", size = 0.5) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 5, family = "sans"),
        axis.text.y = element_text(size = 5, family = "sans"),
        plot.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Correlation Coefficient",
       y = "Count")
ggsave("results/figures/males_cis_05_BAT_inset.pdf", plot = BAT_perm_males_cis_05.inset,
       height = 1.5, width = 2.25, units = "cm")



## determine if any of these categories are ASE outliers
sigsame_ASE_outlier_BAT_males_05 <- sigsame_BAT_males_05 %>% filter(gene %in% ASE_outlier_list) # 4
sigopps_ASE_outlier_BAT_males_05 <- sigopps_BAT_males_05 %>% filter(gene %in% ASE_outlier_list) # 4


### genes that are ASE outliers (NHVT), how many show adaptive and non-adaptive plasticity?
# get list of ASEoutliers
ASE_outlier <- read_delim(here("./data/processed/ASEoutlier_NHVT.txt"))
ASE_outlier_list <- ASE_outlier %>% pull(gene)


plast_outliers_BAT <- cis_merge_EvP_BAT_males_direction %>% filter(gene %in% ASE_outlier_list)

sigsame_BAT_males_ASE_05 = dplyr::filter(plast_outliers_BAT, direction == 1 & padj.e < 0.05 & padj.p < 0.05) # 1
sigopps_BAT_males_ASE_05 = dplyr::filter(plast_outliers_BAT, direction == 0 & padj.e < 0.05 & padj.p < 0.05) # 3
sigboth_BAT_males_ASE_05 = dplyr::filter(plast_outliers_BAT, padj.e < 0.05 & padj.p < 0.05) # 8



















#### BAT Females
### New York vs Brazil (in the warm -- evolved differences)
# filter to get mean of 10 reads across all samples
DE_base_evol_BAT_females <- res_warm_females_BAT_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_evol_BAT_females) # 14,261 genes

### Brazil plasticity (cold vs warm):
# filter to get mean of 10 reads across all samples
DE_base_BZ_plast_BAT_females <- res_BZ_females_BAT_CvW %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_plast_BAT_females) # 14,261 genes

### merge evol DE to plast Brazil DE in BAT
merge_EvP_BAT_females <- dplyr::inner_join(DE_base_evol_BAT_females, DE_base_BZ_plast_BAT_females, by = "gene") %>%
  dplyr::select("gene",
         "lfc.e" = "log2FoldChange.x",
         "padj.e" = "padj.x",
         "lfc.p" = "log2FoldChange.y",
         "padj.p" = "padj.y")
#nrow(merge_EvP_BAT_females) # 14,261 genes


## genes with cis-effects in BAT (cold)

BAT_cold_female_metadata <- read_table(here("data/processed/BAT.FEMALE.COLD.categories.forplot.txt"))

cis_BAT_cold_female = dplyr::filter(BAT_cold_female_metadata, category == "CIS_ONLY") # 365 genes

# save a list of cis-genes
merge_cis_BAT_cold_list <- cis_BAT_cold_female %>% pull(gene)

### filter female genes to only include cis-genes
cis_merge_EvP_BAT_females <- merge_EvP_BAT_females %>% filter(gene %in% merge_cis_BAT_cold_list) # 365 genes


## add a column indicating whether evolved and plastic expression are in the same or opposite direction
cis_merge_EvP_BAT_females_direction <-  cis_merge_EvP_BAT_females %>%
  mutate(direction = case_when((sign(lfc.e)) == (sign(lfc.p)) ~ 1, # same direction (e.g., higher in NY warm and BZ cold)
                               (sign(lfc.e)) != (sign(lfc.p)) ~ 0)) # opposite direction

## pull out categories of genes that we are interested in
## (FDR < 0.05)
sigevol_BAT_females_cis_05 = dplyr::filter(cis_merge_EvP_BAT_females_direction, padj.e < 0.05 & padj.p >= 0.05) # 244 genes
sigplas_BAT_females_cis_05 = dplyr::filter(cis_merge_EvP_BAT_females_direction, padj.p < 0.05 & padj.e >= 0.05) # 18 genes
sigsame_BAT_females_cis_05 = dplyr::filter(cis_merge_EvP_BAT_females_direction, direction == 1 & padj.e < 0.05 & padj.p < 0.05) # 21 genes
sigopps_BAT_females_cis_05 = dplyr::filter(cis_merge_EvP_BAT_females_direction, direction == 0 & padj.e < 0.05 & padj.p < 0.05) # 31 genes
sigboth_BAT_females_cis_05 = dplyr::filter(cis_merge_EvP_BAT_females_direction, padj.e < 0.05 & padj.p < 0.05) # 52 genes


### Correlation Test
##### liver
# check normality of liver data
cis_BAT.model_females <- lm(lfc.e ~ lfc.p, data = sigboth_BAT_females_cis_05)
#check_model(cis_BAT.model_females) # data are not normally distributed - use Spearman correlation

r.obs_BAT_females_cis_05 <- cor(x = sigboth_BAT_females_cis_05$lfc.p, y = sigboth_BAT_females_cis_05$lfc.e, method = "spearman")
P.cor_BAT_females_cis_05 <- cor.test(x = sigboth_BAT_females_cis_05$lfc.p, y = sigboth_BAT_females_cis_05$lfc.e, method = "spearman")$p.value

# permutation function
cor.perm_BAT_females_cis <- function (data = data, nperm = 10000)
{
  random.subsample_BAT_females <- dplyr::slice_sample(data, n=52) # randomly subsample the number of sigboth
  r.per_BAT_females <- sapply (1:nperm, FUN = function (i) cor (x = random.subsample_BAT_females$lfc.p, y = sample (random.subsample_BAT_females$lfc.e)))
  r.per_BAT_females <- c(r.per_BAT_females, r.obs_BAT_females_cis_05)
  hist (r.per_BAT_females, breaks = 100, xlab = "Correlation coefficient", main=NULL,
        col = "white")
  abline (v = r.obs_BAT_females_cis_05, col = 'red', lty = "dotted")
  P.per_BAT_females <- sum (abs (r.per_BAT_females) >= abs (r.obs_BAT_females_cis_05))/(nperm + 1) 
  return (list (r.obs_BAT_females_cis_05 = r.obs_BAT_females_cis_05, P.cor_BAT_females_cis_05 = P.cor_BAT_females_cis_05,
                P.per_BAT_females = P.per_BAT_females, r.per_BAT_females = r.per_BAT_females))
}

BAT_perm_females_cis_05 <- cor.perm_BAT_females_cis(data = cis_merge_EvP_BAT_females) # test using all of the genes
BAT_p_perm_females_cis <- BAT_perm_females_cis_05$P.per_BAT_females # P = 0.172 (barely)

# make BAT_perm a df so you can plot hist in ggplot
df_BAT_perm_males_cis_full_05 <- data.frame(BAT_perm_males_cis_full_05)
# make histogram
BAT_perm_males_cis_full_05.inset <- df_BAT_perm_males_cis_full_05 %>%
  ggplot(aes(x=r.per_BAT_males)) +
  geom_histogram(fill = "white", color = "gray", bins = 100, binwidth = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 400)) +
  scale_x_continuous(expand = c(0,0), limits = c(NA, 0.45)) +
  geom_vline(xintercept = r.obs_BAT_males_cis_full_05, color = "#E09832", linetype = "solid", size = 0.5) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 5, family = "sans"),
        axis.text.y = element_text(size = 5, family = "sans"),
        plot.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Correlation Coefficient",
       y = "Count")
ggsave("results/figures/males_cis_full_05_BAT_inset.pdf", plot = BAT_perm_males_cis_full_05.inset,
       height = 1.5, width = 2.25, units = "cm")



### how many show adaptive and non-adaptive plasticity?
# get list of ASEoutliers

sigsame_ASE_outlier_BAT_males_05 <- sigsame_BAT_males_05 %>% filter(gene %in% ASE_outlier_list) # 4
sigopps_ASE_outlier_BAT_males_05 <- sigopps_BAT_males_05 %>% filter(gene %in% ASE_outlier_list) # 4










##############################################################
# Statistical Analysis - Binomial Test
##############################################################

## males

adap_BAT_males <- sigboth_BAT_males_full_cis_05 %>%
  dplyr::select(gene, direction) %>%
  mutate(adap = ifelse(direction > 0, TRUE, FALSE))

adap_liver_males_cis <- sigboth_liver_males_full_cis_05 %>%
  dplyr::select(gene, direction) %>%
  mutate(adap = ifelse(direction > 0, TRUE, FALSE))

# BAT
dim(adap_BAT_males)[1] # get number of genes in this category (92)
sum(adap_BAT_males$adap) # the number of genes where expression is higher in Brazil (51)
btest_BAT_males = binom.test(x=sum(adap_BAT_males$adap), n = dim(adap_BAT_males)[1])
# P < 0.001

# liver
dim(adap_liver_males_cis)[1] # get number of genes in this category (64)
sum(adap_liver_males_cis$adap) # the number of genes where expression is higher in Brazil (36)
btest_liver_males_cis = binom.test(x=sum(adap_liver_males_cis$adap), n = dim(adap_liver_males_cis)[1])
# P = 0.5625


## females

adap_BAT_females <- sigboth_BAT_females %>%
  dplyr::select(gene, direction) %>%
  mutate(adap = ifelse(direction > 0, TRUE, FALSE))

adap_liver_females <- sigboth_liver_females %>%
  dplyr::select(gene, direction) %>%
  mutate(adap = ifelse(direction > 0, TRUE, FALSE))

# BAT
dim(adap_BAT_females)[1] # get number of genes in this category (475)
sum(adap_BAT_females$adap) # the number of genes where expression is higher in Brazil (302)
btest_BAT_females = binom.test(x=sum(adap_BAT_females$adap), n = dim(adap_BAT_females)[1])
# P < 0.001

# liver
dim(adap_liver_females)[1] # get number of genes in this category (198)
sum(adap_liver_females$adap) # the number of genes where expression is higher in Brazil (142)
btest_liver_females = binom.test(x=sum(adap_liver_females$adap), n = dim(adap_liver_females)[1])
# P < 0.001


