#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots patterns of adaptive (same direction) and non-adaptive (opposite direction) plasticity
# in parental samples. Specifically, we ask if the plastic response of Brazil mice goes in the same direction
# as the evolved response of New York mice at room temperature.
# This script generates Figures 2C (males) and S3D (females) in BallingerMack_2022.

# Main Result: across both tissues and sexes, Brazil mice exhibit patterns of adaptive plasticity

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

set.seed(1991118)

source("./code/postprocess_RNAseq/get_DE_parents_analysis.R") # where DESeq data frames are generated

##############################################################
# Male liver expression (DE in warm (NY vs BZ), and plasticity in Brazil (BZ warm vs cold))
##############################################################

### New York vs Brazil (in the warm -- evolved differences)
# filter to get mean of 10 reads across all samples
DE_base_evol_liver_males <- res_warm_males_liver_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_evol_liver_males) # 14,028 genes

### Brazil plasticity (cold vs warm):
# filter to get mean of 10 reads across all samples
DE_base_BZ_plast_liver_males <- res_BZ_males_liver_CvW %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_plast_liver_males)   # 14,028 genes

### merge evol DE to plast Brazil DE in liver
merge_EvP_liver_males <- dplyr::full_join(DE_base_evol_liver_males, DE_base_BZ_plast_liver_males, by = "gene") %>%
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
sigevol_liver_males = dplyr::filter(merge_EvP_liver_males_direction, padj.e < 0.05 & padj.p >= 0.05) # 3,864 genes
sigplas_liver_males = dplyr::filter(merge_EvP_liver_males_direction, padj.p < 0.05 & padj.e >= 0.05) # 410 genes
sigsame_liver_males = dplyr::filter(merge_EvP_liver_males_direction, direction == 1 & padj.e < 0.05 & padj.p < 0.05) # 296 genes
sigopps_liver_males = dplyr::filter(merge_EvP_liver_males_direction, direction == 0 & padj.e < 0.05 & padj.p < 0.05) # 111 genes
sigboth_liver_males = dplyr::filter(merge_EvP_liver_males_direction, padj.e < 0.05 & padj.p < 0.05) # 407 genes

### plot adaptive and non-adaptive plasticity
liver_males.plot <- merge_EvP_liver_males_direction %>%
  ggplot(aes(x = lfc.p, y = lfc.e)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
  geom_point(data = sigsame_liver_males, aes(x = lfc.p, y = lfc.e), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_point(data = sigopps_liver_males, aes(x = lfc.p, y = lfc.e), fill = "#000000", color = "000000", stroke = 0.175, size = 1.5, shape = 21, alpha = 0.8, show.legend = TRUE) +
  scale_x_continuous(limits = c(-6.25,4.15), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-5,7.75), expand = c(0.01,0.075)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans")) +
  labs(x = "Log<sub>2</sub> Fold Change in BZ (cold vs warm)",
       y = "Log<sub>2</sub> Fold Change in Warm (NY vs BZ)")
#3BAF85
ggsave("results/figures/males_adap_plast_liver.pdf", plot = liver_males.plot, height = 2, width = 2.1)


### Correlation Test

# check normality of liver data
liver.model_males <- lm(lfc.e ~ lfc.p, data = sigboth_liver_males)
#check_model(liver.model_males) # data are not normally distributed - use Spearman correlation

r.obs_liver_males <- cor(x = sigboth_liver_males$lfc.p, y = sigboth_liver_males$lfc.e, method = "spearman")
P.cor_liver_males <- cor.test(x = sigboth_liver_males$lfc.p, y = sigboth_liver_males$lfc.e, method = "spearman")$p.value

# permutation function
cor.perm_liver_males <- function (data = data, nperm = 10000)
{
  random.subsample_liver_males <- dplyr::slice_sample(data, n=407) # randomly subsample the number of sigboth
  r.per_liver_males <- sapply (1:nperm, FUN = function (i) cor (x = random.subsample_liver_males$lfc.p, y = sample (random.subsample_liver_males$lfc.e)))
  r.per_liver_males <- c(r.per_liver_males, r.obs_liver_males)
  hist (r.per_liver_males, breaks = 100, xlab = "Correlation coefficient", main=NULL,
        col = "white")
  abline (v = r.obs_liver_males, col = 'red', lty = "dotted")
  P.per_liver_males <- sum (abs (r.per_liver_males) >= abs (r.obs_liver_males))/(nperm + 1) 
  return (list (r.obs_liver_males = r.obs_liver_males, P.cor_liver_males = P.cor_liver_males,
                P.per_liver_males = P.per_liver_males, r.per_liver_males = r.per_liver_males))
}

liver_perm_males <- cor.perm_liver_males(data = merge_EvP_liver_males) # test using all of the genes
liver_p_perm_males <- liver_perm_males$P.per_liver_males # P < 0.05

# make liver_perm a df so you can plot hist in ggplot
df_liver_perm_males <- data.frame(liver_perm_males)
# make histogram
liver_perm_males.inset <- df_liver_perm_males %>%
  ggplot(aes(x=r.per_liver_males)) +
  geom_histogram(fill = "white", color = "gray", bins = 100, binwidth = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 685)) +
  scale_x_continuous(expand = c(0,0), limits = c(NA, 0.45)) +
  geom_vline(xintercept = r.obs_liver_males, color = "#E09832", linetype = "solid", size = 0.5) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 5, family = "sans"),
        axis.text.y = element_text(size = 5, family = "sans"),
        plot.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Correlation Coefficient",
       y = "Count")
ggsave("results/figures/males_adap_plast_liver_inset.pdf", plot = liver_perm_males.inset,
       height = 1.5, width = 2.25, units = "cm")



##############################################################
# Male BAT expression (DE in warm (NY vs BZ), and plasticity in Brazil (BZ warm vs cold))
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
merge_EvP_BAT_males <- dplyr::full_join(DE_base_evol_BAT_males, DE_base_BZ_plast_BAT_males, by = "gene") %>%
  dplyr::select("gene",
         "lfc.e" = "log2FoldChange.x",
         "padj.e" = "padj.x",
         "lfc.p" = "log2FoldChange.y",
         "padj.p" = "padj.y")
#nrow(merge_EvP_BAT_males) # 14,176 genes

## add a column indicating whether evolved and plastic expression are in the same or opposite direction
merge_EvP_BAT_males_direction <-  merge_EvP_BAT_males %>%
  mutate(direction = case_when((sign(lfc.e)) == (sign(lfc.p)) ~ 1, # same direction
                               (sign(lfc.e)) != (sign(lfc.p)) ~ 0)) # opposite direction

## pull out categories of genes that we are interested in
sigevol_BAT_males = dplyr::filter(merge_EvP_BAT_males_direction, padj.e < 0.05 & padj.p >= 0.05) # 3,215 genes
sigplas_BAT_males = dplyr::filter(merge_EvP_BAT_males_direction, padj.p < 0.05 & padj.e >= 0.05) # 836 genes
sigsame_BAT_males = dplyr::filter(merge_EvP_BAT_males_direction, direction == 1 & padj.e < 0.05 & padj.p < 0.05) # 295 genes
sigopps_BAT_males = dplyr::filter(merge_EvP_BAT_males_direction, direction == 0 & padj.e < 0.05 & padj.p < 0.05) # 203 genes
sigboth_BAT_males = dplyr::filter(merge_EvP_BAT_males_direction, padj.e < 0.05 & padj.p < 0.05) # 498 genes

### plot adaptive and non-adaptive plasticity
BAT_males.plot <- merge_EvP_BAT_males_direction %>%
  ggplot(aes(x = lfc.p, y = lfc.e)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
  geom_point(data = sigsame_BAT_males, aes(x = lfc.p, y = lfc.e), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_point(data = sigopps_BAT_males, aes(x = lfc.p, y = lfc.e), fill = "#000000", color = "000000", stroke = 0.175, size = 1.5, shape = 21, alpha = 0.8, show.legend = TRUE) +
  scale_x_continuous(limits = c(-6.25,4.15), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-5,7.75), expand = c(0.01,0.075)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans")) +
  labs(x = "Log<sub>2</sub> Fold Change in BZ (cold vs warm)",
       y = "Log<sub>2</sub> Fold Change in Warm (NY vs BZ)")
#3BAF85
ggsave("results/figures/males_adap_plast_BAT.pdf", plot = BAT_males.plot, height = 2, width = 2.1)


## Correlation Test

# check normality of BAT data
BAT.model_males <- lm(lfc.e ~ lfc.p, data = sigboth_BAT_males)
#check_model(BAT.model_males) # data are not normally distributed - use Spearman correlation

r.obs_BAT_males <- cor(x = sigboth_BAT_males$lfc.p, y = sigboth_BAT_males$lfc.e, method = "spearman")
P.cor_BAT_males <- cor.test(x = sigboth_BAT_males$lfc.p, y = sigboth_BAT_males$lfc.e, method = "spearman")$p.value

# permutation function
cor.perm_BAT_males <- function (data = data, nperm = 10000)
{
  random.subsample_BAT_males <- dplyr::slice_sample(data, n=498) # randomly sub sample number estimated from sigboth
  r.per_BAT_males <- sapply (1:nperm, FUN = function (i) cor (x = random.subsample_BAT_males$lfc.p, y = sample (random.subsample_BAT_males$lfc.e)))
  r.per_BAT_males <- c(r.per_BAT_males, r.obs_BAT_males)
  hist (r.per_BAT_males, breaks = 100, xlab = "Correlation coefficient", main=NULL,
        col = "white")
  abline (v = r.obs_BAT_males, col = "#FCECAA", lty = "solid")
  P.per_BAT_males <- sum (abs (r.per_BAT_males) >= abs (r.obs_BAT_males))/(nperm + 1) 
  return (list (r.obs_BAT_males = r.obs_BAT_males, P.cor_BAT_males = P.cor_BAT_males,
                P.per_BAT_males = P.per_BAT_males, r.per_BAT_males = r.per_BAT_males))
}

BAT_perm_males <- cor.perm_BAT_males(data = merge_EvP_BAT_males) # test using all of the genes
BAT_perm_pvalue <- BAT_perm_males$P.per # P < 0.05

# make BAT_perm a df so you can plot hist in ggplot
df_BAT_perm_males <- data.frame(BAT_perm_males)
# make histogram
BAT_perm_males.inset <- df_BAT_perm_males %>%
  ggplot(aes(x=r.per_BAT_males)) +
  geom_histogram(fill = "white", color = "gray", bins = 100, binwidth = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 450), breaks = seq(0,400,125)) +
  geom_vline(xintercept = r.obs_BAT_males, color = "#E09832", linetype = "solid", size = 0.5) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 5, family = "sans"),
        plot.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Correlation Coefficient",
       y = "Count")
ggsave("results/figures/males_adap_plast_BAT_inset.pdf", plot = BAT_perm_males.inset,
       height = 1.5, width = 2.25, units = "cm")



##############################################################
# Female liver expression (DE in warm (NY vs BZ), and plasticity in Brazil (BZ warm vs cold))
##############################################################

### New York vs Brazil (in the warm -- evolved differences)
# filter to get mean of 10 reads across all samples
DE_base_evol_liver_females <- res_warm_females_liver_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_evol_liver_females) # 14,163 genes


### Brazil plasticity (cold vs warm):
# filter to get mean of 10 reads across all samples
DE_base_BZ_plast_liver_females <- res_BZ_females_liver_CvW %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_plast_liver_females) # 14,163 genes


### merge evol DE to plast Brazil DE in liver
merge_EvP_liver_females <- dplyr::full_join(DE_base_evol_liver_females, DE_base_BZ_plast_liver_females, by = "gene") %>%
  dplyr::select("gene",
         "lfc.e" = "log2FoldChange.x",
         "padj.e" = "padj.x",
         "lfc.p" = "log2FoldChange.y",
         "padj.p" = "padj.y")
#nrow(merge_EvP_liver_females) 14,163 genes


## add a column indicating whether evolved and plastic expression are in the same or opposite direction
merge_EvP_liver_females_direction <-  merge_EvP_liver_females %>%
  mutate(direction = case_when((sign(lfc.e)) == (sign(lfc.p)) ~ 1, # same direction (e.g., higher in NY warm and BZ cold)
                               (sign(lfc.e)) != (sign(lfc.p)) ~ 0)) # opposite direction

## pull out categories of genes that we are interested in
sigevol_liver_females = dplyr::filter(merge_EvP_liver_females_direction, padj.e < 0.05 & padj.p >= 0.05) # 3,015 genes
sigplas_liver_females = dplyr::filter(merge_EvP_liver_females_direction, padj.p < 0.05 & padj.e >= 0.05) # 406 genes
sigsame_liver_females = dplyr::filter(merge_EvP_liver_females_direction, direction == 1 & padj.e < 0.05 & padj.p < 0.05) # 142 genes
sigopps_liver_females = dplyr::filter(merge_EvP_liver_females_direction, direction == 0 & padj.e < 0.05 & padj.p < 0.05) # 56 genes
sigboth_liver_females = dplyr::filter(merge_EvP_liver_females_direction, padj.e < 0.05 & padj.p < 0.05) # 198 genes

### plot adaptive and non-adaptive plasticity
liver_females.plot <- merge_EvP_liver_females_direction %>%
  ggplot(aes(x = lfc.p, y = lfc.e)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
  geom_point(data = sigsame_liver_females, aes(x = lfc.p, y = lfc.e), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_point(data = sigopps_liver_females, aes(x = lfc.p, y = lfc.e), fill = "#000000", color = "000000", stroke = 0.175, size = 1.5, shape = 21, alpha = 0.8, show.legend = TRUE) +
  scale_x_continuous(limits = c(-6.25,4.15), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-5,7.75), expand = c(0.01,0.075)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans")) +
  labs(x = "Log<sub>2</sub> Fold Change in BZ (cold vs warm)",
       y = "Log<sub>2</sub> Fold Change in Warm (NY vs BZ)")
#3BAF85
ggsave("results/figures/females_adap_plast_liver.pdf", plot = liver_females.plot, height = 2, width = 2.1)


### Correlation Test

# check normality of liver data
liver.model_females <- lm(lfc.e ~ lfc.p, data = sigboth_liver_females)
#check_model(liver.model_females) # data are not normally distributed - use Spearman correlation

r.obs_liver_females <- cor(x = sigboth_liver_females$lfc.p, y = sigboth_liver_females$lfc.e, method = "spearman")
P.cor_liver_females <- cor.test(x = sigboth_liver_females$lfc.p, y = sigboth_liver_females$lfc.e, method = "spearman")$p.value


cor.perm_liver_females <- function (data = data, nperm = 10000)
{
  random.subsample_liver_females <- dplyr::slice_sample(data, n=198) # randomly subsample the number of sigboth
  r.per_liver_females <- sapply (1:nperm, FUN = function (i) cor (x = random.subsample_liver_females$lfc.p, y = sample (random.subsample_liver_females$lfc.e)))
  r.per_liver_females <- c(r.per_liver_females, r.obs_liver_females)
  hist (r.per_liver_females, breaks = 100, xlab = "Correlation coefficient", main=NULL,
        col = "white")
  abline (v = r.obs_liver_females, col = 'red', lty = "dotted")
  P.per_liver_females <- sum (abs (r.per_liver_females) >= abs (r.obs_liver_females))/(nperm + 1) 
  return (list (r.obs_liver_females = r.obs_liver_females, P.cor_liver_females = P.cor_liver_females,
                P.per_liver_females = P.per_liver_females, r.per_liver_females = r.per_liver_females))
}

liver_perm_females <- cor.perm_liver_females(data = merge_EvP_liver_females) # test using all of the genes
liver_p_perm_females <- liver_perm_females$P.per_liver_females # P < 0.05

# make liver_perm a df so you can plot hist in ggplot
df_liver_perm_females <- data.frame(liver_perm_females)
# make histogram
liver_perm_females.inset <- df_liver_perm_females %>%
  ggplot(aes(x=r.per_liver_females)) +
  geom_histogram(fill = "white", color = "gray", bins = 100, binwidth = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 500), breaks = seq(0,500,150)) +
  scale_x_continuous(expand = c(0,0), limits = c(NA, 0.425)) +
  geom_vline(xintercept = r.obs_liver_females, color = "#E09832", linetype = "solid", size = 0.5) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 5, family = "sans"),
        axis.text.y = element_text(size = 5, family = "sans"),
        plot.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Correlation Coefficient",
       y = "Count")
ggsave("results/figures/females_adap_plast_liver_inset.pdf", plot = liver_perm_females.inset,
       height = 1.5, width = 2.25, units = "cm")



##############################################################
# Female BAT expression (DE in warm (NY vs BZ), and plasticity in Brazil (BZ warm vs cold))
##############################################################

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
merge_EvP_BAT_females <- dplyr::full_join(DE_base_evol_BAT_females, DE_base_BZ_plast_BAT_females, by = "gene") %>%
  dplyr::select("gene",
         "lfc.e" = "log2FoldChange.x",
         "padj.e" = "padj.x",
         "lfc.p" = "log2FoldChange.y",
         "padj.p" = "padj.y")
#nrow(merge_EvP_BAT_females) # 14,261 genes

## add a column indicating whether evolved and plastic expression are in the same or opposite direction
merge_EvP_BAT_females_direction <-  merge_EvP_BAT_females %>%
  mutate(direction = case_when((sign(lfc.e)) == (sign(lfc.p)) ~ 1, # same direction
                               (sign(lfc.e)) != (sign(lfc.p)) ~ 0)) # opposite direction

## pull out categories of genes that we are interested in
sigevol_BAT_females = dplyr::filter(merge_EvP_BAT_females_direction, padj.e < 0.05 & padj.p >= 0.05) # 2,577 genes
sigplas_BAT_females = dplyr::filter(merge_EvP_BAT_females_direction, padj.p < 0.05 & padj.e >= 0.05) # 1,160 genes
sigsame_BAT_females = dplyr::filter(merge_EvP_BAT_females_direction, direction == 1 & padj.e < 0.05 & padj.p < 0.05) # 302 genes
sigopps_BAT_females = dplyr::filter(merge_EvP_BAT_females_direction, direction == 0 & padj.e < 0.05 & padj.p < 0.05) # 173 genes
sigboth_BAT_females = dplyr::filter(merge_EvP_BAT_females_direction, padj.e < 0.05 & padj.p < 0.05) # 475 genes

### plot adaptive and non-adaptive plasticity
BAT_females.plot <- merge_EvP_BAT_females_direction %>%
  ggplot(aes(x = lfc.p, y = lfc.e)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
  geom_point(data = sigsame_BAT_females, aes(x = lfc.p, y = lfc.e), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_point(data = sigopps_BAT_females, aes(x = lfc.p, y = lfc.e), fill = "#000000", color = "000000", stroke = 0.175, size = 1.5, shape = 21, alpha = 0.8, show.legend = TRUE) +
  scale_x_continuous(limits = c(-5,15), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-9,22), expand = c(0.01,0.075), breaks = c(-9,0,9,20)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans")) +
  labs(x = "Log<sub>2</sub> Fold Change in BZ (cold vs warm)",
       y = "Log<sub>2</sub> Fold Change in Warm (NY vs BZ)")
#3BAF85
ggsave("results/figures/females_adap_plast_BAT.pdf", plot = BAT_females.plot, height = 2, width = 2.1)


## Correlation Test

# check normality of BAT data
BAT.model_females <- lm(lfc.e ~ lfc.p, data = sigboth_BAT_females)
#check_model(BAT.model_females) # data are not normally distributed - use Spearman correlation

r.obs_BAT_females <- cor(x = sigboth_BAT_females$lfc.p, y = sigboth_BAT_females$lfc.e, method = "spearman")
P.cor_BAT_females <- cor.test(x = sigboth_BAT_females$lfc.p, y = sigboth_BAT_females$lfc.e, method = "spearman")$p.value

# permutation function
cor.perm_BAT_females <- function (data = data, nperm = 10000)
{
  random.subsample_BAT_females <- dplyr::slice_sample(data, n=475) # randomly sub sample number estimated from sigboth
  r.per_BAT_females <- sapply (1:nperm, FUN = function (i) cor (x = random.subsample_BAT_females$lfc.p, y = sample (random.subsample_BAT_females$lfc.e)))
  r.per_BAT_females <- c(r.per_BAT_females, r.obs_BAT_females)
  hist (r.per_BAT_females, breaks = 100, xlab = "Correlation coefficient", main=NULL,
        col = "white")
  abline (v = r.obs_BAT_females, col = "#FCECAA", lty = "solid")
  P.per_BAT_females <- sum (abs (r.per_BAT_females) >= abs (r.obs_BAT_females))/(nperm + 1) 
  return (list (r.obs_BAT_females = r.obs_BAT_females, P.cor_BAT_females = P.cor_BAT_females,
                P.per_BAT_females = P.per_BAT_females, r.per_BAT_females = r.per_BAT_females))
}

BAT_perm_females <- cor.perm_BAT_females(data = merge_EvP_BAT_females) # test using all of the genes
BAT_perm_p_females <- BAT_perm_females$P.per # P < 0.05

# make BAT_perm a df so you can plot hist in ggplot
df_BAT_perm_females <- data.frame(BAT_perm_females)
# make histogram
BAT_perm_females.inset <- df_BAT_perm_females %>%
  ggplot(aes(x=r.per_BAT_females)) +
  geom_histogram(fill = "white", color = "gray", bins = 100, binwidth = NULL) +
  scale_x_continuous(expand = c(0,0), limits = c(NA, 0.45)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 800), breaks = seq(0,750,250)) +
  geom_vline(xintercept = r.obs_BAT_females, color = "#E09832", linetype = "solid", size = 0.5) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 5, family = "sans"),
        plot.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Correlation Coefficient",
       y = "Count")
ggsave("results/figures/females_adap_plast_BAT_inset.pdf", plot = BAT_perm_females.inset,
       height = 1.5, width = 2.25, units = "cm")



##############################################################
# Statistical Analysis - Binomial Test
##############################################################

## males

adap_BAT_males <- sigboth_BAT_males %>%
  dplyr::select(gene, direction) %>%
  mutate(adap = ifelse(direction > 0, TRUE, FALSE))

adap_liver_males <- sigboth_liver_males %>%
  dplyr::select(gene, direction) %>%
  mutate(adap = ifelse(direction > 0, TRUE, FALSE))

# BAT
dim(adap_BAT_males)[1] # get number of genes in this category (498)
sum(adap_BAT_males$adap) # the number of genes where expression is higher in Brazil (295)
btest_BAT_males = binom.test(x=sum(adap_BAT_males$adap), n = dim(adap_BAT_males)[1])
# P < 0.001

# liver
dim(adap_liver_males)[1] # get number of genes in this category (407)
sum(adap_liver_males$adap) # the number of genes where expression is higher in Brazil (296)
btest_liver_males = binom.test(x=sum(adap_liver_males$adap), n = dim(adap_liver_males)[1])
# P < 0.001


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


