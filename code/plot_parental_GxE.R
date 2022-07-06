#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots differential expression in parental samples, with specific focus on patterns of GxE.
# This script generates Figures 2B (males) and S3C (females) in BallingerMack_2022.

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clear R's environment
library(tidyverse)
library(ggtext)
library(DESeq2)

set.seed(19910118)

source("./code/get_DE_parents_analysis.R") # where DESeq data frames are generated

##############################################################
# Male expression patterns
##############################################################

### liver

## filter to get mean of 10 reads across all samples
DE_base_NY_males_liver <- res_NY_males_liver_WvC %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_NY_males_liver) 14,028 genes
#DE_base_NY_males_liver %>% filter(padj < 0.05) %>% nrow() # 211 genes

DE_base_BZ_males_liver <- res_BZ_males_liver_WvC %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_males_liver) 14,028 genes
#DE_base_BZ_males_liver %>% filter(padj < 0.05) %>% nrow() # 817 genes

## merge datasets
merge_plast_males_liver <- dplyr::full_join(DE_base_BZ_males_liver, DE_base_NY_males_liver, by = "gene") %>%
  dplyr::select(
    "gene",
    "lfc.BZ" = "log2FoldChange.x",
    "lfcSE.BZ" = "lfcSE.x",
    "padj.BZ" = "padj.x",
    "lfc.NY" = "log2FoldChange.y",
    "lfcSE.NY" = "lfcSE.y",
    "padj.NY" = "padj.y",
    )
#nrow(merge_plast_males_liver) 14,028 genes
#merge_plast_males_liver %>% filter(padj.BZ < 0.05 | padj.NY < 0.05) %>% nrow() # 922 genes
# roughly 7% of all liver genes are plastic in males

## add a column indicating whether warm and cold groups are in the same direction
merge_plast_males_liver_direction <-  merge_plast_males_liver %>%
  mutate(direction = case_when((sign(lfc.BZ)) == (sign(lfc.NY)) ~ 1, # same direction
                               (sign(lfc.BZ)) != (sign(lfc.NY)) ~ 0)) # opposite direction

## shrink plotting window by assigning large LFC to 5's
merge_plast_males_liver_direction$lfc.NY[merge_plast_males_liver_direction$lfc.NY>5]<-5
merge_plast_males_liver_direction$lfc.BZ[merge_plast_males_liver_direction$lfc.BZ>5]<-5
merge_plast_males_liver_direction$lfc.NY[merge_plast_males_liver_direction$lfc.NY<(-5)]<-(-5)
merge_plast_males_liver_direction$lfc.BZ[merge_plast_males_liver_direction$lfc.BZ<(-5)]<-(-5)

## categorize genes by direction and/or significance
sigNY_liver = dplyr::filter(merge_plast_males_liver_direction, padj.NY < 0.05 & padj.BZ >= 0.05) # pop-specific (NY) response = 105 genes
sigBZ_liver = dplyr::filter(merge_plast_males_liver_direction, padj.BZ < 0.05 & padj.NY >= 0.05) # pop-specific (BZ) response = 711 genes
sigcons_liver = dplyr::filter(merge_plast_males_liver_direction, direction == 1 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, and in same directions (this category is hard to
                                                                                                                # disentangle GxE from G + E) = 106 genes
sigopps_liver = dplyr::filter(merge_plast_males_liver_direction, direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, but in opposite directions = 0 genes
nonsig_liver = dplyr::filter(merge_plast_males_liver_direction, padj.NY >= 0.05 & padj.BZ >= 0.05) # 13,101 genes

## collate all genes that show evidence of GxE
# first, genes identified as GxE by DESEQ (very stringent)
DESEQ_males_liver_gxe <- res_males_liver_GxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10, padj < 0.05)
# save a list of GxE genes defined by DESEQ
DESEQ_males_liver_gxe_list <- DESEQ_males_liver_gxe %>% pull(gene)

# second, genes I identify as GxE from categorization above + DESEQ genes
liver_males_plast_GxE <- merge_plast_males_liver_direction %>%
  filter((padj.NY < 0.05 & padj.BZ >= 0.05) |
           (padj.BZ < 0.05 & padj.NY >= 0.05) |
           (direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) |
           (gene %in% DESEQ_males_liver_gxe_list))
# save a list of all GxE genes
liver_males_plast_GxE_list <- liver_males_plast_GxE %>% pull(gene)
# 817 GxE genes in liver

## double check that DESEQ GxE genes are in my list
#dplyr::anti_join(DESEQ_males_liver_gxe, liver_males_plast_GxE, by = "gene")
# yes

## plot
liver_males_plasticity <- merge_plast_males_liver_direction %>%
  ggplot(aes(x = lfc.NY, y = lfc.BZ)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
  geom_abline(slope = 1, linetype = "dotted", size = 0.25) +
  geom_point(data = nonsig_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "lightgray", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigNY_liver, aes(x = lfc.NY, y = lfc.BZ), fill =  "#1683B5", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigBZ_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigopps_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "#000000", color = "000000", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigcons_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "#B22222", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  scale_x_continuous(limits = c(-5,5), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-5,5),expand = c(0.01,0.075)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans"),
        axis.title = element_blank())
  #labs(x = "Log<sub>2</sub> Fold Change in NY<br>(warm vs cold)",
       #y = "Log<sub>2</sub> Fold Change in BZ<br>(warm vs cold)")
ggsave("results/figures/males_DE_plasticity_liver.pdf", plot = liver_males_plasticity, height = 2, width = 2.1)



### BAT

## filter to get mean of 10 reads across all samples
DE_base_NY_males_BAT <- res_NY_males_BAT_WvC %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_NY_males_BAT) 14,176 genes
#DE_base_NY_males_BAT %>% filter(padj < 0.05) %>% nrow() # 883 genes

DE_base_BZ_males_BAT <- res_BZ_males_BAT_WvC %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_males_BAT) 14,176 genes
#DE_base_BZ_males_BAT %>% filter(padj < 0.05) %>% nrow() # 1,334 genes

## merge datasets
merge_plast_males_BAT <- dplyr::full_join(DE_base_BZ_males_BAT, DE_base_NY_males_BAT, by = "gene") %>%
  dplyr::select(
    "gene",
    "lfc.BZ" = "log2FoldChange.x",
    "lfcSE.BZ" = "lfcSE.x",
    "padj.BZ" = "padj.x",
    "lfc.NY" = "log2FoldChange.y",
    "lfcSE.NY" = "lfcSE.y",
    "padj.NY" = "padj.y",
  )
#nrow(merge_plast_males_BAT) 14,176 genes
#merge_plast_males_BAT %>% filter(padj.BZ < 0.05 | padj.NY < 0.05) %>% nrow() # 1,819 genes
# roughly 13% of all BAT genes are plastic in males

## add a column indicating whether warm and cold groups are in the same direction
merge_plast_males_BAT_direction <-  merge_plast_males_BAT %>%
  mutate(direction = case_when((sign(lfc.BZ)) == (sign(lfc.NY)) ~ 1, # same direction
                               (sign(lfc.BZ)) != (sign(lfc.NY)) ~ 0)) # opposite direction

## shrink plotting window by assigning large LFC to 5's
merge_plast_males_BAT_direction$lfc.NY[merge_plast_males_BAT_direction$lfc.NY>5]<-5
merge_plast_males_BAT_direction$lfc.BZ[merge_plast_males_BAT_direction$lfc.BZ>5]<-5
merge_plast_males_BAT_direction$lfc.NY[merge_plast_males_BAT_direction$lfc.NY<(-5)]<-(-5)
merge_plast_males_BAT_direction$lfc.BZ[merge_plast_males_BAT_direction$lfc.BZ<(-5)]<-(-5)

## categorize genes by direction and/or significance
sigNY_BAT = dplyr::filter(merge_plast_males_BAT_direction, padj.NY < 0.05 & padj.BZ >= 0.05) # pop-specific (NY) response = 485 genes
sigBZ_BAT = dplyr::filter(merge_plast_males_BAT_direction, padj.BZ < 0.05 & padj.NY >= 0.05) # pop-specific (BZ) response = 936 genes
sigcons_BAT = dplyr::filter(merge_plast_males_BAT_direction, direction == 1 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, and in same directions (this category is hard to
                                                                                                          # disentangle GxE from G + E) = 396 genes
sigopps_BAT = dplyr::filter(merge_plast_males_BAT_direction, direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, but in opposite directions = 2 genes
nonsig_BAT = dplyr::filter(merge_plast_males_BAT_direction, padj.NY >= 0.05 & padj.BZ >= 0.05) # 12,338 genes
sigboth_BAT = dplyr::filter(merge_plast_males_BAT_direction, padj.NY < 0.05 & padj.BZ < 0.05) # should be sum of sigcons and sigopps

## collate all genes that show evidence of GxE
# first, genes identified as GxE by DESEQ (very stringent)
DESEQ_males_BAT_gxe <- res_males_BAT_GxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10, padj < 0.05)
# save a list of GxE genes defined by DESEQ
DESEQ_males_BAT_gxe_list <- DESEQ_males_BAT_gxe %>% pull(gene)

# second, genes I identify as GxE from categorization above + DESEQ genes
BAT_males_plast_GxE <- merge_plast_males_BAT_direction %>%
  filter((padj.NY < 0.05 & padj.BZ >= 0.05) |
           (padj.BZ < 0.05 & padj.NY >= 0.05) |
           (direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) |
           (gene %in% DESEQ_males_BAT_gxe_list))
# save a list of all GxE genes
BAT_males_plast_GxE_list <- BAT_males_plast_GxE %>% pull(gene)
# 1,425 genes
## double check that DESEQ GxE genes are in my list
#dplyr::anti_join(DESEQ_males_BAT_gxe, BAT_males_plast_GxE, by = "gene")
# yes

## plot
BAT_males_plasticity <- merge_plast_males_BAT_direction %>%
  ggplot(aes(x = lfc.NY, y = lfc.BZ)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
  geom_abline(slope = 1, linetype = "dotted", size = 0.25) +
  geom_point(data = nonsig_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "lightgray", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigNY_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "#1683B5", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigBZ_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigopps_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "#000000", color = "000000", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigcons_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "#B22222", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  scale_x_continuous(limits = c(-5,5), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-5,5), expand = c(0.01,0.075)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans")) +
  labs(x = "Log<sub>2</sub> Fold Change in NY<br>(warm vs cold)",
       y = "Log<sub>2</sub> Fold Change in BZ<br>(warm vs cold)")
ggsave("results/figures/males_DE_plasticity_BAT.pdf", plot = BAT_males_plasticity, height = 2, width = 2.1)




##############################################################
# Female expression patterns
##############################################################

### liver

## filter to get mean of 10 reads across all samples
DE_base_NY_females_liver <- res_NY_females_liver_WvC %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_NY_females_liver) 14,163 genes
#DE_base_NY_females_liver %>% filter(padj < 0.05) %>% nrow() # 151 genes

DE_base_BZ_females_liver <- res_BZ_females_liver_WvC %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_females_liver) 14,163 genes
#DE_base_BZ_females_liver %>% filter(padj < 0.05) %>% nrow() # 604 genes

## merge datasets
merge_plast_females_liver <- dplyr::full_join(DE_base_BZ_females_liver, DE_base_NY_females_liver, by = "gene") %>%
  dplyr::select(
    "gene",
    "lfc.BZ" = "log2FoldChange.x",
    "lfcSE.BZ" = "lfcSE.x",
    "padj.BZ" = "padj.x",
    "lfc.NY" = "log2FoldChange.y",
    "lfcSE.NY" = "lfcSE.y",
    "padj.NY" = "padj.y",
  )
#nrow(merge_plast_females_liver) 14,163 genes
#merge_plast_females_liver %>% filter(padj.BZ < 0.05 | padj.NY < 0.05) %>% nrow() # 702 genes
# roughly 5% of all liver genes are plastic in females

## add a column indicating whether warm and cold groups are in the same direction
merge_plast_females_liver_direction <-  merge_plast_females_liver %>%
  mutate(direction = case_when((sign(lfc.BZ)) == (sign(lfc.NY)) ~ 1, # same direction
                               (sign(lfc.BZ)) != (sign(lfc.NY)) ~ 0)) # opposite direction

## shrink plotting window by assigning large LFC to 5's
merge_plast_females_liver_direction$lfc.NY[merge_plast_females_liver_direction$lfc.NY>5]<-5
merge_plast_females_liver_direction$lfc.BZ[merge_plast_females_liver_direction$lfc.BZ>5]<-5
merge_plast_females_liver_direction$lfc.NY[merge_plast_females_liver_direction$lfc.NY<(-5)]<-(-5)
merge_plast_females_liver_direction$lfc.BZ[merge_plast_females_liver_direction$lfc.BZ<(-5)]<-(-5)

## categorize genes by direction and/or significance
sigNY_liver_f = dplyr::filter(merge_plast_females_liver_direction, padj.NY < 0.05 & padj.BZ >= 0.05) # pop-specific (NY) response = 98 genes
sigBZ_liver_f = dplyr::filter(merge_plast_females_liver_direction, padj.BZ < 0.05 & padj.NY >= 0.05) # pop-specific (BZ) response = 511 genes
sigcons_liver_f = dplyr::filter(merge_plast_females_liver_direction, direction == 1 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, and in same directions (this category is hard to
# disentangle GxE from G + E) = 53 genes
sigopps_liver_f = dplyr::filter(merge_plast_females_liver_direction, direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, but in opposite directions = 0 genes
nonsig_liver_f = dplyr::filter(merge_plast_females_liver_direction, padj.NY >= 0.05 & padj.BZ >= 0.05) # 13,459 genes

## collate all genes that show evidence of GxE
# first, genes identified as GxE by DESEQ (very stringent)
DESEQ_females_liver_gxe <- res_females_liver_GxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10, padj < 0.05) # 0 GxE genes
# save a list of GxE genes defined by DESEQ
DESEQ_females_liver_gxe_list <- DESEQ_females_liver_gxe %>% pull(gene)

# second, genes I identify as GxE from categorization above + DESEQ genes
liver_females_plast_GxE <- merge_plast_females_liver_direction %>%
  filter((padj.NY < 0.05 & padj.BZ >= 0.05) |
           (padj.BZ < 0.05 & padj.NY >= 0.05) |
           (direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) |
           (gene %in% DESEQ_females_liver_gxe_list))
# save a list of all GxE genes
liver_females_plast_GxE_list <- liver_females_plast_GxE %>% pull(gene)
# 649 GxE genes in female liver

## double check that DESEQ GxE genes are in my list
#dplyr::anti_join(DESEQ_females_liver_gxe, liver_females_plast_GxE, by = "gene")
# yes

## plot
liver_females_plasticity <- merge_plast_females_liver_direction %>%
  ggplot(aes(x = lfc.NY, y = lfc.BZ)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
  geom_abline(slope = 1, linetype = "dotted", size = 0.25) +
  geom_point(data = nonsig_liver_f, aes(x = lfc.NY, y = lfc.BZ), fill = "lightgray", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigNY_liver_f, aes(x = lfc.NY, y = lfc.BZ), fill =  "#1683B5", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigBZ_liver_f, aes(x = lfc.NY, y = lfc.BZ), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigopps_liver_f, aes(x = lfc.NY, y = lfc.BZ), fill = "#000000", color = "000000", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigcons_liver_f, aes(x = lfc.NY, y = lfc.BZ), fill = "#B22222", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  scale_x_continuous(limits = c(-5,5), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-5,5),expand = c(0.01,0.075)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans"),
        axis.title = element_blank())
#labs(x = "Log<sub>2</sub> Fold Change in NY<br>(warm vs cold)",
#y = "Log<sub>2</sub> Fold Change in BZ<br>(warm vs cold)")
ggsave("results/figures/females_DE_plasticity_liver.pdf", plot = liver_females_plasticity, height = 2, width = 2.1)



### BAT

## filter to get mean of 10 reads across all samples
DE_base_NY_females_BAT <- res_NY_females_BAT_WvC %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_NY_females_BAT) 14,261 genes
#DE_base_NY_females_BAT %>% filter(padj < 0.05) %>% nrow() # 594 genes

DE_base_BZ_females_BAT <- res_BZ_females_BAT_WvC %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_females_BAT) 14,261 genes
#DE_base_BZ_females_BAT %>% filter(padj < 0.05) %>% nrow() # 1,635 genes

## merge datasets
merge_plast_females_BAT <- dplyr::full_join(DE_base_BZ_females_BAT, DE_base_NY_females_BAT, by = "gene") %>%
  dplyr::select(
    "gene",
    "lfc.BZ" = "log2FoldChange.x",
    "lfcSE.BZ" = "lfcSE.x",
    "padj.BZ" = "padj.x",
    "lfc.NY" = "log2FoldChange.y",
    "lfcSE.NY" = "lfcSE.y",
    "padj.NY" = "padj.y",
  )
#nrow(merge_plast_females_BAT) 14,261 genes
#merge_plast_females_BAT %>% filter(padj.BZ < 0.05 | padj.NY < 0.05) %>% nrow() # 1,919 genes
# roughly 13% of all BAT genes are plastic in females

## add a column indicating whether warm and cold groups are in the same direction
merge_plast_females_BAT_direction <-  merge_plast_females_BAT %>%
  mutate(direction = case_when((sign(lfc.BZ)) == (sign(lfc.NY)) ~ 1, # same direction
                               (sign(lfc.BZ)) != (sign(lfc.NY)) ~ 0)) # opposite direction

## shrink plotting window by assigning large LFC to 5's
merge_plast_females_BAT_direction$lfc.NY[merge_plast_females_BAT_direction$lfc.NY>5]<-5
merge_plast_females_BAT_direction$lfc.BZ[merge_plast_females_BAT_direction$lfc.BZ>5]<-5
merge_plast_females_BAT_direction$lfc.NY[merge_plast_females_BAT_direction$lfc.NY<(-5)]<-(-5)
merge_plast_females_BAT_direction$lfc.BZ[merge_plast_females_BAT_direction$lfc.BZ<(-5)]<-(-5)

## categorize genes by direction and/or significance
sigNY_BAT_f = dplyr::filter(merge_plast_females_BAT_direction, padj.NY < 0.05 & padj.BZ >= 0.05) # pop-specific (NY) response = 284 genes
sigBZ_BAT_f = dplyr::filter(merge_plast_females_BAT_direction, padj.BZ < 0.05 & padj.NY >= 0.05) # pop-specific (BZ) response = 1,325 genes
sigcons_BAT_f = dplyr::filter(merge_plast_females_BAT_direction, direction == 1 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, and in same directions (this category is hard to
# disentangle GxE from G + E) = 280 genes
sigopps_BAT_f = dplyr::filter(merge_plast_females_BAT_direction, direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, but in opposite directions = 30 genes
nonsig_BAT_f = dplyr::filter(merge_plast_females_BAT_direction, padj.NY >= 0.05 & padj.BZ >= 0.05) # 12,337 genes
sigboth_BAT_f = dplyr::filter(merge_plast_females_BAT_direction, padj.NY < 0.05 & padj.BZ < 0.05) # should be sum of sigcons and sigopps

## collate all genes that show evidence of GxE
# first, genes identified as GxE by DESEQ (very stringent)
DESEQ_females_BAT_gxe <- res_females_BAT_GxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10, padj < 0.05) # 92 genes
# save a list of GxE genes defined by DESEQ
DESEQ_females_BAT_gxe_list <- DESEQ_females_BAT_gxe %>% pull(gene)

# second, genes I identify as GxE from categorization above + DESEQ genes
BAT_females_plast_GxE <- merge_plast_females_BAT_direction %>%
  filter((padj.NY < 0.05 & padj.BZ >= 0.05) |
           (padj.BZ < 0.05 & padj.NY >= 0.05) |
           (direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) |
           (gene %in% DESEQ_females_BAT_gxe_list))
# save a list of all GxE genes
BAT_females_plast_GxE_list <- BAT_females_plast_GxE %>% pull(gene)
# 1,650 genes
## double check that DESEQ GxE genes are in my list
#dplyr::anti_join(DESEQ_females_BAT_gxe, BAT_females_plast_GxE, by = "gene")
# yes

## plot
BAT_females_plasticity <- merge_plast_females_BAT_direction %>%
  ggplot(aes(x = lfc.NY, y = lfc.BZ)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
  geom_abline(slope = 1, linetype = "dotted", size = 0.25) +
  geom_point(data = nonsig_BAT_f, aes(x = lfc.NY, y = lfc.BZ), fill = "lightgray", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigNY_BAT_f, aes(x = lfc.NY, y = lfc.BZ), fill = "#1683B5", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigBZ_BAT_f, aes(x = lfc.NY, y = lfc.BZ), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigopps_BAT_f, aes(x = lfc.NY, y = lfc.BZ), fill = "#000000", color = "000000", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigcons_BAT_f, aes(x = lfc.NY, y = lfc.BZ), fill = "#B22222", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  scale_x_continuous(limits = c(-5,5), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-5,5), expand = c(0.01,0.075)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans")) +
  labs(x = "Log<sub>2</sub> Fold Change in NY<br>(warm vs cold)",
       y = "Log<sub>2</sub> Fold Change in BZ<br>(warm vs cold)")
ggsave("results/figures/females_DE_plasticity_BAT.pdf", plot = BAT_females_plasticity, height = 2, width = 2.1)

