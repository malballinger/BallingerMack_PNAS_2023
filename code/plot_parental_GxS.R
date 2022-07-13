#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots differential expression in parental samples, with specific focus on patterns of GxS.
# This script generates Figures S4 in BallingerMack_2022.

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
# Warm temperature expression patterns
##############################################################

### liver

## filter to get mean of 10 reads across all samples
DE_base_NY_warm_liver_MvF <- res_NY_liver_warm_MvF %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_NY_warm_liver_MvF) # 14,153 genes
#DE_base_NY_warm_liver_MvF %>% filter(padj < 0.05) %>% nrow() # 894 genes

DE_base_BZ_warm_liver_MvF <- res_BZ_liver_warm_MvF %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_warm_liver_MvF) # 14,153 genes
#DE_base_BZ_warm_liver_MvF %>% filter(padj < 0.05) %>% nrow() # 405 genes

## merge datasets
merge_sex_warm_liver <- dplyr::full_join(DE_base_BZ_warm_liver_MvF, DE_base_NY_warm_liver_MvF, by = "gene") %>%
  dplyr::select(
    "gene",
    "lfc.BZ" = "log2FoldChange.x",
    "lfcSE.BZ" = "lfcSE.x",
    "padj.BZ" = "padj.x",
    "lfc.NY" = "log2FoldChange.y",
    "lfcSE.NY" = "lfcSE.y",
    "padj.NY" = "padj.y",
    )
#nrow(merge_sex_warm_liver) 14,153 genes
#merge_sex_warm_liver %>% filter(padj.BZ < 0.05 | padj.NY < 0.05) %>% nrow() # 1,102 genes
# roughly 8% of all liver genes show sexual dimorphism in the warm

## add a column indicating whether BZ and NY groups are in the same direction
merge_sex_warm_liver_direction <-  merge_sex_warm_liver %>%
  mutate(direction = case_when((sign(lfc.BZ)) == (sign(lfc.NY)) ~ 1, # same direction
                               (sign(lfc.BZ)) != (sign(lfc.NY)) ~ 0)) # opposite direction

## shrink plotting window by assigning large LFC to 5's
merge_sex_warm_liver_direction$lfc.NY[merge_sex_warm_liver_direction$lfc.NY>5]<-5
merge_sex_warm_liver_direction$lfc.BZ[merge_sex_warm_liver_direction$lfc.BZ>5]<-5
merge_sex_warm_liver_direction$lfc.NY[merge_sex_warm_liver_direction$lfc.NY<(-5)]<-(-5)
merge_sex_warm_liver_direction$lfc.BZ[merge_sex_warm_liver_direction$lfc.BZ<(-5)]<-(-5)

## categorize genes by direction and/or significance
sigNY_sex_warm_liver = dplyr::filter(merge_sex_warm_liver_direction, padj.NY < 0.05 & padj.BZ >= 0.05) # genotype-specific (NY) response = 697 genes
sigBZ_sex_warm_liver = dplyr::filter(merge_sex_warm_liver_direction, padj.BZ < 0.05 & padj.NY >= 0.05) # genotype-specific (BZ) response = 208 genes
sigcons_sex_warm_liver = dplyr::filter(merge_sex_warm_liver_direction, direction == 1 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, and in same directions (this category is hard to
                                                                                                                # disentangle GxS from G + S) = 197 genes
sigopps_sex_warm_liver = dplyr::filter(merge_sex_warm_liver_direction, direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, but in opposite directions = 0 genes
nonsig_sex_warm_liver = dplyr::filter(merge_sex_warm_liver_direction, padj.NY >= 0.05 & padj.BZ >= 0.05) # 13,046 genes

## collate all genes that show evidence of GxS
# first, genes identified as GxS by DESEQ (very stringent)
DESEQ_sex_liver_warm_gxs <- res_sex_liver_warm_GxS %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10, padj < 0.05)
# save a list of GxS genes defined by DESEQ
DESEQ_sex_liver_warm_gxs_list <- DESEQ_sex_liver_warm_gxs %>% pull(gene)

# second, genes I identify as GxS from categorization above + DESEQ genes
sex_warm_liver_GxS <- merge_sex_warm_liver_direction %>%
  filter((padj.NY < 0.05 & padj.BZ >= 0.05) |
           (padj.BZ < 0.05 & padj.NY >= 0.05) |
           (direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) |
           (gene %in% DESEQ_sex_liver_warm_gxs_list))
# save a list of all GxE genes
sex_warm_liver_GxS_list <- sex_warm_liver_GxS %>% pull(gene)
# 908 GxS genes in warm liver

## double check that DESEQ GxS genes are in my list
#dplyr::anti_join(DESEQ_sex_liver_warm_gxs, sex_warm_liver_GxS, by = "gene")
# yes

## plot
sex_warm_liver_GxS.plot <- merge_sex_warm_liver_direction %>%
  ggplot(aes(x = lfc.NY, y = lfc.BZ)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
  geom_abline(slope = 1, linetype = "dotted", size = 0.25) +
  geom_point(data = nonsig_sex_warm_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "lightgray", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigNY_sex_warm_liver, aes(x = lfc.NY, y = lfc.BZ), fill =  "#1683B5", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigBZ_sex_warm_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigopps_sex_warm_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "#000000", color = "000000", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigcons_sex_warm_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "#B22222", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  scale_x_continuous(limits = c(-5,5), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-5,5),expand = c(0.01,0.075)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans"),
        axis.title = element_blank())
  #labs(x = "Log<sub>2</sub> Fold Change in NY<br>(warm vs cold)",
       #y = "Log<sub>2</sub> Fold Change in BZ<br>(warm vs cold)")
ggsave("results/figures/sex_DE_warm_liver.pdf", plot = sex_warm_liver_GxS.plot, height = 2, width = 2.1)



### BAT

## filter to get mean of 10 reads across all samples
DE_base_NY_warm_BAT_MvF <- res_NY_BAT_warm_MvF %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_NY_warm_BAT_MvF) # 14,162 genes
#DE_base_NY_warm_BAT_MvF %>% filter(padj < 0.05) %>% nrow() # 15 genes

DE_base_BZ_warm_BAT_MvF <- res_BZ_BAT_warm_MvF %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_warm_BAT_MvF) # 14,162 genes
#DE_base_BZ_warm_BAT_MvF %>% filter(padj < 0.05) %>% nrow() # 22 genes

## merge datasets
merge_sex_warm_BAT <- dplyr::full_join(DE_base_BZ_warm_BAT_MvF, DE_base_NY_warm_BAT_MvF, by = "gene") %>%
  dplyr::select(
    "gene",
    "lfc.BZ" = "log2FoldChange.x",
    "lfcSE.BZ" = "lfcSE.x",
    "padj.BZ" = "padj.x",
    "lfc.NY" = "log2FoldChange.y",
    "lfcSE.NY" = "lfcSE.y",
    "padj.NY" = "padj.y",
  )
#nrow(merge_sex_warm_BAT) 14,162 genes
#merge_sex_warm_BAT %>% filter(padj.BZ < 0.05 | padj.NY < 0.05) %>% nrow() # 25 genes
# less than 1% of all BAT genes show sexual dimorphism in the warm

## add a column indicating whether BZ and NY groups are in the same direction
merge_sex_warm_BAT_direction <-  merge_sex_warm_BAT %>%
  mutate(direction = case_when((sign(lfc.BZ)) == (sign(lfc.NY)) ~ 1, # same direction
                               (sign(lfc.BZ)) != (sign(lfc.NY)) ~ 0)) # opposite direction

## shrink plotting window by assigning large LFC to 5's
merge_sex_warm_BAT_direction$lfc.NY[merge_sex_warm_BAT_direction$lfc.NY>5]<-5
merge_sex_warm_BAT_direction$lfc.BZ[merge_sex_warm_BAT_direction$lfc.BZ>5]<-5
merge_sex_warm_BAT_direction$lfc.NY[merge_sex_warm_BAT_direction$lfc.NY<(-5)]<-(-5)
merge_sex_warm_BAT_direction$lfc.BZ[merge_sex_warm_BAT_direction$lfc.BZ<(-5)]<-(-5)

## categorize genes by direction and/or significance
sigNY_sex_warm_BAT = dplyr::filter(merge_sex_warm_BAT_direction, padj.NY < 0.05 & padj.BZ >= 0.05) # genotype-specific (NY) response = 3 genes
sigBZ_sex_warm_BAT = dplyr::filter(merge_sex_warm_BAT_direction, padj.BZ < 0.05 & padj.NY >= 0.05) # genotype-specific (BZ) response = 10 genes
sigcons_sex_warm_BAT = dplyr::filter(merge_sex_warm_BAT_direction, direction == 1 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, and in same directions (this category is hard to
                                                                                                                      # disentangle GxS from G + S) = 12 genes
sigopps_sex_warm_BAT = dplyr::filter(merge_sex_warm_BAT_direction, direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, but in opposite directions = 0 genes
nonsig_sex_warm_BAT = dplyr::filter(merge_sex_warm_BAT_direction, padj.NY >= 0.05 & padj.BZ >= 0.05) # 14,065 genes

## collate all genes that show evidence of GxS
# first, genes identified as GxS by DESEQ (very stringent)
DESEQ_sex_BAT_warm_gxs <- res_sex_BAT_warm_GxS %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10, padj < 0.05)
# 0 genes

# second, genes I identify as GxS from categorization above + DESEQ genes
sex_warm_BAT_GxS <- merge_sex_warm_BAT_direction %>%
  filter((padj.NY < 0.05 & padj.BZ >= 0.05) |
           (padj.BZ < 0.05 & padj.NY >= 0.05) |
           (direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05))
# save a list of all GxE genes
sex_warm_BAT_GxS_list <- sex_warm_BAT_GxS %>% pull(gene)
# 13 GxS genes in warm BAT

## plot
sex_warm_BAT_GxS.plot <- merge_sex_warm_BAT_direction %>%
  ggplot(aes(x = lfc.NY, y = lfc.BZ)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
  geom_abline(slope = 1, linetype = "dotted", size = 0.25) +
  geom_point(data = nonsig_sex_warm_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "lightgray", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigNY_sex_warm_BAT, aes(x = lfc.NY, y = lfc.BZ), fill =  "#1683B5", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigBZ_sex_warm_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigopps_sex_warm_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "#000000", color = "000000", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigcons_sex_warm_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "#B22222", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  scale_x_continuous(limits = c(-5,5), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-5,5),expand = c(0.01,0.075)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans"),
        axis.title = element_blank())
#labs(x = "Log<sub>2</sub> Fold Change in NY<br>(warm vs cold)",
#y = "Log<sub>2</sub> Fold Change in BZ<br>(warm vs cold)")
ggsave("results/figures/sex_DE_warm_liver.pdf", plot = sex_warm_liver_GxS.plot, height = 2, width = 2.1)




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

