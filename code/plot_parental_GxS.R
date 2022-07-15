#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots differential expression in parental samples, with specific focus on patterns of GxS.
# This script generates Figure S4 in BallingerMack_2022.

# Main Result: BAT harbors very little GxS

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
  #labs(x = "Log<sub>2</sub> Fold Change in NY<br>(male vs female)",
       #y = "Log<sub>2</sub> Fold Change in BZ<br>(male vs female)")
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
#labs(x = "Log<sub>2</sub> Fold Change in NY<br>(male vs female)",
#y = "Log<sub>2</sub> Fold Change in BZ<br>(male vs female)")
ggsave("results/figures/sex_DE_warm_BAT.pdf", plot = sex_warm_BAT_GxS.plot, height = 2, width = 2.1)




##############################################################
# Cold temperature expression patterns
##############################################################

### liver

## filter to get mean of 10 reads across all samples
DE_base_NY_cold_liver_MvF <- res_NY_liver_cold_MvF %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_NY_cold_liver_MvF) # 14,124 genes
#DE_base_NY_cold_liver_MvF %>% filter(padj < 0.05) %>% nrow() # 1,012 genes

DE_base_BZ_cold_liver_MvF <- res_BZ_liver_cold_MvF %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_cold_liver_MvF) # 14,124 genes
#DE_base_BZ_cold_liver_MvF %>% filter(padj < 0.05) %>% nrow() # 678 genes

## merge datasets
merge_sex_cold_liver <- dplyr::full_join(DE_base_BZ_cold_liver_MvF, DE_base_NY_cold_liver_MvF, by = "gene") %>%
  dplyr::select(
    "gene",
    "lfc.BZ" = "log2FoldChange.x",
    "lfcSE.BZ" = "lfcSE.x",
    "padj.BZ" = "padj.x",
    "lfc.NY" = "log2FoldChange.y",
    "lfcSE.NY" = "lfcSE.y",
    "padj.NY" = "padj.y",
  )
#nrow(merge_sex_cold_liver) 14,124 genes
#merge_sex_cold_liver %>% filter(padj.BZ < 0.05 | padj.NY < 0.05) %>% nrow() # 1,306 genes
# roughly 9% of all liver genes show sexual dimorphism in the cold

## add a column indicating whether BZ and NY groups are in the same direction
merge_sex_cold_liver_direction <-  merge_sex_cold_liver %>%
  mutate(direction = case_when((sign(lfc.BZ)) == (sign(lfc.NY)) ~ 1, # same direction
                               (sign(lfc.BZ)) != (sign(lfc.NY)) ~ 0)) # opposite direction

## shrink plotting window by assigning large LFC to 5's
merge_sex_cold_liver_direction$lfc.NY[merge_sex_cold_liver_direction$lfc.NY>5]<-5
merge_sex_cold_liver_direction$lfc.BZ[merge_sex_cold_liver_direction$lfc.BZ>5]<-5
merge_sex_cold_liver_direction$lfc.NY[merge_sex_cold_liver_direction$lfc.NY<(-5)]<-(-5)
merge_sex_cold_liver_direction$lfc.BZ[merge_sex_cold_liver_direction$lfc.BZ<(-5)]<-(-5)

## categorize genes by direction and/or significance
sigNY_sex_cold_liver = dplyr::filter(merge_sex_cold_liver_direction, padj.NY < 0.05 & padj.BZ >= 0.05) # genotype-specific (NY) response = 628 genes
sigBZ_sex_cold_liver = dplyr::filter(merge_sex_cold_liver_direction, padj.BZ < 0.05 & padj.NY >= 0.05) # genotype-specific (BZ) response = 294 genes
sigcons_sex_cold_liver = dplyr::filter(merge_sex_cold_liver_direction, direction == 1 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, and in same directions (this category is hard to
                                                                                                                          # disentangle GxS from G + S) = 382 genes
sigopps_sex_cold_liver = dplyr::filter(merge_sex_cold_liver_direction, direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, but in opposite directions = 2 genes
nonsig_sex_cold_liver = dplyr::filter(merge_sex_cold_liver_direction, padj.NY >= 0.05 & padj.BZ >= 0.05) # 12,817 genes

## collate all genes that show evidence of GxS
# first, genes identified as GxS by DESEQ (very stringent)
DESEQ_sex_liver_cold_gxs <- res_sex_liver_cold_GxS %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10, padj < 0.05)
# save a list of GxS genes defined by DESEQ
DESEQ_sex_liver_cold_gxs_list <- DESEQ_sex_liver_cold_gxs %>% pull(gene)

# second, genes I identify as GxS from categorization above + DESEQ genes
sex_cold_liver_GxS <- merge_sex_cold_liver_direction %>%
  filter((padj.NY < 0.05 & padj.BZ >= 0.05) |
           (padj.BZ < 0.05 & padj.NY >= 0.05) |
           (direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) |
           (gene %in% DESEQ_sex_liver_cold_gxs_list))
# save a list of all GxE genes
sex_cold_liver_GxS_list <- sex_cold_liver_GxS %>% pull(gene)
# 947 GxS genes in warm liver

## double check that DESEQ GxS genes are in my list
#dplyr::anti_join(DESEQ_sex_liver_warm_gxs, sex_warm_liver_GxS, by = "gene")
# yes

## plot
sex_cold_liver_GxS.plot <- merge_sex_cold_liver_direction %>%
  ggplot(aes(x = lfc.NY, y = lfc.BZ)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
  geom_abline(slope = 1, linetype = "dotted", size = 0.25) +
  geom_point(data = nonsig_sex_cold_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "lightgray", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigNY_sex_cold_liver, aes(x = lfc.NY, y = lfc.BZ), fill =  "#1683B5", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigBZ_sex_cold_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigopps_sex_cold_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "#000000", color = "000000", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigcons_sex_cold_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "#B22222", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  scale_x_continuous(limits = c(-5,5), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-5,5),expand = c(0.01,0.075)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans"),
        axis.title = element_blank())
#labs(x = "Log<sub>2</sub> Fold Change in NY<br>(male vs female)",
#y = "Log<sub>2</sub> Fold Change in BZ<br>(male vs female)")
ggsave("results/figures/sex_DE_cold_liver.pdf", plot = sex_cold_liver_GxS.plot, height = 2, width = 2.1)



### BAT

## filter to get mean of 10 reads across all samples
DE_base_NY_cold_BAT_MvF <- res_NY_BAT_cold_MvF %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_NY_cold_BAT_MvF) # 14,261 genes
#DE_base_NY_cold_BAT_MvF %>% filter(padj < 0.05) %>% nrow() # 16 genes

DE_base_BZ_cold_BAT_MvF <- res_BZ_BAT_cold_MvF %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_cold_BAT_MvF) # 14,261 genes
#DE_base_BZ_cold_BAT_MvF %>% filter(padj < 0.05) %>% nrow() # 24 genes

## merge datasets
merge_sex_cold_BAT <- dplyr::full_join(DE_base_BZ_cold_BAT_MvF, DE_base_NY_cold_BAT_MvF, by = "gene") %>%
  dplyr::select(
    "gene",
    "lfc.BZ" = "log2FoldChange.x",
    "lfcSE.BZ" = "lfcSE.x",
    "padj.BZ" = "padj.x",
    "lfc.NY" = "log2FoldChange.y",
    "lfcSE.NY" = "lfcSE.y",
    "padj.NY" = "padj.y",
  )
#nrow(merge_sex_cold_BAT) 14,261 genes
#merge_sex_cold_BAT %>% filter(padj.BZ < 0.05 | padj.NY < 0.05) %>% nrow() # 29 genes
# less than 1% of all BAT genes show sexual dimorphism in the cold

## add a column indicating whether BZ and NY groups are in the same direction
merge_sex_cold_BAT_direction <-  merge_sex_cold_BAT %>%
  mutate(direction = case_when((sign(lfc.BZ)) == (sign(lfc.NY)) ~ 1, # same direction
                               (sign(lfc.BZ)) != (sign(lfc.NY)) ~ 0)) # opposite direction

## shrink plotting window by assigning large LFC to 5's
merge_sex_cold_BAT_direction$lfc.NY[merge_sex_cold_BAT_direction$lfc.NY>5]<-5
merge_sex_cold_BAT_direction$lfc.BZ[merge_sex_cold_BAT_direction$lfc.BZ>5]<-5
merge_sex_cold_BAT_direction$lfc.NY[merge_sex_cold_BAT_direction$lfc.NY<(-5)]<-(-5)
merge_sex_cold_BAT_direction$lfc.BZ[merge_sex_cold_BAT_direction$lfc.BZ<(-5)]<-(-5)

## categorize genes by direction and/or significance
sigNY_sex_cold_BAT = dplyr::filter(merge_sex_cold_BAT_direction, padj.NY < 0.05 & padj.BZ >= 0.05) # genotype-specific (NY) response = 5 genes
sigBZ_sex_cold_BAT = dplyr::filter(merge_sex_cold_BAT_direction, padj.BZ < 0.05 & padj.NY >= 0.05) # genotype-specific (BZ) response = 13 genes
sigcons_sex_cold_BAT = dplyr::filter(merge_sex_cold_BAT_direction, direction == 1 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, and in same directions (this category is hard to
                                                                                                                      # disentangle GxS from G + S) = 11 genes
sigopps_sex_cold_BAT = dplyr::filter(merge_sex_cold_BAT_direction, direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, but in opposite directions = 0 genes
nonsig_sex_cold_BAT = dplyr::filter(merge_sex_cold_BAT_direction, padj.NY >= 0.05 & padj.BZ >= 0.05) # 14,231 genes

## collate all genes that show evidence of GxS
# first, genes identified as GxS by DESEQ (very stringent)
DESEQ_sex_BAT_cold_gxs <- res_sex_BAT_cold_GxS %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10, padj < 0.05)
# 0 genes

# second, genes I identify as GxS from categorization above + DESEQ genes
sex_cold_BAT_GxS <- merge_sex_cold_BAT_direction %>%
  filter((padj.NY < 0.05 & padj.BZ >= 0.05) |
           (padj.BZ < 0.05 & padj.NY >= 0.05) |
           (direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05))
# save a list of all GxE genes
sex_cold_BAT_GxS_list <- sex_cold_BAT_GxS %>% pull(gene)
# 18 GxS genes in warm BAT

## plot
sex_cold_BAT_GxS.plot <- merge_sex_cold_BAT_direction %>%
  ggplot(aes(x = lfc.NY, y = lfc.BZ)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
  geom_abline(slope = 1, linetype = "dotted", size = 0.25) +
  geom_point(data = nonsig_sex_cold_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "lightgray", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigNY_sex_cold_BAT, aes(x = lfc.NY, y = lfc.BZ), fill =  "#1683B5", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigBZ_sex_cold_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigopps_sex_cold_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "#000000", color = "000000", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = sigcons_sex_cold_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "#B22222", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = FALSE) +
  scale_x_continuous(limits = c(-5,5), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-5,5),expand = c(0.01,0.075)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans"),
        axis.title = element_blank())
#labs(x = "Log<sub>2</sub> Fold Change in NY<br>(male vs female)",
#y = "Log<sub>2</sub> Fold Change in BZ<br>(male vs female)")
ggsave("results/figures/sex_DE_cold_liver.pdf", plot = sex_cold_liver_GxS.plot, height = 2, width = 2.1)

