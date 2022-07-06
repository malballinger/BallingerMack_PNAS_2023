#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots differential expression in parental samples, with
# specific focus on patterns of GxE.
# This script generates Figures X and SX in BallingerMack_PNAS_2022.

##############################################################
# Required packages
##############################################################

rm(list = ls()) # clear R's environment
library(tidyverse)
library(here)
library(ggtext)
library(cowplot)
library(DESeq2)

set.seed(19910118)

source("./code/get_DE_males_parents.R") # where DESeq data frames are generated

##############################################################
# Male expression patterns (divergence)
##############################################################

### liver

## filter to get mean of 10 reads across all samples
DEsigbase_NY_males_liver <- res_NY_males_liver_WvC %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#14,028 genes

DEsigbase_BZ_males_liver <- res_BZ_males_liver_WvC %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#14,028 genes

## merge datasets
merge_plast_males_liver <- dplyr::inner_join(DEsigbase_BZ_males_liver, DEsigbase_NY_males_liver, by = "gene") %>%
  select(
    "gene",
    "lfc.BZ" = "log2FoldChange.x",
    "lfcSE.BZ" = "lfcSE.x",
    "padj.BZ" = "padj.x",
    "lfc.NY" = "log2FoldChange.y",
    "lfcSE.NY" = "lfcSE.y",
    "padj.NY" = "padj.y",
    )
# 14,028 genes

## add a column indicating whether warm and cold groups are in the same direction
merge_plast_males_liver_direction <-  merge_plast_males_liver %>%
  mutate(direction = case_when((sign(lfc.BZ)) == (sign(lfc.NY)) ~ 1, # same direction
                               (sign(lfc.BZ)) != (sign(lfc.NY)) ~ 0)) # opposite direction

## categorize genes by direction and/or significance
sigNY_liver = dplyr::filter(merge_plast_males_liver_direction, padj.NY < 0.05 & padj.BZ >= 0.05) # pop-specific (NY) response = 485 genes
sigBZ_liver = dplyr::filter(merge_plast_males_liver_direction, padj.BZ < 0.05 & padj.NY >= 0.05) # pop-specific (BZ) response = 936 genes
sigcons_liver = dplyr::filter(merge_plast_males_liver_direction, direction == 1 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, and in same directions (this category is hard to
                                                                                                                # dissentangle GxE from G + E) = 396 genes
sigopps_liver = dplyr::filter(merge_plast_males_liver_direction, direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, but in opposite directions = 2 genes
nonsig_liver = dplyr::filter(merge_plast_males_liver_direction, padj.NY >= 0.05 & padj.BZ >= 0.05) # 12,338 genes

## collate all genes that show evidence of GxE
# first, genes identified as GxE by DESEQ (very stringent)
DESEQ_liver_gxe <- res_males_liver_GxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10,
         padj < 0.05)
# save a list of GxE genes defined by DESEQ
DESEQ_liver_gxe_list <- DESEQ_liver_gxe %>% pull(gene)

# second, genes I identify as GxE from categorization above + DESEQ genes
liver_males_plast_GxE <- merge_plast_males_liver_direction %>%
  filter((padj.NY < 0.05 & padj.BZ >= 0.05) |
           (padj.BZ < 0.05 & padj.NY >= 0.05) |
           (direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) |
           (gene %in% DESEQ_liver_gxe_list))
# save a list of all GxE genes
liver_males_plast_GxE_list <- liver_males_plast_GxE %>% pull(gene)
# 1,429 GxE genes in liver

## double check that DESEQ GxE genes are in my list
#dplyr::anti_join(DESEQ_liver_gxe, liver_males_plast_GxE, by = "gene")
# yes

## shrink plotting window by assigning large LFC to 5's
merge_plast_males_liver_direction$lfc.NY[merge_plast_males_liver_direction$lfc.NY>5]<-5
merge_plast_males_liver_direction$lfc.BZ[merge_plast_males_liver_direction$lfc.BZ>5]<-5
merge_plast_males_liver_direction$lfc.NY[merge_plast_males_liver_direction$lfc.NY<(-5)]<-(-5)
merge_plast_males_liver_direction$lfc.BZ[merge_plast_males_liver_direction$lfc.BZ<(-5)]<-(-5)

## plot
liver_males_plasticity <- merge_plast_males_liver_direction %>%
  ggplot(aes(x = lfc.NY, y = lfc.BZ)) +
  geom_point(data = nonsig_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "lightgray", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_point(data = sigNY_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "#388BFF", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_point(data = sigBZ_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_point(data = sigopps_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "#000000", color = "000000", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_point(data = sigcons_liver, aes(x = lfc.NY, y = lfc.BZ), fill = "#B22222", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_abline(slope = 1, linetype = "dotted") +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  coord_cartesian(xlim = c(-5,5),
                  ylim = c(-5,5)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank())



### BAT

## filter to get mean of 10 reads across all samples
DEsigbase_NY_males_BAT <- res_NY_males_BAT_WvC %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#14,176 genes

DEsigbase_BZ_males_BAT <- res_BZ_males_BAT_WvC %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#14,176 genes

## merge datasets
merge_plast_males_BAT <- dplyr::inner_join(DEsigbase_BZ_males_BAT, DEsigbase_NY_males_BAT, by = "gene") %>%
  select(
    "gene",
    "lfc.BZ" = "log2FoldChange.x",
    "lfcSE.BZ" = "lfcSE.x",
    "padj.BZ" = "padj.x",
    "lfc.NY" = "log2FoldChange.y",
    "lfcSE.NY" = "lfcSE.y",
    "padj.NY" = "padj.y",
  )
# 14,176 genes

## add a column indicating whether warm and cold groups are in the same direction
merge_plast_males_BAT_direction <-  merge_plast_males_BAT %>%
  mutate(direction = case_when((sign(lfc.BZ)) == (sign(lfc.NY)) ~ 1, # same direction
                               (sign(lfc.BZ)) != (sign(lfc.NY)) ~ 0)) # opposite direction

## categorize genes by direction and/or significance
sigNY_BAT = dplyr::filter(merge_plast_males_BAT_direction, padj.NY < 0.05 & padj.BZ >= 0.05) # pop-specific (NY) response
sigBZ_BAT = dplyr::filter(merge_plast_males_BAT_direction, padj.BZ < 0.05 & padj.NY >= 0.05) # pop-specific (BZ) response
sigcons_BAT = dplyr::filter(merge_plast_males_BAT_direction, direction == 1 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, and in same directions (this category is hard to
                                                                                                          # dissentangle GxE from G + E)
sigopps_BAT = dplyr::filter(merge_plast_males_BAT_direction, direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) # both pops show DE, but in opposite directions
nonsig_BAT = dplyr::filter(merge_plast_males_BAT_direction, padj.NY >= 0.05 & padj.BZ >= 0.05)

## collate all genes that show evidence of GxE
# first, genes identified as GxE by DESEQ (very stringent)
DESEQ_BAT_gxe <- res_males_BAT_GxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10, padj < 0.05)
# save a list of GxE genes defined by DESEQ
DESEQ_BAT_gxe_list <- DESEQ_BAT_gxe %>% pull(gene)

# second, genes I identify as GxE from categorization above + DESEQ genes
BAT_males_plast_GxE <- merge_plast_males_BAT_direction %>%
  filter((padj.NY < 0.05 & padj.BZ >= 0.05) |
           (padj.BZ < 0.05 & padj.NY >= 0.05) |
           (direction == 0 & padj.NY < 0.05 & padj.BZ < 0.05) |
           (gene %in% DESEQ_BAT_gxe_list))
# save a list of all GxE genes
BAT_males_plast_GxE_list <- BAT_males_plast_GxE %>% pull(gene)

## double check that DESEQ GxE genes are in my list
#dplyr::anti_join(DESEQ_BAT_gxe, BAT_males_plast_GxE, by = "gene")
# yes

## shrink plotting window by assigning large LFC to 5's
merge_plast_males_BAT_direction$lfc.NY[merge_plast_males_BAT_direction$lfc.NY>5]<-5
merge_plast_males_BAT_direction$lfc.BZ[merge_plast_males_BAT_direction$lfc.BZ>5]<-5
merge_plast_males_BAT_direction$lfc.NY[merge_plast_males_BAT_direction$lfc.NY<(-5)]<-(-5)
merge_plast_males_BAT_direction$lfc.BZ[merge_plast_males_BAT_direction$lfc.BZ<(-5)]<-(-5)

## plot
BAT_males_plasticity <- merge_plast_males_BAT_direction %>%
  ggplot(aes(x = lfc.NY, y = lfc.BZ)) +
  geom_point(data = nonsig_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "lightgray", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_point(data = sigNY_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "#388BFF", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_point(data = sigBZ_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "#E09832", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_point(data = sigopps_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "#000000", color = "000000", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_point(data = sigcons_BAT, aes(x = lfc.NY, y = lfc.BZ), fill = "#B22222", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 0.8, show.legend = TRUE) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_abline(slope = 1, linetype = "dotted") +
  scale_x_continuous(expand = c(0.01,0.01)) +
  scale_y_continuous(expand = c(0.01,0.01)) +
  coord_cartesian(ylim = c(-5,5)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid = element_blank())

cowplot::plot_grid(BAT_males_plasticity, liver_males_plasticity,
                        ncol = 2, nrow = 1)
ggsave("results/figures/males_DE_plasticity.png", height = 3, width = 5)









### how to add axes titles manually in cowplot:

# t <- add_sub(plot, expression("Log"[2]*~"Fold Change (BZ vs NY) in Warm"), hjust = 0)
# ggdraw(t)
# p <- add_sub(t, expression("Log"[2]*~"Fold Change (BZ vs NY) in Cold"), 0, 1.5, angle=90)
# ggdraw(p)
