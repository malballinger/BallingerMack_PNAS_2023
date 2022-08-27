#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots cis and trans patterns for each tissue and temperature separately,
# and generates Figures 3B and 3C in BallingerMack_PNAS_2022.

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clear R's environment
library(tidyverse)
library(here)
library(ggtext)
library(ggalluvial)
library(cowplot)
library(DESeq2)

set.seed(19910118)

##############################################################
# Liver cis, trans (both warm and cold)
##############################################################

### warm

liver_warm_male_metadata <- read_delim(here("./data/processed/Liver.MALE.WARM.categories.forplot.txt")) %>%
  select("gene" = "...1",
         "lfc.parental" = "log2FoldChange...2",
         "padj.parental" = "padj...3",
         "lfc.F1" = "log2FoldChange...4",
         "padj.F1" = "padj...5",
         "padj.PvF1" = "padj...6",
         "category" = "Amb&Conserved")

## shrink plotting window by assigning large LFC to 5's
liver_warm_male_metadata$lfc.parental[liver_warm_male_metadata$lfc.parental>3.5]<-3.5
liver_warm_male_metadata$lfc.F1[liver_warm_male_metadata$lfc.F1>3.5]<-3.5
liver_warm_male_metadata$lfc.parental[liver_warm_male_metadata$lfc.parental<(-3.5)]<-(-3.5)
liver_warm_male_metadata$lfc.F1[liver_warm_male_metadata$lfc.F1<(-3.5)]<-(-3.5)

## categorize genes by direction and/or significance
cis_warm_liver = dplyr::filter(liver_warm_male_metadata, category == "CIS_ONLY") # 431 genes
trans_warm_liver = dplyr::filter(liver_warm_male_metadata, category == "TRANS_ONLY") # 244 genes
cistrans_warm_liver = dplyr::filter(liver_warm_male_metadata, category == "Cis&Trans") # 358 genes
nonsig_warm_liver = dplyr::filter(liver_warm_male_metadata, category == "Amb&Conserved") # 4,864 genes

## plot
liver_warm_cistrans <- liver_warm_male_metadata %>%
  ggplot(aes(x = lfc.parental, y = lfc.F1)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_abline(slope = 1, linetype = "dotted") +
  geom_point(data = nonsig_warm_liver, aes(x = lfc.parental, y = lfc.F1), fill = "darkgray", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = cistrans_warm_liver, aes(x = lfc.parental, y = lfc.F1), fill =  "#29AF7F", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = trans_warm_liver, aes(x = lfc.parental, y = lfc.F1), fill = "#3B0F6F", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = cis_warm_liver, aes(x = lfc.parental, y = lfc.F1), fill = "#FC4E07", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  scale_x_continuous(limits = c(-3.5,3.5), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(-3.5,3.5), expand = c(0.01,0.01)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text (size = 10, family = "sans"),
        panel.grid = element_blank())
  # labs(x = "Log<sub>2</sub> Fold Change (NY parent / BZ parent)",
  #      y = "Log<sub>2</sub> Fold Change (NY allele / BZ allele)")

ggsave("results/figures/liver_warm_cistrans.pdf", plot = liver_warm_cistrans, height = 2.75, width = 2.75)


### cold

liver_cold_male_metadata <- read_delim(here("./data/processed/Liver.MALE.COLD.categories.forplot.txt")) %>%
  select("gene" = "...1",
         "lfc.parental" = "log2FoldChange...2",
         "padj.parental" = "padj...3",
         "lfc.F1" = "log2FoldChange...4",
         "padj.F1" = "padj...5",
         "padj.PvF1" = "liverpvals_coldmale_adj",
         "category" = "Amb&Conserved")

## shrink plotting window by assigning large LFC to 5's
liver_cold_male_metadata$lfc.parental[liver_cold_male_metadata$lfc.parental>3.5]<-3.5
liver_cold_male_metadata$lfc.F1[liver_cold_male_metadata$lfc.F1>3.5]<-3.5
liver_cold_male_metadata$lfc.parental[liver_cold_male_metadata$lfc.parental<(-3.5)]<-(-3.5)
liver_cold_male_metadata$lfc.F1[liver_cold_male_metadata$lfc.F1<(-3.5)]<-(-3.5)

## categorize genes by direction and/or significance
cis_cold_liver = dplyr::filter(liver_cold_male_metadata, category == "CIS_ONLY") # 450 genes
trans_cold_liver = dplyr::filter(liver_cold_male_metadata, category == "TRANS_ONLY") # 198 genes
cistrans_cold_liver = dplyr::filter(liver_cold_male_metadata, category == "Cis&Trans") # 343 genes
nonsig_cold_liver = dplyr::filter(liver_cold_male_metadata, category == "Amb&Conserved") # 4,979 genes

## plot
liver_cold_cistrans <- liver_cold_male_metadata %>%
  ggplot(aes(x = lfc.parental, y = lfc.F1)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_abline(slope = 1, linetype = "dotted") +
  geom_point(data = nonsig_cold_liver, aes(x = lfc.parental, y = lfc.F1), fill = "darkgray", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = cistrans_cold_liver, aes(x = lfc.parental, y = lfc.F1), fill =  "#29AF7F", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = trans_cold_liver, aes(x = lfc.parental, y = lfc.F1), fill = "#3B0F6F", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = cis_cold_liver, aes(x = lfc.parental, y = lfc.F1), fill = "#FC4E07", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  scale_x_continuous(limits = c(-3.5,3.5), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(-3.5,3.5), expand = c(0.01,0.01)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text (size = 10, family = "sans"),
        panel.grid = element_blank())

ggsave("results/figures/liver_cold_cistrans.pdf", plot = liver_cold_cistrans, height = 2.75, width = 2.75)


### aluvial plot for liver

WC_liver <- dplyr::full_join(liver_warm_male_metadata, liver_cold_male_metadata, by = "gene") %>%
  select("gene",
         "category_warm" = "category.x",
         "category_cold" = "category.y") %>%
  mutate(across(category_warm, replace_na, "Amb&Conserved"),
         across(category_cold, replace_na, "Amb&Conserved"))

WC_liver$both=paste(WC_liver$category_warm, WC_liver$category_cold, sep="_")

ord <- c("Amb&Conserved", "TRANS_ONLY", "Cis&Trans", "CIS_ONLY") # define order of stratum (for plotting purposes)

WC_liver_format <- WC_liver %>%
  group_by(category_warm, category_cold, both) %>%
  summarize(n = n(), .groups = "drop") %>% # counts the number of genes in the "both" categories
  mutate(freq = n / sum(n)) %>% # calculates the frequency of each category
  filter(both != "Amb&Conserved_Amb&Conserved") %>% # removes the category that shows Amb&Cons across both environments (largest ## of genes)
  mutate(category_warm = factor(category_warm, levels = ord))

cp1 <- c("#FC4E07", "#29AF7F", "#3B0F6F", "darkgray", "#FC4E07", "#29AF7F", "#3B0F6F", "darkgray") # true colors
cp2 <- c("#FC7945", "#4AB08C", "#482870", "gray", "#FC7945", "#4AB08C", "#482870", "gray") # saturated true colors

WC_liver_alluvial <- WC_liver_format %>%
  ggplot(aes(y = n, axis1 = category_warm, axis2 = category_cold)) +
  geom_alluvium(aes(fill= category_warm, color = category_warm), size = 0.7, alpha = 0.6, curve_type = "cubic", width = 1/4, show.legend = FALSE,
                aes.bind = "flows") +
  geom_stratum(aes(order = both), color = cp2, fill = cp2, alpha = 0.99, width = 1/4, show.legend = FALSE) +
  scale_fill_manual(breaks = c("Amb&Conserved", "Cis&Trans", "TRANS_ONLY", "CIS_ONLY"),
                    values = c("darkgray", "#29AF7F", "#3B0F6F", "#FC4E07"),
                    labels = c("Conserved", "Cis + Trans", "Trans", "Cis")) +
  scale_color_manual(breaks = c("Amb&Conserved", "Cis&Trans", "TRANS_ONLY", "CIS_ONLY"),
                    values = c("darkgray", "#29AF7F", "#3B0F6F", "#FC4E07"),
                    labels = c("Conserved", "Cis + Trans", "Trans", "Cis")) +
  scale_x_continuous(breaks = 1:2, labels = c("Warm", "Cold"), expand = c(0.025,0.025)) +
  scale_y_continuous(limits = c(0,1400), breaks = seq(0, 1400, 400),
                     expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text (size = 10, family = "sans"))

ggsave("results/figures/Liver_cistrans_alluvial.pdf", plot = WC_liver_alluvial, height = 2.6, width = 2)



##############################################################
# BAT cis, trans (Dboth warm and cold)
##############################################################

### warm

BAT_warm_male_metadata <- read_delim(here("./data/processed/BAT.MALE.WARM.categories.forplot.txt")) %>%
  select("gene" = "...1",
         "lfc.parental" = "log2FoldChange...2",
         "padj.parental" = "padj...3",
         "lfc.F1" = "log2FoldChange...4",
         "padj.F1" = "padj...5",
         "padj.PvF1" = "batpadj",
         "category" = "Amb&Conserved")

## shrink plotting window by assigning large LFC to 5's
BAT_warm_male_metadata$lfc.parental[BAT_warm_male_metadata$lfc.parental>3.5]<-3.5
BAT_warm_male_metadata$lfc.F1[BAT_warm_male_metadata$lfc.F1>3.5]<-3.5
BAT_warm_male_metadata$lfc.parental[BAT_warm_male_metadata$lfc.parental<(-3.5)]<-(-3.5)
BAT_warm_male_metadata$lfc.F1[BAT_warm_male_metadata$lfc.F1<(-3.5)]<-(-3.5)

## categorize genes by direction and/or significance
cis_warm_BAT = dplyr::filter(BAT_warm_male_metadata, category == "CIS_ONLY") # 452 genes
trans_warm_BAT = dplyr::filter(BAT_warm_male_metadata, category == "TRANS_ONLY") # 242 genes
cistrans_warm_BAT = dplyr::filter(BAT_warm_male_metadata, category == "Cis&Trans") # 378 genes
nonsig_warm_BAT = dplyr::filter(BAT_warm_male_metadata, category == "Amb&Conserved") # 4,825 genes

## plot
BAT_warm_cistrans <- BAT_warm_male_metadata %>%
  ggplot(aes(x = lfc.parental, y = lfc.F1)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_abline(slope = 1, linetype = "dotted") +
  geom_point(data = nonsig_warm_BAT, aes(x = lfc.parental, y = lfc.F1), fill = "darkgray", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = cistrans_warm_BAT, aes(x = lfc.parental, y = lfc.F1), fill =  "#29AF7F", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = trans_warm_BAT, aes(x = lfc.parental, y = lfc.F1), fill = "#3B0F6F", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = cis_warm_BAT, aes(x = lfc.parental, y = lfc.F1), fill = "#FC4E07", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  scale_x_continuous(limits = c(-3.5,3.5), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(-3.5,3.5), expand = c(0.01,0.01)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text (size = 10, family = "sans"),
        panel.grid = element_blank())

ggsave("results/figures/BAT_warm_cistrans.pdf", plot = BAT_warm_cistrans, height = 2.75, width = 2.75)


### cold

BAT_cold_male_metadata <- read_delim(here("./data/processed/BAT.MALE.COLD.categories.forplot.txt")) %>%
  select("gene" = "...1",
         "lfc.parental" = "log2FoldChange...2",
         "padj.parental" = "padj...3",
         "lfc.F1" = "log2FoldChange...4",
         "padj.F1" = "padj...5",
         "padj.PvF1" = "batpvals_coldmale_adj",
         "category" = "Amb&Conserved")

## shrink plotting window by assigning large LFC to 5's
BAT_cold_male_metadata$lfc.parental[BAT_cold_male_metadata$lfc.parental>3.5]<-3.5
BAT_cold_male_metadata$lfc.F1[BAT_cold_male_metadata$lfc.F1>3.5]<-3.5
BAT_cold_male_metadata$lfc.parental[BAT_cold_male_metadata$lfc.parental<(-3.5)]<-(-3.5)
BAT_cold_male_metadata$lfc.F1[BAT_cold_male_metadata$lfc.F1<(-3.5)]<-(-3.5)

## categorize genes by direction and/or significance
cis_cold_BAT = dplyr::filter(BAT_cold_male_metadata, category == "CIS_ONLY") # 478 genes
trans_cold_BAT = dplyr::filter(BAT_cold_male_metadata, category == "TRANS_ONLY") # 170 genes
cistrans_cold_BAT = dplyr::filter(BAT_cold_male_metadata, category == "Cis&Trans") # 312 genes
nonsig_cold_BAT = dplyr::filter(BAT_cold_male_metadata, category == "Amb&Conserved") # 5,010 genes

## plot
BAT_cold_cistrans <- BAT_cold_male_metadata %>%
  ggplot(aes(x = lfc.parental, y = lfc.F1)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_abline(slope = 1, linetype = "dotted") +
  geom_point(data = nonsig_cold_BAT, aes(x = lfc.parental, y = lfc.F1), fill = "darkgray", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = cistrans_cold_BAT, aes(x = lfc.parental, y = lfc.F1), fill =  "#29AF7F", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = trans_cold_BAT, aes(x = lfc.parental, y = lfc.F1), fill = "#3B0F6F", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  geom_point(data = cis_cold_BAT, aes(x = lfc.parental, y = lfc.F1), fill = "#FC4E07", color = "grey100", stroke = 0.3, size = 2.25, shape = 21, alpha = 0.8, show.legend = FALSE) +
  scale_x_continuous(limits = c(-3.5,3.5), expand = c(0.01,0.01)) +
  scale_y_continuous(limits = c(-3.5,3.5), expand = c(0.01,0.01)) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text (size = 10, family = "sans"),
        panel.grid = element_blank())

ggsave("results/figures/BAT_cold_cistrans.pdf", plot = BAT_cold_cistrans, height = 2.75, width = 2.75)


### aluvial plot for BAT

WC_BAT <- dplyr::full_join(BAT_warm_male_metadata, BAT_cold_male_metadata, by = "gene") %>%
  select("gene",
         "category_warm" = "category.x",
         "category_cold" = "category.y") %>%
  mutate(across(category_warm, replace_na, "Amb&Conserved"), # genes that are not expressed are placed in Amb&Conserved category
         across(category_cold, replace_na, "Amb&Conserved"))

WC_BAT$both=paste(WC_BAT$category_warm, WC_BAT$category_cold, sep="_")

ord <- c("Amb&Conserved", "TRANS_ONLY", "Cis&Trans", "CIS_ONLY") # define order of stratum (for plotting purposes)

WC_BAT_format <- WC_BAT %>%
  group_by(category_warm, category_cold, both) %>%
  summarize(n = n(), .groups = "drop") %>%
  mutate(freq = n / sum(n)) %>%
  filter(both != "Amb&Conserved_Amb&Conserved") %>%
  mutate(category_warm = factor(category_warm, levels = ord))

# orange, green, purple
cp1 <- c("#FC4E07", "#29AF7F", "#3B0F6F", "darkgray", "#FC4E07", "#29AF7F", "#3B0F6F", "darkgray") # true colors
cp2 <- c("#FC7945", "#4AB08C", "#482870", "gray", "#FC7945", "#4AB08C", "#482870", "gray") # saturated true colors
cp3 <- c("lightgrey", "lightgrey", "#3B0F6F", "lightgrey", "lightgrey", "lightgrey", "#3B0F6F", "lightgrey")
cp4 <- c("darkgrey", "darkgrey", "#3B0F6F", "darkgrey", "darkgrey", "darkgrey", "#3B0F6F", "darkgrey")
cp5 <- c("#FC4E07", "lightgrey", "lightgrey", "lightgrey", "#FC4E07", "lightgrey", "lightgrey", "lightgrey")
cp6 <- c("#FC4E07", "darkgrey", "darkgrey", "darkgrey", "#FC4E07", "darkgrey", "darkgrey", "darkgrey")

WC_BAT_alluvial <- WC_BAT_format %>%
  ggplot(aes(y = n, axis1 = category_warm, axis2 = category_cold)) +
  geom_alluvium(aes(fill= category_warm, color = category_warm), size = 0.8, alpha = 0.6, curve_type = "cubic", width = 1/4, show.legend = FALSE,
                aes.bind = "flows") +
  geom_stratum(aes(order = both), color = cp2, fill = cp2, alpha = 0.99, width = 1/4, show.legend = FALSE) +
  scale_fill_manual(breaks = c("Amb&Conserved", "Cis&Trans", "TRANS_ONLY", "CIS_ONLY"),
                    values = c("darkgray", "#29AF7F", "#3B0F6F", "#FC4E07"),
                    labels = c("Conserved", "Cis + Trans", "Trans", "Cis")) +
  scale_color_manual(breaks = c("Amb&Conserved", "Cis&Trans", "TRANS_ONLY", "CIS_ONLY"),
                     values = c("darkgray", "#29AF7F", "#3B0F6F", "#FC4E07"),
                     labels = c("Conserved", "Cis + Trans", "Trans", "Cis")) +
  scale_x_continuous(breaks = 1:2, labels = c("Warm", "Cold"), expand = c(0.025,0.025)) +
  scale_y_continuous(limits = c(0,1400), breaks = seq(0, 1400, 400),
                     expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text (size = 10, family = "sans"))

ggsave("results/figures/BAT_cistrans_alluvial.pdf", plot = WC_BAT_alluvial, height = 2.6, width = 2)

WC_BAT_alluvial_transhi <- WC_BAT_format %>%
  ggplot(aes(y = n, axis1 = category_warm, axis2 = category_cold)) +
  geom_alluvium(aes(fill= category_warm, color = category_warm), size = 0.8, alpha = 0.6, curve_type = "cubic", width = 1/4, show.legend = FALSE,
                aes.bind = "flows") +
  geom_stratum(aes(order = both), color = cp4, fill = cp3, alpha = 0.99, width = 1/4, show.legend = FALSE) +
  #geom_text(stat = "stratum", aes(label = after_stat(stratum)),reverse = TRUE) +
  scale_fill_manual(breaks = c("Amb&Conserved", "Cis&Trans", "TRANS_ONLY", "CIS_ONLY"),
                    values = c("darkgray", "darkgray", "#3B0F6F", "darkgray"),
                    labels = c("Conserved", "Cis + Trans", "Trans", "Cis")) +
  scale_color_manual(breaks = c("Amb&Conserved", "Cis&Trans", "TRANS_ONLY", "CIS_ONLY"),
                     values = c("darkgray", "darkgray", "#3B0F6F", "darkgray"),
                     labels = c("Conserved", "Cis + Trans", "Trans", "Cis")) +
  scale_x_continuous(breaks = 1:2, labels = c("warm", "cold"), expand = c(0.025,0.025)) +
  scale_y_continuous(limits = c(0,1400), breaks = seq(0, 1400, 400),
                     expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12)) +
  labs(x = NULL,
       y = "Number of genes")

WC_BAT_alluvial_cishi <- WC_BAT_format %>%
  ggplot(aes(y = n, axis1 = category_warm, axis2 = category_cold)) +
  geom_alluvium(aes(fill= category_warm, color = category_warm), size = 0.8, alpha = 0.6, curve_type = "cubic", width = 1/4, show.legend = FALSE,
                aes.bind = "flows") +
  geom_stratum(aes(order = both), color = cp6, fill = cp5, alpha = 0.99, width = 1/4, show.legend = FALSE) +
  #geom_text(stat = "stratum", aes(label = after_stat(stratum)),reverse = TRUE) +
  scale_fill_manual(breaks = c("Amb&Conserved", "Cis&Trans", "TRANS_ONLY", "CIS_ONLY"),
                    values = c("darkgray", "darkgray", "darkgray", "#FC4E07"),
                    labels = c("Conserved", "Cis + Trans", "Trans", "Cis")) +
  scale_color_manual(breaks = c("Amb&Conserved", "Cis&Trans", "TRANS_ONLY", "CIS_ONLY"),
                     values = c("darkgray", "darkgray", "darkgray", "#FC4E07"),
                     labels = c("Conserved", "Cis + Trans", "Trans", "Cis")) +
  scale_x_continuous(breaks = 1:2, labels = c("warm", "cold"), expand = c(0.025,0.025)) +
  scale_y_continuous(limits = c(0,1400), breaks = seq(0, 1400, 400),
                     expand = c(0,0)) +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text = element_markdown(size = 12),
        axis.title.y = element_markdown(size = 12)) +
  labs(x = NULL,
       y = "Number of genes")

# 
# fig3 <- plot_grid(BAT_warm_cistrans, BAT_cold_cistrans, WC_BAT_alluvial,
#           liver_warm_cistrans, liver_cold_cistrans, WC_liver_alluvial,
#           nrow = 2, ncol = 3, rel_widths = c(2,2,1.5), rel_heights = c(1,1))
# 
# ggsave("results/figures/fig3_test.pdf", plot = fig3, height = 7, width = 5)
