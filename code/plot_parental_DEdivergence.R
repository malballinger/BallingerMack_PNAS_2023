#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots differential expression in parental samples, with specific focus on patterns of divergence.
# This script generates Figures 1E (males) and S3A (females) in BallingerMack_2022.

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
# Male expression patterns (divergence)
##############################################################

### liver

DE_base_warm_males_liver <- res_warm_males_liver_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_warm_males_liver) 14,028 genes
#DE_base_warm_males_liver %>% filter(padj < 0.05) %>% nrow() # 4,271 genes

DE_base_cold_males_liver <- res_cold_males_liver_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_cold_males_liver) 14,028 genes
#DE_base_cold_males_liver %>% filter(padj < 0.05) %>% nrow() # 3,622 genes

## merge datasets
merge_div_males_liver <- dplyr::full_join(DE_base_warm_males_liver, DE_base_cold_males_liver, by = c("gene")) %>%
  dplyr::select(
    "gene",
    "lfc.w" = "log2FoldChange.x",
    "lfcSE.w" = "lfcSE.x",
    "padj.w" = "padj.x",
    "lfc.c" = "log2FoldChange.y",
    "lfcSE.c" = "lfcSE.y",
    "padj.c" = "padj.y"
    )
#nrow(merge_div_males_liver) 14,028 genes
#merge_div_males_liver %>% filter(padj.w < 0.05 | padj.c < 0.05) %>% nrow() # 5,236 genes
# roughly 37% of all liver genes are DE between NY and BZ in males

## add a column indicating whether warm and cold groups are in the same direction
merge_div_males_liver_direction <-  merge_div_males_liver %>%
  mutate(direction = case_when((sign(lfc.w)) == (sign(lfc.c)) ~ 1, # same direction
                               (sign(lfc.w)) != (sign(lfc.c)) ~ 0)) # opposite direction

# shrink plotting window by assigning large LFC to 10's
merge_div_males_liver_direction$lfc.w[merge_div_males_liver_direction$lfc.w>10]<-10
merge_div_males_liver_direction$lfc.c[merge_div_males_liver_direction$lfc.c>10]<-10
merge_div_males_liver_direction$lfc.w[merge_div_males_liver_direction$lfc.w<(-10)]<-(-10)
merge_div_males_liver_direction$lfc.c[merge_div_males_liver_direction$lfc.c<(-10)]<-(-10)

# categorize genes by direction and/or significance
sigcold_liver = dplyr::filter(merge_div_males_liver_direction, padj.c < 0.05 & padj.w >= 0.05) # 965 genes
sigconserve_liver = dplyr::filter(merge_div_males_liver_direction, direction == 1 & padj.w < 0.05 & padj.c < 0.05) # 2,657 genes
sigwarm_liver = dplyr::filter(merge_div_males_liver_direction, padj.w < 0.05 & padj.c >= 0.05) # 1,614 genes
sigopps_liver = dplyr::filter(merge_div_males_liver_direction, direction == 0 & padj.w < 0.05 & padj.c < 0.05)
nonsig_liver = dplyr::filter(merge_div_males_liver_direction, padj.w >= 0.05 & padj.c >= 0.05)

# plot
liver_males_divergence <- merge_div_males_liver_direction %>%
  ggplot(aes(x = lfc.w, y = lfc.c)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_abline(slope = 1, linetype = "dotted") +
  geom_point(data = sigconserve_liver, aes(x = lfc.w, y = lfc.c), fill = "#44B379", color = "grey100", stroke = 0.4, size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = nonsig_liver, aes(x = lfc.w, y = lfc.c), fill = "darkgray", color = "grey100", stroke = 0.4, size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = sigwarm_liver, aes(x = lfc.w, y = lfc.c), fill = "#E8E430", color = "grey100", stroke = 0.4, size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = sigcold_liver, aes(x = lfc.w, y = lfc.c), fill = "#452D72", color = "grey100", stroke = 0.4, size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = sigopps_liver, aes(x = lfc.w, y = lfc.c), fill = "#000000", color = "000000", size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
  scale_x_continuous(limits = c(-10,10), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-10,10), expand = c(0.01,0.075)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text (size = 12.5, family = "sans"),
        axis.title = element_blank())
  #       axis.title.y = element_markdown(vjust = 1, size = 14, family = "sans"),
  #       axis.title.x = element_markdown(vjust = -1.5, size = 14, family = "sans")) +
  # labs(x = "Log<sub>2</sub> Fold Change in warm<br>(NY vs BZ)",
  #      y = "Log<sub>2</sub> Fold Change in cold<br>(NY vs BZ)")
ggsave("results/figures/males_DE_divergence_liver.pdf", plot = liver_males_divergence, height = 3.25, width = 3.5)


### BAT

DE_base_warm_males_BAT <- res_warm_males_BAT_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_warm_males_BAT) 14,176 genes
#DE_base_warm_males_BAT %>% filter(padj < 0.05) %>% nrow() # 3,713 genes

DE_base_cold_males_BAT <- res_cold_males_BAT_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_cold_males_BAT) 14,176 genes
#DE_base_cold_males_BAT %>% filter(padj < 0.05) %>% nrow() # 3,662 genes

## merge datasets
merge_div_males_BAT <- dplyr::full_join(DE_base_warm_males_BAT, DE_base_cold_males_BAT, by = c("gene")) %>%
  dplyr::select(
    "gene",
    "lfc.w" = "log2FoldChange.x",
    "lfcSE.w" = "lfcSE.x",
    "padj.w" = "padj.x",
    "lfc.c" = "log2FoldChange.y",
    "lfcSE.c" = "lfcSE.y",
    "padj.c" = "padj.y"
  )
# 14,176 genes
#merge_div_males_BAT %>% filter(padj.w < 0.05 | padj.c < 0.05) %>% nrow() # 4,853 genes
# roughly 34% of all BAT genes are DE between NY and BZ in males

## add a column indicating whether warm and cold groups are in the same direction
merge_div_males_BAT_direction <-  merge_div_males_BAT %>%
  mutate(direction = case_when((sign(lfc.w)) == (sign(lfc.c)) ~ 1, # same direction
                               (sign(lfc.w)) != (sign(lfc.c)) ~ 0)) # opposite direction


# shrink plotting window by assigning large LFC to 10's
merge_div_males_BAT_direction$lfc.w[merge_div_males_BAT_direction$lfc.w>10]<-10
merge_div_males_BAT_direction$lfc.c[merge_div_males_BAT_direction$lfc.c>10]<-10
merge_div_males_BAT_direction$lfc.w[merge_div_males_BAT_direction$lfc.w<(-10)]<-(-10)
merge_div_males_BAT_direction$lfc.c[merge_div_males_BAT_direction$lfc.c<(-10)]<-(-10)

# categorize genes by direction and/or significance
sigcold_BAT = dplyr::filter(merge_div_males_BAT_direction, padj.c < 0.05 & padj.w >= 0.05) # 1,140 genes
sigconserve_BAT = dplyr::filter(merge_div_males_BAT_direction, direction == 1 & padj.w < 0.05 & padj.c < 0.05) # 2,515 genes
sigwarm_BAT = dplyr::filter(merge_div_males_BAT_direction, padj.w < 0.05 & padj.c >= 0.05) # 1,191 genes
sigopps_BAT = dplyr::filter(merge_div_males_BAT_direction, direction == 0 & padj.w < 0.05 & padj.c < 0.05)
nonsig_BAT = dplyr::filter(merge_div_males_BAT_direction, padj.w >= 0.05 & padj.c >= 0.05)

# plot
BAT_males_divergence <- merge_div_males_BAT_direction %>%
  ggplot(aes(x = lfc.w, y = lfc.c)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_abline(slope = 1, linetype = "dotted") +
  geom_point(data = sigconserve_BAT, aes(x = lfc.w, y = lfc.c), fill = "#44B379", color = "grey100", stroke = 0.4, size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = nonsig_BAT, aes(x = lfc.w, y = lfc.c), fill = "darkgray", color = "grey100", stroke = 0.4, size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = sigwarm_BAT, aes(x = lfc.w, y = lfc.c), fill = "#E8E430", color = "grey100", stroke = 0.4, size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = sigcold_BAT, aes(x = lfc.w, y = lfc.c), fill = "#452D72", color = "grey100", stroke = 0.4, size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = sigopps_BAT, aes(x = lfc.w, y = lfc.c), fill = "#000000", color = "000000", size = 3.5, shape = 21, alpha = 1, show.legend = FALSE) +
  scale_x_continuous(limits = c(-10,10), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-10,10), expand = c(0.01,0.075)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text (size = 12.5, family = "sans"),
        axis.title = element_blank())
#       axis.title.y = element_markdown(vjust = 1, size = 14, family = "sans"),
#       axis.title.x = element_markdown(vjust = -1.5, size = 14, family = "sans")) +
# labs(x = "Log<sub>2</sub> Fold Change in warm<br>(NY vs BZ)",
#      y = "Log<sub>2</sub> Fold Change in cold<br>(NY vs BZ)")
ggsave("results/figures/males_DE_divergence_BAT.pdf", plot = BAT_males_divergence, height = 3.25, width = 3.5)




##############################################################
# Female expression patterns (divergence)
##############################################################

### liver

DE_base_warm_females_liver <- res_warm_females_liver_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_warm_females_liver) 14,163 genes
#DE_base_warm_females_liver %>% filter(padj < 0.05) %>% nrow() # 3,213 genes

DE_base_cold_females_liver <- res_cold_females_liver_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_cold_females_liver) 14,163 genes
#DE_base_cold_females_liver %>% filter(padj < 0.05) %>% nrow() # 3,959 genes

## merge datasets
merge_div_females_liver <- dplyr::full_join(DE_base_warm_females_liver, DE_base_cold_females_liver, by = c("gene")) %>%
  dplyr::select(
    "gene",
    "lfc.w" = "log2FoldChange.x",
    "lfcSE.w" = "lfcSE.x",
    "padj.w" = "padj.x",
    "lfc.c" = "log2FoldChange.y",
    "lfcSE.c" = "lfcSE.y",
    "padj.c" = "padj.y"
  )
#nrow(merge_div_females_liver) 14,163 genes
#merge_div_females_liver %>% filter(padj.w < 0.05 | padj.c < 0.05) %>% nrow() # 4,781 genes
# roughly 34% of all liver genes are DE between NY and BZ in females

## add a column indicating whether warm and cold groups are in the same direction
merge_div_females_liver_direction <-  merge_div_females_liver %>%
  mutate(direction = case_when((sign(lfc.w)) == (sign(lfc.c)) ~ 1, # same direction
                               (sign(lfc.w)) != (sign(lfc.c)) ~ 0)) # opposite direction

# shrink plotting window by assigning large LFC to 10's
merge_div_females_liver_direction$lfc.w[merge_div_females_liver_direction$lfc.w>10]<-10
merge_div_females_liver_direction$lfc.c[merge_div_females_liver_direction$lfc.c>10]<-10
merge_div_females_liver_direction$lfc.w[merge_div_females_liver_direction$lfc.w<(-10)]<-(-10)
merge_div_females_liver_direction$lfc.c[merge_div_females_liver_direction$lfc.c<(-10)]<-(-10)

# categorize genes by direction and/or significance
sigcold_liver_f = dplyr::filter(merge_div_females_liver_direction, padj.c < 0.05 & padj.w >= 0.05) # 1,568 genes
sigconserve_liver_f = dplyr::filter(merge_div_females_liver_direction, direction == 1 & padj.w < 0.05 & padj.c < 0.05) # 2,391 genes
sigwarm_liver_f = dplyr::filter(merge_div_females_liver_direction, padj.w < 0.05 & padj.c >= 0.05) # 822 genes
sigopps_liver_f = dplyr::filter(merge_div_females_liver_direction, direction == 0 & padj.w < 0.05 & padj.c < 0.05)
nonsig_liver_f = dplyr::filter(merge_div_females_liver_direction, padj.w >= 0.05 & padj.c >= 0.05)

# plot
liver_females_divergence <- merge_div_females_liver_direction %>%
  ggplot(aes(x = lfc.w, y = lfc.c)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
  geom_abline(slope = 1, linetype = "dotted", size = 0.25) +
  geom_point(data = sigconserve_liver_f, aes(x = lfc.w, y = lfc.c), fill = "#44B379", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = nonsig_liver_f, aes(x = lfc.w, y = lfc.c), fill = "darkgray", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = sigwarm_liver_f, aes(x = lfc.w, y = lfc.c), fill = "#E8E430", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = sigcold_liver_f, aes(x = lfc.w, y = lfc.c), fill = "#452D72", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = sigopps_liver_f, aes(x = lfc.w, y = lfc.c), fill = "#000000", color = "000000", size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
  scale_x_continuous(limits = c(-10,10), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-10,10), expand = c(0.01,0.075)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans"),
        axis.title = element_blank())
#       axis.title.y = element_markdown(vjust = 1, size = 14, family = "sans"),
#       axis.title.x = element_markdown(vjust = -1.5, size = 14, family = "sans")) +
# labs(x = "Log<sub>2</sub> Fold Change in warm<br>(NY vs BZ)",
#      y = "Log<sub>2</sub> Fold Change in cold<br>(NY vs BZ)")
ggsave("results/figures/females_DE_divergence_liver.pdf", plot = liver_females_divergence, height = 2, width = 2.1)


### BAT

DE_base_warm_females_BAT <- res_warm_females_BAT_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_warm_females_BAT) 14,261 genes
#DE_base_warm_females_BAT %>% filter(padj < 0.05) %>% nrow() # 3,052 genes

DE_base_cold_females_BAT <- res_cold_females_BAT_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_cold_females_BAT) 14,261 genes
#DE_base_cold_females_BAT %>% filter(padj < 0.05) %>% nrow() # 3,992 genes

## merge datasets
merge_div_females_BAT <- dplyr::full_join(DE_base_warm_females_BAT, DE_base_cold_females_BAT, by = c("gene")) %>%
  dplyr::select(
    "gene",
    "lfc.w" = "log2FoldChange.x",
    "lfcSE.w" = "lfcSE.x",
    "padj.w" = "padj.x",
    "lfc.c" = "log2FoldChange.y",
    "lfcSE.c" = "lfcSE.y",
    "padj.c" = "padj.y"
  )
#nrow(merge_div_females_BAT) 14,261 genes
#merge_div_females_BAT %>% filter(padj.w < 0.05 | padj.c < 0.05) %>% nrow() # 4,834 genes
# roughly 34% of all BAT genes are DE between NY and BZ in females

## add a column indicating whether warm and cold groups are in the same direction
merge_div_females_BAT_direction <-  merge_div_females_BAT %>%
  mutate(direction = case_when((sign(lfc.w)) == (sign(lfc.c)) ~ 1, # same direction
                               (sign(lfc.w)) != (sign(lfc.c)) ~ 0)) # opposite direction


# shrink plotting window by assigning large LFC to 10's
merge_div_females_BAT_direction$lfc.w[merge_div_females_BAT_direction$lfc.w>10]<-10
merge_div_females_BAT_direction$lfc.c[merge_div_females_BAT_direction$lfc.c>10]<-10
merge_div_females_BAT_direction$lfc.w[merge_div_females_BAT_direction$lfc.w<(-10)]<-(-10)
merge_div_females_BAT_direction$lfc.c[merge_div_females_BAT_direction$lfc.c<(-10)]<-(-10)

# categorize genes by direction and/or significance
sigcold_BAT_f = dplyr::filter(merge_div_females_BAT_direction, padj.c < 0.05 & padj.w >= 0.05) # 1,782 genes
sigconserve_BAT_f = dplyr::filter(merge_div_females_BAT_direction, direction == 1 & padj.w < 0.05 & padj.c < 0.05) # 2,171 genes
sigwarm_BAT_f = dplyr::filter(merge_div_females_BAT_direction, padj.w < 0.05 & padj.c >= 0.05) # 842 genes
sigopps_BAT_f = dplyr::filter(merge_div_females_BAT_direction, direction == 0 & padj.w < 0.05 & padj.c < 0.05) # 39 genes
nonsig_BAT_f = dplyr::filter(merge_div_females_BAT_direction, padj.w >= 0.05 & padj.c >= 0.05)

# plot
BAT_females_divergence <- merge_div_females_BAT_direction %>%
  ggplot(aes(x = lfc.w, y = lfc.c)) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.25) +
  geom_hline(yintercept = 0, linetype = "dotted", size = 0.25) +
  geom_abline(slope = 1, linetype = "dotted", size = 0.25) +
  geom_point(data = sigconserve_BAT_f, aes(x = lfc.w, y = lfc.c), fill = "#44B379", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = nonsig_BAT_f, aes(x = lfc.w, y = lfc.c), fill = "darkgray", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = sigwarm_BAT_f, aes(x = lfc.w, y = lfc.c), fill = "#E8E430", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = sigcold_BAT_f, aes(x = lfc.w, y = lfc.c), fill = "#452D72", color = "grey100", stroke = 0.25, size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
  geom_point(data = sigopps_BAT_f, aes(x = lfc.w, y = lfc.c), fill = "#000000", color = "000000", size = 2, shape = 21, alpha = 1, show.legend = FALSE) +
  scale_x_continuous(limits = c(-10,10), expand = c(0.01,0.075)) +
  scale_y_continuous(limits = c(-10,10), expand = c(0.01,0.075)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 7.5, family = "sans"),
        axis.title = element_blank())
#       axis.title.y = element_markdown(vjust = 1, size = 14, family = "sans"),
#       axis.title.x = element_markdown(vjust = -1.5, size = 14, family = "sans")) +
# labs(x = "Log<sub>2</sub> Fold Change in warm<br>(NY vs BZ)",
#      y = "Log<sub>2</sub> Fold Change in cold<br>(NY vs BZ)")
ggsave("results/figures/females_DE_divergence_BAT.pdf", plot = BAT_females_divergence, height = 2, width = 2.1)
