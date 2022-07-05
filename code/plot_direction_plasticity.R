#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots patterns of adaptive (same direction) and non-adaptive (opposite direction) plasticity
# in male parental samples. Specifically,we ask if the plastic response of Brazil mice goes in the same direction as the evolved response
# of New York mice at room temperature.
# This script generates Figure 2C in BallingerMack_2022.

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

source("./code/get_DE_males_parents.R") # where DESeq datasets are generated

##############################################################
# Liver expression (DE in warm (NY vs BZ), and plasticity in Brazil (BZ warm vs cold))
##############################################################

### New York vs Brazil (in the warm -- evolved differences)
# filter to get mean of 10 reads across all samples
DE_base_evol_liver_males <- res_warm_males_liver_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_evol_liver)   # 14,028 genes


### Brazil plasticity (cold vs warm):
# filter to get mean of 10 reads across all samples
DE_base_BZ_plast_liver_males <- res_BZ_males_liver_CvW %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_plast_liver)   # 14,028 genes


### merge evol DE to plast Brazil DE in liver
merge_EvP_liver_males <- dplyr::full_join(DE_base_evol_liver_males, DE_base_BZ_plast_liver_males, by = "gene") %>%
  select("gene",
         "lfc.e" = "log2FoldChange.x",
         "padj.e" = "padj.x",
         "lfc.p" = "log2FoldChange.y",
         "padj.p" = "padj.y")
# 14,028 genes


## add a column indicating whether warm and cold groups are in the same direction
merge_EvP_liver_males_direction <-  merge_EvP_liver_males %>%
  mutate(direction = case_when((sign(lfc.e)) == (sign(lfc.p)) ~ 1, # same direction
                               (sign(lfc.e)) != (sign(lfc.p)) ~ 0)) # opposite direction

## shrink plotting window by assigning large LFC to 5's
# merge_EvP_liver_direction$lfc.p[merge_EvP_liver_direction$lfc.p>5]<-5
# merge_EvP_liver_direction$lfc.e[merge_EvP_liver_direction$lfc.e>5]<-5
# merge_EvP_liver_direction$lfc.p[merge_EvP_liver_direction$lfc.p<(-5)]<-(-5)
# merge_EvP_liver_direction$lfc.e[merge_EvP_liver_direction$lfc.e<(-5)]<-(-5)

## pull out categories of genes that we are interested in
sigevol_liver_males = dplyr::filter(merge_EvP_liver_males_direction, padj.e < 0.05 & padj.p >= 0.05) # 3,864 genes
sigplas_liver_males = dplyr::filter(merge_EvP_liver_males_direction, padj.p < 0.05 & padj.e >= 0.05) # 410 genes
sigsame_liver_males = dplyr::filter(merge_EvP_liver_males_direction, direction == 1 & padj.e < 0.05 & padj.p < 0.05) # 296 genes
sigopps_liver_males = dplyr::filter(merge_EvP_liver_males_direction, direction == 0 & padj.e < 0.05 & padj.p < 0.05) # 111 genes
sigboth_liver_males = dplyr::filter(merge_EvP_liver_males_direction, padj.e < 0.05 & padj.p < 0.05) # 407 genes

### plot adaptive and non-adaptive plasticity
all_BZ_liver_males.plot <-
  ggplot(data = merge_EvP_liver_males_direction, aes(x = lfc.p, y = lfc.e)) +
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
ggsave("results/figures/males_adap_plast_liver.pdf", plot = all_BZ_liver_males.plot, height = 2, width = 2.1)


##############################################################
# Correlation Test - Liver
##############################################################

# check normality of liver data
liver.model_males <- lm(lfc.e ~ lfc.p,
                     data = sigboth_liver_males)
# check_model(liver.model_males) # data are not normally distributed - use Spearman correlation

r.obs_liver_males <- cor(x = sigboth_liver_males$lfc.p, y = sigboth_liver_males$lfc.e, method = "spearman")
P.cor_liver_males <- cor.test(x = sigboth_liver_males$lfc.p, y = sigboth_liver_males$lfc.e, method = "spearman")$p.value


cor.perm <- function (data = data, nperm = 10000)
{
  random_subsample <- dplyr::slice_sample(data, n=407) # randomly subsample the number of sigboth
  r.per <- sapply (1:nperm, FUN = function (i) cor (x = random_subsample$lfc.p, y = sample (random_subsample$lfc.e)))
  r.per <- c(r.per, r.obs_liver_males)
  hist (r.per, breaks = 100, xlab = "Correlation coefficient", main=NULL,
        col = "white")
  abline (v = r.obs_liver_males, col = 'red', lty = "dotted")
  P.per <- sum (abs (r.per) >= abs (r.obs_liver_males))/(nperm + 1) 
  return (list (r.obs_liver_males = r.obs_liver_males, P.cor_liver_males = P.cor_liver_males, P.per = P.per, r.per = r.per))
}

liver_perm_males <- cor.perm(data = merge_EvP_liver_males)
liver_p_per <- liver_perm_males$P.per # P < 0.05

# make liver_perm a df so you can plot hist in ggplot
df_liver_perm_males <- data.frame(liver_perm_males)
# make histogram
liver_perm_males.inset <-
  ggplot(df_liver_perm_males, aes(x=r.per)) +
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
       y = "Frequency")
ggsave("results/figures/males_adap_plast_liver_inset.pdf", plot = liver_perm_males.inset,
       height = 1.5, width = 2.25, units = "cm")

# join liver adaptive/nonadaptive plot with permutation model
# liver_nonadap <-
#   ggdraw() +
#   draw_plot(liver.plot) +
#   draw_plot(liver_perm.inset, x= 0.12, y = 0.15, width = 0.28, height = 0.2)



##############################################################
# BAT expression (DE in warm (NY vs BZ), and plasticity in Brazil (BZ warm vs cold))   --- REP 1
##############################################################

### New York vs Brazil (in the warm -- evolved differences)
# filter to get mean of 10 reads across all samples
DE_base_evol_BAT_males <- res_warm_males_BAT_NYvsBZ %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_evol_BAT_males)   # 14,176 genes


### Brazil plasticity (cold vs warm):
# filter to get mean of 10 reads across all samples
DE_base_BZ_plast_BAT_males <- res_BZ_males_BAT_CvW %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10)
#nrow(DE_base_BZ_plast_BAT_males)   # 14,176 genes


### merge evol DE to plast Brazil DE in BAT
merge_EvP_BAT_males <- dplyr::full_join(DE_base_evol_BAT_males, DE_base_BZ_plast_BAT_males, by = "gene") %>%
  select("gene",
         "lfc.e" = "log2FoldChange.x",
         "padj.e" = "padj.x",
         "lfc.p" = "log2FoldChange.y",
         "padj.p" = "padj.y")
# 14,176 genes


## add a column indicating whether warm and cold groups are in the same direction
merge_EvP_BAT_males_direction <-  merge_EvP_BAT_males %>%
  mutate(direction = case_when((sign(lfc.e)) == (sign(lfc.p)) ~ 1, # same direction
                               (sign(lfc.e)) != (sign(lfc.p)) ~ 0)) # opposite direction

## shrink plotting window by assigning large LFC to 5's
# merge_EvP_liver_direction$lfc.p[merge_EvP_liver_direction$lfc.p>5]<-5
# merge_EvP_liver_direction$lfc.e[merge_EvP_liver_direction$lfc.e>5]<-5
# merge_EvP_liver_direction$lfc.p[merge_EvP_liver_direction$lfc.p<(-5)]<-(-5)
# merge_EvP_liver_direction$lfc.e[merge_EvP_liver_direction$lfc.e<(-5)]<-(-5)

## pull out categories of genes that we are interested in
sigevol_BAT_males = dplyr::filter(merge_EvP_BAT_males_direction, padj.e < 0.05 & padj.p >= 0.05) # 3,215 genes
sigplas_BAT_males = dplyr::filter(merge_EvP_BAT_males_direction, padj.p < 0.05 & padj.e >= 0.05) # 836 genes
sigsame_BAT_males = dplyr::filter(merge_EvP_BAT_males_direction, direction == 1 & padj.e < 0.05 & padj.p < 0.05) # 295 genes
sigopps_BAT_males = dplyr::filter(merge_EvP_BAT_males_direction, direction == 0 & padj.e < 0.05 & padj.p < 0.05) # 203 genes
sigboth_BAT_males = dplyr::filter(merge_EvP_BAT_males_direction, padj.e < 0.05 & padj.p < 0.05) # 498 genes


### plot adaptive and non-adaptive plasticity
BAT.plot <-
  ggplot(data = merge_EvP_BAT_males_direction, aes(x = lfc.p, y = lfc.e)) +
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
ggsave("results/figures/males_adap_plast_BAT.pdf", plot = BAT.plot, height = 2, width = 2.1)


##############################################################
# Correlation Test - BAT
##############################################################

# check normality of liver data
BAT.model_males <- lm(lfc.e ~ lfc.p,
                  data = sigboth_BAT_males)
# check_model(BAT.model_males) # data are not normally distributed - use Spearman correlation

r.obs_BAT_males <- cor(x = sigboth_BAT_males$lfc.p, y = sigboth_BAT_males$lfc.e, method = "spearman")
P.cor_BAT_males <- cor.test(x = sigboth_BAT_males$lfc.p, y = sigboth_BAT_males$lfc.e, method = "spearman")$p.value

cor.perm.BAT <- function (data = data, nperm = 10000)
{
  random_subsample <- dplyr::slice_sample(data, n=498) # randomly subsample number estimated from sigboth
  r.per <- sapply (1:nperm, FUN = function (i) cor (x = random_subsample$lfc.p, y = sample (random_subsample$lfc.e)))
  r.per <- c(r.per, r.obs_BAT_males)
  hist (r.per, breaks = 100, xlab = "Correlation coefficient", main=NULL,
        col = "white")
  abline (v = r.obs_BAT_males, col = "#FCECAA", lty = "solid")
  P.per <- sum (abs (r.per) >= abs (r.obs_BAT_males))/(nperm + 1) 
  return (list (r.obs_BAT_males = r.obs_BAT_males, P.cor_BAT_males = P.cor_BAT_males, P.per = P.per, r.per = r.per))
}

BAT_perm_males <- cor.perm.BAT(data = merge_EvP_BAT_males)
BAT_perm_pvalue <- BAT_perm_males$P.per # P < 0.05

# make BAT_perm a df so you can plot hist in ggplot
df_BAT_perm_males <- as_tibble(BAT_perm_males)
# make histogram
BAT_perm_males.inset <-
  ggplot(df_BAT_perm_males, aes(x=r.per)) +
  geom_histogram(fill = "white", color = "gray", bins = 100, binwidth = NULL) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 450), breaks = seq(0,400,125)) +
  geom_vline(xintercept = r.obs_BAT_males, color = "#E09832", linetype = "solid", size = 0.5) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 5, family = "sans"),
        plot.background = element_blank(),
        panel.grid = element_blank()) +
  labs(x = "Correlation Coefficient",
       y = "Frequency")
ggsave("results/figures/males_adap_plast_BAT_inset.pdf", plot = BAT_perm_males.inset,
       height = 1.5, width = 2.25, units = "cm")
# 
# # join BAT adaptive/nonadaptive plot with permutation model
# BAT_nonadap <-
#   ggdraw() +
#   draw_plot(BAT.plot) +
#   draw_plot(BAT_perm.inset, x= 0.12, y = 0.78, width = 0.28, height = 0.2)








#### Bar plots of gene categories ####
# 
# level_order <- c(1,0)
# 
# # BAT
# format_BAT <- sigboth_BAT %>%
#   select(gene, direction) %>%
#   mutate(direction = as_factor(direction))
# 
# BAT_bar <-
#   ggplot(data = format_BAT, aes(x = factor(direction, level = level_order), fill = direction)) +
#   geom_bar(color = "black", show.legend = FALSE, size = 1) +
#   scale_fill_manual(values = c("white", "darkgray")) +
#   scale_x_discrete(labels = c("Adaptive<br>plasticity", "Non-adaptive<br>plasticity")) +
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 300)) +
#   #theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_markdown(),
#         plot.title = element_markdown(size = 8, margin = margin(b = 0, l = 200))) +
#   labs(x = NULL,
#        y = "Number of genes",
#        title = "*P* < 0.001")
# 
# 
# # liver
# format_liver <- sigboth_liver %>%
#   select(gene, direction) %>%
#   mutate(direction = as_factor(direction))
# 
# liver_bar <-
#   ggplot(data = format_liver, aes(x = factor(direction, level = level_order), fill = direction)) +
#   geom_bar(color = "black", show.legend = FALSE, size = 1) +
#   scale_fill_manual(values = c("white", "darkgray")) +
#   scale_x_discrete(labels = c("Adaptive<br>plasticity", "Non-adaptive<br>plasticity")) +
#   scale_y_continuous(expand = c(0, 0), limits = c(0, 300)) +
#   #theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_markdown(),
#         plot.title = element_markdown(size = 8, margin = margin(b = 0, l = 200))) +
#   labs(x = NULL,
#        y = "Number of genes",
#        title = "*P* < 0.001")
# 
# 
# ### Binomial Test
# 
# no.adap_BAT <- sigboth_BAT %>%
#   select(gene, direction) %>%
#   mutate(adap = ifelse(direction > 0, TRUE, FALSE))
# 
# no.adap_liver <- sigboth_liver %>%
#   select(gene, direction) %>%
#   mutate(adap = ifelse(direction > 0, TRUE, FALSE))
# 
# # BAT
# dim(no.adap_BAT)[1] # get number of genes in this category
# sum(no.adap_BAT$adap) # the number of genes where expression is higher in Brazil
# btest = binom.test(x=sum(no.adap_BAT$adap), n = dim(no.adap_BAT)[1])
# # P < 0.001
# 
# # liver
# dim(no.adap_liver)[1] # get number of genes in this category
# sum(no.adap_liver$adap) # the number of genes where expression is higher in Brazil
# btest = binom.test(x=sum(no.adap_liver$adap), n = dim(no.adap_liver)[1])
# # P < 0.001
# 
# 

