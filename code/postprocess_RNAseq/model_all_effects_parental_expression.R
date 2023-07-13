#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script runs a fully parameterized model on parental gene expression.
# Results are included in the SI Appendix (Table S5) of BallingerMack_2023.

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clear R's environment
library(DESeq2)
library(tidyverse)
library(here)

##############################################################
# Import data
##############################################################

Counts_Metadata <- read_table(here("./data/raw/ReadCounts/all_parents_counts.txt"), col_names = TRUE, col_types = NULL)

SampleInfo <- read_delim(here("./data/raw/ReadCounts/all_parents_sample_info_fullmodel.txt"))
colnames(SampleInfo)[1] <- "sampleID"


##############################################################
# Convert tibbles into matrices
##############################################################

tibb_to_mtrx <- function(tibb)
{
  
  genes <- pull(tibb, gene) # extract the column of all the genes
  matrix <- as.matrix(tibb[,-1]) # create matrix of tibble by omitting first column (which has the genes)
  rownames(matrix) <- genes # put genes back into the now matrix (but the genes column does not have a header)
  
  return(matrix)
  
}

all_cts_mtrx <- tibb_to_mtrx(tibb = Counts_Metadata) # SampleInfo


##############################################################
# Construct DESeq datasets
##############################################################

dds_full <- DESeqDataSetFromMatrix(countData = all_cts_mtrx,
                                   colData = SampleInfo,
                                   design = ~ condition + sex + population + temperature + 
                                            condition:sex + condition:population + condition:temperature + population:sex + sex:temperature + population:temperature +
                                            condition:population:sex + condition:sex:temperature + condition:population:temperature + sex:population:temperature +
                                            condition:sex:population:temperature)


# Get the model matrix
mod_mat_full <- model.matrix(design(dds_full), colData(dds_full))

# Define coefficient vectors for each condition
# T
liver <- colMeans(mod_mat_full[dds_full$condition == "LIVER", ])
bat <- colMeans(mod_mat_full[dds_full$condition == "BAT", ])
# S
males <- colMeans(mod_mat_full[dds_full$sex == "M", ])
females <- colMeans(mod_mat_full[dds_full$sex == "F", ])
# G
NewYork <- colMeans(mod_mat_full[dds_full$population == "NewYork", ])
Brazil <- colMeans(mod_mat_full[dds_full$population == "Brazil", ])
# E
warm <- colMeans(mod_mat_full[dds_full$temperature == "Warm", ])
cold <- colMeans(mod_mat_full[dds_full$temperature == "Cold", ])
# TxS
liv_M <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$sex == "M", ])
liv_F <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$sex == "F", ])
bat_M <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$sex == "M", ])
bat_F <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$sex == "F", ])
# TxG
liv_NY <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$population == "NewYork", ])
liv_BZ <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$population == "Brazil", ])
bat_NY <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$population == "NewYork", ])
bat_BZ <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$population == "Brazil", ])
# TxE
liv_W <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$temperature == "Warm", ])
liv_C <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$temperature == "Cold", ])
bat_W <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$temperature == "Warm", ])
bat_C <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$temperature == "Cold", ])
# GxS
NY_M <- colMeans(mod_mat_full[dds_full$population == "NewYork" & dds_full$sex == "M", ])
NY_F <- colMeans(mod_mat_full[dds_full$population == "NewYork" & dds_full$sex == "F", ])
BZ_M <- colMeans(mod_mat_full[dds_full$population == "Brazil" & dds_full$sex == "M", ])
BZ_F <- colMeans(mod_mat_full[dds_full$population == "Brazil" & dds_full$sex == "F", ])
# SxE
M_W <- colMeans(mod_mat_full[dds_full$sex == "M" & dds_full$temperature == "Warm", ])
M_C <- colMeans(mod_mat_full[dds_full$sex == "M" & dds_full$temperature == "Cold", ])
F_W <- colMeans(mod_mat_full[dds_full$sex == "F" & dds_full$temperature == "Warm", ])
F_C <- colMeans(mod_mat_full[dds_full$sex == "F" & dds_full$temperature == "Cold", ])
# GxE
NY_W <- colMeans(mod_mat_full[dds_full$population == "NewYork" & dds_full$temperature == "Warm", ])
NY_C <- colMeans(mod_mat_full[dds_full$population == "NewYork" & dds_full$temperature == "Cold", ])
BZ_W <- colMeans(mod_mat_full[dds_full$population == "Brazil" & dds_full$temperature == "Warm", ])
BZ_C <- colMeans(mod_mat_full[dds_full$population == "Brazil" & dds_full$temperature == "Cold", ])
# TxGxS
liv_NY_M <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$population == "NewYork" & dds_full$sex == "M", ])
liv_NY_F <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$population == "NewYork" & dds_full$sex == "F", ])
liv_BZ_M <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$population == "Brazil" & dds_full$sex == "M", ])
liv_BZ_F <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$population == "Brazil" & dds_full$sex == "F", ])
bat_NY_M <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$population == "NewYork" & dds_full$sex == "M", ])
bat_NY_F <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$population == "NewYork" & dds_full$sex == "F", ])
bat_BZ_M <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$population == "Brazil" & dds_full$sex == "M", ])
bat_BZ_F <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$population == "Brazil" & dds_full$sex == "F", ])
# TxSxE
liv_M_W <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$sex == "M" & dds_full$temperature == "Warm", ])
liv_M_C <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$sex == "M" & dds_full$temperature == "Cold", ])
liv_F_W <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$sex == "F" & dds_full$temperature == "Warm", ])
liv_F_C <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$sex == "F" & dds_full$temperature == "Cold", ])
bat_M_W <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$sex == "M" & dds_full$temperature == "Warm", ])
bat_M_C <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$sex == "M" & dds_full$temperature == "Cold", ])
bat_F_W <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$sex == "F" & dds_full$temperature == "Warm", ])
bat_F_C <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$sex == "F" & dds_full$temperature == "Cold", ])
# TxGxE
liv_NY_W <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$population == "NewYork" & dds_full$temperature == "W", ])
liv_NY_C <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$population == "NewYork" & dds_full$temperature == "C", ])
liv_BZ_W <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$population == "Brazil" & dds_full$temperature == "W", ])
liv_BZ_C <- colMeans(mod_mat_full[dds_full$condition == "LIVER" & dds_full$population == "Brazil" & dds_full$temperature == "C", ])
bat_NY_W <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$population == "NewYork" & dds_full$temperature == "W", ])
bat_NY_C <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$population == "NewYork" & dds_full$temperature == "C", ])
bat_BZ_W <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$population == "Brazil" & dds_full$temperature == "W", ])
bat_BZ_C <- colMeans(mod_mat_full[dds_full$condition == "BAT" & dds_full$population == "Brazil" & dds_full$temperature == "C", ])
# SxGxE
M_NY_W <- colMeans(mod_mat_full[dds_full$sex == "M" & dds_full$population == "NewYork" & dds_full$temperature == "W", ])
M_NY_C <- colMeans(mod_mat_full[dds_full$sex == "M" & dds_full$population == "NewYork" & dds_full$temperature == "C", ])
M_BZ_W <- colMeans(mod_mat_full[dds_full$sex == "M" & dds_full$population == "Brazil" & dds_full$temperature == "W", ])
M_BZ_C <- colMeans(mod_mat_full[dds_full$sex == "M" & dds_full$population == "Brazil" & dds_full$temperature == "C", ])
F_NY_W <- colMeans(mod_mat_full[dds_full$sex == "F" & dds_full$population == "NewYork" & dds_full$temperature == "W", ])
F_NY_C <- colMeans(mod_mat_full[dds_full$sex == "F" & dds_full$population == "NewYork" & dds_full$temperature == "C", ])
F_BZ_W <- colMeans(mod_mat_full[dds_full$sex == "F" & dds_full$population == "Brazil" & dds_full$temperature == "W", ])
F_BZ_C <- colMeans(mod_mat_full[dds_full$sex == "F" & dds_full$population == "Brazil" & dds_full$temperature == "C", ])


##############################################################
# Run DESeq on full & reduced model model
##############################################################

dds_full_deseq <- DESeq(dds_full, test = "Wald")

ddsClean <- dds_full_deseq[which(mcols(dds_full_deseq)$betaConv),] # omit the 24 rows that didn't converge as these are typically genes with very small counts and little power


# T
res_full_T <- results(ddsClean, contrast = liver - bat, pAdjustMethod = "BH")
sig_genes_full_T <- res_full_T %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)

# S
res_full_S <- results(ddsClean, contrast = males - females, pAdjustMethod = "BH")
sig_genes_full_S <- res_full_S %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)

# G
res_full_G <- results(ddsClean, contrast = NewYork - Brazil, pAdjustMethod = "BH")
sig_genes_full_G <- res_full_G %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)

# E
res_full_E <- results(ddsClean, contrast = warm - cold, pAdjustMethod = "BH")
sig_genes_full_E <- res_full_E %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)

# TxS
res_full_TxS <- results(ddsClean, contrast = (liv_M - liv_F) - (bat_M - bat_F), pAdjustMethod = "BH")
sig_genes_full_TxS <- res_full_TxS %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)

# TxG
res_full_TxG <- results(ddsClean, contrast = (liv_NY - liv_BZ) - (bat_NY - bat_BZ), pAdjustMethod = "BH")
sig_genes_full_TxG <- res_full_TxG %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)

# TxE
res_full_TxE <- results(ddsClean, contrast = (liv_W - liv_C) - (bat_W - bat_C), pAdjustMethod = "BH")
sig_genes_full_TxE <- res_full_TxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)

# GxS
res_full_GxS <- results(ddsClean, contrast = (NY_M - NY_F) - (BZ_M - BZ_F), pAdjustMethod = "BH")
sig_genes_full_GxS <- res_full_GxS %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)

# SxE
res_full_SxE <- results(ddsClean, contrast = (M_W - M_C) - (F_W - F_C), pAdjustMethod = "BH")
sig_genes_full_SxE <- res_full_SxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)

# GxE
res_full_GxE <- results(ddsClean, contrast = (NY_W - NY_C) - (BZ_W - BZ_C), pAdjustMethod = "BH")
sig_genes_full_GxE <- res_full_GxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)

# TxGxS
res_full_TxGxS <- results(ddsClean, contrast = ((liv_NY_M - liv_NY_F) - (liv_BZ_M - liv_BZ_F)) - ((bat_NY_M - bat_NY_F) - (bat_BZ_M - bat_BZ_F)), pAdjustMethod = "BH")
sig_genes_full_TxGxS <- res_full_TxGxS %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)

# TxSxE
res_full_TxSxE <- results(ddsClean, contrast = ((liv_M_W - liv_M_C) - (liv_F_W - liv_F_C)) - ((bat_M_W - bat_M_C) - (bat_F_W - bat_F_C)), pAdjustMethod = "BH")
sig_genes_full_TxSxE <- res_full_TxSxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)

# TxGxE
res_full_TxGxE <- results(ddsClean, name = "conditionLIVER.populationNewYork.temperatureWarm", pAdjustMethod = "BH")
sig_genes_full_TxGxE <- res_full_TxGxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)

# SxGxE
res_full_SxGxE <- results(ddsClean, name = "sexM.populationNewYork.temperatureWarm", pAdjustMethod = "BH")
sig_genes_full_SxGxE <- res_full_SxGxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)

# TxSxGxE
res_full_TxSxGxE <- results(ddsClean, name = "conditionLIVER.sexM.populationNewYork.temperatureWarm", pAdjustMethod = "BH")
sig_genes_full_TxSxGxE <- res_full_TxSxGxE %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)


# Run DESeq using reduced model to test if 4-way interaction has any significant explanatory power
# Likelihood ratio test for genotype by treatment by tissue by sex interaction:
dds_reduced <- DESeq(dds_full, test = "LRT",
                     reduced = ~ condition + sex + population + temperature +
                       condition:sex + condition:population + condition:temperature + population:sex + sex:temperature + population:temperature +
                       condition:population:sex + condition:sex:temperature + condition:population:temperature + sex:population:temperature)
# 24 rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT

# 4-way interaction
res_reduced <- results(dds_reduced, pAdjustMethod = "BH")
sig_genes_4way_interaction <- res_reduced %>%
  data.frame() %>%
  rownames_to_column(var = "gene") %>%
  as_tibble() %>%
  filter(baseMean >= 10) %>%
  filter(padj < 0.05)
