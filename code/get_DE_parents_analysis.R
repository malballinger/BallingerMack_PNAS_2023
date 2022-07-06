#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script quantifies differential expression in the parental samples.
# This script generates contrasts for downstream analyses in BallingerMack_PNAS_2022.


##############################################################
# Required packages
##############################################################

rm(list = ls()) # clear R's environment
library(tidyverse)
library(here)
library(DESeq2)

set.seed(19910118)

source("./code/clean_parentalReadCounts.R") # where matrices are generated

##############################################################
# Construct DESeq datasets
##############################################################

### males -- LIVER
dds_males_liver <- DESeqDataSetFromMatrix(countData = males_liver_cts_mtrx,
                                       colData = sampleinfo_males_liver,
                                       design = ~ population + temperature + population:temperature)
# get the model matrix
mod_mat_males_liver <- model.matrix(design(dds_males_liver), colData(dds_males_liver))

# Define coefficient vectors for each condition
NewYork_warm_males_liver <- colMeans(mod_mat_males_liver[dds_males_liver$population == "NewYork" & dds_males_liver$temperature == "Warm", ])
NewYork_cold_males_liver <- colMeans(mod_mat_males_liver[dds_males_liver$population == "NewYork" & dds_males_liver$temperature == "Cold", ])
Brazil_warm_males_liver <- colMeans(mod_mat_males_liver[dds_males_liver$population == "Brazil" & dds_males_liver$temperature == "Warm", ])
Brazil_cold_males_liver <- colMeans(mod_mat_males_liver[dds_males_liver$population == "Brazil" & dds_males_liver$temperature == "Cold", ])
Warm_males_liver <- colMeans(mod_mat_males_liver[dds_males_liver$temperature == "Warm", ])
Cold_males_liver <- colMeans(mod_mat_males_liver[dds_males_liver$temperature == "Cold", ])
NewYork_males_liver <- colMeans(mod_mat_males_liver[dds_males_liver$population == "NewYork", ])
Brazil_males_liver <- colMeans(mod_mat_males_liver[dds_males_liver$population == "Brazil", ])


### males -- BAT
dds_males_BAT <- DESeqDataSetFromMatrix(countData = males_BAT_cts_mtrx,
                                          colData = sampleinfo_males_BAT,
                                          design = ~ population + temperature + population:temperature)
# get the model matrix
mod_mat_males_BAT <- model.matrix(design(dds_males_BAT), colData(dds_males_BAT))

# Define coefficient vectors for each condition
NewYork_warm_males_BAT <- colMeans(mod_mat_males_BAT[dds_males_BAT$population == "NewYork" & dds_males_BAT$temperature == "Warm", ])
NewYork_cold_males_BAT <- colMeans(mod_mat_males_BAT[dds_males_BAT$population == "NewYork" & dds_males_BAT$temperature == "Cold", ])
Brazil_warm_males_BAT <- colMeans(mod_mat_males_BAT[dds_males_BAT$population == "Brazil" & dds_males_BAT$temperature == "Warm", ])
Brazil_cold_males_BAT <- colMeans(mod_mat_males_BAT[dds_males_BAT$population == "Brazil" & dds_males_BAT$temperature == "Cold", ])
Warm_males_BAT <- colMeans(mod_mat_males_BAT[dds_males_BAT$temperature == "Warm", ])
Cold_males_BAT <- colMeans(mod_mat_males_BAT[dds_males_BAT$temperature == "Cold", ])
NewYork_males_BAT <- colMeans(mod_mat_males_BAT[dds_males_BAT$population == "NewYork", ])
Brazil_males_BAT <- colMeans(mod_mat_males_BAT[dds_males_BAT$population == "Brazil", ])


### females -- LIVER
dds_females_liver <- DESeqDataSetFromMatrix(countData = females_liver_cts_mtrx,
                                          colData = sampleinfo_females_liver,
                                          design = ~ population + temperature + population:temperature)
# get the model matrix
mod_mat_females_liver <- model.matrix(design(dds_females_liver), colData(dds_females_liver))

# Define coefficient vectors for each condition
NewYork_warm_females_liver <- colMeans(mod_mat_females_liver[dds_females_liver$population == "NewYork" & dds_females_liver$temperature == "Warm", ])
NewYork_cold_females_liver <- colMeans(mod_mat_females_liver[dds_females_liver$population == "NewYork" & dds_females_liver$temperature == "Cold", ])
Brazil_warm_females_liver <- colMeans(mod_mat_females_liver[dds_females_liver$population == "Brazil" & dds_females_liver$temperature == "Warm", ])
Brazil_cold_females_liver <- colMeans(mod_mat_females_liver[dds_females_liver$population == "Brazil" & dds_females_liver$temperature == "Cold", ])
Warm_females_liver <- colMeans(mod_mat_females_liver[dds_females_liver$temperature == "Warm", ])
Cold_females_liver <- colMeans(mod_mat_females_liver[dds_females_liver$temperature == "Cold", ])
NewYork_females_liver <- colMeans(mod_mat_females_liver[dds_females_liver$population == "NewYork", ])
Brazil_females_liver <- colMeans(mod_mat_females_liver[dds_females_liver$population == "Brazil", ])


### females -- BAT
dds_females_BAT <- DESeqDataSetFromMatrix(countData = females_BAT_cts_mtrx,
                                        colData = sampleinfo_females_BAT,
                                        design = ~ population + temperature + population:temperature)
# get the model matrix
mod_mat_females_BAT <- model.matrix(design(dds_females_BAT), colData(dds_females_BAT))

# Define coefficient vectors for each condition
NewYork_warm_females_BAT <- colMeans(mod_mat_females_BAT[dds_females_BAT$population == "NewYork" & dds_females_BAT$temperature == "Warm", ])
NewYork_cold_females_BAT <- colMeans(mod_mat_females_BAT[dds_females_BAT$population == "NewYork" & dds_females_BAT$temperature == "Cold", ])
Brazil_warm_females_BAT <- colMeans(mod_mat_females_BAT[dds_females_BAT$population == "Brazil" & dds_females_BAT$temperature == "Warm", ])
Brazil_cold_females_BAT <- colMeans(mod_mat_females_BAT[dds_females_BAT$population == "Brazil" & dds_females_BAT$temperature == "Cold", ])
Warm_females_BAT <- colMeans(mod_mat_females_BAT[dds_females_BAT$temperature == "Warm", ])
Cold_females_BAT <- colMeans(mod_mat_females_BAT[dds_females_BAT$temperature == "Cold", ])
NewYork_females_BAT <- colMeans(mod_mat_females_BAT[dds_females_BAT$population == "NewYork", ])
Brazil_females_BAT <- colMeans(mod_mat_females_BAT[dds_females_BAT$population == "Brazil", ])

##############################################################
# Run DESeq
##############################################################

### males -- LIVER
deseq_males_liver <- DESeq(dds_males_liver, test = "Wald")
#saveRDS(deseq_males_liver, file = "./data/processed/DESeq_males_Liver.rda", compress = TRUE)
#res_males_liver <- results(deseq_males_liver)

# New York vs Brazil (in the warm (evolved diffs)):
res_warm_males_liver_NYvsBZ <- results(deseq_males_liver, contrast = NewYork_warm_males_liver - Brazil_warm_males_liver,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in BZ warm
#write.csv(res_warm_males_liver_NYvsBZ, file = "data/processed/DE_warm_males_liver_NYvsBZ.csv")

# Brazil vs NewYork (in the warm (evolved diffs)):
res_warm_males_liver_BZvsNY <- results(deseq_males_liver, contrast = Brazil_warm_males_liver - NewYork_warm_males_liver,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in NY warm
#write.csv(res_warm_males_liver_BZvsNY, file = "data/processed/DE_warm_males_liver_BZvsNY.csv")

# New York vs Brazil (in the cold):
res_cold_males_liver_NYvsBZ <- results(deseq_males_liver, contrast = NewYork_cold_males_liver - Brazil_cold_males_liver,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in BZ cold
# Brazil vs NewYork (in the cold):
res_cold_males_liver_BZvsNY <- results(deseq_males_liver, contrast = Brazil_cold_males_liver - NewYork_cold_males_liver,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in NY cold
# Warm vs Cold (for New York):
res_NY_males_liver_WvC <- results(deseq_males_liver, contrast = NewYork_warm_males_liver - NewYork_cold_males_liver,
                                  pAdjustMethod = "BH") # (-)LFC == higher expression in NY cold
# Warm vs Cold (for Brazil):
res_BZ_males_liver_WvC <- results(deseq_males_liver, contrast = Brazil_warm_males_liver - Brazil_cold_males_liver,
                                  pAdjustMethod = "BH") # (-)LFC == higher expression in BZ cold
# Cold vs Warm (for New York):
res_NY_males_liver_CvW <- results(deseq_males_liver, contrast = NewYork_cold_males_liver - NewYork_warm_males_liver,
                                  pAdjustMethod = "BH") # (-)LFC == higher expression in NY warm
# Cold vs Warm (for Brazil):
res_BZ_males_liver_CvW <- results(deseq_males_liver, contrast = Brazil_cold_males_liver - Brazil_warm_males_liver,
                                  pAdjustMethod = "BH") # (-)LFC == higher expression in BZ warm
# Effect of genotype ("G")
res_males_liver_G <- results(deseq_males_liver, contrast = NewYork_males_liver - Brazil_males_liver,
                             pAdjustMethod = "BH") # (-)LFC == higher expression in BZ males
# Effect of environment ("E")
res_males_liver_E <- results(deseq_males_liver, contrast = Warm_males_liver - Cold_males_liver,
                             pAdjustMethod = "BH") # (-)LFC == higher expression in cold males
# DESEQ's definition of genotype-by-environment ("GxE")
res_males_liver_GxE <- results(deseq_males_liver, contrast = (NewYork_warm_males_liver - NewYork_cold_males_liver) -
                                                              (Brazil_warm_males_liver - Brazil_cold_males_liver),
                               pAdjustMethod = "BH")

### males -- BAT
deseq_males_BAT <- DESeq(dds_males_BAT, test = "Wald")
#res_males_BAT <- results(deseq_males_BAT)

# New York vs Brazil (in the warm (evolved diffs)):
res_warm_males_BAT_NYvsBZ <- results(deseq_males_BAT, contrast = NewYork_warm_males_BAT - Brazil_warm_males_BAT,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in BZ warm
# Brazil vs NewYork (in the warm (evolved diffs)):
res_warm_males_BAT_BZvsNY <- results(deseq_males_BAT, contrast = Brazil_warm_males_BAT - NewYork_warm_males_BAT,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in NY warm
# New York vs Brazil (in the cold):
res_cold_males_BAT_NYvsBZ <- results(deseq_males_BAT, contrast = NewYork_cold_males_BAT - Brazil_cold_males_BAT,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in BZ cold
# Brazil vs NewYork (in the cold):
res_cold_males_BAT_BZvsNY <- results(deseq_males_BAT, contrast = Brazil_cold_males_BAT - NewYork_cold_males_BAT,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in NY cold
# Warm vs Cold (for New York):
res_NY_males_BAT_WvC <- results(deseq_males_BAT, contrast = NewYork_warm_males_BAT - NewYork_cold_males_BAT,
                                  pAdjustMethod = "BH") # (-)LFC == higher expression in NY cold
# Warm vs Cold (for Brazil):
res_BZ_males_BAT_WvC <- results(deseq_males_BAT, contrast = Brazil_warm_males_BAT - Brazil_cold_males_BAT,
                                  pAdjustMethod = "BH") # (-)LFC == higher expression in BZ cold
# Cold vs Warm (for New York):
res_NY_males_BAT_CvW <- results(deseq_males_BAT, contrast = NewYork_cold_males_BAT - NewYork_warm_males_BAT,
                                  pAdjustMethod = "BH") # (-)LFC == higher expression in NY warm
# Cold vs Warm (for Brazil):
res_BZ_males_BAT_CvW <- results(deseq_males_BAT, contrast = Brazil_cold_males_BAT - Brazil_warm_males_BAT,
                                  pAdjustMethod = "BH") # (-)LFC == higher expression in BZ warm
# Effect of genotype ("G")
res_males_BAT_G <- results(deseq_males_BAT, contrast = NewYork_males_BAT - Brazil_males_BAT,
                             pAdjustMethod = "BH") # (-)LFC == higher expression in BZ males
# Effect of environment ("E")
res_males_BAT_E <- results(deseq_males_BAT, contrast = Warm_males_BAT - Cold_males_BAT,
                             pAdjustMethod = "BH") # (-)LFC == higher expression in cold males
# DESEQ's definition of genotype-by-environment ("GxE")
res_males_BAT_GxE <- results(deseq_males_BAT, contrast = (NewYork_warm_males_BAT - NewYork_cold_males_BAT) -
                                                          (Brazil_warm_males_BAT - Brazil_cold_males_BAT),
                             pAdjustMethod = "BH")

### females -- LIVER
deseq_females_liver <- DESeq(dds_females_liver, test = "Wald")
#res_females_liver <- results(deseq_females_liver)

# New York vs Brazil (in the warm (evolved diffs)):
res_warm_females_liver_NYvsBZ <- results(deseq_females_liver, contrast = NewYork_warm_females_liver - Brazil_warm_females_liver,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in BZ warm
# Brazil vs NewYork (in the warm (evolved diffs)):
res_warm_females_liver_BZvsNY <- results(deseq_females_liver, contrast = Brazil_warm_females_liver - NewYork_warm_females_liver,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in NY warm
# New York vs Brazil (in the cold):
res_cold_females_liver_NYvsBZ <- results(deseq_females_liver, contrast = NewYork_cold_females_liver - Brazil_cold_females_liver,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in BZ cold
# Brazil vs NewYork (in the cold):
res_cold_females_liver_BZvsNY <- results(deseq_females_liver, contrast = Brazil_cold_females_liver - NewYork_cold_females_liver,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in NY cold
# Warm vs Cold (for New York):
res_NY_females_liver_WvC <- results(deseq_females_liver, contrast = NewYork_warm_females_liver - NewYork_cold_females_liver,
                                  pAdjustMethod = "BH") # (-)LFC == higher expression in NY cold
# Warm vs Cold (for Brazil):
res_BZ_females_liver_WvC <- results(deseq_females_liver, contrast = Brazil_warm_females_liver - Brazil_cold_females_liver,
                                  pAdjustMethod = "BH") # (-)LFC == higher expression in BZ cold
# Cold vs Warm (for New York):
res_NY_females_liver_CvW <- results(deseq_females_liver, contrast = NewYork_cold_females_liver - NewYork_warm_females_liver,
                                  pAdjustMethod = "BH") # (-)LFC == higher expression in NY warm
# Cold vs Warm (for Brazil):
res_BZ_females_liver_CvW <- results(deseq_females_liver, contrast = Brazil_cold_females_liver - Brazil_warm_females_liver,
                                  pAdjustMethod = "BH") # (-)LFC == higher expression in BZ warm
# Effect of genotype ("G")
res_females_liver_G <- results(deseq_females_liver, contrast = NewYork_females_liver - Brazil_females_liver,
                             pAdjustMethod = "BH") # (-)LFC == higher expression in BZ females
# Effect of environment ("E")
res_females_liver_E <- results(deseq_females_liver, contrast = Warm_females_liver - Cold_females_liver,
                             pAdjustMethod = "BH") # (-)LFC == higher expression in cold females

### females -- BAT
deseq_females_BAT <- DESeq(dds_females_BAT, test = "Wald")
#res_females_BAT <- results(deseq_females_BAT)

# New York vs Brazil (in the warm (evolved diffs)):
res_warm_females_BAT_NYvsBZ <- results(deseq_females_BAT, contrast = NewYork_warm_females_BAT - Brazil_warm_females_BAT,
                                     pAdjustMethod = "BH") # (-)LFC == higher expression in BZ warm
# Brazil vs NewYork (in the warm (evolved diffs)):
res_warm_females_BAT_BZvsNY <- results(deseq_females_BAT, contrast = Brazil_warm_females_BAT - NewYork_warm_females_BAT,
                                     pAdjustMethod = "BH") # (-)LFC == higher expression in NY warm
# New York vs Brazil (in the cold):
res_cold_females_BAT_NYvsBZ <- results(deseq_females_BAT, contrast = NewYork_cold_females_BAT - Brazil_cold_females_BAT,
                                     pAdjustMethod = "BH") # (-)LFC == higher expression in BZ cold
# Brazil vs NewYork (in the cold):
res_cold_females_BAT_BZvsNY <- results(deseq_females_BAT, contrast = Brazil_cold_females_BAT - NewYork_cold_females_BAT,
                                     pAdjustMethod = "BH") # (-)LFC == higher expression in NY cold
# Warm vs Cold (for New York):
res_NY_females_BAT_WvC <- results(deseq_females_BAT, contrast = NewYork_warm_females_BAT - NewYork_cold_females_BAT,
                                pAdjustMethod = "BH") # (-)LFC == higher expression in NY cold
# Warm vs Cold (for Brazil):
res_BZ_females_BAT_WvC <- results(deseq_females_BAT, contrast = Brazil_warm_females_BAT - Brazil_cold_females_BAT,
                                pAdjustMethod = "BH") # (-)LFC == higher expression in BZ cold
# Cold vs Warm (for New York):
res_NY_females_BAT_CvW <- results(deseq_females_BAT, contrast = NewYork_cold_females_BAT - NewYork_warm_females_BAT,
                                pAdjustMethod = "BH") # (-)LFC == higher expression in NY warm
# Cold vs Warm (for Brazil):
res_BZ_females_BAT_CvW <- results(deseq_females_BAT, contrast = Brazil_cold_females_BAT - Brazil_warm_females_BAT,
                                pAdjustMethod = "BH") # (-)LFC == higher expression in BZ warm
# Effect of genotype ("G")
res_females_BAT_G <- results(deseq_females_BAT, contrast = NewYork_females_BAT - Brazil_females_BAT,
                           pAdjustMethod = "BH") # (-)LFC == higher expression in BZ females
# Effect of environment ("E")
res_females_BAT_E <- results(deseq_females_BAT, contrast = Warm_females_BAT - Cold_females_BAT,
                           pAdjustMethod = "BH") # (-)LFC == higher expression in cold females






















# res_males_liver_E %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% View()
ex_gene <- plotCounts(deseq_males_liver, gene = "ENSMUSG00000027474", intgroup = c("temperature", "population"), returnData = TRUE)
level_order <- c("Warm", "Cold")
ggplot(data=ex_gene, aes(x = factor(temperature, level = level_order), y=count, color = population)) + #, group = population)) +
  geom_boxplot(position = "identity", fill = NA) + labs(x = "Temp", y = "Normalized Counts")
level_order2 <- c("Brazil", "NewYork")
ggplot(data=ex_gene, aes(x = factor(population, level = level_order2), y=count, color = temperature)) + #, group = population)) +
  geom_boxplot(position = "identity", fill = NA) + labs(x = "Population", y = "Normalized Counts")

