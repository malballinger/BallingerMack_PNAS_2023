#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script quantifies differential expression in the parental samples and
# generates contrasts for downstream analyses in BallingerMack_2022.

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clear R's environment
library(tidyverse)
library(DESeq2)

set.seed(19910118)

source("./code/postprocess_RNAseq/clean_parentalReadCounts.R") # where matrices are generated

##############################################################
# Construct DESeq datasets
##############################################################

#### G x E contrasts

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


#### G x S contrasts

### LIVER -- warm
dds_sex_liver_warm <- DESeqDataSetFromMatrix(countData = warm_liver_cts_mtrx,
                                             colData = sampleinfo_liver_warm,
                                             design = ~ population + sex + population:sex)
# get the model matrix
mod_mat_sex_liver_warm <- model.matrix(design(dds_sex_liver_warm), colData(dds_sex_liver_warm))

# Define coefficient vectors for each category
NewYork_males_liver_warm <- colMeans(mod_mat_sex_liver_warm[dds_sex_liver_warm$population == "NewYork" & dds_sex_liver_warm$sex == "M", ])
NewYork_females_liver_warm <- colMeans(mod_mat_sex_liver_warm[dds_sex_liver_warm$population == "NewYork" & dds_sex_liver_warm$sex == "F", ])
Brazil_males_liver_warm <- colMeans(mod_mat_sex_liver_warm[dds_sex_liver_warm$population == "Brazil" & dds_sex_liver_warm$sex == "M", ])
Brazil_females_liver_warm <- colMeans(mod_mat_sex_liver_warm[dds_sex_liver_warm$population == "Brazil" & dds_sex_liver_warm$sex == "F", ])
Males_liver_warm <- colMeans(mod_mat_sex_liver_warm[dds_sex_liver_warm$sex == "M", ])
Females_liver_warm <- colMeans(mod_mat_sex_liver_warm[dds_sex_liver_warm$sex == "F", ])
NewYork_liver_warm <- colMeans(mod_mat_sex_liver_warm[dds_sex_liver_warm$population == "NewYork", ])
Brazil_liver_warm <- colMeans(mod_mat_sex_liver_warm[dds_sex_liver_warm$population == "Brazil", ])


### LIVER -- cold
dds_sex_liver_cold <- DESeqDataSetFromMatrix(countData = cold_liver_cts_mtrx,
                                             colData = sampleinfo_liver_cold,
                                             design = ~ population + sex + population:sex)
# get the model matrix
mod_mat_sex_liver_cold <- model.matrix(design(dds_sex_liver_cold), colData(dds_sex_liver_cold))

# Define coefficient vectors for each category
NewYork_males_liver_cold <- colMeans(mod_mat_sex_liver_cold[dds_sex_liver_cold$population == "NewYork" & dds_sex_liver_cold$sex == "M", ])
NewYork_females_liver_cold <- colMeans(mod_mat_sex_liver_cold[dds_sex_liver_cold$population == "NewYork" & dds_sex_liver_cold$sex == "F", ])
Brazil_males_liver_cold <- colMeans(mod_mat_sex_liver_cold[dds_sex_liver_cold$population == "Brazil" & dds_sex_liver_cold$sex == "M", ])
Brazil_females_liver_cold <- colMeans(mod_mat_sex_liver_cold[dds_sex_liver_cold$population == "Brazil" & dds_sex_liver_cold$sex == "F", ])
Males_liver_cold <- colMeans(mod_mat_sex_liver_cold[dds_sex_liver_cold$sex == "M", ])
Females_liver_cold <- colMeans(mod_mat_sex_liver_cold[dds_sex_liver_cold$sex == "F", ])
NewYork_liver_cold <- colMeans(mod_mat_sex_liver_cold[dds_sex_liver_cold$population == "NewYork", ])
Brazil_liver_cold <- colMeans(mod_mat_sex_liver_cold[dds_sex_liver_cold$population == "Brazil", ])


### BAT -- warm
dds_sex_BAT_warm <- DESeqDataSetFromMatrix(countData = warm_BAT_cts_mtrx,
                                           colData = sampleinfo_BAT_warm,
                                           design = ~ population + sex + population:sex)
# get the model matrix
mod_mat_sex_BAT_warm <- model.matrix(design(dds_sex_BAT_warm), colData(dds_sex_BAT_warm))

# Define coefficient vectors for each category
NewYork_males_BAT_warm <- colMeans(mod_mat_sex_BAT_warm[dds_sex_BAT_warm$population == "NewYork" & dds_sex_BAT_warm$sex == "M", ])
NewYork_females_BAT_warm <- colMeans(mod_mat_sex_BAT_warm[dds_sex_BAT_warm$population == "NewYork" & dds_sex_BAT_warm$sex == "F", ])
Brazil_males_BAT_warm <- colMeans(mod_mat_sex_BAT_warm[dds_sex_BAT_warm$population == "Brazil" & dds_sex_BAT_warm$sex == "M", ])
Brazil_females_BAT_warm <- colMeans(mod_mat_sex_BAT_warm[dds_sex_BAT_warm$population == "Brazil" & dds_sex_BAT_warm$sex == "F", ])
Males_BAT_warm <- colMeans(mod_mat_sex_BAT_warm[dds_sex_BAT_warm$sex == "M", ])
Females_BAT_warm <- colMeans(mod_mat_sex_BAT_warm[dds_sex_BAT_warm$sex == "F", ])
NewYork_BAT_warm <- colMeans(mod_mat_sex_BAT_warm[dds_sex_BAT_warm$population == "NewYork", ])
Brazil_BAT_warm <- colMeans(mod_mat_sex_BAT_warm[dds_sex_BAT_warm$population == "Brazil", ])


### BAT -- cold
dds_sex_BAT_cold <- DESeqDataSetFromMatrix(countData = cold_BAT_cts_mtrx,
                                           colData = sampleinfo_BAT_cold,
                                           design = ~ population + sex + population:sex)
# get the model matrix
mod_mat_sex_BAT_cold <- model.matrix(design(dds_sex_BAT_cold), colData(dds_sex_BAT_cold))

# Define coefficient vectors for each category
NewYork_males_BAT_cold <- colMeans(mod_mat_sex_BAT_cold[dds_sex_BAT_cold$population == "NewYork" & dds_sex_BAT_cold$sex == "M", ])
NewYork_females_BAT_cold <- colMeans(mod_mat_sex_BAT_cold[dds_sex_BAT_cold$population == "NewYork" & dds_sex_BAT_cold$sex == "F", ])
Brazil_males_BAT_cold <- colMeans(mod_mat_sex_BAT_cold[dds_sex_BAT_cold$population == "Brazil" & dds_sex_BAT_cold$sex == "M", ])
Brazil_females_BAT_cold <- colMeans(mod_mat_sex_BAT_cold[dds_sex_BAT_cold$population == "Brazil" & dds_sex_BAT_cold$sex == "F", ])
Males_BAT_cold <- colMeans(mod_mat_sex_BAT_cold[dds_sex_BAT_cold$sex == "M", ])
Females_BAT_cold <- colMeans(mod_mat_sex_BAT_cold[dds_sex_BAT_cold$sex == "F", ])
NewYork_BAT_cold <- colMeans(mod_mat_sex_BAT_cold[dds_sex_BAT_cold$population == "NewYork", ])
Brazil_BAT_cold <- colMeans(mod_mat_sex_BAT_cold[dds_sex_BAT_cold$population == "Brazil", ])


##############################################################
# Run DESeq
##############################################################

### males -- LIVER
deseq_males_liver <- DESeq(dds_males_liver, test = "Wald")
#saveRDS(deseq_males_liver, file = "./data/processed/DESeq_males_Liver.rda", compress = TRUE)
#res_males_liver <- results(deseq_males_liver)

# New York vs Brazil (in the warm):
res_warm_males_liver_NYvsBZ <- results(deseq_males_liver, contrast = NewYork_warm_males_liver - Brazil_warm_males_liver,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in BZ warm
#write.csv(res_warm_males_liver_NYvsBZ, file = "data/processed/DE_warm_males_liver_NYvsBZ.csv")

# Brazil vs New York (in the warm:
res_warm_males_liver_BZvsNY <- results(deseq_males_liver, contrast = Brazil_warm_males_liver - NewYork_warm_males_liver,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in NY warm
#write.csv(res_warm_males_liver_BZvsNY, file = "data/processed/DE_warm_males_liver_BZvsNY.csv")

# New York vs Brazil (in the cold):
res_cold_males_liver_NYvsBZ <- results(deseq_males_liver, contrast = NewYork_cold_males_liver - Brazil_cold_males_liver,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in BZ cold
# Brazil vs New York (in the cold):
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

# New York vs Brazil (in the warm):
res_warm_males_BAT_NYvsBZ <- results(deseq_males_BAT, contrast = NewYork_warm_males_BAT - Brazil_warm_males_BAT,
                                     pAdjustMethod = "BH") # (-)LFC == higher expression in BZ warm
# Brazil vs New York (in the warm):
res_warm_males_BAT_BZvsNY <- results(deseq_males_BAT, contrast = Brazil_warm_males_BAT - NewYork_warm_males_BAT,
                                     pAdjustMethod = "BH") # (-)LFC == higher expression in NY warm
# New York vs Brazil (in the cold):
res_cold_males_BAT_NYvsBZ <- results(deseq_males_BAT, contrast = NewYork_cold_males_BAT - Brazil_cold_males_BAT,
                                     pAdjustMethod = "BH") # (-)LFC == higher expression in BZ cold
# Brazil vs New York (in the cold):
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

# New York vs Brazil (in the warm):
res_warm_females_liver_NYvsBZ <- results(deseq_females_liver, contrast = NewYork_warm_females_liver - Brazil_warm_females_liver,
                                         pAdjustMethod = "BH") # (-)LFC == higher expression in BZ warm
# Brazil vs New York (in the warm):
res_warm_females_liver_BZvsNY <- results(deseq_females_liver, contrast = Brazil_warm_females_liver - NewYork_warm_females_liver,
                                         pAdjustMethod = "BH") # (-)LFC == higher expression in NY warm
# New York vs Brazil (in the cold):
res_cold_females_liver_NYvsBZ <- results(deseq_females_liver, contrast = NewYork_cold_females_liver - Brazil_cold_females_liver,
                                         pAdjustMethod = "BH") # (-)LFC == higher expression in BZ cold
# Brazil vs New York (in the cold):
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
# DESEQ's definition of genotype-by-environment ("GxE")
res_females_liver_GxE <- results(deseq_females_liver, contrast = (NewYork_warm_females_liver - NewYork_cold_females_liver) -
                                                                  (Brazil_warm_females_liver - Brazil_cold_females_liver),
                                 pAdjustMethod = "BH")


### females -- BAT
deseq_females_BAT <- DESeq(dds_females_BAT, test = "Wald")
#res_females_BAT <- results(deseq_females_BAT) %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10) # xxxx genes

# New York vs Brazil (in the warm):
res_warm_females_BAT_NYvsBZ <- results(deseq_females_BAT, contrast = NewYork_warm_females_BAT - Brazil_warm_females_BAT,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in BZ warm
# Brazil vs New York (in the warm):
res_warm_females_BAT_BZvsNY <- results(deseq_females_BAT, contrast = Brazil_warm_females_BAT - NewYork_warm_females_BAT,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in NY warm
# New York vs Brazil (in the cold):
res_cold_females_BAT_NYvsBZ <- results(deseq_females_BAT, contrast = NewYork_cold_females_BAT - Brazil_cold_females_BAT,
                                       pAdjustMethod = "BH") # (-)LFC == higher expression in BZ cold
# Brazil vs New York (in the cold):
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
# DESEQ's definition of genotype-by-environment ("GxE")
res_females_BAT_GxE <- results(deseq_females_BAT, contrast = (NewYork_warm_females_BAT - NewYork_cold_females_BAT) -
                                                              (Brazil_warm_females_BAT - Brazil_cold_females_BAT),
                               pAdjustMethod = "BH")


### LIVER -- warm
deseq_sex_liver_warm <- DESeq(dds_sex_liver_warm, test = "Wald")
#no_genes_tested_sex_liver_warm <- results(deseq_sex_liver_warm) %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10) # 14,153 genes

# New York vs Brazil (males):
res_males_liver_warm_NYvBZ <- results(deseq_sex_liver_warm, contrast = NewYork_males_liver_warm - Brazil_males_liver_warm,
                                      pAdjustMethod = "BH") # (-)LFC == higher expression in BZ males
#no_DEsigbase_genes_males_liver_warm_NYvBZ <- res_males_liver_warm_NYvBZ %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10, padj < 0.05) # 4,422 genes

# New York vs Brazil (females):
res_females_liver_warm_NYvBZ <- results(deseq_sex_liver_warm, contrast = NewYork_females_liver_warm - Brazil_females_liver_warm,
                                        pAdjustMethod = "BH") # (-)LFC == higher expression in BZ females
#no_DEsigbase_genes_females_liver_warm_NYvBZ <- res_females_liver_warm_NYvBZ %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10, padj < 0.05) # 3,524 genes

# Males vs Females (New York):
res_NY_liver_warm_MvF <- results(deseq_sex_liver_warm, contrast = NewYork_males_liver_warm - NewYork_females_liver_warm,
                                 pAdjustMethod = "BH") # (-)LFC == higher expression in NY females
# Males vs Females (Brazil):
res_BZ_liver_warm_MvF <- results(deseq_sex_liver_warm, contrast = Brazil_males_liver_warm - Brazil_females_liver_warm,
                                 pAdjustMethod = "BH") # (-)LFC == higher expression in BZ females
# Effect of genotype ("G")
res_sex_liver_warm_G <- results(deseq_sex_liver_warm, contrast = NewYork_liver_warm - Brazil_liver_warm,
                                pAdjustMethod = "BH") # (-)LFC == higher expression in BZ
# Effect of environment ("S")
res_sex_liver_warm_S <- results(deseq_sex_liver_warm, contrast = Males_liver_warm - Females_liver_warm,
                                pAdjustMethod = "BH") # (-)LFC == higher expression in Females
# DESEQ's definition of genotype-by-sex ("GxS")
res_sex_liver_warm_GxS <- results(deseq_sex_liver_warm, contrast = (NewYork_males_liver_warm - NewYork_females_liver_warm) -
                                                                    (Brazil_males_liver_warm - Brazil_females_liver_warm),
                                  pAdjustMethod = "BH")
#no_DEsigbase_genes_sex_liver_warm_GxS <- res_sex_liver_warm_GxS %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10, padj < 0.05) # 22 genes


### LIVER -- cold
deseq_sex_liver_cold <- DESeq(dds_sex_liver_cold, test = "Wald")
#no_genes_tested_sex_liver_cold <- results(deseq_sex_liver_cold) %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10) # 14,124 genes

# New York vs Brazil (males):
res_males_liver_cold_NYvBZ <- results(deseq_sex_liver_cold, contrast = NewYork_males_liver_cold - Brazil_males_liver_cold,
                                      pAdjustMethod = "BH") # (-)LFC == higher expression in BZ males
#no_DEsigbase_genes_males_liver_cold_NYvBZ <- res_males_liver_cold_NYvBZ %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10, padj < 0.05) # 3,305 genes

# New York vs Brazil (females):
res_females_liver_cold_NYvBZ <- results(deseq_sex_liver_cold, contrast = NewYork_females_liver_cold - Brazil_females_liver_cold,
                                        pAdjustMethod = "BH") # (-)LFC == higher expression in BZ females
#no_DEsigbase_genes_females_liver_cold_NYvBZ <- res_females_liver_cold_NYvBZ %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10, padj < 0.05) # 3,796 genes

# Males vs Females (New York):
res_NY_liver_cold_MvF <- results(deseq_sex_liver_cold, contrast = NewYork_males_liver_cold - NewYork_females_liver_cold,
                                 pAdjustMethod = "BH") # (-)LFC == higher expression in NY females
# Males vs Females (Brazil):
res_BZ_liver_cold_MvF <- results(deseq_sex_liver_cold, contrast = Brazil_males_liver_cold - Brazil_females_liver_cold,
                                 pAdjustMethod = "BH") # (-)LFC == higher expression in BZ females
# Effect of genotype ("G")
res_sex_liver_cold_G <- results(deseq_sex_liver_cold, contrast = NewYork_liver_cold - Brazil_liver_cold,
                                pAdjustMethod = "BH") # (-)LFC == higher expression in BZ
# Effect of environment ("S")
res_sex_liver_cold_S <- results(deseq_sex_liver_cold, contrast = Males_liver_cold - Females_liver_cold,
                                pAdjustMethod = "BH") # (-)LFC == higher expression in Females
# DESEQ's definition of genotype-by-sex ("GxS")
res_sex_liver_cold_GxS <- results(deseq_sex_liver_cold, contrast = (NewYork_males_liver_cold - NewYork_females_liver_cold) -
                                                                    (Brazil_males_liver_cold - Brazil_females_liver_cold),
                                  pAdjustMethod = "BH")
#no_DEsigbase_genes_sex_liver_cold_GxS <- res_sex_liver_cold_GxS %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10, padj < 0.05) # 75 genes


### BAT -- warm
deseq_sex_BAT_warm <- DESeq(dds_sex_BAT_warm, test = "Wald")
#no_genes_tested_sex_BAT_warm <- results(deseq_sex_BAT_warm) %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10) # 14,162 genes

# New York vs Brazil (males):
res_males_BAT_warm_NYvBZ <- results(deseq_sex_BAT_warm, contrast = NewYork_males_BAT_warm - Brazil_males_BAT_warm,
                                    pAdjustMethod = "BH") # (-)LFC == higher expression in BZ males
#no_DEsigbase_genes_males_BAT_warm_NYvBZ <- res_males_BAT_warm_NYvBZ %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10, padj < 0.05) # 3,995 genes

# New York vs Brazil (females):
res_females_BAT_warm_NYvBZ <- results(deseq_sex_BAT_warm, contrast = NewYork_females_BAT_warm - Brazil_females_BAT_warm,
                                      pAdjustMethod = "BH") # (-)LFC == higher expression in BZ females
#no_DEsigbase_genes_females_BAT_warm_NYvBZ <- res_females_BAT_warm_NYvBZ %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10, padj < 0.05) # 3,544 genes

# Males vs Females (New York):
res_NY_BAT_warm_MvF <- results(deseq_sex_BAT_warm, contrast = NewYork_males_BAT_warm - NewYork_females_BAT_warm,
                               pAdjustMethod = "BH") # (-)LFC == higher expression in NY females
# Males vs Females (Brazil):
res_BZ_BAT_warm_MvF <- results(deseq_sex_BAT_warm, contrast = Brazil_males_BAT_warm - Brazil_females_BAT_warm,
                               pAdjustMethod = "BH") # (-)LFC == higher expression in BZ females
# Effect of genotype ("G")
res_sex_BAT_warm_G <- results(deseq_sex_BAT_warm, contrast = NewYork_BAT_warm - Brazil_BAT_warm,
                              pAdjustMethod = "BH") # (-)LFC == higher expression in BZ
# Effect of environment ("S")
res_sex_BAT_warm_S <- results(deseq_sex_BAT_warm, contrast = Males_BAT_warm - Females_BAT_warm,
                              pAdjustMethod = "BH") # (-)LFC == higher expression in Females
# DESEQ's definition of genotype-by-sex ("GxS")
res_sex_BAT_warm_GxS <- results(deseq_sex_BAT_warm, contrast = (NewYork_males_BAT_warm - NewYork_females_BAT_warm) -
                                                                (Brazil_males_BAT_warm - Brazil_females_BAT_warm),
                                pAdjustMethod = "BH")
#no_DEsigbase_genes_sex_BAT_warm_GxS <- res_sex_BAT_warm_GxS %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10, padj < 0.05) # 0 genes


### BAT -- cold
deseq_sex_BAT_cold <- DESeq(dds_sex_BAT_cold, test = "Wald")
#no_genes_tested_sex_BAT_cold <- results(deseq_sex_BAT_cold) %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10) # 14,261 genes

# New York vs Brazil (males):
res_males_BAT_cold_NYvBZ <- results(deseq_sex_BAT_cold, contrast = NewYork_males_BAT_cold - Brazil_males_BAT_cold,
                                    pAdjustMethod = "BH") # (-)LFC == higher expression in BZ males
#no_DEsigbase_genes_males_BAT_cold_NYvBZ <- res_males_BAT_cold_NYvBZ %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10, padj < 0.05) # 3,272 genes

# New York vs Brazil (females):
res_females_BAT_cold_NYvBZ <- results(deseq_sex_BAT_cold, contrast = NewYork_females_BAT_cold - Brazil_females_BAT_cold,
                                      pAdjustMethod = "BH") # (-)LFC == higher expression in BZ females
#no_DEsigbase_genes_females_BAT_cold_NYvBZ <- res_females_BAT_cold_NYvBZ %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10, padj < 0.05) # 3,884 genes

# Males vs Females (New York):
res_NY_BAT_cold_MvF <- results(deseq_sex_BAT_cold, contrast = NewYork_males_BAT_cold - NewYork_females_BAT_cold,
                               pAdjustMethod = "BH") # (-)LFC == higher expression in NY females
# Males vs Females (Brazil):
res_BZ_BAT_cold_MvF <- results(deseq_sex_BAT_cold, contrast = Brazil_males_BAT_cold - Brazil_females_BAT_cold,
                               pAdjustMethod = "BH") # (-)LFC == higher expression in BZ females
# Effect of genotype ("G")
res_sex_BAT_cold_G <- results(deseq_sex_liver_cold, contrast = NewYork_liver_cold - Brazil_liver_cold,
                                pAdjustMethod = "BH") # (-)LFC == higher expression in BZ
# Effect of environment ("S")
res_sex_BAT_cold_S <- results(deseq_sex_BAT_cold, contrast = Males_BAT_cold - Females_BAT_cold,
                              pAdjustMethod = "BH") # (-)LFC == higher expression in Females
# DESEQ's definition of genotype-by-sex ("GxS")
res_sex_BAT_cold_GxS <- results(deseq_sex_BAT_cold, contrast = (NewYork_males_BAT_cold - NewYork_females_BAT_cold) -
                                                                (Brazil_males_BAT_cold - Brazil_females_BAT_cold),
                                pAdjustMethod = "BH")
#no_DEsigbase_genes_sex_BAT_cold_GxS <- res_sex_BAT_cold_GxS %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% filter(baseMean >= 10, padj < 0.05) # 0 genes






















## plot genes of interest
# res_males_liver_E %>% data.frame() %>% rownames_to_column(var = "gene") %>% as_tibble() %>% View()
# ex_gene <- plotCounts(deseq_males_liver, gene = "ENSMUSG00000021928", intgroup = c("temperature", "population"), returnData = TRUE)
# level_order <- c("Warm", "Cold")
# ggplot(data=ex_gene, aes(x = factor(temperature, level = level_order), y=count, color = population)) + #, group = population)) +
#   geom_boxplot(position = "identity", fill = NA) + labs(x = "Temp", y = "Normalized Counts")
# level_order2 <- c("Brazil", "NewYork")
# ggplot(data=ex_gene, aes(x = factor(population, level = level_order2), y=count, color = temperature)) + #, group = population)) +
#   geom_boxplot(position = "identity", fill = NA) + labs(x = "Population", y = "Normalized Counts")

