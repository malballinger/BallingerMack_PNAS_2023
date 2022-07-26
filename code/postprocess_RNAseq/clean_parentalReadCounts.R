#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script filters and subsets the counts metadata files, and
# generates matrices that are used for downstream analyses in BallingerMack_2022.

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clear R's environment
library(tidyverse)
library(here)

##############################################################
# Import data
##############################################################

Counts_Metadata <- read_table(here("./data/raw/ReadCounts/all_parents_counts.txt"), col_names = TRUE, col_types = NULL)

SampleInfo <- read_delim(here("./data/raw/ReadCounts/all_parents_sample_info.txt"))
colnames(SampleInfo)[1] <- "sampleID"

##############################################################
# Subset data
##############################################################

## liver only
counts_all_liver <- Counts_Metadata %>%
  dplyr::select(gene, starts_with("MBSD_0"))

all_liver_sampleinfo <- SampleInfo %>%
  filter(condition == "LIVER")

## BAT only
counts_all_BAT <- Counts_Metadata %>%
  dplyr::select(gene, starts_with("MBSD_1"))

all_BAT_sampleinfo <- SampleInfo %>%
  filter(condition == "BAT")


## warm only (divergence)
counts_warm <- Counts_Metadata %>% # males listed first (BZ, then NY); females listed in second half
  dplyr::select(gene, starts_with(c("MBSD_002", "MBSD_101", "MBSD_003", "MBSD_102", "MBSD_009", "MBSD_108", # start of BZ warm males
                             "MBSD_012", "MBSD_111", "MBSD_013", "MBSD_112", "MBSD_022", "MBSD_121",
                             "MBSD_055", "MBSD_154", "MBSD_061", "MBSD_160", "MBSD_065", "MBSD_164", # start of NY warm males
                             "MBSD_072", "MBSD_171", "MBSD_075", "MBSD_174", "MBSD_077", "MBSD_176",
                             "MBSD_004", "MBSD_103", "MBSD_006", "MBSD_105", "MBSD_008", "MBSD_107", # start of BZ warm females
                             "MBSD_010", "MBSD_109", "MBSD_017", "MBSD_116", "MBSD_019", "MBSD_118",
                             "MBSD_054", "MBSD_153", "MBSD_056", "MBSD_155", "MBSD_062", "MBSD_161", # start of NY warm females
                             "MBSD_068", "MBSD_167", "MBSD_073", "MBSD_172", "MBSD_076", "MBSD_175")))

sampleinfo_warm <- SampleInfo %>%
  filter(temperature == "Warm") %>% # males listed first (BZ, then NY); females listed in second half
  arrange(match(sampleID, c("MBSD_002.counts", "MBSD_101.counts", "MBSD_003.counts", "MBSD_102.counts", "MBSD_009.counts", "MBSD_108.counts", # start of BZ males
                            "MBSD_012.counts", "MBSD_111.counts", "MBSD_013.counts", "MBSD_112.counts", "MBSD_022.counts", "MBSD_121.counts",
                            "MBSD_055.counts", "MBSD_154.counts", "MBSD_061.counts", "MBSD_160.counts", "MBSD_065.counts", "MBSD_164.counts", # start of NY males
                            "MBSD_072.counts", "MBSD_171.counts", "MBSD_075.counts", "MBSD_174.counts", "MBSD_077.counts", "MBSD_176.counts",
                            "MBSD_004.counts", "MBSD_103.counts", "MBSD_006.counts", "MBSD_105.counts", "MBSD_008.counts", "MBSD_107.counts", # start of BZ females
                            "MBSD_010.counts", "MBSD_109.counts", "MBSD_017.counts", "MBSD_116.counts", "MBSD_019.counts", "MBSD_118.counts",
                            "MBSD_054.counts", "MBSD_153.counts", "MBSD_056.counts", "MBSD_155.counts", "MBSD_062.counts", "MBSD_161.counts", # start of NY females
                            "MBSD_068.counts", "MBSD_167.counts", "MBSD_073.counts", "MBSD_172.counts", "MBSD_076.counts", "MBSD_175.counts")))

# warm, liver only
counts_liver_warm <- counts_warm %>%
  dplyr::select(gene, starts_with("MBSD_0"))

sampleinfo_liver_warm <- sampleinfo_warm %>%
  filter(condition == "LIVER") %>% # males listed first (BZ, then NY); females listed in second half
  arrange(match(sampleID, c("MBSD_002.counts", "MBSD_003.counts", "MBSD_009.counts", "MBSD_012.counts", "MBSD_013.counts", "MBSD_022.counts", # start of males
                            "MBSD_055.counts", "MBSD_061.counts", "MBSD_065.counts", "MBSD_072.counts", "MBSD_075.counts", "MBSD_077.counts",
                            "MBSD_004.counts", "MBSD_006.counts", "MBSD_008.counts", "MBSD_010.counts", "MBSD_017.counts", "MBSD_019.counts", # start of females
                            "MBSD_054.counts", "MBSD_056.counts", "MBSD_062.counts", "MBSD_068.counts", "MBSD_073.counts", "MBSD_076.counts")))

# warm, BAT only
counts_BAT_warm <- counts_warm %>%
  dplyr::select(gene, starts_with("MBSD_1"))

sampleinfo_BAT_warm <- sampleinfo_warm %>%
  filter(condition == "BAT") %>% # males listed first (BZ, then NY); females listed in second half
  arrange(match(sampleID, c("MBSD_101.counts", "MBSD_102.counts", "MBSD_108.counts", "MBSD_111.counts", "MBSD_112.counts", "MBSD_121.counts", # start of males
                            "MBSD_154.counts", "MBSD_160.counts", "MBSD_164.counts", "MBSD_171.counts", "MBSD_174.counts", "MBSD_176.counts",
                            "MBSD_103.counts", "MBSD_105.counts", "MBSD_107.counts", "MBSD_109.counts", "MBSD_116.counts", "MBSD_118.counts", # start of females
                            "MBSD_153.counts", "MBSD_155.counts", "MBSD_161.counts", "MBSD_167.counts", "MBSD_172.counts", "MBSD_175.counts")))


## cold only
counts_cold <- Counts_Metadata %>% # males listed first (BZ, then NY); females listed in second half
  dplyr::select(gene, starts_with(c("MBSD_005", "MBSD_104", "MBSD_007", "MBSD_106", "MBSD_011", "MBSD_110", # start of BZ males
                             "MBSD_015", "MBSD_114", "MBSD_020", "MBSD_119", "MBSD_024", "MBSD_123",
                             "MBSD_057", "MBSD_156", "MBSD_059", "MBSD_158", "MBSD_063", "MBSD_162", # start of NY males
                             "MBSD_067", "MBSD_166", "MBSD_069", "MBSD_168", "MBSD_074", "MBSD_173",
                             "MBSD_001", "MBSD_100", "MBSD_014", "MBSD_113", "MBSD_016", "MBSD_115", # start of BZ females
                             "MBSD_018", "MBSD_117", "MBSD_021", "MBSD_120", "MBSD_023", "MBSD_122",
                             "MBSD_058", "MBSD_157", "MBSD_060", "MBSD_159", "MBSD_064", "MBSD_163", # start of NY females
                             "MBSD_066", "MBSD_165", "MBSD_070", "MBSD_169", "MBSD_071", "MBSD_170")))

sampleinfo_cold <- SampleInfo %>%
  filter(temperature == "Cold") %>% # males listed first (BZ, then NY); females listed in second half
  arrange(match(sampleID, c("MBSD_005.counts", "MBSD_104.counts", "MBSD_007.counts", "MBSD_106.counts", "MBSD_011.counts", "MBSD_110.counts", # start of BZ males
                            "MBSD_015.counts", "MBSD_114.counts", "MBSD_020.counts", "MBSD_119.counts", "MBSD_024.counts", "MBSD_123.counts",
                            "MBSD_057.counts", "MBSD_156.counts", "MBSD_059.counts", "MBSD_158.counts", "MBSD_063.counts", "MBSD_162.counts", # start of NY males
                            "MBSD_067.counts", "MBSD_166.counts", "MBSD_069.counts", "MBSD_168.counts", "MBSD_074.counts", "MBSD_173.counts",
                            "MBSD_001.counts", "MBSD_100.counts", "MBSD_014.counts", "MBSD_113.counts", "MBSD_016.counts", "MBSD_115.counts", # start of BZ females
                            "MBSD_018.counts", "MBSD_117.counts", "MBSD_021.counts", "MBSD_120.counts", "MBSD_023.counts", "MBSD_122.counts",
                            "MBSD_058.counts", "MBSD_157.counts", "MBSD_060.counts", "MBSD_159.counts", "MBSD_064.counts", "MBSD_163.counts", # start of NY females
                            "MBSD_066.counts", "MBSD_165.counts", "MBSD_070.counts", "MBSD_169.counts", "MBSD_071.counts", "MBSD_170.counts")))

# cold, liver only
counts_liver_cold <- counts_cold %>%
  dplyr::select(gene, starts_with("MBSD_0"))

sampleinfo_liver_cold <- sampleinfo_cold %>%
  filter(condition == "LIVER") %>% # males listed first (BZ, then NY); females listed in second half
  arrange(match(sampleID, c("MBSD_005.counts", "MBSD_007.counts", "MBSD_011.counts", "MBSD_015.counts", "MBSD_020.counts", "MBSD_024.counts", # start of males
                            "MBSD_057.counts", "MBSD_059.counts", "MBSD_063.counts", "MBSD_067.counts", "MBSD_069.counts", "MBSD_074.counts",
                            "MBSD_001.counts", "MBSD_014.counts", "MBSD_016.counts", "MBSD_018.counts", "MBSD_021.counts", "MBSD_023.counts", # start of females
                            "MBSD_058.counts", "MBSD_060.counts", "MBSD_064.counts", "MBSD_066.counts", "MBSD_070.counts", "MBSD_071.counts")))

# cold, BAT only
counts_BAT_cold <- counts_cold %>%
  dplyr::select(gene, starts_with("MBSD_1"))

sampleinfo_BAT_cold <- sampleinfo_cold %>%
  filter(condition == "BAT") %>% # males listed first (BZ, then NY); females listed in second half
  arrange(match(sampleID, c("MBSD_104.counts", "MBSD_106.counts", "MBSD_110.counts", "MBSD_114.counts", "MBSD_119.counts", "MBSD_123.counts", # start of males
                            "MBSD_156.counts", "MBSD_158.counts", "MBSD_162.counts", "MBSD_166.counts", "MBSD_168.counts", "MBSD_173.counts",
                            "MBSD_100.counts", "MBSD_113.counts", "MBSD_115.counts", "MBSD_117.counts", "MBSD_120.counts", "MBSD_122.counts", # start of females
                            "MBSD_157.counts", "MBSD_159.counts", "MBSD_163.counts", "MBSD_165.counts", "MBSD_169.counts", "MBSD_170.counts")))


## males only
counts_males <- Counts_Metadata %>% # BZ listed first; NY listed in second half
  dplyr::select(gene, starts_with(c("MBSD_002", "MBSD_101", "MBSD_003", "MBSD_102", "MBSD_005", "MBSD_104", # start of BZ males
                             "MBSD_007", "MBSD_106", "MBSD_009", "MBSD_108", "MBSD_011", "MBSD_110",
                             "MBSD_012", "MBSD_111", "MBSD_013", "MBSD_112", "MBSD_015", "MBSD_114",
                             "MBSD_020", "MBSD_119", "MBSD_022", "MBSD_121", "MBSD_024", "MBSD_123",
                             "MBSD_055", "MBSD_154", "MBSD_057", "MBSD_156", "MBSD_059", "MBSD_158", # start of NY males
                             "MBSD_061", "MBSD_160", "MBSD_063", "MBSD_162", "MBSD_065", "MBSD_164",
                             "MBSD_067", "MBSD_166", "MBSD_069", "MBSD_168", "MBSD_072", "MBSD_171",
                             "MBSD_074", "MBSD_173", "MBSD_075", "MBSD_174", "MBSD_077", "MBSD_176")))

sampleinfo_males <- SampleInfo %>%
  filter(sex == "M") %>%
  arrange(match(sampleID, c("MBSD_002.counts", "MBSD_101.counts", "MBSD_003.counts", "MBSD_102.counts", "MBSD_005.counts", "MBSD_104.counts", # start of BZ males
                            "MBSD_007.counts", "MBSD_106.counts", "MBSD_009.counts", "MBSD_108.counts", "MBSD_011.counts", "MBSD_110.counts",
                            "MBSD_012.counts", "MBSD_111.counts", "MBSD_013.counts", "MBSD_112.counts", "MBSD_015.counts", "MBSD_114.counts",
                            "MBSD_020.counts", "MBSD_119.counts", "MBSD_022.counts", "MBSD_121.counts", "MBSD_024.counts", "MBSD_123.counts",
                            "MBSD_055.counts", "MBSD_154.counts", "MBSD_057.counts", "MBSD_156.counts", "MBSD_059.counts", "MBSD_158.counts", # start of NY males
                            "MBSD_061.counts", "MBSD_160.counts", "MBSD_063.counts", "MBSD_162.counts", "MBSD_065.counts", "MBSD_164.counts",
                            "MBSD_067.counts", "MBSD_166.counts", "MBSD_069.counts", "MBSD_168.counts", "MBSD_072.counts", "MBSD_171.counts",
                            "MBSD_074.counts", "MBSD_173.counts", "MBSD_075.counts", "MBSD_174.counts", "MBSD_077.counts", "MBSD_176.counts")))

# males, liver only
counts_males_liver <- counts_males %>%
  dplyr::select(gene, starts_with("MBSD_0"))

sampleinfo_males_liver <- sampleinfo_males %>%
  filter(condition == "LIVER") %>% # BZ listed first; NY listed in second half
  arrange(match(sampleID, c("MBSD_002.counts", "MBSD_003.counts", "MBSD_005.counts", "MBSD_007.counts", "MBSD_009.counts", "MBSD_011.counts", # start of BZ males
                            "MBSD_012.counts", "MBSD_013.counts", "MBSD_015.counts", "MBSD_020.counts", "MBSD_022.counts", "MBSD_024.counts",
                            "MBSD_055.counts", "MBSD_057.counts", "MBSD_059.counts", "MBSD_061.counts", "MBSD_063.counts", "MBSD_065.counts", # start of NY males
                            "MBSD_067.counts", "MBSD_069.counts", "MBSD_072.counts", "MBSD_074.counts", "MBSD_075.counts", "MBSD_077.counts")))

# males, BAT only
counts_males_BAT <- counts_males %>%
  dplyr::select(gene, starts_with("MBSD_1"))

sampleinfo_males_BAT <- sampleinfo_males %>%
  filter(condition == "BAT") %>% # BZ listed first; NY listed in second half
  arrange(match(sampleID, c("MBSD_101.counts", "MBSD_102.counts", "MBSD_104.counts", "MBSD_106.counts", "MBSD_108.counts", "MBSD_110.counts", # start of BZ males
                            "MBSD_111.counts", "MBSD_112.counts", "MBSD_114.counts", "MBSD_119.counts", "MBSD_121.counts", "MBSD_123.counts",
                            "MBSD_154.counts", "MBSD_156.counts", "MBSD_158.counts", "MBSD_160.counts", "MBSD_162.counts", "MBSD_164.counts", # start of NY males
                            "MBSD_166.counts", "MBSD_168.counts", "MBSD_171.counts", "MBSD_173.counts", "MBSD_174.counts", "MBSD_176.counts")))


# warm males, both tissues
counts_males_warm <- Counts_Metadata %>% # BZ listed first; NY listed in second half
  dplyr::select(gene, starts_with(c("MBSD_002", "MBSD_101", "MBSD_003", "MBSD_102", "MBSD_009", "MBSD_108", # start of warm BZ males
                             "MBSD_012", "MBSD_111", "MBSD_013", "MBSD_112", "MBSD_022", "MBSD_121",
                             "MBSD_055", "MBSD_154", "MBSD_061", "MBSD_160", "MBSD_065", "MBSD_164", # start of warm NY males
                             "MBSD_072", "MBSD_171", "MBSD_075", "MBSD_174", "MBSD_077", "MBSD_176")))


sampleinfo_males_warm <- SampleInfo %>%
  filter(temperature == "Warm") %>%
  filter(sex == "M") %>% # BZ listed first; NY listed in second half
  arrange(match(sampleID, c("MBSD_002.counts", "MBSD_101.counts", "MBSD_003.counts", "MBSD_102.counts", "MBSD_009.counts", "MBSD_108.counts", # start of warm BZ males
                            "MBSD_012.counts", "MBSD_111.counts", "MBSD_013.counts", "MBSD_112.counts", "MBSD_022.counts", "MBSD_121.counts",
                            "MBSD_055.counts", "MBSD_154.counts", "MBSD_061.counts", "MBSD_160.counts", "MBSD_065.counts", "MBSD_164.counts", # start of warm NY males
                            "MBSD_072.counts", "MBSD_171.counts", "MBSD_075.counts", "MBSD_174.counts", "MBSD_077.counts", "MBSD_176.counts")))

# warm males, liver only
counts_males_warm_liver <- counts_males_warm %>%
  dplyr::select(gene, starts_with("MBSD_0"))

sampleinfo_males_warm_liver <- sampleinfo_males_warm %>%
  filter(condition == "LIVER") %>% # BZ listed first, then NY
  arrange(match(sampleID, c("MBSD_002.counts", "MBSD_003.counts", "MBSD_009.counts", "MBSD_012.counts", "MBSD_013.counts", "MBSD_022.counts",
                            "MBSD_055.counts", "MBSD_061.counts", "MBSD_065.counts", "MBSD_072.counts", "MBSD_075.counts", "MBSD_077.counts")))

# warm males, BAT only
counts_males_warm_BAT <- counts_males_warm %>%
  dplyr::select(gene, starts_with("MBSD_1"))

sampleinfo_males_warm_BAT <- sampleinfo_males_warm %>%
  filter(condition == "BAT") %>% # BZ listed first, then NY
  arrange(match(sampleID, c("MBSD_101.counts", "MBSD_102.counts", "MBSD_108.counts", "MBSD_111.counts", "MBSD_112.counts", "MBSD_121.counts",
                            "MBSD_154.counts", "MBSD_160.counts", "MBSD_164.counts", "MBSD_171.counts", "MBSD_174.counts", "MBSD_176.counts")))

# cold males, both tissues
counts_males_cold <- Counts_Metadata %>% # BZ listed first; NY listed in second half
  dplyr::select(gene, starts_with(c("MBSD_005", "MBSD_104", "MBSD_007", "MBSD_106", "MBSD_011", "MBSD_110", # start of cold BZ males
                             "MBSD_015", "MBSD_114", "MBSD_020", "MBSD_119", "MBSD_024", "MBSD_123",
                             "MBSD_057", "MBSD_156", "MBSD_059", "MBSD_158", "MBSD_063", "MBSD_162", # start of cold NY males
                             "MBSD_067", "MBSD_166", "MBSD_069", "MBSD_168", "MBSD_074", "MBSD_173")))

sampleinfo_males_cold <- SampleInfo %>%
  filter(temperature == "Cold") %>%
  filter(sex == "M") %>% # BZ listed first; NY listed in second half
  arrange(match(sampleID, c("MBSD_005.counts", "MBSD_104.counts", "MBSD_007.counts", "MBSD_106.counts", "MBSD_011.counts", "MBSD_110.counts", # start of cold BZ males
                            "MBSD_015.counts", "MBSD_114.counts", "MBSD_020.counts", "MBSD_119.counts", "MBSD_024.counts", "MBSD_123.counts",
                            "MBSD_057.counts", "MBSD_156.counts", "MBSD_059.counts", "MBSD_158.counts", "MBSD_063.counts", "MBSD_162.counts", # start of cold NY males
                            "MBSD_067.counts", "MBSD_166.counts", "MBSD_069.counts", "MBSD_168.counts", "MBSD_074.counts", "MBSD_173.counts")))

# cold males, liver only
counts_males_cold_liver <- counts_males_cold %>%
  dplyr::select(gene, starts_with("MBSD_0"))

sampleinfo_males_cold_liver <- sampleinfo_males_cold %>%
  filter(condition == "LIVER") %>% # BZ listed first, then NY
  arrange(match(sampleID, c("MBSD_005.counts", "MBSD_007.counts", "MBSD_011.counts", "MBSD_015.counts", "MBSD_020.counts", "MBSD_024.counts",
                            "MBSD_057.counts", "MBSD_059.counts", "MBSD_063.counts", "MBSD_067.counts", "MBSD_069.counts", "MBSD_074.counts")))

# cold males, BAT only
counts_males_cold_BAT <- counts_males_cold %>%
  dplyr::select(gene, starts_with("MBSD_1"))

sampleinfo_males_cold_BAT <- sampleinfo_males_cold %>%
  filter(condition == "BAT") %>% # BZ listed first, then NY
  arrange(match(sampleID, c("MBSD_104.counts", "MBSD_106.counts", "MBSD_110.counts", "MBSD_114.counts", "MBSD_119.counts", "MBSD_123.counts",
                            "MBSD_156.counts", "MBSD_158.counts", "MBSD_162.counts", "MBSD_166.counts", "MBSD_168.counts", "MBSD_173.counts")))


## females only, both tissues
counts_females <- Counts_Metadata %>% # BZ listed first (warm, then cold); NY listed in second half
  dplyr::select(gene, starts_with(c("MBSD_004", "MBSD_103", "MBSD_006", "MBSD_105", "MBSD_008", "MBSD_107", # start of BZ warm females
                             "MBSD_010", "MBSD_109", "MBSD_017", "MBSD_116", "MBSD_019", "MBSD_118",
                             "MBSD_001", "MBSD_100", "MBSD_014", "MBSD_113", "MBSD_016", "MBSD_115", # start of BZ cold females
                             "MBSD_018", "MBSD_117", "MBSD_021", "MBSD_120", "MBSD_023", "MBSD_122",
                             "MBSD_054", "MBSD_153", "MBSD_056", "MBSD_155", "MBSD_062", "MBSD_161", # start of NY warm females
                             "MBSD_068", "MBSD_167", "MBSD_073", "MBSD_172", "MBSD_076", "MBSD_175",
                             "MBSD_058", "MBSD_157", "MBSD_060", "MBSD_159", "MBSD_064", "MBSD_163", # start of NY cold females
                             "MBSD_066", "MBSD_165", "MBSD_070", "MBSD_169", "MBSD_071", "MBSD_170")))


sampleinfo_females <- SampleInfo %>%
  filter(sex == "F") %>% # BZ listed first (warm, then cold); NY listed in second half
  arrange(match(sampleID, c("MBSD_004.counts", "MBSD_103.counts", "MBSD_006.counts", "MBSD_105.counts", "MBSD_008.counts", "MBSD_107.counts", # start of BZ warm females
                            "MBSD_010.counts", "MBSD_109.counts", "MBSD_017.counts", "MBSD_116.counts", "MBSD_019.counts", "MBSD_118.counts",
                            "MBSD_001.counts", "MBSD_100.counts", "MBSD_014.counts", "MBSD_113.counts", "MBSD_016.counts", "MBSD_115.counts", # start of BZ cold females
                            "MBSD_018.counts", "MBSD_117.counts", "MBSD_021.counts", "MBSD_120.counts", "MBSD_023.counts", "MBSD_122.counts",
                            "MBSD_054.counts", "MBSD_153.counts", "MBSD_056.counts", "MBSD_155.counts", "MBSD_062.counts", "MBSD_161.counts", # start of NY warm females
                            "MBSD_068.counts", "MBSD_167.counts", "MBSD_073.counts", "MBSD_172.counts", "MBSD_076.counts", "MBSD_175.counts",
                            "MBSD_058.counts", "MBSD_157.counts", "MBSD_060.counts", "MBSD_159.counts", "MBSD_064.counts", "MBSD_163.counts", # start of NY cold females
                            "MBSD_066.counts", "MBSD_165.counts", "MBSD_070.counts", "MBSD_169.counts", "MBSD_071.counts", "MBSD_170.counts")))

# females, liver only
counts_females_liver <- counts_females %>%
  dplyr::select(gene, starts_with("MBSD_0"))

sampleinfo_females_liver <- sampleinfo_females %>%
  filter(condition == "LIVER") %>% # BZ listed first (warm, then cold); NY listed in second half
  arrange(match(sampleID, c("MBSD_004.counts", "MBSD_006.counts", "MBSD_008.counts", "MBSD_010.counts", "MBSD_017.counts", "MBSD_019.counts", # BZ warm females
                            "MBSD_001.counts", "MBSD_014.counts", "MBSD_016.counts", "MBSD_018.counts", "MBSD_021.counts", "MBSD_023.counts", # BZ cold females
                            "MBSD_054.counts", "MBSD_056.counts", "MBSD_062.counts", "MBSD_068.counts", "MBSD_073.counts", "MBSD_076.counts", # NY warm females
                            "MBSD_058.counts", "MBSD_060.counts", "MBSD_064.counts", "MBSD_066.counts", "MBSD_070.counts", "MBSD_071.counts"))) # NY cold females

# females, BAT only
counts_females_BAT <- counts_females %>%
  dplyr::select(gene, starts_with("MBSD_1"))

sampleinfo_females_BAT <- sampleinfo_females %>%
  filter(condition == "BAT") %>% # BZ listed first (warm, then cold); NY listed in second half
  arrange(match(sampleID, c("MBSD_103.counts", "MBSD_105.counts", "MBSD_107.counts", "MBSD_109.counts", "MBSD_116.counts", "MBSD_118.counts", # BZ warm females
                            "MBSD_100.counts", "MBSD_113.counts", "MBSD_115.counts", "MBSD_117.counts", "MBSD_120.counts", "MBSD_122.counts", # BZ cold females
                            "MBSD_153.counts", "MBSD_155.counts", "MBSD_161.counts", "MBSD_167.counts", "MBSD_172.counts", "MBSD_175.counts", # NY warm females
                            "MBSD_157.counts", "MBSD_159.counts", "MBSD_163.counts", "MBSD_165.counts", "MBSD_169.counts", "MBSD_170.counts"))) # NY cold females


# warm females, both tissues
counts_females_warm <- Counts_Metadata %>% # BZ listed first, then NY
  dplyr::select(gene, starts_with(c("MBSD_004", "MBSD_103", "MBSD_006", "MBSD_105", "MBSD_008", "MBSD_107", # start of BZ warm females
                             "MBSD_010", "MBSD_109", "MBSD_017", "MBSD_116", "MBSD_019", "MBSD_118",
                             "MBSD_054", "MBSD_153", "MBSD_056", "MBSD_155", "MBSD_062", "MBSD_161", # start of NY warm females
                             "MBSD_068", "MBSD_167", "MBSD_073", "MBSD_172", "MBSD_076", "MBSD_175")))

sampleinfo_females_warm <- SampleInfo %>%
  filter(temperature == "Warm") %>%
  filter(sex == "F") %>% # BZ listed first, then NY
  arrange(match(sampleID, c("MBSD_004.counts", "MBSD_103.counts", "MBSD_006.counts", "MBSD_105.counts", "MBSD_008.counts", "MBSD_107.counts", # start of BZ warm females
                            "MBSD_010.counts", "MBSD_109.counts", "MBSD_017.counts", "MBSD_116.counts", "MBSD_019.counts", "MBSD_118.counts",
                            "MBSD_054.counts", "MBSD_153.counts", "MBSD_056.counts", "MBSD_155.counts", "MBSD_062.counts", "MBSD_161.counts", # start of NY warm females
                            "MBSD_068.counts", "MBSD_167.counts", "MBSD_073.counts", "MBSD_172.counts", "MBSD_076.counts", "MBSD_175.counts")))

# warm females, liver only
counts_females_warm_liver <- counts_females_warm %>%
  dplyr::select(gene, starts_with("MBSD_0"))

sampleinfo_females_warm_liver <- sampleinfo_females_warm %>%
  filter(condition == "LIVER") %>% # BZ listed first, then NY
  arrange(match(sampleID, c("MBSD_004.counts", "MBSD_006.counts", "MBSD_008.counts", "MBSD_010.counts", "MBSD_017.counts", "MBSD_019.counts",
                            "MBSD_054.counts", "MBSD_056.counts", "MBSD_062.counts", "MBSD_068.counts", "MBSD_073.counts", "MBSD_076.counts")))

# warm females, BAT only
counts_females_warm_BAT <- counts_females_warm %>%
  dplyr::select(gene, starts_with("MBSD_1"))

sampleinfo_females_warm_BAT <- sampleinfo_females_warm %>%
  filter(condition == "BAT") %>% # BZ listed first, then NY
  arrange(match(sampleID, c("MBSD_103.counts", "MBSD_105.counts", "MBSD_107.counts", "MBSD_109.counts", "MBSD_116.counts", "MBSD_118.counts",
                            "MBSD_153.counts", "MBSD_155.counts", "MBSD_161.counts", "MBSD_167.counts", "MBSD_172.counts", "MBSD_175.counts")))

# cold females, both tissues
counts_females_cold <- Counts_Metadata %>% # BZ listed first, then NY
  dplyr::select(gene, starts_with(c("MBSD_001", "MBSD_100", "MBSD_014", "MBSD_113", "MBSD_016", "MBSD_115", # start of BZ cold females
                             "MBSD_018", "MBSD_117", "MBSD_021", "MBSD_120", "MBSD_023", "MBSD_122",
                             "MBSD_058", "MBSD_157", "MBSD_060", "MBSD_159", "MBSD_064", "MBSD_163", # start of NY cold females
                             "MBSD_066", "MBSD_165", "MBSD_070", "MBSD_169", "MBSD_071", "MBSD_170")))

sampleinfo_females_cold <- SampleInfo %>%
  filter(temperature == "Cold") %>%
  filter(sex == "F") %>% # BZ listed first, then NY
  arrange(match(sampleID, c("MBSD_001.counts", "MBSD_100.counts", "MBSD_014.counts", "MBSD_113.counts", "MBSD_016.counts", "MBSD_115.counts", # start of BZ cold females
                            "MBSD_018.counts", "MBSD_117.counts", "MBSD_021.counts", "MBSD_120.counts", "MBSD_023.counts", "MBSD_122.counts",
                            "MBSD_058.counts", "MBSD_157.counts", "MBSD_060.counts", "MBSD_159.counts", "MBSD_064.counts", "MBSD_163.counts", # start of NY cold females
                            "MBSD_066.counts", "MBSD_165.counts", "MBSD_070.counts", "MBSD_169.counts", "MBSD_071.counts", "MBSD_170.counts")))

# cold females, liver only
counts_females_cold_liver <- counts_females_cold %>%
  dplyr::select(gene, starts_with("MBSD_0"))

sampleinfo_females_cold_liver <- sampleinfo_females_cold %>%
  filter(condition == "LIVER") %>% # BZ listed first, then NY
  arrange(match(sampleID, c("MBSD_001.counts", "MBSD_014.counts", "MBSD_016.counts", "MBSD_018.counts", "MBSD_021.counts", "MBSD_023.counts",
                            "MBSD_058.counts", "MBSD_060.counts", "MBSD_064.counts", "MBSD_066.counts", "MBSD_070.counts", "MBSD_071.counts")))

# cold females, BAT only
counts_females_cold_BAT <- counts_females_cold %>%
  dplyr::select(gene, starts_with("MBSD_1"))

sampleinfo_females_cold_BAT <- sampleinfo_females_cold %>%
  filter(condition == "BAT") %>% # BZ listed first, then NY
  arrange(match(sampleID, c("MBSD_100.counts", "MBSD_113.counts", "MBSD_115.counts", "MBSD_117.counts", "MBSD_120.counts", "MBSD_122.counts",
                            "MBSD_157.counts", "MBSD_159.counts", "MBSD_163.counts", "MBSD_165.counts", "MBSD_169.counts", "MBSD_170.counts")))



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
all_liver_cts_mtrx <- tibb_to_mtrx(tibb = counts_all_liver) # all_liver_sampleinfo
all_BAT_cts_mtrx <- tibb_to_mtrx(tibb = counts_all_BAT) # all_BAT_sampleinfo
warm_cts_mtrx <- tibb_to_mtrx(tibb = counts_warm) # sampleinfo_warm
warm_liver_cts_mtrx <- tibb_to_mtrx(tibb = counts_liver_warm) # sampleinfo_liver_warm
warm_BAT_cts_mtrx <- tibb_to_mtrx(tibb = counts_BAT_warm) # sampleinfo_BAT_warm
cold_cts_mtrx <- tibb_to_mtrx(tibb = counts_cold) # sampleinfo_cold
cold_liver_cts_mtrx <- tibb_to_mtrx(tibb = counts_liver_cold) # sampleinfo_liver_cold
cold_BAT_cts_mtrx <- tibb_to_mtrx(tibb = counts_BAT_cold) # sampleinfo_BAT_cold

males_cts_mtrx <- tibb_to_mtrx(tibb = counts_males) # sampleinfo_males
males_liver_cts_mtrx <- tibb_to_mtrx(tibb = counts_males_liver) # sampleinfo_males_liver
males_BAT_cts_mtrx <- tibb_to_mtrx(tibb = counts_males_BAT) # sampleinfo_males_BAT
males_warm_cts_mtrx <- tibb_to_mtrx(tibb = counts_males_warm) # sampleinfo_males_warm
males_warm_liver_cts_mtrx <- tibb_to_mtrx(tibb = counts_males_warm_liver) # sampleinfo_males_warm_liver
males_warm_BAT_cts_mtrx <- tibb_to_mtrx(tibb = counts_males_warm_BAT) # sampleinfo_males_warm_BAT
males_cold_cts_mtrx <- tibb_to_mtrx(tibb = counts_males_cold) # sampleinfo_males_cold
males_cold_liver_cts_mtrx <- tibb_to_mtrx(tibb = counts_males_cold_liver) # sampleinfo_males_cold_liver
males_cold_BAT_cts_mtrx <- tibb_to_mtrx(tibb = counts_males_cold_BAT) # sampleinfo_males_cold_BAT

females_cts_mtrx <- tibb_to_mtrx(tibb = counts_females) # sampleinfo_females
females_liver_cts_mtrx <- tibb_to_mtrx(tibb = counts_females_liver) # sampleinfo_females_liver
females_BAT_cts_mtrx <- tibb_to_mtrx(tibb = counts_females_BAT) # sampleinfo_females_BAT
females_warm_cts_mtrx <- tibb_to_mtrx(tibb = counts_females_warm) # sampleinfo_females_warm
females_warm_liver_cts_mtrx <- tibb_to_mtrx(tibb = counts_females_warm_liver) # sampleinfo_females_warm_liver
females_warm_BAT_cts_mtrx <- tibb_to_mtrx(tibb = counts_females_warm_BAT) # sampleinfo_females_warm_BAT
females_cold_cts_mtrx <- tibb_to_mtrx(tibb = counts_females_cold) # sampleinfo_females_cold
females_cold_liver_cts_mtrx <- tibb_to_mtrx(tibb = counts_females_cold_liver) # sampleinfo_females_cold_liver
females_cold_BAT_cts_mtrx <- tibb_to_mtrx(tibb = counts_females_cold_BAT) # sampleinfo_females_cold_BAT

