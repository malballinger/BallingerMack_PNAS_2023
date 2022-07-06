#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script cleans the metadata files: data/raw/post_dissection_metadata_RAW_2021-02-18.csv &
# data/raw/thermal_phenotypes.csv.
# The cleaned datasets are joined together and used for downstream analyses in BallingerMack_2022.

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clears R's environment
library(tidyverse)
library(here)

##############################################################
# Import & clean data
##############################################################

PostDissectionMetaData <- read_csv(here("data/raw/post_dissection_metadata_RAW_2021-02-18.csv")) %>%
  select(-BMI_kg_m2, -WeaningBodyWeight_g, -WeaningTailLength_mm, -StartExpBodyWeight_g, -StartExpTailLength_mm, -BodyLength_mm, -Age_weeks) %>%
  select(-last_col()) %>%
  rename("DOD" = "DateMeasured") %>%
  filter(Population == "BRAZIL" | Population == "NEW_YORK") %>% # only keep parental populations (remove F1 hybrids)
  filter(Line == "193x255" | Line == "19x13") %>% # only want MANA and SARA
  filter(Environment == "RT") %>% # only want warm-temperature mice
  mutate(Sex = fct_recode(Sex, "Female" = "F", "Male" = "M")) %>% # spells out males and females
  mutate(Environment = fct_recode(Environment, "Warm" = "RT")) %>% # spells out warm
  select(-Line) # don't need this since only using one line per pop

ThermalConductanceData <- read_csv(here("data/raw/thermal_phenotypes.csv")) %>%
  select(-Short_ID) %>%
  rename("Mouse_ID" = "Full_ID") %>%
  mutate(Sex = fct_recode(Sex, "Female" = "F", "Male" = "M"), # spells out males and females
         Environment = fct_recode(Environment, "Warm" = "RT")) # spells out warm
  
# join data frames together
full_phenotypes <- dplyr::full_join(PostDissectionMetaData, ThermalConductanceData) 

# save new dataset in processed folder
write.csv(full_phenotypes, file = "data/processed/all_phenotypes.csv", row.names = TRUE)

