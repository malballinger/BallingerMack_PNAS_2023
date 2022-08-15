#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script models both relative and residual thermal conductance and extremity length from
# Brazil (MANA) and New York (SARA) house mice.
# Data were cleaned using the script ./clean_phenotypes.R.
# Statistical models are used in downstream analyses of BallingerMack_2022.

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clear R's environment
library(tidyverse)
library(here)
library(lmerTest)
library(lme4)
library(emmeans)
library(car)
library(nlme)
library(effectsize)
library(report)

set.seed(19910118)

##############################################################
# Import data
##############################################################

phenotype_data <- read_csv(here("data/processed/all_phenotypes.csv"))

##############################################################
# Phenotypes 'relative' to body mass and residuals
##############################################################

### Thermal Conductance differences

## Note: it isn't well established that conductance co-varies with body size so we only
## analyze and report raw conductance values

ggplot(data = phenotype_data, aes(x = Population, y = AVG_conductance, color = Population)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(show.legend = FALSE)

TC_mod <- lmer(AVG_conductance ~ Population + Sex + (1|Generation),
               data = phenotype_data)
car::Anova(TC_mod, type = "III")
effectsize::omega_squared(TC_mod)



### Body Mass differences
# filter data to only include age-controlled mice (from experiment of Ballinger and Nachman, 2022)

age_ctrl_data <- phenotype_data %>%
  filter(Age_days < 100)

ggplot(data = age_ctrl_data, aes(x = Population, y = BodyWeight_g, color = Population)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(show.legend = FALSE)

BW_mod <- lmer(BodyWeight_g ~ Population + Sex + (1|Generation),
               data = age_ctrl_data)
car::Anova(BW_mod, type = "III")
effectsize::omega_squared(BW_mod)



### Tail Length differences

ggplot(data = age_ctrl_data, aes(x = Population, y = FinalTailLength_mm, color = Population)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(show.legend = FALSE)

TL_mod <- lmer(FinalTailLength_mm ~ BodyWeight_g +  Population + Sex + (1|Generation),
                   data = age_ctrl_data)
car::Anova(TL_mod, type = "III")
effectsize::omega_squared(TL_mod)

# body mass is a significant main effect of tail length
# let's see if relative tail length is better

age_ctrl_data <- age_ctrl_data %>%
  mutate(rel_TL = FinalTailLength_mm / BodyWeight_g)

ggplot(data = age_ctrl_data, aes(x = Population, y = rel_TL, color = Population)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(show.legend = FALSE)

rel_TL_mod <- lmer(rel_TL ~  Population + Sex + (1|Generation),
               data = age_ctrl_data)
car::Anova(rel_TL_mod, type = "III")
effectsize::omega_squared(rel_TL_mod)

# let's see how relative compares to residuals

residsTLBW <- lm(FinalTailLength_mm ~ BodyWeight_g, data = age_ctrl_data,
                 na.action = na.exclude)
age_ctrl_data$resids_TLBW <- resid(residsTLBW)

ggplot(data = age_ctrl_data, aes(x = Population, y = resids_TLBW, color = Population)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(show.legend = FALSE)

resids_TL_mod <- lmer(resids_TLBW ~ Population + Sex + (1|Generation),
                      data = age_ctrl_data)
car::Anova(resids_TL_mod, type = "III")
effectsize::omega_squared(resids_TL_mod)

# resids and relative give similar results



### Ear Length differences

ggplot(data = age_ctrl_data, aes(x = Population, y = EarLength_mm, color = Population)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(show.legend = FALSE)

EL_mod <- lmer(EarLength_mm ~ BodyWeight_g +  Population + Sex + (1|Generation),
               data = age_ctrl_data)
car::Anova(EL_mod, type = "III")
effectsize::omega_squared(EL_mod)

# body mass is a significant main effect of ear length
# let's see if relative ear length is better

age_ctrl_data <- age_ctrl_data %>%
  mutate(rel_EL = EarLength_mm / BodyWeight_g)

ggplot(data = age_ctrl_data, aes(x = Population, y = rel_EL, color = Population)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(show.legend = FALSE)

rel_EL_mod <- lmer(rel_EL ~  Population + Sex + (1|Generation),
                   data = age_ctrl_data)
car::Anova(rel_EL_mod, type = "III")
effectsize::omega_squared(rel_EL_mod)

# let's see how relative compares to residuals

residsELBW <- lm(EarLength_mm ~ BodyWeight_g, data = age_ctrl_data,
                 na.action = na.exclude)
age_ctrl_data$resids_ELBW <- resid(residsELBW)

ggplot(data = age_ctrl_data, aes(x = Population, y = resids_ELBW, color = Population)) +
  geom_boxplot(show.legend = FALSE) +
  geom_jitter(show.legend = FALSE)

resids_EL_mod <- lmer(resids_ELBW ~ Population + Sex + (1|Generation),
                      data = age_ctrl_data)
car::Anova(resids_EL_mod, type = "III")
effectsize::omega_squared(resids_EL_mod)

# resids and relative give roughly similar results



## Overall, relative and residuals give us similar results (line differences)



