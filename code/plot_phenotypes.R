#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script plots differences in morphology and thermal pelage conductance between
# Brazil (MANA) and New York (SARA) mice.
# The generated plots are used in Figure 1B in BallingerMack_2022.

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clears R's environment
library(tidyverse)
library(here)
library(glue)
library(ggtext)

#calc_element("axis.text.x", theme_bw()) # determine default colors

##############################################################
# Import data
##############################################################

phenotype_data <- read_csv(here("data/processed/all_phenotypes.csv")) %>%
  mutate(rel_conductance = AVG_conductance / BodyWeight_g,
         rel_TL = FinalTailLength_mm / BodyWeight_g,
         rel_EL = EarLength_mm / BodyWeight_g)

## Get values from statistical analyses

# Refer to ./model_phenotypes.R for model comparisons and statistical analyses

BW_deets <- ("&#42;line<br>\
              &#42;sex")

TC_deets <- ("&#42;line")

EL_deets <- ("&#42;line<br>\
              &#42;sex")

TL_deets <- ("&#42;line<br>\
              &#42;sex")

##############################################################
# Plot data
##############################################################

## thermal conductance

TC_plot <- phenotype_data %>%
  ggplot(aes(x=Population, y=AVG_conductance)) +
  geom_boxplot(aes(fill = Population), position = position_dodge(width=0.25), alpha=0.8, size = 0.5, outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(position=position_jitter(width = 0.15, seed=19910118), alpha = 0.15, show.legend = FALSE) +
  scale_fill_manual(values=c("#E09832", "#1683B5"),
                    breaks=c("BRAZIL", "NEW_YORK"),
                    labels=c("Brazil", "New York")) +
  scale_x_discrete(breaks=c("BRAZIL", "NEW_YORK"),
                   labels=c("Brazil", "New York")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_markdown(vjust = 1, size = 13, family = "sans"),
        axis.text = element_text(size = 12.5, family = "sans"),
        axis.text.x = element_text(color = "black"),
        plot.tag = element_markdown(family = "sans", size = 10.5, face = "italic", hjust = 1, color = "grey30"),
        plot.tag.position = c(0.97,0.95)) +
  labs(x = "",
       y = "Pelage Conductance",
       tag = TC_deets)
ggsave("results/figures/pelage_cond.pdf", plot = TC_plot, height = 2.75, width = 2.8, units = "in")


### for body mass and extremity length, filter for age
age_ctrl_data <- phenotype_data %>% filter(Age_days < 100)

## body mass

BW_plot <- age_ctrl_data %>%
  ggplot(aes(x=Population, y=BodyWeight_g)) +
  geom_boxplot(aes(fill = Population), position = position_dodge(width=0.25), alpha=0.8, size = 0.5, outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(position=position_jitter(width = 0.15, seed=19910118), alpha = 0.15, show.legend = FALSE) +
  scale_fill_manual(values=c("#E09832", "#1683B5"),
                    breaks=c("BRAZIL", "NEW_YORK"),
                    labels=c("Brazil", "New York")) +
  scale_x_discrete(breaks=c("BRAZIL", "NEW_YORK"),
                   labels=c("Brazil", "New York")) +
  scale_y_continuous(breaks = seq(from=9, to=20, by=2.5), labels = seq(from=9, to=20, by=2.5), limits = c(9,20)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 1, size = 13, family = "sans"),
        axis.text = element_text(size = 12.5, family = "sans"),
        axis.text.x = element_text(color = "black"),
        plot.tag = element_markdown(family = "sans", size = 10.5, face = "italic", hjust = 1, color = "grey30"),
        plot.tag.position = c(0.97,0.93)) +
  labs(x = "",
       y = "Body Mass (g)",
       tag = BW_deets)
ggsave("results/figures/BW_diff.pdf", plot = BW_plot, height = 2.75, width = 2.8, units = "in")

# plot tail length

rel_TL_plot <- age_ctrl_data %>%
  ggplot(aes(x=Population, y=rel_TL)) +
  geom_boxplot(aes(fill = Population), position = position_dodge(width=0.25), alpha=0.8, size = 0.5, outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(position=position_jitter(width = 0.15, seed=19910118), alpha = 0.15, show.legend = FALSE) +
  scale_fill_manual(values=c("#E09832", "#1683B5"),
                    breaks=c("BRAZIL", "NEW_YORK"),
                    labels=c("Brazil", "New York")) +
  scale_x_discrete(breaks=c("BRAZIL", "NEW_YORK"),
                   labels=c("Brazil", "New York")) +
  scale_y_continuous(breaks = seq(from=4, to=8.5, by=1.5), labels = seq(from=4, to=8.5, by=1.5), limits = c(4,8.5)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 1, size = 13, family = "sans"),
        axis.text = element_text(size = 12.5, family = "sans"),
        axis.text.x = element_text(color = "black"),
        plot.tag = element_markdown(family = "sans", size = 10.5, face = "italic", hjust = 1, color = "grey30"),
        plot.tag.position = c(0.97,0.93)) +
  labs(x = "",
       y = "Relative Tail Length",
       tag = TL_deets)
ggsave("results/figures/rel_TL.pdf", plot = rel_TL_plot, height = 2.75, width = 2.8, units = "in")


# plot ear length

rel_EL_plot <- age_ctrl_data %>%
  ggplot(aes(x=Population, y=rel_EL)) +
  geom_boxplot(aes(fill = Population), position = position_dodge(width=0.25), alpha=0.8, size = 0.5, outlier.shape = NA, show.legend = FALSE) +
  geom_jitter(position=position_jitter(width = 0.15, seed=19910118), alpha = 0.15, show.legend = FALSE) +
  scale_fill_manual(values=c("#E09832", "#1683B5"),
                    breaks=c("BRAZIL", "NEW_YORK"),
                    labels=c("Brazil", "New York")) +
  scale_x_discrete(breaks=c("BRAZIL", "NEW_YORK"),
                   labels=c("Brazil", "New York")) +
  scale_y_continuous(breaks = seq(from=0.75, to=1.6, by=0.25), labels = seq(from=0.75, to=1.6, by=0.25), limits = c(0.75,1.6)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(vjust = 1, size = 13, family = "sans"),
        axis.text = element_text(size = 12.5, family = "sans"),
        axis.text.x = element_text(color = "black"),
        plot.tag = element_markdown(family = "sans", size = 10.5, face = "italic", hjust = 1, color = "grey30"),
        plot.tag.position = c(0.97,0.93)) +
  labs(x = "",
       y = "Relative Ear Length",
       tag = EL_deets)
ggsave("results/figures/rel_EL.pdf", plot = rel_EL_plot, height = 2.75, width = 2.8, units = "in")

