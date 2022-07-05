#!/usr/bin/env Rscript --vanilla

##############################################################

# Author: Mallory A. Ballinger

# This script extracts climate data (mean temperature) from WorldClim dataset.
# This script generates Figure 1A in BallingerMack_2022.

##############################################################
# Required packages
##############################################################

#rm(list = ls()) # clears R's environment
library(tidyverse)
library(raster)
library(sf)
library(rgdal)
library(ggtext)

##############################################################
# Download data
##############################################################

# download monthly average temperature data
tmean_data <- getData(name = "worldclim", var = "tmean", res = 2.5)

# convert temperature values to Celsius
gain(tmean_data) <- 0.1

# calculate mean of the monthly average temperatures
tmean_mean <- mean(tmean_data)

# crop where we want data from
c1 = crop(tmean_mean, extent(-130, -30, -40, 55))

# transform into data frame
c1_df = as.data.frame(c1, xy = TRUE, na.rm = TRUE)

#asp <- tmaptools::get_asp_ratio(c1) # returns 1.043626

##############################################################
# Plot
##############################################################

AMtemp_plot <-
  ggplot() +
  geom_hline(yintercept = 0, linetype = 1, color = "grey95") +
  geom_raster(data = c1_df,
              aes(x=x, y=y, fill=layer)) +
  coord_quickmap(xlim = c(-130,-25),
                 ylim = c(-40,55),
                 clip = "off",
                 expand = FALSE) +
  scale_x_continuous(limits = c(-130,-30),
                     breaks = seq(-130, -30, 20),
                     labels = c("130\u00B0W", "110\u00B0W", "90\u00B0W", "70\u00B0W", "50\u00B0W", "30\u00B0W")) +
  scale_y_continuous(limits = c(-40,55),
                     breaks = seq(-40, 55, 20),
                     labels = c("40\u00B0S", "20\u00B0S", "0\u00B0", "20\u00B0N", "40\u00B0N")) +
  scale_fill_viridis_c(limits = c(-10,30),
                       breaks = seq(-5, 25, 7.5),
                       labels = c("-5\u00B0C", "2.5\u00B0C", "10\u00B0C", "17.5\u00B0C", "25\u00B0C"),
                       #option = "magman",
                       ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        aspect.ratio = 1/1,
        legend.position = c(0.13,0.2),
        legend.title = element_text(size = 10, family = "sans", color = "grey30"),
        legend.spacing.x = unit(1.0, 'mm'),
        legend.spacing.y = unit(0.5, "mm"),
        legend.text = element_text(size = 10, family = "sans", color = "grey30"),
        axis.text = element_text (size = 12.5, family = "sans"),
        axis.title.y = element_text(vjust = 1, size = 14, family = "sans"),
        axis.title.x = element_text(vjust = -1, size = 14, family = "sans")) +
  labs(x = "Longitude",
       y = "Latitude",
       fill = "Annual Mean\nTemperature")

ggsave("results/figures/AMtemp_plot.pdf", plot = AMtemp_plot, height = 5.25, width = 5.25)

