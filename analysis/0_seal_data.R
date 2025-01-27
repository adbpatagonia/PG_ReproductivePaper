# 2025-01-27
# script to plot the seal data for the
# Environmental effects on harp seal reproduction in the NW Atlantic paper

# libraries ----
library(data.table)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(GGally)
library(viridis)
library(xlsx)
library(scales)

ggplot2::theme_set(theme_light())

# data -----
## population numbers -----

seal.pop <- read.xlsx(paste0(here::here(), "/data/seal/Population_Size_1952-2024.xlsx"), sheetName = "Population_Size_1952-2024")

p.seal.pop <- ggplot(seal.pop, aes(x = Year, y = N)) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.4) +
  geom_line() +
  scale_y_continuous(labels = scales::comma) +
  NULL
