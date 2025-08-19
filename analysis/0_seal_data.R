# 2025-08-19
# read population numbers and biological rates for the
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

## biological rates -----
biolrates <- fread( file = paste0(here::here(), "/data/seal/BiologicalRates.csv"))

# wrangle data ----
## population numbers -----
seal.pop <- seal.pop %>%
  rename(cohortyear = Year)
## merge ----
seal.data <- left_join(seal.pop, biolrates, by = 'cohortyear')


# plots -----
## population numbers -----
p.seal.pop <- ggplot(seal.pop, aes(x = cohortyear, y = N)) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.4) +
  geom_line() +
  xlab("Year") +
  scale_y_continuous(labels = scales::comma) +
  NULL


ggplot(seal.data, aes(x = cohortyear, y = abrate)) +
  # geom_smooth(span = 0.3) +
  geom_point() +
  geom_line(lty=2)

ggplot(seal.data, aes(x = cohortyear, y = pregrate)) +
  # geom_smooth(span = 0.3) +
  geom_point() +
  geom_line(lty=2) +
  NULL
