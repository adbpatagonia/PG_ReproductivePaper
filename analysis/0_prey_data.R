# 205-06-25
# script to plot the prey data for the
# Environmental effects on harp seal reproduction in the NW Atlantic paper

# data sent by MKA to GBS
# email in data/prey/Gmail - Fwd_ fish data.pdf
# data in "data/prey/Data_for_Garry.xlsx"

# libraries ----
library(readxl)
library(janitor)
library(data.table)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(GGally)
library(viridis)

ggplot2::theme_set(theme_light())

# read data -----
prey_fall <- read_excel(path = paste0(here::here(), "/data/prey/Data_for_Garry.xlsx"), sheet = "bmass & abund 2J3KL Fall RV ") %>%
  janitor::clean_names() %>%
  data.table()


capelin_spring <- read_excel(path = paste0(here::here(), "/data/prey/Data_for_Garry.xlsx"), sheet = "3L Spring Acoustic Survey") %>%
  janitor::clean_names() %>%
  data.table()

# wrangle data -----
# remove assigned or interpolated capelin biomass values
capelin_spring[!is.na(x3), capelin_kt := NA ]

capelin_spring[, log_capelin := log(capelin_kt)]
prey_fall[, log_biomass := log(biomass_tonnes)]


# plots -----
p.prey.fall <- ggplot(prey_fall, aes(x = year, y = biomass_tonnes)) +
  geom_point() +
  geom_line(lty = 2) +
  facet_wrap(species ~. , scales = "free_y")

p.prey.fall.log <- ggplot(prey_fall, aes(x = year, y = log_biomass)) +
  geom_point() +
  geom_line(lty = 2) +
  facet_wrap(species ~. , scales = "free_y") +
  ylab("ln(biomass (tonnes))")

p.prey.spring <- ggplot(capelin_spring, aes(x = year, y = capelin_kt)) +
  geom_point() +
  geom_line(lty = 2)

p.prey.spring.log <- ggplot(capelin_spring, aes(x = year, y = log_capelin)) +
  geom_point() +
  geom_line(lty = 2) +
  ylab("ln(Capelin (kilotonnes))")


