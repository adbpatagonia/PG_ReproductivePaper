# This script takes the dataset used to calculate biological rate
# created in script 0_seal_biologicalrate.R
# and bootstraps the pregnancy and abortion rates to obtain CIs
# by resampling ID.Sex, within cohortyear and boot iteration
s1 <- Sys.time()
# libraries -----
library(data.table)
library(tidyverse)

# data -----
dat.biolrates <- fread(paste0(here::here(), "/data/seal/CleanDatasetForBiologicalRates.csv"))


# parameters -----
nboot <- 10000
ci <- 95

# bootstrap ----
dat.br <- dat.biolrates[, .(ID.Sex, cohortyear, maturity, pregnancy, EP)]
dat.br <- dat.br[order(dat.br$ID.Sex),]

## resample ID Sex ----
ranids <- replicate(expr = dat.br[, .(ID.Sex.boot = sample(x = ID.Sex,
                                                           replace = TRUE,
                                                           prob = NULL)),
                                  keyby = .(cohortyear)],
                    n = nboot,
                    simplify = FALSE)

# transform to data.table
ranidsex <- rbindlist(ranids, idcol = "nboot")
ranidsex[, cohortyear := NULL]

## get maturity and pregnancy info for resampled ID.Sex ----
setkey(ranidsex, ID.Sex.boot)
setkey(dat.br, ID.Sex)
ran.dat.br <- merge(ranidsex, dat.br, by.x = 'ID.Sex.boot', by.y = 'ID.Sex', all.x = TRUE)


## obtain totals for maturity, pregnancy and EP of resampled data ----
bootmat <- ran.dat.br[, .N, by = c('cohortyear', 'nboot','maturity')] %>%
  filter(maturity == 1) %>%
  select(-maturity) %>%
  dplyr::rename(n.mat = N)
bootpreg <- ran.dat.br[, .N, by = c('cohortyear', 'nboot','pregnancy')] %>%
  filter(pregnancy == 1) %>%
  select(-pregnancy) %>%
  dplyr::rename(n.preg = N)
bootep <- ran.dat.br[, .N, by = c('cohortyear', 'nboot','EP')] %>%
  filter(EP == 1) %>%
  select(-EP) %>%
  dplyr::rename(n.ep = N) %>%
  right_join(
    expand.grid(cohortyear = unique(dat.br$cohortyear),
                nboot = 1:nboot)
  ) %>%
  mutate(n.ep = ifelse(is.na(n.ep), 0, n.ep)) %>%
  arrange(cohortyear, nboot)

bootdat <- merge(bootmat, bootpreg) %>%
  left_join(bootep) %>%
  data.table

## calculate pregnancy rate ----
# pregnancy rate = No. of pregnant females/No. of mature females
bootdat[, pregrate := n.preg/n.mat]

##  calculate abortion rate ----
bootdat[, totpreg := n.preg + n.ep]
bootdat[, abrate := n.ep/totpreg]

# obtain quantiles ----
bootCI <- bootdat %>%
  group_by(cohortyear) %>%
  reframe(preglb = quantile(pregrate, probs = (1-ci/100)/2,          na.rm = TRUE),
          pregub = quantile(pregrate, probs = ci/100 + (1-ci/100)/2, na.rm = TRUE),
          ablb =   quantile(abrate,   probs = (1-ci/100)/2,          na.rm = TRUE),
          abub =   quantile(abrate,   probs = ci/100 + (1-ci/100)/2, na.rm = TRUE)) %>%
  data.table()

# output -----
fwrite(x = bootCI,
       na = NA,
       file = paste0(here::here(), "/output/BiologicalRates_bootstrap.csv"))

s2 <- Sys.time()
s2-s1
rm(list = ls())
