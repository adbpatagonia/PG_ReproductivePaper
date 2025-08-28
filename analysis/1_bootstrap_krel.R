# This script takes the dataset used to calculate krel
# created in script 0_seal_biologicalrate.R
# and bootstraps the krel to obtain CIs
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
dat.mod.krel <- dat.biolrates[Code.Female.Maturity %in% c(8, 2, 4, 1)] %>%
  mutate(cyear =  as.factor(cohortyear)) %>%
    select(ID.Sex, cohortyear, cyear, Code.Female.Maturity,  krel) %>%
  filter(!is.na(krel))

dat.mod.krel <- dat.mod.krel[order(dat.mod.krel$ID.Sex),]

## resample ID Sex ----
ranids <- replicate(expr = dat.mod.krel[, .(ID.Sex.boot = sample(x = ID.Sex,
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
setkey(dat.mod.krel, ID.Sex)
ran.dat.mod.krel <- merge(ranidsex, dat.mod.krel, by.x = 'ID.Sex.boot', by.y = 'ID.Sex', all.x = TRUE)



# obtain quantiles ----
bootCI <- ran.dat.mod.krel %>%
  group_by(cohortyear, Code.Female.Maturity) %>%
  reframe(krelmean = mean(krel),
          krelmedian = median(krel),
          krellb = quantile(krel, probs = (1-ci/100)/2,          na.rm = TRUE),
          krelub = quantile(krel, probs = ci/100 + (1-ci/100)/2, na.rm = TRUE)) %>%
  data.table()



# output -----
fwrite(x = bootCI,
       na = NA,
       file = paste0(here::here(), "/output/RelativeCondition_bootstrap.csv"))

s2 <- Sys.time()
s2-s1
rm(list = ls())
