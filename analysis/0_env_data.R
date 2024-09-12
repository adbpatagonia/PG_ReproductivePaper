# 2024-09-10
# script to plot the environmental data for the
# Environmental effects on harp seal reproduction in the NW Atlantic paper

# libraries ----
library(data.table)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(GGally)
library(viridis)

ggplot2::theme_set(theme_light())

# functions -----
source(paste0(here::here(),"/R/HighstatLibV7.R"))
source(paste0(here::here(),"/R/standardize_x.r"))


# data ----
## NAO ----
nao <- fread(paste0(here::here(), "/data/environment/NAO/nao_station_djfm.txt"),  fill=TRUE)

nao <- nao[,1:2]
names(nao) <- c('year', 'winterNAO')




## AMO ----
### Kaplan ----
#### unsmoothed ----
amo.kaplan.unsm <- fread(paste0(here::here(),
                                "/data/environment/AMO/amon.us.long.data.txt"),
                         skip = 1,
                         nrows = 168)
names(amo.kaplan.unsm) <- c('year', 'jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')
names(amo.kaplan.unsm) <- c('year', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12')
amo.kaplan.unsm <- amo.kaplan.unsm %>%
  pivot_longer(!year, names_to = "month", values_to = "amo.unsm") %>%
  mutate(month = as.numeric(month)) %>%
  data.table()


#### smoothed ----
amo.kaplan.sm <- fread(paste0(here::here(),
                              "/data/environment/AMO/amon.sm.long.data"),
                       skip = 1,
                       nrows = 168)
names(amo.kaplan.sm) <- c('year', 'jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec')
names(amo.kaplan.sm) <- c('year', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12')
amo.kaplan.sm <- amo.kaplan.sm %>%
  pivot_longer(!year, names_to = "month", values_to = "amo.sm") %>%
  mutate(month = as.numeric(month)) %>%
  data.table()

#### merge Kaplan ----
amo.kaplan <- merge(amo.kaplan.unsm, amo.kaplan.sm)
amo.kaplan[, tim := year + month/12]
amo.kaplan[amo.sm == -99.990, amo.sm := NA]
amo.kaplan[amo.unsm == -99.990, amo.unsm := NA]

# let's make sure we can reproduce the smoothed time series
amo.kaplan[, amo.sm.calc := frollmean(amo.unsm,
                                      n = 121,
                                      fill=NA,
                                      algo=c("exact"),
                                      align=c("center"),
                                      na.rm = FALSE,
                                      hasNA = TRUE, adaptive=FALSE)]

# lines are completely superimposed - we have reproduced the smoothed time series
ggplot(amo.kaplan, aes(tim, amo.unsm)) +
  geom_line() +
  geom_line(aes(y = amo.sm), linewidth = 2) +
  geom_line(aes(y = amo.sm.calc),
            col = 'red',
            linewidth = 2, linetype = 2) +
  ylab("AMO")

### annualize AMO ----
# follow the same strategy as for NAO and AO, i.e. get the JFM mean
annual.amo.kaplan <- amo.kaplan[month < 4] %>%
  group_by(year) %>%
  reframe(winter.amo.sm = mean(amo.sm),
          winter.amo.unsm = mean(amo.unsm)) %>%
  data.table()

annual.amo.kaplan[, winter.amo.sm.calc := frollmean(winter.amo.unsm,
                                      n = 10,
                                      fill=NA,
                                      algo=c("exact"),
                                      align=c("center"),
                                      na.rm = FALSE,
                                      hasNA = TRUE, adaptive=FALSE)]


# lines are completely superimposed - we have reproduced the smoothed time series
ggplot(annual.amo.kaplan, aes(year, winter.amo.unsm)) +
  geom_line() +
  geom_line(aes(y = winter.amo.sm), linewidth = 2) +
  geom_line(aes(y = winter.amo.sm.calc),
            col = 'red',
            linewidth = 2, linetype = 2) +
  ylab("AMO")

### NOAA ----
amo.noaa <- fread(paste0(here::here(),
                         "/data/environment/AMO/ersst.v5.amo.dat.txt"),
                  skip = 1)

names(amo.noaa) <- c('year', 'month', 'ssta')
amo.noaa[, tim := year + month/12]

amo.noaa[, ssta.sm := frollmean(ssta,
                                n = 121,
                                fill=NA,
                                algo=c("exact"),
                                align=c("center"),
                                na.rm = FALSE,
                                hasNA = FALSE, adaptive=FALSE)]


## AO ----
ao <- fread(paste0(here::here(), "/data/environment/AO/monthly.ao.index.b50.current.ascii.txt"))

names(ao) <- c('year', 'month', 'AO')
ao[, tim := year + month/12]

### smoothed ----
# using the same smoother as for AMO
ao[, ao.sm := frollmean(AO,
                          n = 121,
                          fill=NA,
                          algo=c("exact"),
                          align=c("center"),
                          na.rm = FALSE,
                          hasNA = FALSE,
                          adaptive=FALSE)]

### seasonal mean ----
# Zhang et al 2021 present the JFM mean
ao.seasonal <- ao %>%
  filter(month < 4) %>%
  group_by(year) %>%
  reframe(ao.seasonal = mean(AO)) %>%
  data.table()

### annual mean ----
# Zhang et al 2021 present the JFM mean
ao.annual <- ao %>%
  # filter(month < 4) %>%
  group_by(year) %>%
  reframe(ao.annual = mean(AO)) %>%
  data.table()

ao.both <- merge(ao.seasonal, ao.annual) %>%
  pivot_longer(!year,  names_to = "index", values_to = "ao") %>%
  group_by(index) %>%
  mutate(ao_std = standardize_x(ao))


## ice area cover ----
ice <- fread(paste0(here::here(), "/data/environment/IceCoverage/ecoast_sdtt_1969_2024_0129_0129.csv"),
             skip = 9)

names(ice) <- c('week', 'status', 'perc_interpolated', 'total_concentration', 'no_data',
                'old_ice', 'first_year_ice', 'young_ice', 'new_ice', 'average_concentration',
                'median_concentration')
ice <- ice[1:56]

ice[, year := as.numeric(substr(week, start = 1, stop = 4))]

## NLCI ----
nlci <- fread(paste0(here::here(), "/data/environment/NLCI/version_2023/NL_climate_index.csv"))
names(nlci) <- c('year', 'NLCI')


## merge ----
env.dat <- merge(ao.seasonal, nlci,
                 by = 'year')

env.dat <- merge(env.dat, ice[,.(year, first_year_ice)],
                 by = 'year')
env.dat <- merge(env.dat, nao,
                 by = 'year')
env.dat <- merge(env.dat, annual.amo.kaplan[,.(year, winter.amo.sm)],
                 by = 'year')

# correlations -----
vifs <-  corvif(env.dat %>% select(-year)) %>%
  as.data.table(keep.rownames = TRUE) %>%
  arrange(GVIF)
names(vifs) <- c('variable', 'VIF')

vifs_nonlci <-  corvif(env.dat %>% select(-year, -NLCI)) %>%
  as.data.table(keep.rownames = TRUE) %>%
  arrange(GVIF)
names(vifs_nonlci) <- c('variable', 'VIF')

vifs.reduced <-  corvif(env.dat %>% select(-year, -NLCI, -winterNAO)) %>%
  as.data.table(keep.rownames = TRUE) %>%
  arrange(GVIF)
names(vifs.reduced) <- c('variable', 'VIF')

# plots ------
p.ice <- ggplot(ice, aes(year, first_year_ice*100)) + geom_line() + ylab("Percent Ice coverage in January 29")
p.nao <- ggplot(nao[year > 1949], aes(year, winterNAO)) + geom_line()
p.nlci <- ggplot(nlci, aes(year, NLCI)) + geom_line()
p.amo.kaplan <-
  ggplot(amo.kaplan, aes(tim, amo.unsm)) +
  geom_line() +
  geom_line(aes(y = amo.sm), linewidth = 2) +
  ylab("Kaplan AMO")

p.amo.noaa <-
  ggplot(amo.noaa, aes(tim, ssta)) +
  geom_line() +
  geom_line(aes(y = ssta.sm), linewidth = 2) +
  ylab("NOAA AMO")

p.amo.kaplan.annual <-
  ggplot(annual.amo.kaplan, aes(year, winter.amo.unsm)) +
  geom_line() +
  geom_line(aes(y = winter.amo.sm), linewidth = 2) +
  geom_line(aes(y = winter.amo.sm.calc),
            col = 'red',
            linewidth = 2, linetype = 1) +
  ylab("Kaplan Winter (JFM) AMO")

p.ao <- ggplot(ao, aes(tim, AO)) + geom_line()
p.ao.sm <- ggplot(ao, aes(tim, ao.sm)) +
  geom_line() +
  ylab("121-month smoothed AO")

p.ao.seasonal <- ggplot(ao.seasonal[year > 1969], aes(year, ao.seasonal)) +
  geom_line() +
  ylab("Mean AOI from January to March")

p.ao.seasonal.bars <- ggplot(ao.annual[year > 1981], aes(year, ao.annual)) +
  geom_bar(stat = "identity") +
  ylab("Annual Mean AO")

p.ao.compare <- ggplot(ao.both, aes(year, ao_std, color = index)) +
  geom_line() +
  ylab("Standardized AOI") +
  scale_colour_manual(values =  plasma(5)[c(2,4)]) +
  scale_x_continuous(breaks = seq(1950, 2030, 10), minor_breaks = seq(1955, 2025, 10)) +
  theme(legend.title = element_blank(),
        legend.position = 'top')
p.env.corrs <- ggpairs(env.dat %>% select(-year))
