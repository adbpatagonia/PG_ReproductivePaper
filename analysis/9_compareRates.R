# functions -----
source(paste0(here::here(), "/R/not_in.r"))

# data -----
source( paste0(here::here(), "/analysis/0_seal_data.R"))

# read ICES data

compare.rates.ICES <-
  bind_cols(
    fread('assets/Stenson_etal_ICES.csv') %>%
      rename(cohortyear = cohyear) %>%
      select(cohortyear,
             n.mature.ICES = mature,
             n.pregnant.ICES = pregnant,
             fecrate, abrate.ICES = abrate,
             n.EP.ICES = EP) %>%
      filter(!is.na(n.mature.ICES)) %>%
      mutate(n.EP.ICES = ifelse(is.na(n.EP.ICES), 0, n.EP.ICES)) ,
    biolrates %>%
      select(cohortyear,
             n.mature.now = n.mature,
             n.pregnant.now = n.pregnant,
             pregrate,abrate.now = abrate,
             n.EP.now = n.EP) %>%
      filter(!is.na(n.mature.now)) %>%
      filter(cohortyear > 1978) %>%
      filter(cohortyear < 2015) %>%
      select(-cohortyear)) %>%
  select(cohortyear,
         n.mature.ICES, n.mature.now,
         n.pregnant.ICES, n.pregnant.now, fecrate, pregrate,
         n.EP.ICES, n.EP.now,
         abrate.ICES, abrate.now) %>%
  mutate(diff.mature = n.mature.ICES -  n.mature.now,
         diff.pregnant = n.pregnant.ICES - n.pregnant.now,
         diff.EP = n.EP.ICES - n.EP.now) %>%
  data.table()


# pregnancy rate = No. of pregnant females/No. of mature females
#  No. of pregnant females = pregrate * No. of mature females

# abortion rate = No. of abortions/(No. of abortions + No. of viable pregnancies)
# No. of abortions = (abortion rate * No. of viable pregnancies) / (1 - No. of abortions)
compare.rates <-
  bind_cols(
    fread('assets/fecdata_Stenson_etal_2020.csv') %>%
      na.omit() %>%
      mutate(n.pregnant.resdoc = as.integer(n.mature.resdoc * pregrate.resdoc)) %>%
      mutate(n.EP.resdoc = as.integer((abrate.resdoc * n.pregnant.resdoc)/(1 - abrate.resdoc))),
    biolrates %>%
      select(cohortyear,
             n.mature.now = n.mature,
             n.pregnant.now = n.pregnant,
             pregrate.now= pregrate,abrate.now = abrate,
             n.EP.now = n.EP) %>%
      filter(!is.na(n.mature.now)) %>%
      filter(cohortyear > 1953) %>%
      filter(cohortyear < 2020) %>%
      select(-cohortyear)) %>%
  select(cohortyear,
         n.mature.resdoc, n.mature.now,
         n.pregnant.resdoc, n.pregnant.now, pregrate.resdoc, pregrate.now,
         n.EP.resdoc, n.EP.now,
         abrate.resdoc, abrate.now) %>%
  data.table()


compare.rates <- compare.rates[cohortyear < 1979,
              (c('n.pregnant.now',
                 'n.pregnant.resdoc',
                 'n.mature.resdoc',
                 'n.mature.now',
                 'n.EP.now',
                 'n.EP.resdoc')) := NA] %>%

  mutate(diff.mature = n.mature.resdoc -  n.mature.now,
         diff.pregnant = n.pregnant.resdoc - n.pregnant.now,
         diff.EP = n.EP.resdoc - n.EP.now) %>%
  fwrite(., file = 'output/CompareDatasets.csv')

# remove years when everything is identical
yr_remove <- pull(compare.rates[diff.mature == 0 & diff.pregnant == 0 & diff.EP == 0, .(cohortyear)])
compare.rates <- compare.rates[cohortyear %!in% yr_remove]

# remove old data
compare.rates <- compare.rates[cohortyear > 1978]
# compare datasets ----
dat.biolrates <- fread(paste0(here::here(), "/data/seal/CleanDatasetForBiologicalRates.csv"))
dat.ICES <- fread(paste0(here::here(), "/assets/fecdata_Stenson_etal_2015.csv"),
                  skip = 3) %>%
  rename(cohortyear = cohyear,
         EP = EarlyPupper) %>%
  select(-V14, -V15, -V16, -V17, -V18)


## EP ----
ids_ep_ices <- paste0(pull(dat.ICES[EP == 1, .(ID)] ),"F")
ids_nonep_ices <- paste0(pull(dat.ICES[EP == 0, .(ID)] ),"F")


ids_ep_now <- pull(dat.biolrates[EP == 1, .(ID.Sex)] )
ids_ep_now <- substr(ids_ep_now, 1, 8)
ids_nonep_now <- pull(dat.biolrates[EP == 1, .(ID.Sex)] )
ids_nonep_now <- substr(ids_nonep_now, 1, 8)

dat.biolrates[ID.Sex %in% ids_ep_ices & EP == 0,.(Month, Day)] %>%
  unique() %>%
  arrange_all()

dat.biolrates[ID.Sex %in% ids_nonep_ices & EP == 1]

dat.ICES[ID %in% ids_nonep_now & EP == 1]
dat.ICES[ID %in% ids_ep_now & EP == 0]




dat.ICES[ID.Sex %in% ids_ep_ices & Month == 1 & Day == 28]


compare.rates %>%
  mutate(diff.pregrate = pregrate.resdoc - pregrate.now) %>%
  filter(abs(diff.pregrate )> 0.05)

compare.rates %>%
  mutate(diff.abrate = abrate.resdoc - abrate.now) %>%
  filter(abs(diff.abrate )> 0.05)
