# 2025-01-27
# script to plot the seal data for the
# Environmental effects on harp seal reproduction in the NW Atlantic paper


# Notes on Code.Female.Maturity
# Teams chat with Shelley

# From Stenson et al (ICES)
# Thus, we assumed that females with a CA diameter greater than 12 mm represent a
# pregnancy from the current year that has been terminated within the previous 30 d.
# Similarly, seals with an active pregnancy have CLs measuring .13 mm across.
# Therefore, seals that lacked a developing foetus, but had a CL of ≥13 mm or CA ≥12 mm,
# a rugose uterus and a large uterine horn, were assumed to have pupped recently
# and of those, we assumed early puppers if they were collected before feb 20

# Code.Female.Maturity 8 meet the criteria set out for CA and uterus
# Code.Female.Maturity 18 meet some criteria but not others OR are close but not quite
# Therefore:
# 1. Do not consider Code.Female.Maturity 18 as abortions
# 2. Consider Code.Female.Maturity 8 + collection before Feb 20 = abortion


# libraries ----
library(data.table)
library(ggplot2)
library(plotly)
library(tidyverse)
library(lubridate)
library(GGally)
library(viridis)
library(openxlsx)
library(scales)
library(gghalves)
library(ggbeeswarm)
library(kableExtra)

ggplot2::theme_set(theme_light())

# functions -----
source(paste0(here::here(), "/R/mround.R"))
source(paste0(here::here(), "/R/plot_lw.R"))
source(paste0(here::here(), "/R/krel.R"))
source(paste0(here::here(), "/R/not_in.r"))

# read data -----
## population numbers -----
seal.pop <- openxlsx::read.xlsx(xlsxFile = paste0(here::here(), "/data/seal/Population_Size_1952-2024.xlsx"),
                                sheet = "Population_Size_1952-2024") %>%
  data.table()

## ovary ----
ovary <- openxlsx::read.xlsx(xlsxFile = paste0(here::here(), "/data-raw/seal/qry_HarpOvaryDataForABuren_July2025.xlsx"),
                             sheet = "qry_HarpOvaryDataForABuren_July") %>%
  data.table()

## morphagerepro ----
morphagerepro <- openxlsx::read.xlsx(xlsxFile = paste0(here::here(), "/data-raw/seal/qry_CollectionMorphAgeReproDataforABuren_July2025_v2.xlsx"),
                                     sheet = "qry_CollectionMorphAgeReproData") %>%
  data.table()

## biol rates old ----
# we don't have individual data for these years
# these values were taken from the file assets/fecdata.xlsx
# the file was used to produce the estimates for the Stenson et al ICES paper - stored here: "D:\Buren_files\MEGA\papersAle\harps_fecundity\data\fecdata.xlsx"
biolrates.old <- read.csv(paste0(here::here(), "/data/seal/BiologicalRates_1954-1978.csv")) %>%
  data.table()

# parameters ----
# early pupper dates ----
last.ep.date <- 51

# wrangle data -----
## biol rates old ----
# replace NA in n.EP to zero
biolrates.old$n.EP <- as.integer(biolrates.old$n.EP)
biolrates.old[is.na(n.EP), n.EP := 0]

## morphagerepro ----
### cohort year ----
# From Shelley Lang:  base your cohort year on the actual collection year
# with the adjustment for month of year
# Cohort year is Sept. 1 to Aug. 31
morphagerepro[, cohortyear := as.integer(substr(ID.Sex, 1, 4))]
morphagerepro[, cohortyear := Year]
morphagerepro[Month > 8 & Month <13, cohortyear := Year + 1]
# there is one seal 20130096F with 9999 as month -
# Season in the DB shows up as W - therefore leave as is


morphagerepro$date <- (with(morphagerepro, as.Date(paste(formatC(Day, width=2, flag="0"),'/',formatC(Month, width=2, flag="0"),'/',Year,sep=''),format='%d/%m/%Y')))
morphagerepro$doy <- yday(morphagerepro$date)
### remove dead on beach ----
# Screen out ID Sex 20180054F (MarineMammal ID 64849) – seal was found dead on a beach
morphagerepro <- morphagerepro[ID.Sex != "20180054F"]

### Age ----
morphagerepro[Final.Cohort.Age == 98]
morphagerepro[Final.Cohort.Age == 98, Final.Cohort.Age := 99]

### remove fetal weight ----
# remove the weight of the foetus for pregnant females
# I obtained monthly mean foetal weights from the PG_Growth work
# Ask Shelley to send updated values
# ADB and SL looked at potential outliers of foetus weights. script 9_PotentialFoetusOutliers.r
# We came to these 3 conclusions:
# 1. ID.Sex = 20110194M. This should be dropped
# 2. cohort years 2013 and 2014. The weight for the foetuses are expressed in lbs. Transfor them
# Note that the weights of the larger seals for those cohort years are expressed in kgs
# 3. ID.Sex = 20169906M. This is the heaviest foetus collected in Dec, but "but it's also a long seal (male) so I think it might be fine (data matches the original entry)"
# Will keep as is
fwd <- 2.9
fwj <- 6.3
fwf <- 8.3

# drop 2011 foetus
morphagerepro <- morphagerepro[ID.Sex != "20110194M"]

# transform 2013 & 2014 foetus weight to kg
morphagerepro[Final.Cohort.Age == 91 & cohortyear %in% c(2013,2014), Body.Weight := Body.Weight * 0.453592]

fwd <- mean(morphagerepro[Final.Cohort.Age == 91 & Month == 12, Body.Weight])
fwj <- mean(morphagerepro[Final.Cohort.Age == 91 & Month == 1,  Body.Weight])
fwf <- mean(morphagerepro[Final.Cohort.Age == 91 & Month == 2,  Body.Weight])
morphagerepro$Body.Weight <- ifelse(morphagerepro$Code.Female.Maturity %in% 2 & morphagerepro$Month == 12,
                                    morphagerepro$Body.Weight - fwd,
                                    morphagerepro$Body.Weight)
morphagerepro$Body.Weight <- ifelse(morphagerepro$Code.Female.Maturity %in% 2 & morphagerepro$Month == 1,
                                    morphagerepro$Body.Weight - fwj,
                                    morphagerepro$Body.Weight)
morphagerepro$Body.Weight <- ifelse(morphagerepro$Code.Female.Maturity %in% 2 & morphagerepro$Month == 2,
                                    morphagerepro$Body.Weight - fwf,
                                    morphagerepro$Body.Weight)



### remove foetus, stillborn & starvling ----
ids_pelremov <- pull(morphagerepro[Code.PelageType %in% 91:93, .(ID.Sex)])
ids_ageremov <- pull(morphagerepro[Final.Cohort.Age %in% 91:93, .(ID.Sex)])

morphagerepro <- morphagerepro[ID.Sex %!in% ids_pelremov]
morphagerepro <- morphagerepro[ID.Sex %!in% ids_ageremov]

ovary <- ovary[ID.Sex %!in% ids_pelremov]
ovary <- ovary[ID.Sex %!in% ids_ageremov]

### remove age zero (YOY) ----
ids_yoy <-  pull(morphagerepro[Final.Cohort.Age == 0, .(ID.Sex)])
morphagerepro <- morphagerepro[ID.Sex %!in% ids_yoy]
ovary <- ovary[ID.Sex %!in% ids_yoy]

# age = 80: incomplete tooth but seal is at least 20 years old
morphagerepro[Final.Cohort.Age == 80]

## no duplicates
which(duplicated(morphagerepro$ID.Sex))
morphagerepro[anyDuplicated(morphagerepro$ID.Sex)]

### length -----
morphagerepro$Body.Length %>% max(na.rm = TRUE)
morphagerepro[Body.Length > 200]
morphagerepro$Body.Length %>% min(na.rm = TRUE)
morphagerepro[Body.Length ==31]


# from Shelley
# assets/DataDecisions_ADB_2025-08_SL.docx
morphagerepro[ID.Sex %in% c("20181046F", "20113315F", "20192983F"), Body.Length  := NA]
morphagerepro[ID.Sex == "20224484F", Body.Length  := 174]

# from P Goulet
# 20161235M the girth and length have been inverted, should be girth=130  and length=163

morphagerepro[ID.Sex == "20161235M"]
morphagerepro[ID.Sex == "20161235M", Body.Length  := 163]
morphagerepro[ID.Sex == "20161235M",  Maximum.Girth := 130]


### check inconsistencies with previous instructions from Rachel -----
# edits included in file README_Buren.txt, commit a415c479c87771fccf382a41c95f1ec712a1e696
# - remove for now. needs to be checked against datasheets.
# 2025: update from Shelley L. THe first 4 are excluded in the flags
# Exclude.From.Age  and Eclude.From.Repro
# the last seal should not be included
removefornow <- c('20161194F', '20161214F', '20161251F', '20171691F', '20181055F')
# insert into Do not use table and do not use for anything
remove <- c('20113802F', '20130095F')
#  change to age 99
changeage99 <- c('20100066F','20123810F','20140236F','20150057F','20103066F')
# change to age 1
changeage1 <- c('20130046F', '20130056F', '20161352F')
# change to age 20
changeage20 <- ('20161096F')

### changeage99
morphagerepro[ID.Sex %in% changeage99,.(ID.Sex, Final.Cohort.Age)]

### changeage1
# these should be immature
morphagerepro[ID.Sex %in% changeage1  ,.(ID.Sex, Code.Female.Maturity)]
# no data for these in tbl_ovary
ovary[ID.Sex %in% morphagerepro[ID.Sex %in% changeage1,.(ID.Sex)]]

# no morph data and no ovary data - I will change code.female.maturity to 1, until I am told otherwise
ids99 <- morphagerepro[ID.Sex %in% changeage99 & Code.Female.Maturity > 1, .(ID.Sex)] %>% pull()
morphagerepro[ID.Sex %in% ids99, Code.Female.Maturity := 1]
morphagerepro[ID.Sex %in% ids99, Female.Maturity := "Immature"]

### changeage20
morphagerepro[ID.Sex %in% changeage20 ,.(ID.Sex, Code.Female.Maturity, Female.Maturity)]
# no data for these in tbl_ovary
ovary[ID.Sex %in% morphagerepro[ID.Sex %in% changeage20 ,.(ID.Sex)]]

### remove
morphagerepro[ID.Sex %in% remove,.(ID.Sex, Code.Female.Maturity, Female.Maturity)]
ovary[ID.Sex %in% remove]

# I will remove until I am told otherwise
ovary <- ovary[!ID.Sex %in% remove]
morphagerepro <- morphagerepro[!ID.Sex %in% remove]


### exclude from repro -----
ids_exc_repro <- morphagerepro[Exclude.From.Repro == TRUE, .(ID.Sex)] %>% pull()

# 20224400F – Exclude from ALL analyses
morphagerepro <- morphagerepro[ID.Sex != "20224400F"]

# set female maturity to NA
morphagerepro[Exclude.From.Repro == TRUE, Code.Female.Maturity := NA]
morphagerepro[Exclude.From.Repro == TRUE, Female.Maturity := NA]

### exclude from age -----
morphagerepro[Exclude.From.Age == TRUE]
# set Final.Cohort.Age to NA
morphagerepro[Exclude.From.Age == TRUE, Final.Cohort.Age := NA]



### Age ----
# ages from Bonnie S in emails to GBS - Dec 1 and 4, 2023
# age[idsex == "20171709F", age := 25]
# age[idsex == "20082059F", age := 99]
# age[idsex == "20103001M", age := 27]
# These have NA in Final.Cohort.Age.
# Decision: I will use age from Bonnie until I am told otherwise
morphagerepro[ID.Sex == "20171709F", Final.Cohort.Age := 25]
morphagerepro[ID.Sex == "20082059F", Final.Cohort.Age := 99]
morphagerepro[ID.Sex == "20103001M", Final.Cohort.Age := 27]

### weight ----
# this is from the work on growth
# at least 3 seals have the wrong weight,
ids_w <-  c("19940928F", "20042978F", "20062138F")
# I received feedback from Shelley Lang - different for each seal
# assets/DataDecisions_ADB_2025-08_SL.docx
# morphagerepro[ID.Sex %in% ids_w, Body.Weight := NA]
morphagerepro[ID.Sex == "19940928F", Final.Cohort.Age := 91]
morphagerepro[ID.Sex == "20042978F", Final.Cohort.Age := 0]
morphagerepro[ID.Sex == "20062138F", Final.Cohort.Age := 0]

morphagerepro[ID.Sex == "20062110M", Final.Cohort.Age := 0]
morphagerepro[ID.Sex == "20062167M", Final.Cohort.Age := 0]

# from Shelley
# assets/DataDecisions_ADB_2025-08_SL.docx
morphagerepro[ID.Sex %in% c( "20192983F", "20215080F"  ), Body.Weight  := NA]
morphagerepro[ID.Sex == "20183420M", Body.Weight  := 152]


### Senescence ----
# OK
femcode6 <- c('20091859F','20111392F','20130006F','20140355F','20170015F','20171718F')
morphagerepro[ID.Sex %in% femcode6]

### Remove OU ----
# OK - variable Exclude.from.repro == TRUE
removeOU <- c('20130099F','20160009F','20171726F')
morphagerepro[ID.Sex %in% removeOU]




### female maturity -----
#### NA Female Maturity -----
# large number of females with NA
morphagerepro[is.na(Code.Female.Maturity) & Sex == "F"]

ids_narepro <- morphagerepro[is.na(Code.Female.Maturity) & Sex == "F", .(ID.Sex)] %>% pull()
ovary[ID.Sex %in% ids_narepro,. (ID.Sex, Code.Female.Maturity, Female.Maturity)]

#### potential inconsistencies given age -----
# 2025-08-22: I received feedback from Shelley L
# assets/DataDecisions_ADB_2025-08_SL_V2.docx
# I implemented those decisions below
morphagerepro[Final.Cohort.Age < 4 & Code.Female.Maturity > 1,
              .(ID.Sex, Final.Cohort.Age, Code.Female.Maturity, Female.Maturity, Code.PelageType,
                Body.Length, Body.Weight)] %>%
  kable(.) %>%
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE)

ids_young <- pull(morphagerepro[Final.Cohort.Age < 4 & Code.Female.Maturity > 1,
                                .(ID.Sex)] )

ovary[ID.Sex %in% ids_young,
      .(ID.Sex, Alb.Ov.No.of.CL, Alb.Ov.No.of.CA,
        Lut.Ov.No.of.corpora.lutea, Lut.Ov.No.of.CA,
        Implanted.Embryo, Code.Female.Maturity)] %>%
  kable(.) %>%
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE)

# 19790413F OK. Code 3
morphagerepro[ID.Sex == '19790413F', .(Code.Female.Maturity)];ovary[ID.Sex == '19790413F', .(Code.Female.Maturity)]

# 19791108F . Exclude from repro
ovary <- ovary[ID.Sex != '19791108F']
morphagerepro[ID.Sex == '19791108F', Exclude.From.Repro := TRUE]
morphagerepro[ID.Sex == '19791108F', Code.Female.Maturity := NA]
morphagerepro[ID.Sex == '19791108F', Female.Maturity := NA]

# 19791269F & 19791284F . OK. Code 3
morphagerepro[ID.Sex == '19791269F', .(Code.Female.Maturity)];ovary[ID.Sex == '19791269F', .(Code.Female.Maturity)]
morphagerepro[ID.Sex == '19791284F', .(Code.Female.Maturity)];ovary[ID.Sex == '19791284F', .(Code.Female.Maturity)]

# 19792001F. It should be 8. OK. . Maturity Code was changed from 4 to 8 in 2010 (NW)
morphagerepro[ID.Sex == '19792001F', .(Code.Female.Maturity)]

# 19831148F & 19831228F & 19841901F & 19873080F  . 3 OK
morphagerepro[ID.Sex == '19831148F', .(Code.Female.Maturity)];ovary[ID.Sex == '19831148F', .(Code.Female.Maturity)]
morphagerepro[ID.Sex == '19831228F', .(Code.Female.Maturity)];ovary[ID.Sex == '19831228F', .(Code.Female.Maturity)]
morphagerepro[ID.Sex == '19841901F', .(Code.Female.Maturity)];ovary[ID.Sex == '19841901F', .(Code.Female.Maturity)]
morphagerepro[ID.Sex == '19873080F', .(Code.Female.Maturity)];ovary[ID.Sex == '19873080F', .(Code.Female.Maturity)]

# Exclude From All
# 19880346F, 19880354F, 19880427F and 19880433F.
# 19922202F
# 19932551F
# 19940841F
# 19980231F
ids_exc_all <- c('19880346F', '19880354F', '19880427F' , '19880433F',
                 '19922202F', '19932551F', '19940841F', '19980231F')
morphagerepro <- morphagerepro[ID.Sex %!in% ids_exc_all]
ovary <- ovary[ID.Sex %!in% ids_exc_all]

# 19981951F. OK. Code 3
morphagerepro[ID.Sex == '19981951F', .(Code.Female.Maturity)];ovary[ID.Sex == '19981951F', .(Code.Female.Maturity)]

# 20003503F. Exclude from Age
morphagerepro[ID.Sex == '19981951F', Final.Cohort.Age := NA]
morphagerepro[ID.Sex == '19981951F', Exclude.From.Age := TRUE]

# 20161204F . Recode Maturity as 1
morphagerepro[ID.Sex == '20161204F', Code.Female.Maturity := 1]
ovary[ID.Sex == '20161204F', Code.Female.Maturity := 1]
morphagerepro[ID.Sex == '20161204F', Female.Maturity := 'Immature']
ovary[ID.Sex == '20161204F', Female.Maturity := 'Immature']

# 20193682F. Exclude from Repro
morphagerepro[ID.Sex == '20193682F', Exclude.From.Repro := TRUE]
morphagerepro[ID.Sex == '20193682F', Code.Female.Maturity := NA]
morphagerepro[ID.Sex == '20193682F', Female.Maturity := NA]

ovary <- ovary[ID.Sex != '20193682F']

# ids_young3 <- pull(morphagerepro[Final.Cohort.Age < 3 & Code.Female.Maturity > 1,
#                                  .(ID.Sex)] )
# morphagerepro[ID.Sex %in% ids_young3, Code.Female.Maturity := 1]
# morphagerepro[ID.Sex %in% ids_young3, Female.Maturity := "Immature"]
#
# ovary[ID.Sex %in% ids_young3, Code.Female.Maturity := 1]
# ovary[ID.Sex %in% ids_young3, Female.Maturity := "Immature"]

##### Femmat unknown ----
# age = 1, thus, set as immature
ids_fem0 <- pull(morphagerepro[Final.Cohort.Age < 3 & Code.Female.Maturity == 0,
                               .(ID.Sex)] )
morphagerepro[ID.Sex %in% ids_fem0]
morphagerepro[ID.Sex %in% ids_fem0, Code.Female.Maturity := 1]
morphagerepro[ID.Sex %in% ids_fem0, Female.Maturity := "Immature"]

ovary[ID.Sex %in% ids_fem0, Code.Female.Maturity := 1]
ovary[ID.Sex %in% ids_fem0, Female.Maturity := "Immature"]

#### newborns to beaters ----
morphagerepro[is.na(Female.Maturity)]
# set Immature for newborns to beaters
# simply making sure that they are all small
morphagerepro[Code.PelageType < 7,.(Body.Weight)] %>% max(., na.rm = TRUE)
morphagerepro[Code.PelageType < 7, Female.Maturity := "Immature"]
morphagerepro[Code.PelageType < 7, Code.Female.Maturity := 1]

# these were detected using the plotly below
morphagerepro[ID.Sex %in% c('20100080F', '20060444F', '20002258F', '19791358F', '20060440F')]
# most had Code.Pelage.Type = 7, will leave as is
# this seal is extremely small and Code.Pelage.Type=99 (blank)
# change female maturity
morphagerepro[ID.Sex %in% c('20100080F'), Female.Maturity := "Immature"]
morphagerepro[ID.Sex %in% c('20100080F'), Code.Female.Maturity := 1]

# find others with blank pelage type that are potentially immature
# these look fine
morphagerepro[Code.PelageType == 99 & Code.Female.Maturity > 1]

####  consider only beater and older ----
ids_youngremove <-  pull(morphagerepro[Code.PelageType %in% c(0:5), .(ID.Sex)])
morphagerepro <- morphagerepro[ID.Sex %!in% ids_youngremove]
ovary <- ovary[ID.Sex %!in% ids_youngremove]


### Fem mat 18 ----
# Ask shelley what the 18 is, and why these do not have morph data
morphagerepro[Code.Female.Maturity==18, .(ID.Sex, Code.Female.Maturity, Female.Maturity, Body.Length, Body.Weight, date)] %>%
  arrange(date) %>%
  kable(.) %>%
  kable_styling(bootstrap_options = c("striped", "hover"),
                full_width = FALSE)

## ovary ----
### remove NA female maturity -----
# samples collected were not O&U
ovary[is.na(Code.Female.Maturity)]
ovary <- ovary[!is.na(Code.Female.Maturity)]

ovary[Code.Female.Maturity == 0]

morphagerepro[ID.Sex == '20220020F']

# 20224400F – Exclude from ALL analyses
ovary <- ovary[ID.Sex != "20224400F"]

### remove dead on beach ----
# Screen out ID Sex 20180054F (MarineMammal ID 64849) – seal was found dead on a beach
# "PG_ReproductivePaper\data\seal\ScreenOutSeal.pdf"
ovary <- ovary[ID.Sex != "20180054F"]

### exclude from repro -----
ovary <- ovary[Exclude.From.Repro == FALSE]
ovary <- ovary[!ID.Sex %in% ids_exc_repro]

## no duplicates
which(duplicated(ovary$ID.Sex))
ovary[anyDuplicated(ovary$ID.Sex)]

### Senescence ----
# OK
ovary[ID.Sex %in% femcode6]

## define maturity codes  ----
# 0=FALSE  --- 1=TRUE
unique(morphagerepro$Code.Female.Maturity)
morphagerepro$maturity <- as.integer(NA)
morphagerepro$maturity <- ifelse(morphagerepro$Code.Female.Maturity == 1, 0,
                                 ifelse(is.na(morphagerepro$Code.Female.Maturity), NA,
                                        ifelse(morphagerepro$Code.Female.Maturity > 1, 1, 10)))

# this looks OK
morphagerepro %>%
  distinct(Code.Female.Maturity, Female.Maturity, maturity) %>%
  arrange_all()

## define pregnancy codes  ------
# 0=FALSE  --- 1=TRUE

# there are 106 seals that had Code.Implanted.Embryo == 99 (blank)
# and at the same time Code.Female.Maturity == 2 (Mature; pregnant (implanted embryo))
# or Code.Female.Maturity == 3 (Mature; pregnant (delay period)). These will not make a difference in the end, as these codes are used for March-August
# I will set Code.Implanted.Embryo=1 (Present) for those seals
ovary[Code.Implanted.Embryo == 99 ] %>%
  group_by(Code.Female.Maturity, Female.Maturity) %>%
  tally()

ovary[Code.Female.Maturity == 8 ] %>%
  group_by(Code.Implanted.Embryo, Implanted.Embryo) %>%
  tally()
ovary[Code.Female.Maturity == 8 & Code.Implanted.Embryo == 99]

ids_impemb <- pull(ovary[Code.Implanted.Embryo == 99 & Code.Female.Maturity %in% c(2, 3), .(ID.Sex)])
ovary[ID.Sex %in% ids_impemb, Code.Implanted.Embryo := 1]
ovary[ID.Sex %in% ids_impemb, Implanted.Embryo := "Present"]

# set blank code implanted embryo to absent for the rest of the maturity stages,
# except codes 9 and 18
ids_nonimpemb <- pull(ovary[Code.Implanted.Embryo == 99 & Code.Female.Maturity %!in% c(9, 18), .(ID.Sex)])
ovary[ID.Sex %in% ids_nonimpemb, Code.Implanted.Embryo := 2]
ovary[ID.Sex %in% ids_nonimpemb, Implanted.Embryo := "Absent"]

ovary[Code.Implanted.Embryo == 99, .(Code.Female.Maturity)] %>% unique() %>% arrange_all()

# we need to define pregnancy stage for all seals, as these will affect pregnancy rate
# there are seals Code.Female.Maturity %!in% c(9, 18) that will have NA
# there are only 3 seals that would count to pregnancy rate
# The 2 seals collected in 1981 were coded as mature, not pregnant
# Here, I will drop them. SL to confirm apporach
right_join(
  morphagerepro[,.(ID.Sex, Year, Month)] ,
  ovary, by = "ID.Sex") %>%
  filter(Code.Implanted.Embryo == 99) %>%
  filter(Code.Female.Maturity %in% c(9, 18)) %>%
  filter(Exclude.From.Repro == FALSE) %>%
  filter(Month %in% c(10:12, 1, 2))

ids_uncpreg <- right_join(
  morphagerepro[,.(ID.Sex, Year, Month)] ,
  ovary, by = "ID.Sex") %>%
  filter(Code.Implanted.Embryo == 99) %>%
  filter(Code.Female.Maturity %in% c(9, 18)) %>%
  filter(Exclude.From.Repro == FALSE) %>%
  filter(Month %in% c(10:12, 1, 2)) %>%
  select(ID.Sex) %>%
  pull()

ovary <- ovary[ID.Sex %!in% ids_uncpreg]
morphagerepro <- morphagerepro[ID.Sex %!in% ids_uncpreg]

# this looks OK - the blank in Implanted.Embryo was caught outside the period of interest Month %in% c(10:12, 1, 2)
ovary %>%
  distinct(Code.Female.Maturity, Code.Implanted.Embryo, Implanted.Embryo) %>%
  # filter(Code.Implanted.Embryo == 99 ) %>%
  arrange_all()

# join ovary data with morphagerepro
morphagerepro <- left_join(morphagerepro,
                           ovary[, .(ID.Sex, Code.Implanted.Embryo, Implanted.Embryo)],
                           by = "ID.Sex")
morphagerepro %>%
  filter(Month %in% c(10:12, 1, 2)) %>%
  filter(Exclude.From.Repro == FALSE) %>%
  filter(Sex == "F") %>%
  distinct(Code.Female.Maturity, Female.Maturity, Code.Implanted.Embryo, Implanted.Embryo) %>%
  arrange_all()


morphagerepro$pregnancy <- as.integer(NA)
morphagerepro[Code.Implanted.Embryo == 1, pregnancy := 1]
morphagerepro[Code.Implanted.Embryo == 2, pregnancy := 0]

morphagerepro %>%
  filter(Month %in% c(10:12, 1, 2)) %>%
  filter(Exclude.From.Repro == FALSE) %>%
  filter(Sex == "F") %>%
  distinct(Code.Female.Maturity, Female.Maturity, maturity, Code.Implanted.Embryo, Implanted.Embryo, pregnancy) %>%
  arrange_all()


## define early puppers  -----
# 0=FALSE  --- 1=TRUE
morphagerepro$EP <- 0
morphagerepro[Code.Female.Maturity == 8 &
                doy < last.ep.date, EP := 1]

# this looks OK
morphagerepro %>%
  filter(Month %in% c(10:12, 1, 2)) %>%
  filter(Exclude.From.Repro == FALSE) %>%
  filter(Sex == "F") %>%
  distinct(Code.Female.Maturity, Female.Maturity, maturity, Code.Implanted.Embryo, Implanted.Embryo, pregnancy, EP) %>%
  arrange_all()

### Code 8, later in the season set as pregnant  ----
# There are records of Code.Female.Maturity = 8 until doy 67
# Consider those between doy 51 and 67 as pregnant
morphagerepro[Code.Female.Maturity == 8] %>%
  distinct(doy, maturity, pregnancy)%>%
  arrange_all()

morphagerepro[Code.Female.Maturity == 8 &
                data.table::between(doy , last.ep.date, 70), pregnancy := 1]

# plots -----
## population numbers -----
p.seal.pop <- ggplot(seal.pop, aes(x = Year, y = N)) +
  geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI), alpha = 0.4) +
  geom_line() +
  scale_y_continuous(labels = scales::comma) +
  NULL

## LW -----

plotdat <- morphagerepro[Sex == "F" &
                           !is.na(Body.Length) &
                           !is.na(Body.Weight) ,
                         .(Body.Length, Body.Weight, ID.Sex, Code.PelageType, cohortyear, Female.Maturity, Code.Female.Maturity)] %>%
  arrange(Body.Length)
names(plotdat) <- c('length', 'weight', 'idsex', 'codepelagetype', 'cohortyear', 'femmat', 'codefemmat')

# consider only beater and older
plotdat[is.na(codepelagetype)]
plotdat <- plotdat[codepelagetype %!in%  c(0:5)]
# plotdat3 <- plotdat[codepelagetype > 5]

# potential outliers
# I received feedback from Shelley Lang - different for each seal
# assets/DataDecisions_ADB_2025-08_SL.docx
id_outs <- c(
  '20113315F',
  '20215080F',
  '20181046F',
  '20224400F',
  '20224484F',
  '20192983F'
  # '19962312F',
  # '20072884F'
)

outs <- plotdat[idsex %in% id_outs]


outs[,.(idsex, length, weight, femmat)]

# plotdat <- plotdat[idsex %!in% id_outs]

p.lw <-   plot_lw(plotdat$length, plotdat$weight, plotdat$idsex)

p.lw <- p.lw +
  geom_point(data = plotdat,
             alpha = 0.6,
             pch = 16,
             aes(x = length,
                 y = weight,
                 label = idsex,
                 color = femmat)) +
  # geom_point(data = outs, aes(x = length,
  #                             y = weight,
  #                             label = idsex),
  #            size = 3, color = "black", fill = "black") +
  theme(legend.position = 'bottom')

ggplotly(p.lw)



plot_lw( plotdat[!idsex %in% id_outs]$length,
         plotdat[!idsex %in% id_outs]$weight,
         plotdat[!idsex %in% id_outs]$idsex) +
  geom_point(data = plotdat[!idsex %in% id_outs],
             alpha = 0.4,
             pch = 16,
             aes(x = length,
                 y = weight,
                 label = idsex,
                 color = as.factor(cohortyear)))


plotdat[!idsex %in% id_outs] %>%
  ggplot(., aes(x = length,
                y = weight)) +
  geom_point(alpha = 0.4,
             pch = 16) +
  geom_smooth()

## Relative condition -----
plotdat <- plotdat %>%
  mutate(krel = krel(length, weight))
# outliers?
# will leave them in
plotdat[!idsex %in% id_outs] %>%
  mutate(krel = krel(length, weight)) %>%
  filter(krel > 1.5)
pot_outs <- plotdat[!idsex %in% id_outs] %>%
  mutate(krel = krel(length, weight)) %>%
  filter(krel > 1.5) %>% pull(idsex)

p <-   plot_lw(plotdat$length, plotdat$weight, plotdat$idsex)
p +
  geom_point(data = plotdat,
             alpha = 0.6,
             pch = 16,
             aes(x = length,
                 y = weight,
                 label = idsex,
                 color = femmat)) +
  geom_point(data = plotdat[idsex %in% pot_outs], aes(x = length,
                                                      y = weight,
                                                      label = idsex),
             size = 3, color = "black", fill = "black") +
  theme(legend.position = 'bottom')
rm(p)

ggplotly(
  plotdat[!idsex %in% id_outs] %>%
    mutate(krel = krel(length, weight)) %>%
    left_join(morphagerepro[,.(idsex = ID.Sex, Female.Maturity, Month )]) %>%
    ggplot(.,aes(x = (cohortyear),
                 y = krel,
                 fill = Female.Maturity)) +
    geom_hline(yintercept = 1) +
    geom_quasirandom(
      # side = 1,
      alpha = .8, pch = 16,
      # corral = "gutter",
      # method = "compactswarm",
      cex = 1,
      # corral.width = 0.4,
      # stat = density,
      aes(color = factor(cohortyear), label = idsex)) +
    labs(x = "Cohort year", y = "Relative condition") +

    facet_wrap(.~Female.Maturity) +
    scale_x_continuous(breaks = seq(1978, 2022, 2), limits = c(1978, 2022)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45),
          panel.grid.minor = element_blank()) +
    geom_smooth()
)



ggplotly(
  plotdat[!idsex %in% id_outs] %>%
    mutate(krel = krel(length, weight)) %>%
    filter(codefemmat %in%  c(2:5, 8,7)) %>%
    left_join(morphagerepro[,.(idsex = ID.Sex,  Month )]) %>%
    transform(femmat = factor(femmat,
                              levels = c(
                                'Mature; pregnant (implanted embryo)',
                                'Mature; lactating',
                                'Mature; post lactation/prior to ovulation',
                                'Mature; recently given birth outside of normal period',
                                'Mature; pregnant (delay period)',
                                'Mature; not-pregnant, parous'
                              ))) %>%
    ggplot(.,aes(x = (cohortyear),
                 y = krel,
                 fill = femmat )) +
    # geom_half_violin(side = "l", alpha = 0.6, trim = FALSE, fill = 'gray60') +
    # geom_half_point(alpha = 0.3) +
    # geom_beeswarm(
    #   # side = 1,
    #   alpha = .8, pch = 16,
    #   corral = "gutter",
    #   method = "compactswarm",
    #   cex = 1,
    #   corral.width = 0.4,
    #   # stat = density,
    #   aes(color = factor(cohortyear))) +
  geom_hline(yintercept = 1) +
    geom_quasirandom(
      # side = 1,
      alpha = .8, pch = 16,
      # corral = "gutter",
      # method = "compactswarm",
      cex = 1,
      # corral.width = 0.4,
      # stat = density,
      aes(color = factor(cohortyear), label = idsex)) +
    labs(x = "Cohort year", y = "Relative condition") +
    facet_grid(femmat ~Month) +
    scale_x_continuous(breaks = seq(1978, 2022, 2), limits = c(1978, 2022)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45),
          panel.grid.minor = element_blank())# +
  # geom_smooth()
)


# bring relative condition to dataset -----

morphagerepro <- left_join(morphagerepro,
                           plotdat[,.(ID.Sex = idsex, krel)])

## age - weight ----
ggplotly(
  ggplot(morphagerepro %>%
           filter(Final.Cohort.Age<70) %>%
           filter(Sex == "F"),
         aes(Final.Cohort.Age, Body.Weight, label = ID.Sex)) +
    geom_point(alpha = 0.4)
)


# biological rates ----

## define dataset to calculate biological rates -----
dat.biolrates <- morphagerepro %>%
  filter(Sex == "F") %>%
  filter(Exclude.From.Repro == FALSE) %>%
  filter(Month %in% c(10:12, 1, 2))


dat.biolrates %>% distinct(Code.Female.Maturity, Female.Maturity) %>% arrange_all()

dat.biolrates %>%
  distinct(Code.Female.Maturity, Female.Maturity, maturity, Code.Implanted.Embryo, Implanted.Embryo, pregnancy, EP) %>%
  arrange_all()

## pregnancy rate ----
# pregnancy rate = No. of pregnant females/No. of mature females
### number of mature and immature females, by cohort year ----
mat.summary <- dat.biolrates[!is.na(maturity), .N, by = c('maturity', 'cohortyear')] %>%
  pivot_wider(names_from = maturity, values_from = N) %>%
  rename(n.mature = `1`,
         n.immature = `0`) %>%
  data.table()

### number of pregnant and non-pregnant females, by cohort year ----
preg.summary <- dat.biolrates[!is.na(pregnancy), .N, by = c('pregnancy', 'cohortyear')] %>%
  pivot_wider(names_from = pregnancy, values_from = N) %>%
  rename(n.pregnant = `1`,
         n.nonpregnant = `0`) %>%
  data.table()

### merge
biolrates <- merge(mat.summary, preg.summary, by = "cohortyear")

## abortion rate ----
# abortion rate = No. of abortions/(No. of abortions + No. of viable pregnancies)
### number of early puppers, by cohort year ----
ep.summary <- dat.biolrates[!is.na(EP), .N, by = c('EP', 'cohortyear')] %>%
  pivot_wider(names_from = EP, values_from = N) %>%
  rename(n.EP = `1`,
         n.nonEP = `0`) %>%
  # drop n.nonEP
  select(-n.nonEP) %>%
  data.table()

# years when n.EP is NA means that there were no early puppers
# replace by zero
ep.summary[is.na(n.EP), n.EP := 0]

### merge
biolrates <- merge(biolrates, ep.summary, by = "cohortyear")

## biol rates + biol rtes.old ---------
biolrates <- rbind( biolrates.old, biolrates)

### calculate pregnancy rate ----
# pregnancy rate = No. of pregnant females/No. of mature females
biolrates[, pregrate := n.pregnant/n.mature]

###  calculate abortion rate ----
biolrates[, totpreg := n.pregnant + n.EP]
biolrates[, abrate := n.EP/totpreg]

biolrates <- left_join(
  data.table(cohortyear = 1950:2022), biolrates)

ggplot(biolrates, aes(x = cohortyear, y = abrate)) +
  # geom_smooth(span = 0.3) +
  geom_point() +
  geom_line(lty=2)

ggplot(biolrates, aes(x = cohortyear, y = pregrate)) +
  # geom_smooth(span = 0.3) +
  geom_point() +
  geom_line(lty=2) +
  NULL

ggplot(dat.biolrates, aes(x = krel, y = EP, colour = as.factor(cohortyear)) )+
  geom_point()

# summarize condition ----
summaryK <- dat.biolrates %>%
  group_by(cohortyear) %>%
  reframe(meanK = mean (krel, na.rm = TRUE),
          sdK = sd (krel, na.rm = TRUE)) %>%
  data.table()

biolrates <- left_join(biolrates, summaryK)

# output ----
fwrite(x = dat.biolrates,
       na = NA,
       file = paste0(here::here(), "/data/seal/CleanDatasetForBiologicalRates.csv"))

fwrite(x = biolrates,
       na = NA,
       file = paste0(here::here(), "/data/seal/BiologicalRates.csv"))
