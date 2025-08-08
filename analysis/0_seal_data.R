# 2025-01-27
# script to plot the seal data for the
# Environmental effects on harp seal reproduction in the NW Atlantic paper

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

# read data -----
## population numbers -----
seal.pop <- openxlsx::read.xlsx(xlsxFile = paste0(here::here(), "/data/seal/Population_Size_1952-2024.xlsx"),
                                sheet = "Population_Size_1952-2024") %>%
  data.table()

## ovary ----
ovary <- openxlsx::read.xlsx(xlsxFile = paste0(here::here(), "/data/seal/qry_HarpOvaryDataForABuren_July2025.xlsx"),
                             sheet = "qry_HarpOvaryDataForABuren_July") %>%
  data.table()

## morphagerepro ----
morphagerepro <- openxlsx::read.xlsx(xlsxFile = paste0(here::here(), "/data/seal/qry_CollectionMorphAgeReproDataforABuren_July2025_v2.xlsx"),
                                     sheet = "qry_CollectionMorphAgeReproData") %>%
  data.table()

# wrangle data -----


## morphagerepro ----
### remove dead on beach ----
# Screen out ID Sex 20180054F (MarineMammal ID 64849) – seal was found dead on a beach
morphagerepro <- morphagerepro[ID.Sex != "20180054F"]

### remove foetus, stillborn & starvling ----
morphagerepro <- morphagerepro[Final.Cohort.Age < 90]

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

# from P Goulet
# 20161235M the girth and length have been inverted, should be girth=130  and length=163

morphagerepro[ID.Sex == "20161235M"]
morphagerepro[ID.Sex == "20161235M"]
morphagerepro[ID.Sex == "20161235M", Body.Length  := 163]
morphagerepro[ID.Sex == "20161235M",  Maximum.Girth := 130]


### check inconsistencies with previous instructions from Rachel -----
# edits included in file README_Buren.txt, commit a415c479c87771fccf382a41c95f1ec712a1e696
# - remove for now. needs to be checked against datasheets.
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
# these should be immature
morphagerepro[ID.Sex %in% changeage99 & Code.Female.Maturity > 1,.(ID.Sex, Code.Female.Maturity)]
# no data for these in tbl_ovary
ovary[ID.Sex %in% morphagerepro[ID.Sex %in% changeage99 & Code.Female.Maturity > 1,.(ID.Sex)]]

### changeage1
# these should be immature
morphagerepro[ID.Sex %in% changeage1  ,.(ID.Sex, Code.Female.Maturity)]
# no data for these in tbl_ovary
ovary[ID.Sex %in% morphagerepro[ID.Sex %in% changeage1,.(ID.Sex)]]

# no morph data and no ovary data - I will change code.female.maturity to 1, until I am told otherwise
morphagerepro[ID.Sex %in% changeage99 & Code.Female.Maturity > 1, Code.Female.Maturity := 1]

### changeage20
morphagerepro[ID.Sex %in% changeage20 ,.(ID.Sex, Code.Female.Maturity, Female.Maturity)]
# no data for these in tbl_ovary
ovary[ID.Sex %in% morphagerepro[ID.Sex %in% changeage20 ,.(ID.Sex)]]

### remove
morphagerepro[ID.Sex %in% remove,.(ID.Sex, Code.Female.Maturity, Female.Maturity)]
ovary[ID.Sex %in% remove]

# I will remove until I am told otherwise
ovary <- ovary[!ID.Sex %in% remove]


### removefornow
morphagerepro[ID.Sex %in% removefornow,.(ID.Sex, Code.Female.Maturity, Female.Maturity)]
ovary[ID.Sex %in% removefornow]

# I will remove until I am told otherwise
morphagerepro <- morphagerepro[!ID.Sex %in% removefornow]
ovary <- ovary[!ID.Sex %in% removefornow]


### exclude from repro -----
ids_exc_repro <- morphagerepro[Exclude.From.Repro == TRUE, .(ID.Sex)] %>% pull()

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
# at least 3 seals have the wrong weight, turn them to NA
ids_w <-  c("19940928F", "20042978F", "20062138F")
morphagerepro[ID.Sex %in% ids_w, Body.Weight := NA]

### Senescence ----
# OK
femcode6 <- c('20091859F','20111392F','20130006F','20140355F','20170015F','20171718F')
morphagerepro[ID.Sex %in% femcode6]

### Remove OU ----
# OK - variable Exclude.from.repro == TRUE
removeOU <- c('20130099F','20160009F','20171726F')
morphagerepro[ID.Sex %in% removeOU]


### remove fetal weight ----
# remove the weight of the foetus for pregnant females
# I obtained monthly mean foetal weights from the PG_Growth work
# Ask Shelley to send updated values
fwd <- 2.9
fwj <- 6.3
fwf <- 8.3
morphagerepro$Body.Weight <- ifelse(morphagerepro$Code.Female.Maturity %in% 2 & morphagerepro$Month == 12,
                                    morphagerepro$Body.Weight - fwd,
                                    morphagerepro$Body.Weight)
morphagerepro$Body.Weight <- ifelse(morphagerepro$Code.Female.Maturity %in% 2 & morphagerepro$Month == 1,
                                    morphagerepro$Body.Weight - fwj,
                                    morphagerepro$Body.Weight)
morphagerepro$Body.Weight <- ifelse(morphagerepro$Code.Female.Maturity %in% 2 & morphagerepro$Month == 2,
                                    morphagerepro$Body.Weight - fwf,
                                    morphagerepro$Body.Weight)

### cohort year ----
morphagerepro[, cohortyear := as.integer(substr(ID.Sex, 1, 4))]

### female maturity -----
#### potential inconsistencies given age -----
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


#### newborns to beaters ----
morphagerepro[is.na(Female.Maturity)]
# set Immature for newborns to beaters
# simply making sure that they are all small
morphagerepro[Code.PelageType < 7,.(Body.Weight)] %>% max(., na.rm = TRUE)
morphagerepro[Code.PelageType < 7, Female.Maturity := "Immature"]
morphagerepro[Code.PelageType < 7, Code.Female.Maturity := 1]



## ovary ----
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
plotdat <- plotdat[codepelagetype > 5]

# potential outliers
id_outs <- c(
  '20113315F',
  '20215080F',
  '20181046F',
  '20224400F',
  '20224484F',
  '20192983F',
  '19962312F',
  '20072884F'
)

outs <- plotdat[idsex %in% id_outs]

p.lw <-   plot_lw(plotdat$length, plotdat$weight, plotdat$idsex)

p.lw <- p.lw +
  geom_point(data = plotdat,
             alpha = 0.6,
             pch = 16,
             aes(x = length,
                 y = weight,
                 label = idsex,
                 color = femmat)) +
  geom_point(data = outs, aes(x = length,
                              y = weight,
                              label = idsex),
             size = 3, color = "black", fill = "black") +
  theme(legend.position = 'bottom')

# ggplotly(p.lw)



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
