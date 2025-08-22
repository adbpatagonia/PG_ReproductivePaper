# Sent the output file to SHelley L
library(data.table)
library(ggplot2)
library(plotly)
library(tidyverse)

morphagerepro <- openxlsx::read.xlsx(xlsxFile = paste0(here::here(), "/data-raw/seal/qry_CollectionMorphAgeReproDataforABuren_July2025_v2.xlsx"),
                                     sheet = "qry_CollectionMorphAgeReproData") %>%
  data.table()


# these are the mean monthly  foetal weights I used for the growth paper
# fwd <- 2.9
# fwj <- 6.3
# fwf <- 8.3
#
# I used these values to identify potential outliers in the foetus data
morphagerepro[, cohortyear := as.integer(substr(ID.Sex, 1, 4))]
morphagerepro[, cohortyear := Year]
morphagerepro %>%
  filter(Final.Cohort.Age == 91 ) %>%
  filter(Month == 12 & Body.Weight > 6 |
           Month == 1 & Body.Weight > 10 |
           Month == 2 & Body.Weight > 12) %>%
  arrange(ID.Sex) #%>%
  # filter(cohortyear %in% c(2013,2014)) %>%
  # mutate(Body.Weight = ifelse(cohortyear %in% c(2013,2014), Body.Weight * 0.453592, Body.Weight))
  # fwrite(., file = 'output/PotentialFoetusOutliers.csv')


# morphagerepro[cohortyear %in% c(2013,2014), Body.Weight:= Body.Weight * 0.453592]


ggplotly(
  morphagerepro %>%
    filter(Final.Cohort.Age == 91 ) %>%
    filter(Month %in% c(1, 2, 12)) %>%
    group_by(Year, Month) %>%
    reframe(meanw = mean(Body.Weight)) %>%
    ggplot(., aes(x = Year, y = meanw)) +
    geom_point() +
    facet_grid(.~Month)
)

# ggplotly(
  morphagerepro %>%
    # filter(Final.Cohort.Age == 91 ) %>%
    # filter(Month %in% c(1, 2, 12)) %>%

    ggplot(., aes(x = cohortyear, y = Body.Weight)) +
    geom_point() +
    geom_point(data = morphagerepro[cohortyear %in% c(2013,2014)], col = 'red') +
    facet_grid(.~Month)
# )
