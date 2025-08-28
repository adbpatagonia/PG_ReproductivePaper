# libraries ------
library(mgcv)
library(gratia)
library(gam.hp)
library(plotly)
library(ggbeeswarm)

# data -----
source( paste0(here::here(), "/analysis/0_seal_data.R"))
dat.biolrates <- fread(paste0(here::here(), "/data/seal/CleanDatasetForBiologicalRates.csv"))

# functions ----
source( paste0(here::here(), "/R/y.transf.betareg.r"))


# pregnancy rate model ----
# change the pregrate in 1979 to the value used in the resdoc
# will circle back to this once I finish comparing rates
seal.data[cohortyear == 1979, pregrate := 0.941]

# change pregrate and abrate in 2012 to the values used in the resdoc
# will circle back to this once I finish comparing rates
seal.data[cohortyear == 2012, pregrate := 0.627]
seal.data[cohortyear == 2012, abrate := 0.036]

# change abrate in 2004 to the values used in the resdoc
# will circle back to this once I finish comparing rates
seal.data[cohortyear == 2004, abrate := 0.238]

# fit model
m.preg <- gam(pregrate ~ s(sealpop) +
                s(abrate) ,
              data = seal.data,
              method = "REML",
              family = betar())


appraise(m.preg)
k.check(m.preg)
concurvity(m.preg)
concurvity(m.preg, full = FALSE)
m.drop <- update(m.preg, . ~ . - s(abrate))

AIC(m.preg, m.drop)
draw(m.preg, residuals = FALSE)
draw(m.preg, residuals = TRUE)
summary(m.preg)

preg.fits <- predict(m.preg,
                     type = "response",
                     se.fit = TRUE) %>%
  data.frame() %>%
  mutate(lci.preg = fit - 1.96 * se.fit,
         uci.preg = fit + 1.96 * se.fit) %>%
  select(-se.fit) %>%
  rename(fit.preg = fit) %>%
  bind_cols(
    seal.data[!is.na(pregrate), .(cohortyear, sealpop, pregrate, abrate)])

ggplot(preg.fits, aes(x = cohortyear, y = pregrate)) +
  geom_point() +
  ylim(0,1) +
  geom_line(aes(y = fit.preg)) +
  geom_ribbon(aes(y = fit.preg, ymin = lci.preg, ymax = uci.preg), alpha = 0.2)

plot.gamhp(gam.hp(m.preg), plot.perc = TRUE)

ggplot(preg.fits, aes(y = fit.preg, x = pregrate)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  ylim(0, 1) +
  xlim(0, 1)
