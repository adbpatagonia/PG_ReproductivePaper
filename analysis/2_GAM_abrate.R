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

# wrnagle data -----
# change abrate in 2012 to the values used in the resdoc
# will circle back to this once I finish comparing rates
seal.data[cohortyear == 2012, abrate := 0.036]

# change abrate in 2004 to the values used in the resdoc
# will circle back to this once I finish comparing rates
seal.data[cohortyear == 2004, abrate := 0.238]



# abortion rate -----
ggplot(biolrates, aes(meanK, abrate)) +
  geom_point()

ggplot(biolrates[!is.na(meanK)], aes(x = cohortyear, y = meanK) )+
  geom_point() +
  # geom_smooth() +
  NULL

ggplot(biolrates[!is.na(meanK)], aes(x = cohortyear, y = abrate) )+
  geom_point() +
  # geom_smooth() +
  NULL

# fit model
# family beta works in the open interval (0,1)
# there are several values of 0 in the dataset. Need to fix

dat.mod <- seal.data[!is.na(abrate), .(cohortyear, abrate, meanK)]
dat.mod <- dat.mod[!is.na(meanK)]

dat.mod[, abrate_adj := y.transf.betareg(abrate)]

m.ab <- gam(abrate_adj ~ s(meanK) ,
            data = dat.mod,
            method = "REML",
            family = betar())


appraise(m.ab)
k.check(m.ab)
concurvity(m.ab)
concurvity(m.ab, full = FALSE)

draw(m.ab, residuals = FALSE)
# the outlier is cohortyear 2012
# the second outlier is cohortyear 2004
draw(m.ab, residuals = TRUE)
summary(m.ab)

ab.fits <- predict(m.ab,
                   type = "response",
                   se.fit = TRUE) %>%
  data.frame() %>%
  mutate(lci.ab = fit - 1.96 * se.fit,
         uci.ab = fit + 1.96 * se.fit) %>%
  select(-se.fit) %>%
  rename(fit.ab = fit) %>%
  bind_cols(
    dat.mod)

ggplot(ab.fits, aes(x = meanK, y = abrate_adj)) +
  # ylim(0,1) +
  xlim(0.9, 1.35) +
  geom_line(aes(y = fit.ab)) +
  geom_ribbon(aes(y = fit.ab, ymin = lci.ab, ymax = uci.ab), alpha = 0.2) +
  geom_point()

# plot.gamhp(gam.hp(m.ab), plot.perc = TRUE)

ggplot(ab.fits, aes(y = fit.ab, x = abrate_adj)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  ylim(0, .5) +
  xlim(0, .5)


# abortions ----
ggplot(dat.biolrates, aes(x = krel, y = EP, colour = as.factor(cohortyear)) )+
  geom_point()



dat.mod <- dat.biolrates[!is.na(krel),.(cohortyear, EP, krel)]
m.abs <- glm(EP ~ krel ,
             data = dat.mod,
             # method = "REML",
             family = binomial())


appraise(m.abs)
draw(m.abs)
summary(m.abs)
plot(m.abs)
k.check(m.abs)


newdat <- data.frame(krel = seq(min(dat.mod$krel), max(dat.mod$krel), length.out = 200))
newdat <- predict(m.abs, newdat, type = "response", se.fit = TRUE) %>%
  data.frame() %>%
  mutate(lci.ab = fit - 1.96 * se.fit,
         uci.ab = fit + 1.96 * se.fit) %>%
  select(-se.fit) %>%
  rename(fit.ab = fit) %>%
  bind_cols(newdat)

abs.fits <- predict(m.abs,
                    type = "response",
                    se.fit = TRUE) %>%
  data.frame() %>%
  mutate(lci.ab = fit - 1.96 * se.fit,
         uci.ab = fit + 1.96 * se.fit) %>%
  select(-se.fit) %>%
  rename(fit.ab = fit) %>%
  bind_cols(
    dat.mod)

ggplot(dat.mod, aes(x = krel, y = EP)) +
  # ylim(0,1) +
  xlim(0.9, 1.35) +
  geom_line(data = newdat, aes(y = fit.ab)) +
  geom_ribbon(data = newdat,aes(y = fit.ab, ymin = lci.ab, ymax = uci.ab), alpha = 0.2) +
  geom_point()

