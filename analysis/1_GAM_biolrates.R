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

# condition -----
sum.relk <- dat.biolrates %>%
  group_by(Code.Female.Maturity, Female.Maturity, cohortyear) %>%
  reframe(meanK = mean(krel, na.rm = TRUE),
          sdK = sd(krel, na.rm = TRUE)) %>%
  data.table()

ggplotly(
  ggplot(dat.biolrates[Code.Female.Maturity %in% c(8, 2, 4, 1)], aes(x = cohortyear, y = krel) )+
    # geom_point() +
    facet_grid(Female.Maturity ~ .) +
    geom_hline(yintercept = 1) +
    geom_quasirandom(
      # side = 1,
      alpha = .3, pch = 16,
      # corral = "gutter",
      # method = "compactswarm",
      cex = 1,
      # corral.width = 0.4,
      # stat = density,
      aes(color = factor(cohortyear), label = ID.Sex)) +
    # stat_smooth() +
    theme(legend.position = 'none') +
    geom_linerange(data = sum.relk[Code.Female.Maturity %in% c(8, 2, 4, 1)],
                   aes(x = cohortyear,
                       color = factor(cohortyear),
                       y = meanK,
                       ymin = meanK - sdK,
                       ymax = meanK + sdK) ) +
    geom_point(data = sum.relk[Code.Female.Maturity %in% c(8, 2, 4, 1)],
               aes(x = cohortyear,
                   color = factor(cohortyear),
                   y = meanK) ) +
    NULL
)


dat.mod.krel <- dat.biolrates[Code.Female.Maturity %in% c(8, 2, 4, 1)] %>%
  mutate(cyear =  as.factor(cohortyear)) %>%
  mutate(Female.Maturity = as.factor(Female.Maturity)) %>%
  select(cohortyear, cyear, Code.Female.Maturity, Female.Maturity, krel) %>%
  filter(!is.na(krel)) %>%
  transform(Female.Maturity = factor(Female.Maturity,
                                     levels = c(
                                     "Mature; pregnant (implanted embryo)",
                                     "Mature; recently given birth outside of normal period",
                                     "Immature",
                                     "Mature; not-pregnant, parous"
                                     )))

dat.mod.krel %>% distinct(Female.Maturity)
dat.mod.krel %>% distinct(cyear)


#' @ADB: need to fit HGAMs as per Pedersen et al
#'source( paste0(here::here(), "/analysis/1a_GAM_biolrates_krel.R"))
m.krel <- gam(krel ~ s(cohortyear) +
                Female.Maturity,
              method = "REML",
              data = dat.mod.krel)
summary(m.krel)
anova(m.krel)
draw(m.krel)
draw(m.krel, residuals = TRUE)
# performance::check_model(m.krel)

plot.gamhp(gam.hp(m.krel), plot.perc = TRUE)
