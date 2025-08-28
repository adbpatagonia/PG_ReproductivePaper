# libraries ------
library(mgcv)
library(gratia)
library(gam.hp)
library(plotly)
library(ggbeeswarm)
library(ggprism)
library(ggtext)

# data -----
source( paste0(here::here(), "/analysis/0_seal_data.R"))
dat.biolrates <- fread(paste0(here::here(), "/data/seal/CleanDatasetForBiologicalRates.csv"))
source(paste0(here::here(), "/analysis/0_env_data.R"))
source(paste0(here::here(), "/analysis/0_prey_data.R"))
bootCI_krel <- fread(paste0(here::here(), "/output/RelativeCondition_bootstrap.csv"))

# functions ----
source( paste0(here::here(), "/R/y.transf.betareg.r"))

# define model data -----
dat.mod.krel <- dat.biolrates[Code.Female.Maturity %in% c(8, 2, 4, 1)] %>%
  mutate(cyear =  as.factor(cohortyear)) %>%
  mutate(Female.Maturity = as.factor(Female.Maturity)) %>%
  select(ID.Sex, cohortyear, cyear, Code.Female.Maturity, Female.Maturity, krel) %>%
  filter(!is.na(krel)) %>%
  transform(Female.Maturity = factor(Female.Maturity,
                                     levels = c(
                                       "Mature; pregnant (implanted embryo)",
                                       "Mature; recently given birth outside of normal period",
                                       "Immature",
                                       "Mature; not-pregnant, parous"
                                     ),
                                     labels = c(
                                       "Mature; pregnant",
                                       "Mature; early pupper",
                                       "Immature",
                                       "Mature; not-pregnant"
                                     )))

dat.mod.krel <- left_join(dat.mod.krel,
                          env.dat %>%
                            rename(cohortyear = year)
) %>%
  data.table()

# EDA -----
sum.relk <- dat.mod.krel %>%
  group_by(Code.Female.Maturity, Female.Maturity, cohortyear) %>%
  reframe(meanK = mean(krel, na.rm = TRUE),
          sdK = sd(krel, na.rm = TRUE)) %>%
  data.table()

sum.relk <- sum.relk %>%
  left_join(bootCI_krel, by = c("cohortyear", "Code.Female.Maturity"))


p.cond.femmat.data <- ggplot(dat.mod.krel, aes(x = cohortyear, y = krel) )+
  facet_grid(Female.Maturity ~ .) +
  geom_hline(yintercept = 1) +
  geom_quasirandom(
    alpha = .3, pch = 16,
    cex = 1,
    aes(color = factor(cohortyear), label = ID.Sex)) +
  geom_linerange(data = sum.relk,
                 aes(x = cohortyear,
                     color = factor(cohortyear),
                     y = meanK,
                     ymin = krellb   ,
                     ymax = krelub) ) +
  geom_point(data = sum.relk,
             aes(x = cohortyear,
                 color = factor(cohortyear),
                 y = meanK) ) +
  ylab("Relative Condition (*K<sub>r</sub>*)") +
  xlab("Cohort year") +
  scale_x_continuous(breaks = seq(1980, 2020, 10),
                     minor_breaks = seq(1985, 2015, 10),
                     guide = guide_prism_minor()) +
  scale_y_continuous(breaks = seq(0.6, 1.6, 0.20),
                     minor_breaks = seq(0.5, 1.7, 0.20),
                     guide = guide_prism_minor()) +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(colour = 'black'),
        axis.title.y = element_markdown()) +
  NULL




ggplotly(
  ggplot(dat.mod.krel, aes(x = NLCI, y = krel) )+
    # geom_point() +
    # facet_grid(Female.Maturity ~ .) +
    geom_hline(yintercept = 1) +
    geom_quasirandom(
      # side = 1,
      alpha = .3, pch = 16,
      # corral = "gutter",
      # method = "compactswarm",
      cex = 1,
      # corral.width = 0.4,
      # stat = density,
      aes(label = ID.Sex)) +
    # stat_smooth() +
    theme(legend.position = 'none') +
    geom_smooth() +
    # geom_linerange(data = sum.relk,
    #                aes(x = cohortyear,
    #                    color = factor(cohortyear),
    #                    y = meanK,
    #                    ymin = meanK - sdK,
    #                    ymax = meanK + sdK) ) +
    # geom_point(data = sum.relk,
    #            aes(x = cohortyear,
    #                color = factor(cohortyear),
    #                y = meanK) ) +
    NULL
)







# model G ----
# A single common (global) smoother for all observations
# two smoothers:
# 1. a Thin plate regression spline of cohortyear, and
# 2. a random effect for Female.Maturity to model reach-specific intercepts
# The random effect smoother (bs="re") that we used for the Female.Maturity factor
# always has a k value equal to the number of levels in the grouping variable
krel_modG <- gam(krel ~ s(cohortyear, k = 30, bs = "tp") +
                   s(Female.Maturity, k = length(levels(dat.mod.krel$Female.Maturity)), bs = "re"),
                 data = dat.mod.krel,
                 method = "REML", family = "gaussian")
## Model check ----
k.check(krel_modG)
appraise(krel_modG)
concurvity(krel_modG)
concurvity(krel_modG, full = FALSE)
# draw(krel_modG)
# visreg::visreg(krel_modG,"cohortyear", by = "Female.Maturity")
# overall concurvity
o_conc <- concrvity(krel_modG)
draw(o_conc)
# pairwise concurvity
p_conc <- concrvity(krel_modG, pairwise = TRUE)
draw(p_conc)

plot.gamhp(gam.hp(krel_modG), plot.perc = TRUE)

# model GS ----
# A single common smoother plus group-level smoothers that have the same wiggliness
# analogue to a GLMM with varying slopes
# two smoothers:
# 1. a Thin plate regression spline of cohortyear, and
# 2.  group-level smoother (reach-specific)

krel_modGS <- gam(krel ~ s(cohortyear,
                           # k = 30,
                           bs = "tp") +
                    s(cohortyear, Female.Maturity, k = 30, bs = "fs", m = 2),
                  data = dat.mod.krel,
                  method = "REML", family = "gaussian")

## Model check ----
k.check(krel_modGS)
appraise(krel_modGS)
concurvity(krel_modGS)
concurvity(krel_modGS, full = FALSE)
# draw(krel_modGS)
visreg::visreg(krel_modGS,"cohortyear", by = "Female.Maturity")
# overall concurvity
o_conc <- concrvity(krel_modGS)
draw(o_conc)
# pairwise concurvity
p_conc <- concrvity(krel_modGS, pairwise = TRUE)
draw(p_conc)

plot.gamhp(gam.hp(krel_modGS), plot.perc = TRUE)

# model GI ----
# A single common smoother plus group-level smoothers that have their own level of wiggliness
# three smoothers:
# 1. a Thin plate regression spline of cohortyear, and
# 2.  group-level smoother (reach-specific). Diffrence from model GS: m = 1
# 3. a random effect for Female.Maturity to model reach-specific intercepts
# The random effect smoother (bs="re") that we used for the Female.Maturity factor
# always has a k value equal to the number of levels in the grouping variable

krel_modGI <- gam(krel ~ s(cohortyear, k = 30, bs = "tp") +
                    s(cohortyear, Female.Maturity, k = 30, bs = "fs", m = 1) +
                    s(Female.Maturity, k = length(levels(dat.mod.krel$Female.Maturity)), bs = "re"),
                  data = dat.mod.krel,
                  method = "REML", family = "gaussian")

## Model check ----
k.check(krel_modGI)
appraise(krel_modGI)
concurvity(krel_modGI)
concurvity(krel_modGI, full = FALSE)
# draw(krel_modGI)
# visreg::visreg(krel_modGI,"cohortyear", by = "Female.Maturity")
# overall concurvity
o_conc <- concrvity(krel_modGI)
draw(o_conc)
# pairwise concurvity
p_conc <- concrvity(krel_modGI, pairwise = TRUE)
draw(p_conc)

plot.gamhp(gam.hp(krel_modGI), plot.perc = TRUE)

# model S ----
# Model S (shared smoothers) is model GS without the global smoother term
# This model assumes all groups have the same smoothness, but that the individual shapes of the smooth terms are not related
# If in a study there are very few data points in each grouping level (relative to the strength of the functional relationship of interest), estimates from model S will typically be much more variable than from model GS
# one smoother:
# 1.  group-level smoother (reach-specific)

krel_modS <- gam(krel ~ s(cohortyear, Female.Maturity, k = 30, bs = "fs", m = 2),
                 data = dat.mod.krel,
                 method = "REML", family = "gaussian")

## Model check ----
k.check(krel_modS)
appraise(krel_modS)
concurvity(krel_modS)
concurvity(krel_modS, full = FALSE)
# draw(krel_modS)
# visreg::visreg(krel_modS,"cohortyear", by = "Female.Maturity")
# overall concurvity
o_conc <- concrvity(krel_modS)
draw(o_conc)
# pairwise concurvity
p_conc <- concrvity(krel_modS, pairwise = TRUE)
draw(p_conc)


# plot.gamhp(gam.hp(krel_modS), plot.perc = TRUE)

# model I ----
# Model I is model GI without the first term
#  group-level smoothers that have their own level of wiggliness
# two smoothers:
# 1.  group-level smoother (reach-specific). Diffrence from model GS: m = 1
# 2. a random effect for Female.Maturity to model reach-specific intercepts
# The random effect smoother (bs="re") that we used for the Female.Maturity factor
# always has a k value equal to the number of levels in the grouping variable
krel_modI <- gam(krel ~ s(cohortyear, Female.Maturity, k = 30, bs = "fs", m = 1) +
                   s(Female.Maturity, k = length(levels(dat.mod.krel$Female.Maturity)), bs = "re"),
                 data = dat.mod.krel,
                 method = "REML", family = "gaussian")

## Model check ----
k.check(krel_modI)
appraise(krel_modI)
concurvity(krel_modI)
concurvity(krel_modI, full = FALSE)
# draw(krel_modI)
# visreg::visreg(krel_modI,"cohortyear", by = "Female.Maturity")
# overall concurvity
o_conc <- concrvity(krel_modI)
draw(o_conc)
# pairwise concurvity
p_conc <- concrvity(krel_modI, pairwise = TRUE)
draw(p_conc)

plot.gamhp(gam.hp(krel_modI), plot.perc = TRUE)

# model selection -----
ms <- AIC(krel_modI,
          krel_modGI,
          krel_modS,
          krel_modGS,
          krel_modG) %>%
  arrange(AIC)
ms$deltaAIC <- ms$AIC - min(ms$AIC)
ms$wi <- exp(-0.5 * ms$deltaAIC)/sum(exp(-0.5 * ms$deltaAIC))
ms$er <- max(ms$wi)/ms$wi

ms.krel <- data.table(ms, keep.rownames = TRUE) %>%
  rename(model = rn)


# plot best model I ----
startyr <- min(dat.mod.krel$cohortyear)
lastyr <- max(dat.mod.krel$cohortyear)
wd <- 0.4
## effects plot ----
p.krel_effects_modI <- draw(krel_modI)
p.krel_effectscohortyear_modI <- p.krel_effects_modI[[1]] +
  scale_x_continuous(    breaks = seq(1980, 2020, 10),
                         minor_breaks = seq(1985, 2025, 10),
                         guide = guide_prism_minor()) +
  geom_hline(yintercept = 0) +
  xlab("Cohort year") +
  scale_x_continuous(breaks = min(dat.mod.krel$cohortyear, na.rm = TRUE):max(dat.mod.krel$cohortyear, na.rm = TRUE)) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        panel.grid = element_blank())

p.krel_effectscohortyear_FemMat_modI <- p.krel_effects_modI[[2]] +
  # scale_x_continuous(breaks = seq(startyr,lastyr,5)) +
  scale_colour_manual(name = '', values = viridis::viridis(3)[-3]) +
  geom_hline(yintercept = 0) +
  xlab("Female Maturity") +
  # ylab("PARTIAL EFFECT") +
  theme(legend.position = "bottom",
        panel.grid = element_blank())

## response plot -----
# setup prediction data
krel_modI_pred <- with(dat.mod.krel,
                       expand.grid(cohortyear = min(cohortyear):max(cohortyear),
                                   Female.Maturity = levels(Female.Maturity)))
# make the prediction, add this and a column of standard errors to the prediction data.frame.
krel_modI_pred <- cbind(krel_modI_pred,
                        predict(krel_modI,
                                krel_modI_pred,
                                se.fit = TRUE,
                                type = "response"))
# plot
p.krel_bestmodel_modI <- p.cond.femmat.data +
  geom_ribbon(aes(ymin = (fit - 2*se.fit),
                  ymax = (fit + 2*se.fit),
                  x = cohortyear
  ),
  color = "transparent",
  data = krel_modI_pred,
  alpha = 0.15,
  inherit.aes = FALSE) +
  geom_line(aes(y = fit),
            col = 'gray30',
            data = krel_modI_pred)

## Female maturity ----
# Get coefficients
coef_summary <- summary(krel_modI)

# Extract the smooths
ranef_vals <- coef(krel_modI)[grep("s\\(Female.Maturity\\)", names(coef(krel_modI)))]
ranef_vals

# estimated effect of maturity class (including both the random effect and the factor-smooth at the chosen cohortyear
newdat <- expand.grid(
  cohortyear = seq(1980, 2022, 1),
  Female.Maturity = levels(dat.mod.krel$Female.Maturity)
)

pred <- predict(krel_modI, newdat, se.fit = TRUE)
pred <- cbind(newdat, fit = pred$fit, se = pred$se.fit)
p.cond.femmat.effect <- ggplot(pred, aes(cohortyear, y = fit, color = factor(Female.Maturity),
                                         fill = factor(Female.Maturity))) +
  geom_hline(yintercept = 1) +
  geom_point(position = position_dodge2(width = 0.5)) +
  geom_line(position = position_dodge2(width = 0.5)) +
  geom_ribbon(aes(ymin = (fit - 2*se),
                  ymax = (fit + 2*se),
                  x = cohortyear),
              alpha = 0.15 ,
              color = "transparent"
  ) +
  xlab("Cohort Year") +
  ylab("Relative Condition (*K<sub>r</sub>*)") +
  scale_colour_manual(values = viridis::viridis(5)[-3]) +
  scale_fill_manual(values = viridis::viridis(5)[-3]) +
  scale_x_continuous(breaks = seq(1980, 2020, 10),
                     minor_breaks = seq(1985, 2015, 10),
                     guide = guide_prism_minor()) +
  scale_y_continuous(breaks = seq(0.8, 1.4, 0.10),
                     minor_breaks = seq(0.85, 1.35, 0.10),
                     guide = guide_prism_minor()) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.title.y = element_markdown()) +
  NULL



merge(pred, ao.seasonal, by.x = 'cohortyear', by.y = 'year', all.x = TRUE) %>%

  ggplot(., aes(cohortyear, y = fit, color = factor(Female.Maturity),
                fill = factor(Female.Maturity))) +
  geom_hline(yintercept = 1) +
  geom_line(aes(y = ao.seasonal + 1), col = 'black') +
  geom_point(position = position_dodge2(width = 0.5)) +
  geom_line(position = position_dodge2(width = 0.5)) +
  geom_ribbon(aes(ymin = (fit - 2*se),
                  ymax = (fit + 2*se),
                  x = cohortyear),
              alpha = 0.15 ,
              color = "transparent"
  ) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  xlab("Cohort Year") +
  scale_colour_manual(values = viridis::viridis(5)[-3]) +
  scale_fill_manual(values = viridis::viridis(5)[-3]) +
  scale_x_continuous(breaks = seq(1980, 2020, 10),
                     minor_breaks = seq(1985, 2015, 10),
                     guide = guide_prism_minor()) +
  scale_y_continuous(breaks = seq(0.8, 1.4, 0.10),
                     minor_breaks = seq(0.85, 1.35, 0.10),
                     guide = guide_prism_minor()) +
  ylab('Relative Condition')




merge(pred, ice, by.x = 'cohortyear', by.y = 'year', all.x = TRUE) %>%

  ggplot(., aes(first_year_ice, y = fit, color = factor(Female.Maturity),
                fill = factor(Female.Maturity))) +
  # xlim(-1,1) +
  # geom_hline(yintercept = 1) +
  # geom_line(aes(y = biomass_tonnes ), col = 'black') +
  geom_point(position = position_dodge2(width = 0.5)) +
  geom_smooth() +
  # geom_line(position = position_dodge2(width = 0.5)) +
  # geom_ribbon(aes(ymin = (fit - 2*se),
  #                 ymax = (fit + 2*se),
  #                 x = cohortyear),
  #             alpha = 0.15 ,
  #             color = "transparent"
  # ) +
  theme(legend.position = 'bottom',
        legend.title = element_blank(),
        panel.grid = element_blank()) +
  # xlab("Cohort Year") +
  scale_colour_manual(values = viridis::viridis(5)[-3]) +
  scale_fill_manual(values = viridis::viridis(5)[-3]) +
  # scale_x_continuous(breaks = seq(1980, 2020, 10),
  #                    minor_breaks = seq(1985, 2015, 10),
  #                    guide = guide_prism_minor()) +
  # scale_y_continuous(breaks = seq(0.8, 1.4, 0.10),
  #                    minor_breaks = seq(0.85, 1.35, 0.10),
  #                    guide = guide_prism_minor()) +
  ylab('Relative Condition')

## deviance partition -----
p.krel.partition <- plot.gamhp(gam.hp(krel_modI), plot.perc = TRUE) +
  scale_y_continuous(    breaks = seq(0, 70, 10),
                         minor_breaks = seq(5, 65, 10),
                         guide = guide_prism_minor()) +
  scale_x_discrete(labels = c("s(Cohort year,\nby(Female Maturity))",
                              "s(Female Maturity),\nrandom effect")) +
  theme(axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"))

## predicted vs observed ----

krel.pred.obs <- krel_modI_pred %>%
  right_join(sum.relk)

p.krel.obs.pred <- ggplot(krel.pred.obs, aes(y = fit, x = meanK, color = Female.Maturity)) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_linerange(aes(xmin = krellb,
                     xmax = krelub),
                 alpha = 0.25) +
  geom_linerange(aes(ymin = fit - 2*se.fit,
                     ymax = fit + 2*se.fit),
                 alpha = 0.25) +
  geom_point() +
  scale_y_continuous( limits = c(0.81, 1.38),
                      breaks = seq(0.8, 1.4, 0.20),
                      minor_breaks = seq(0.9, 1.3, 0.20),
                      guide = guide_prism_minor()) +
  scale_x_continuous( limits = c(0.55, 1.68),
                      breaks = seq(0.6, 1.6, 0.20),
                      minor_breaks = seq(0.5, 1.7, 0.20),
                      guide = guide_prism_minor()) +
  scale_colour_manual(values = viridis::viridis(5)[-3]) +
  ylab("Predicted relative condition (*K<sub>r</sub>*)") +
  xlab("Observed relative condition (*K<sub>r</sub>*)")  +
  theme(panel.grid = element_blank(),
        legend.title = element_blank(),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown(),
        legend.position = 'bottom')

p.krel.obs.pred.facet <- p.krel.obs.pred +
  facet_grid(Female.Maturity ~ .) +
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(colour = 'black'),
        axis.title.y = element_markdown(),
        axis.title.x = element_markdown()) +
  NULL


# output ------
fwrite(x = ms,
       row.names = TRUE,
       file =  paste0(here::here(), "/output/Condition/Condition_ModelSelection_CohYear-FemMat.csv"))

ggsave(plot = p.krel.partition,
       filename = paste0(here::here(), "/output/Condition/Condition_Model_Partition.png"),
       height = 5,
       width = 5)

ggsave(plot = p.cond.femmat.effect,
       filename = paste0(here::here(), "/output/Condition/Condition_Model_CohYear-FemMat_PartialEffect.png"),
       height = 5,
       width = 13)

ggsave(plot = p.cond.femmat.data,
       filename = paste0(here::here(), "/output/Condition/Condition_CohYear-FemMat.png"),
       height = 8,
       width = 13)

ggsave(plot = p.krel_bestmodel_modI,
       filename = paste0(here::here(), "/output/Condition/Condition_Model_CohYear-FemMat_fit.png"),
       height = 8,
       width = 13)

ggsave(plot = p.krel.obs.pred,
       filename = paste0(here::here(), "/output/Condition/Condition_Model_Observed.png"),
       height = 5,
       width = 5)
ggsave(plot = p.krel.obs.pred.facet,
       filename = paste0(here::here(), "/output/Condition/Condition_Model_Observed_FemMat.png"),
       height = 8,
       width = 5)
