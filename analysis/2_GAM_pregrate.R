# libraries ------
library(mgcv)
library(gratia)
library(gam.hp)
library(plotly)
library(ggbeeswarm)
library(ggprism)
library(ggpubr)

# data -----
source( paste0(here::here(), "/analysis/0_seal_data.R"))
dat.biolrates <- fread(paste0(here::here(), "/data/seal/CleanDatasetForBiologicalRates.csv"))
bootCI <- fread(paste0(here::here(), "/output/BiologicalRates_bootstrap.csv"))

# functions ----
source( paste0(here::here(), "/R/y.transf.betareg.r"))


# pregnancy rate model ----
# change the pregrate in 1979 to the value used in the resdoc
# will circle back to this once I finish comparing rates
seal.data[cohortyear == 1979, pregrate := 0.941]

# change pregrate and abrate in 2012 to the values used in the resdoc
# will circle back to this once I finish comparing rates
seal.data[cohortyear == 2012, pregrate := 0.627]

## fit model ----
m.preg <- gam(pregrate ~ s(sealpop) +
                s(abrate) ,
              data = seal.data,
              method = "REML",
              family = betar())

## check model ----
appraise(m.preg)
k.check(m.preg)
concurvity(m.preg)
concurvity(m.preg, full = FALSE)
m.drop <- update(m.preg, . ~ . - s(abrate))

AIC(m.preg, m.drop)
draw(m.preg, residuals = FALSE)
draw(m.preg, residuals = TRUE)
summary(m.preg)

# plot ----
dev.expl <- round(summary(m.preg)$dev.expl, 2) * 100
label.preg.model <-   paste0("Pregnancy rate ~ s(Population size) + s(Abortion rate)\nExplained deviance: ",
                             dev.expl, "%")

## model fit -----
preg.fits <- predict(m.preg,
                     type = "response",
                     se.fit = TRUE) %>%
  data.frame() %>%
  mutate(lci.preg = fit - 1.96 * se.fit,
         uci.preg = fit + 1.96 * se.fit) %>%
  select(-se.fit) %>%
  rename(fit.preg = fit) %>%
  bind_cols(
    seal.data[!is.na(pregrate), .(cohortyear, sealpop, pregrate)])


preg.fits <- preg.fits %>%
  left_join(bootCI[, .(cohortyear, preglb, pregub)])

# preg.fits <- expand.grid(cohortyear = seq(min(preg.fits$cohortyear),
#                                           max(preg.fits$cohortyear),
#                                           1)) %>%
#   left_join(preg.fits)

p.preg <- ggplot(preg.fits, aes(x = cohortyear, y = pregrate)) +
  geom_ribbon(aes(y = fit.preg, ymin = lci.preg, ymax = uci.preg),
              alpha = 0.2) +
  geom_line(aes(y = fit.preg)) +
  geom_point() +
  geom_linerange(aes(ymin = preglb,
                     ymax = pregub),
                 alpha = 0.5, col = 'gray40') +
  scale_x_continuous(    breaks = seq(1950, 2020, 10),
                         minor_breaks = seq(1955, 2025, 10),
                         guide = guide_prism_minor()) +
  scale_y_continuous(    breaks = seq(0, 1, 0.10),
                         minor_breaks = seq(0.05, 0.95, 0.10),
                         guide = guide_prism_minor()) +
  xlab("Cohort year") +
  ylab("Pregnancy rate") +
  annotate("text", x = 1952, y = .3, label = label.preg.model,
           fontface = "plain",
           family = 'sans',
           size = 3.5,
           hjust = 0) +
  theme(panel.grid = element_blank())

## deviance partition -----
p.preg.partition <- plot.gamhp(gam.hp(m.preg), plot.perc = TRUE) +
  scale_y_continuous(    breaks = seq(0, 50, 10),
                         minor_breaks = seq(5, 55, 10),
                         guide = guide_prism_minor()) +
  scale_x_discrete(labels = c("s(Population size)",
                              "s(Abortion rate)")) +
  theme(axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "white"))

## predicted vs observed ----
p.preg.obs.pred <- ggplot(preg.fits, aes(y = fit.preg, x = pregrate)) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  geom_point() +
  geom_linerange(aes(xmin = preglb,
                     xmax = pregub),
                 alpha = 0.5, col = 'gray40') +
  geom_linerange(aes(ymin = lci.preg,
                     ymax = uci.preg),
                 alpha = 0.5, col = 'gray40') +
  scale_x_continuous( limits = c(0.1, 1),
                      breaks = seq(0, 1, 0.10),
                      minor_breaks = seq(0.05, 0.95, 0.10),
                      guide = guide_prism_minor()) +
  scale_y_continuous( limits = c(0.1, 1),
                      breaks = seq(0, 1, 0.10),
                      minor_breaks = seq(0.05, 0.95, 0.10),
                      guide = guide_prism_minor()) +
  ylab("Predicted pregnancy rate") +
  xlab("Observed pregnancy rate")  +
  theme(panel.grid = element_blank())

## Partial effects plot -----
# Extract smooth estimates
sealpop_eff <- smooth_estimates(m.preg, select = "s(sealpop)") %>%
  mutate(term = "Seal population") %>%
  rename(value = sealpop)
abrate_eff  <- smooth_estimates(m.preg, select = "s(abrate)") %>%
  mutate(term = "Abortion rate") %>%
  rename(value = abrate)

# Combine and plot
effs <- rbind(sealpop_eff, abrate_eff) %>%
  data.table()

res <- rbindlist(l = list(
  data.table(m.preg$model %>%
               select(-sealpop) %>%
               rename(value = abrate) %>%
               mutate(term = "Abortion rate"),
             res = m.preg$residuals),
  data.table(m.preg$model %>%
               select(-abrate) %>%
               rename(value = sealpop) %>%
               mutate(term = "Seal population"),
             res = m.preg$residuals)
))

p.preg.partialeffects <- ggplot(effs, aes(x = value )) +
  geom_line( aes(y = .estimate) )+
  geom_ribbon(alpha = 0.2, aes(ymin = .estimate - .se, ymax = .estimate + .se)) +
  geom_point(data = res, aes(y = res)) +
  facet_wrap(~term, scales = "free") +
  labs(y = "Partial effect", x = "Predictor value",
       title = "Partial smooth effects from GAM") +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(colour = 'black'))


p.resids <- draw(m.preg, residuals = TRUE)
p.resid.spop <- p.resids[[1]] +
  xlab("Population size") +
  geom_hline(yintercept = 0) +
  theme(panel.grid = element_blank())

p.resid.spop$layers[[1]]$aes_params$colour <- "gray30"
p.resid.spop$layers[[1]]$aes_params$alpha <- 0.5
p.resid.spop$labels$title <- "s(Population size)"

p.resid.abr <- p.resids[[2]] +
  xlab("Abortion rate") +
  ylab("") +
  geom_hline(yintercept = 0) +
  theme(panel.grid = element_blank())

p.resid.abr$layers[[1]]$aes_params$colour <- "gray30"
p.resid.abr$layers[[1]]$aes_params$alpha <- 0.5
p.resid.abr$labels$title <- "s(Abortion rate)"

p.preg.partialeffects <- ggpubr::ggarrange(p.resid.spop, p.resid.abr,
                  ncol = 2
)


# output ----
ggsave(plot = p.preg.obs.pred,
       filename = paste0(here::here(), "/output/Pregnancy_Model_Observed.png"),
       height = 5,
       width = 5)
ggsave(plot = p.preg,
       filename = paste0(here::here(), "/output/Pregnancy_Model.png"),
       height = 5,
       width = 13)
ggsave(plot = p.preg.partition,
       filename = paste0(here::here(), "/output/Pregnancy_Model_Partition.png"),
       height = 5,
       width = 5)
ggsave(plot = p.preg.partialeffects,
       filename = paste0(here::here(), "/output/Pregnancy_Model_PartialEffects.png"),
       height = 5,
       width = 8)
