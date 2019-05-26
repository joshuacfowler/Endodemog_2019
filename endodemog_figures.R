## Title: Grass endophyte vital rate variance with a bayesian framework
## Purpose: Visualizes output of endophyte effect on year variance from 
## model outputs for Survival, Growth and Flowering models 
## Authors: Joshua and Tom

library(tidyverse)
library(rstan)
library(StanHeaders)
library(shinystan)
library(bayesplot)
library(devtools)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

#  source in raw output from models
## Read in the survival model output for all species
survPOAL <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_POAL_withplot.rds")
survPOSY <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_POSY_withplot.rds")
survLOAR <- readRDS(file = "model_run_MAR7/endodemog_surv_LOAR_withplot.rds")
survELVI <- readRDS(file = "model_run_MAR7/endodemog_surv_ELVI_withplot.rds")
survELRI <- readRDS(file = "model_run_MAR7/endodemog_surv_ELRI_withplot.rds")
survFESU <- readRDS(file = "model_run_MAR7/endodemog_surv_FESU_withplot.rds")
survAGPE <- readRDS(file = "model_run_MAR7/endodemog_surv_AGPE_withplot.rds")


## Read in the flower model output for all species
flwPOAL <- readRDS(file = "model_run_MAR7/endodemog_flw_POAL_withplot.rds")
flwPOSY <- readRDS(file = "model_run_MAR7/endodemog_flw_POSY_withplot.rds")
flwLOAR <- readRDS(file = "model_run_MAR7/endodemog_flw_LOAR_withplot.rds")
flwELVI <- readRDS(file = "model_run_MAR7/endodemog_flw_ELVI_withplot.rds")
flwELRI <- readRDS(file = "model_run_MAR7/endodemog_flw_ELRI_withplot.rds")
flwFESU <- readRDS(file = "model_run_MAR7/endodemog_flw_FESU_withplot.rds")
flwAGPE <- readRDS(file = "model_run_MAR7/endodemog_flw_AGPE_withplot.rds")



## Read in the growth model output for all species
growPOAL <- readRDS(file = "model_run_MAR7/endodemog_grow_POAL_withplot.rds")
growPOSY <- readRDS(file = "model_run_MAR7/endodemog_grow_POSY_withplot.rds")
growLOAR <- readRDS(file = "model_run_MAR7/endodemog_grow_LOAR_withplot.rds")
growELVI <- readRDS(file = "model_run_MAR7/endodemog_grow_ELVI_withplot.rds")
growELRI <- readRDS(file = "model_run_MAR7/endodemog_grow_ELRI_withplot.rds")
growFESU <- readRDS(file = "model_run_MAR7/endodemog_grow_FESU_withplot.rds")
growAGPE <- readRDS(file = "model_run_MAR7/endodemog_grow_AGPE_withplot.rds")

# save posteriors and summaries within dataframes
params = c("beta", "tau_year", "sigma_e", "tau_plot")

post_survPOAL <- as.data.frame(survPOAL, pars = params)
post_survPOSY <- as.data.frame(survPOSY, pars = params)
post_survLOAR <- as.data.frame(survLOAR, pars = params)
post_survELVI <- as.data.frame(survELVI, pars = params)
post_survELRI <- as.data.frame(survELRI, pars = params)
post_survFESU <- as.data.frame(survFESU, pars = params)
post_survAGPE <- as.data.frame(survAGPE, pars = params)

sum_survPOAL <- summary(survPOAL, pars = params)
sum_survPOSY <- summary(survPOSY, pars = params)
sum_survLOAR <- summary(survLOAR, pars = params)
sum_survELVI <- summary(survELVI, pars = params)
sum_survELRI <- summary(survELRI, pars = params)
sum_survFESU <- summary(survFESU, pars = params)
sum_survAGPE <- summary(survAGPE, pars = params)




post_flwPOAL <- as.data.frame(flwPOAL, pars = params)
post_flwPOSY <- as.data.frame(flwPOSY, pars = params)
post_flwLOAR <- as.data.frame(flwLOAR, pars = params)
post_flwELVI <- as.data.frame(flwELVI, pars = params)
post_flwELRI <- as.data.frame(flwELRI, pars = params)
post_flwFESU <- as.data.frame(flwFESU, pars = params)
post_flwAGPE <- as.data.frame(flwAGPE, pars = params)

sum_flwPOAL <- summary(flwPOAL, pars = params)
sum_flwPOSY <- summary(flwPOSY, pars = params)
sum_flwLOAR <- summary(flwLOAR, pars = params)
sum_flwELVI <- summary(flwELVI, pars = params)
sum_flwELRI <- summary(flwELRI, pars = params)
sum_flwFESU <- summary(flwFESU, pars = params)
sum_flwAGPE <- summary(flwAGPE, pars = params)



post_growPOAL <- as.data.frame(growPOAL, pars = params)
post_growPOSY <- as.data.frame(growPOSY, pars = params)
post_growLOAR <- as.data.frame(growLOAR, pars = params)
post_growELVI <- as.data.frame(growELVI, pars = params)
post_growELRI <- as.data.frame(growELRI, pars = params)
post_growFESU <- as.data.frame(growFESU, pars = params)
post_growAGPE <- as.data.frame(growAGPE, pars = params)

sum_growPOAL <- summary(growPOAL, pars = params)
sum_growPOSY <- summary(growPOSY, pars = params)
sum_growLOAR <- summary(growLOAR, pars = params)
sum_growELVI <- summary(growELVI, pars = params)
sum_growELRI <- summary(growELRI, pars = params)
sum_growFESU <- summary(growFESU, pars = params)
sum_growAGPE <- summary(growAGPE, pars = params)

#Looking at different color palettes
# Classic palette BuPu, with 4 colors
purps = brewer.pal(9, "Purples")
bupu = brewer.pal(9,"BuPu")
oran = brewer.pal(9, "Oranges")


brewer.pal(9, "Oranges")
orrd = brewer.pal(9, "OrRd")
# I can add more tones to this palette :
purps = colorRampPalette(purps)(11)
# pie(rep(1, length(purps)), col = purps , main="")
bupu = colorRampPalette(bupu)(11)
# pie(rep(1, length(bupu)), col = bupu , main="")
oran = colorRampPalette(oran)(11)
# pie(rep(1, length(oran)), col = oran , main="")
orrd = colorRampPalette(orrd)(11)
# pie(rep(1, length(orrd)), col = orrd , main="")

yearcolors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
                "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a",
                "#ffff99")
colors2 <- c("#ff7f00","#6a3d9a") #"#ff7f00" = E-, "#e31a1c" = E+

# Variance posterior distributions by species

# survival
AGPE_s <- ggplot() +
  geom_histogram(data = post_survAGPE, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .4, bins = 200) +
  geom_histogram(data = post_survAGPE, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .4, bins = 200) +
  labs(x = expression(σ[year]), y = "Posterior density", title = "Variance Estimate") + theme_classic()

ELRI_s <- ggplot() +
  geom_histogram(data = post_survELRI, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .4, bins = 200) +
  geom_histogram(data = post_survELRI, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .4, bins = 200) +
  labs(x = expression(σ[year]), y = "Posterior density", title = "Variance Estimate") + theme_classic()

ELVI_s <- ggplot() +
  geom_histogram(data = post_survELVI, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .5, bins = 200) +
  geom_histogram(data = post_survELVI, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .5, bins = 200) +
  labs(x = expression(σ[year]), y = "Posterior density", title = "Variance Estimate") + theme_classic()

FESU_s <- ggplot() +
  geom_histogram(data = post_survFESU, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .5, bins = 200) +
  geom_histogram(data = post_survFESU, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .5, bins = 200) +
  labs(x = expression(σ[year]), y = "Posterior density", title = "Variance Estimate") + theme_classic()

LOAR_s <- ggplot() +
  geom_histogram(data = post_survLOAR, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .5, bins = 200) +
  geom_histogram(data = post_survLOAR, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .5, bins = 200) +
  labs(x = expression(σ[year]), y = "Posterior density", title = "Variance Estimate") + theme_classic()

POAL_s <- ggplot() +
  geom_histogram(data = post_survPOAL, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .4, bins = 200) +
  geom_histogram(data = post_survPOAL, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .4, bins = 200) +
  labs(x = expression(σ[year]), y = "Posterior density", title = "Variance Estimate") + theme_classic()


POSY_s <- ggplot() +
  geom_histogram(data = post_survPOSY, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .4, bins = 200) +
  geom_histogram(data = post_survPOSY, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .4, bins = 200) +
  labs(x = expression(σ[year]), y = "Posterior density", title = "Variance Estimate") + theme_classic()


# flowering
AGPE_f <- ggplot() +
  geom_histogram(data = post_flwAGPE, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .4, bins = 200) +
  geom_histogram(data = post_flwAGPE, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .4, bins = 200) +
  labs(x = expression(σ[year]), y = "", title = "AGPE") + theme_classic()

ELRI_f <- ggplot() +
  geom_histogram(data = post_flwELRI, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .4, bins = 200) +
  geom_histogram(data = post_flwELRI, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .4, bins = 200) +
  labs(x = expression(σ[year]), y = "", title = "ELRI") + theme_classic()

ELVI_f <- ggplot() +
  geom_histogram(data = post_flwELVI, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .5, bins = 200) +
  geom_histogram(data = post_flwELVI, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .5, bins = 200) +
  labs(x = expression(σ[year]), y = "", title = "ELVI") + theme_classic()

FESU_f <- ggplot() +
  geom_histogram(data = post_flwFESU, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .5, bins = 200) +
  geom_histogram(data = post_flwFESU, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .5, bins = 200) +
  labs(x = expression(σ[year]), y = "", title = "FESU") + theme_classic()

LOAR_f <- ggplot() +
  geom_histogram(data = post_flwLOAR, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .5, bins = 200) +
  geom_histogram(data = post_flwLOAR, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .5, bins = 200) +
  labs(x = expression(σ[year]), y = "", title = "LOAR") + theme_classic()

POAL_f <- ggplot() +
  geom_histogram(data = post_flwPOAL, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .4, bins = 200) +
  geom_histogram(data = post_flwPOAL, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .4, bins = 200) +
  labs(x = expression(σ[year]), y = "", title = "POAL") + theme_classic()

POSY_f <- ggplot() +
  geom_histogram(data = post_flwPOSY, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .4, bins = 200) +
  geom_histogram(data = post_flwPOSY, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .4, bins = 200) +
  labs(x = expression(σ[year]), y = "", title = "POSY") + theme_classic()


# growth
AGPE_g <- ggplot() +
  geom_histogram(data = post_growAGPE, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .4, bins = 200) +
  geom_histogram(data = post_growAGPE, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .4, bins = 200) +
  labs(x = expression(σ[year]), y = "", title = "AGPE") + theme_classic()

ELRI_g <- ggplot() +
  geom_histogram(data = post_growELRI, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .4, bins = 200) +
  geom_histogram(data = post_growELRI, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .4, bins = 200) +
  labs(x = expression(σ[year]), y = "", title = "ELRI") + theme_classic()

ELVI_g <- ggplot() +
  geom_histogram(data = post_growELVI, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .5, bins = 200) +
  geom_histogram(data = post_growELVI, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .5, bins = 200) +
  labs(x = expression(σ[year]), y = "", title = "ELVI") + theme_classic()

FESU_g <- ggplot() +
  geom_histogram(data = post_growFESU, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .5, bins = 200) +
  geom_histogram(data = post_growFESU, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .5, bins = 200) +
  labs(x = expression(σ[year]), y = "", title = "FESU") + theme_classic()

LOAR_g <- ggplot() +
  geom_histogram(data = post_growLOAR, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .5, bins = 200) +
  geom_histogram(data = post_growLOAR, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .5, bins = 200) +
  labs(x = expression(σ[year]), y = "", title = "LOAR") + theme_classic()

POAL_g <- ggplot() +
  geom_histogram(data = post_growPOAL, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .4, bins = 200) +
  geom_histogram(data = post_growPOAL, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .4, bins = 200) +
  labs(x = expression(σ[year]), y = "", title = "POAL") + theme_classic()

POSY_g <- ggplot() +
  geom_histogram(data = post_growPOSY, aes(`sigma_e[1]`), fill = "#ff7f00", alpha = .4, bins = 200) +
  geom_histogram(data = post_growPOSY, aes(`sigma_e[2]`), fill = "#6a3d9a", alpha = .4, bins = 200) +
  labs(x = expression(σ[year]), y = "", title = "POSY") + theme_classic()



survvar <- grid.arrange(AGPE_s, ELRI_s, ELVI_s, FESU_s, LOAR_s, POAL_s, POSY_s,  ncol= 1)
titlesurv <- annotate_figure(survvar, top = "Survival", left = "")
flwrvar <- grid.arrange(AGPE_f, ELRI_f, ELVI_f, FESU_f, LOAR_f, POAL_f, POSY_f, ncol= 1)
titleflw <- annotate_figure(flwrvar, top = "Flowering", left = "")
growvar <- grid.arrange(AGPE_g, ELRI_g, ELVI_g, FESU_g, LOAR_g, POAL_g, POSY_g,  ncol= 1)
titlegrow <- annotate_figure(growvar, top = "Growth", left = "")
var <- grid.arrange(titlesurv, titleflw, titlegrow, ncol = 3)
titlevar <- annotate_figure(var, top = "Interannual Variance", left = "Posterior Density")

titlevar


# boxplots of variance estimates
AGPE_melt_surv <- post_survAGPE %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "AGPE") %>% 
  melt()
sumAGPE_melt_surv <- rownames_to_column(as.data.frame(sum_survAGPE$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "AGPE")

ELVI_melt_surv <- post_survELVI %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "ELVI") %>% 
  melt()
sumELVI_melt_surv <- rownames_to_column(as.data.frame(sum_survELVI$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "ELVI")

ELRI_melt_surv <- post_survELRI %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "ELRI") %>% 
  melt()
sumELRI_melt_surv <- rownames_to_column(as.data.frame(sum_survELRI$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "ELRI")

FESU_melt_surv <- post_survFESU %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "FESU") %>% 
  melt()
sumFESU_melt_surv <- rownames_to_column(as.data.frame(sum_survFESU$summary)) %>%
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "FESU")

LOAR_melt_surv <- post_survLOAR %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "LOAR") %>% 
 melt()
sumLOAR_melt_surv <- rownames_to_column(as.data.frame(sum_survLOAR$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "LOAR")

POAL_melt_surv <- post_survPOAL %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "POAL") %>% 
  melt()
sumPOAL_melt_surv <- rownames_to_column(as.data.frame(sum_survPOAL$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "POAL")

POSY_melt_surv <- post_survPOSY %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "POSY") %>% 
  melt()
sumPOSY_melt_surv <- rownames_to_column(as.data.frame(sum_survPOSY$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "POSY")

surv_species <- AGPE_melt_surv %>% 
  merge(ELRI_melt_surv, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(ELVI_melt_surv, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(FESU_melt_surv, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(LOAR_melt_surv, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(POAL_melt_surv, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(POSY_melt_surv, by = c("species", "value", "variable"), all = TRUE)

sum_surv_species <- sumAGPE_melt_surv %>% 
  merge(sumELRI_melt_surv, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumELVI_melt_surv, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumFESU_melt_surv, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumLOAR_melt_surv, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumPOAL_melt_surv, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumPOSY_melt_surv, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE)

AGPE_melt_flw <- post_flwAGPE %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "AGPE") %>% 
  melt()
sumAGPE_melt_flw <- rownames_to_column(as.data.frame(sum_flwAGPE$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "AGPE")

ELVI_melt_flw <- post_flwELVI %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "ELVI") %>% 
  melt()
sumELVI_melt_flw <- rownames_to_column(as.data.frame(sum_flwELVI$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "ELVI")

ELRI_melt_flw <- post_flwELRI %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "ELRI") %>% 
  melt()
sumELRI_melt_flw <- rownames_to_column(as.data.frame(sum_flwELRI$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "ELRI")

FESU_melt_flw <- post_flwFESU %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "FESU") %>% 
  melt()
sumFESU_melt_flw <- rownames_to_column(as.data.frame(sum_flwFESU$summary)) %>%
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "FESU")

LOAR_melt_flw <- post_flwLOAR %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "LOAR") %>% 
  melt()
sumLOAR_melt_flw <- rownames_to_column(as.data.frame(sum_flwLOAR$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "LOAR")
POAL_melt_flw <- post_flwPOAL %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "POAL") %>% 
  melt()
sumPOAL_melt_flw <- rownames_to_column(as.data.frame(sum_flwPOAL$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "POAL")
POSY_melt_flw <- post_flwPOSY %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "POSY") %>% 
  melt()
sumPOSY_melt_flw <- rownames_to_column(as.data.frame(sum_flwPOSY$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "POSY")

flw_species <- AGPE_melt_flw %>% 
  merge(ELRI_melt_flw, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(ELVI_melt_flw, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(FESU_melt_flw, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(LOAR_melt_flw, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(POAL_melt_flw, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(POSY_melt_flw, by = c("species", "value", "variable"), all = TRUE)

sum_flw_species <- sumAGPE_melt_flw %>% 
  merge(sumELRI_melt_flw, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumELVI_melt_flw, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumFESU_melt_flw, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumLOAR_melt_flw, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumPOAL_melt_flw, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumPOSY_melt_flw, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE)

AGPE_melt_grow <- post_growAGPE %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "AGPE") %>% 
  melt()
sumAGPE_melt_grow <- rownames_to_column(as.data.frame(sum_growAGPE$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "AGPE")

ELVI_melt_grow <- post_growELVI %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "ELVI") %>% 
  melt()
sumELVI_melt_grow <- rownames_to_column(as.data.frame(sum_growELVI$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "ELVI")

ELRI_melt_grow <- post_growELRI %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "ELRI") %>% 
  melt()
sumELRI_melt_grow <- rownames_to_column(as.data.frame(sum_growELRI$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "ELRI")

FESU_melt_grow <- post_growFESU %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "FESU") %>% 
  melt()
sumFESU_melt_grow <- rownames_to_column(as.data.frame(sum_growFESU$summary)) %>%
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "FESU")

LOAR_melt_grow <- post_growLOAR %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "LOAR") %>% 
  melt()
sumLOAR_melt_grow <- rownames_to_column(as.data.frame(sum_growLOAR$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "LOAR")
POAL_melt_grow <- post_growPOAL %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "POAL") %>% 
  melt()
sumPOAL_melt_grow <- rownames_to_column(as.data.frame(sum_growPOAL$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "POAL")
POSY_melt_grow <- post_growPOSY %>% 
  select("sigma_e[1]", "sigma_e[2]") %>%
  mutate(species = "POSY") %>% 
  melt()
sumPOSY_melt_grow <- rownames_to_column(as.data.frame(sum_growPOSY$summary)) %>% 
  select("rowname", "mean", "2.5%", "97.5%") %>% 
  filter(rowname == "sigma_e[1]" | rowname == "sigma_e[2]") %>% 
  mutate(species = "POSY")

grow_species <- AGPE_melt_grow %>% 
  merge(ELRI_melt_grow, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(ELVI_melt_grow, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(FESU_melt_grow, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(LOAR_melt_grow, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(POAL_melt_grow, by = c("species", "value", "variable"), all = TRUE) %>% 
  merge(POSY_melt_grow, by = c("species", "value", "variable"), all = TRUE)

sum_grow_species <- sumAGPE_melt_grow %>% 
  merge(sumELRI_melt_grow, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumELVI_melt_grow, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumFESU_melt_grow, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumLOAR_melt_grow, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumPOAL_melt_grow, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE) %>% 
  merge(sumPOSY_melt_grow, by = c("rowname", "species", "mean", "2.5%", "97.5%"), all = TRUE)

# options

ggplot() +
  geom_boxplot(data = surv_species, aes(y = value, x = variable, fill = variable)) +
  facet_wrap(~species, nrow = 2) + 
  labs(x = "species", y = expression(σ[year]), title = "POAL Survival Variance") + theme_classic() +
  scale_fill_manual(values = colors2)

ggplot() +
  geom_boxplot(data = surv_species, aes(y = value, fill = variable)) +
  facet_grid(row = vars(species), col = vars(variable)) + 
  labs(x = "species", y = expression(σ[year]), title = "POAL Survival Variance") + theme_classic() +
  scale_fill_manual(values = colors2)

ggplot(data = sum_surv_species) +
  geom_point(aes(y = mean, x = species, color = rowname)) +
  geom_errorbar(aes(x = species, ymin = `2.5%`, ymax = `97.5%`, color = rowname), width = .5) +
  labs(x = "species", y = expression(σ[year]), title = "POAL Survival Variance") + theme_classic() +
  scale_color_manual(values = colors2)

# means with 90%CI
surv <- ggplot(data = sum_surv_species) +
  geom_point(aes(y = mean, x = species, color = rowname), position = position_dodge( width = .7), size = 2) + 
  geom_errorbar(aes(x = species, ymin = `2.5%`, ymax = `97.5%`, color = rowname), width = 0, position = position_dodge(width = .7), lwd = 1, alpha = .5) +
  labs(x = "", y = expression(σ[year]), title = "Survival") + theme_classic() +
  scale_color_manual(values = colors2) 

flw <- ggplot(data = sum_flw_species) +
  geom_point(aes(y = mean, x = species, color = rowname), position = position_dodge( width = .7), size = 2) + 
  geom_errorbar(aes(x = species, ymin = `2.5%`, ymax = `97.5%`, color = rowname), width = 0, position = position_dodge(width = .7), lwd = 1, alpha = .5) +
  labs(x = "", y = expression(σ[year]), title = "Flowering") + theme_classic() +
  scale_color_manual(values = colors2) 

grow <- ggplot(data = sum_grow_species) +
  geom_point(aes(y = mean, x = species, color = rowname), position = position_dodge( width = .7), size = 2) + 
  geom_errorbar(aes(x = species, ymin = `2.5%`, ymax = `97.5%`, color = rowname), width = 0, position = position_dodge(width = .7), lwd = 1, alpha = .5) +
  labs(x = "", y = expression(σ[year]), title = "Growth") + theme_classic() +
  scale_color_manual(values = colors2) 

mCIplot <- grid.arrange(surv, grow, nrow= 2)
titlemCIplot <- annotate_figure(mCIplot, top = "", bottom = "Species", left = "")

titlemCIplot




###### plots of Survival model fits#######

######## POAL model fits
# actual data points for probability of survival
POALsurv_bin0 <- POAL_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  filter(endo == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())

POALsurv_bin1 <- POAL_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  filter(endo == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())
POALsurv_bin <- POAL_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())


POALsurv_binmean <- POAL_data %>% 
  select(surv_t1, logsize_t, endo) %>%
  mutate(endo = recode(endo, "0" = "E-", "1" = "E+", minus = "E-", plus = "E+")) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(POAL_data$logsize_t)
xdummy <- seq(min(POAL_data$logsize_t), max(POAL_data$logsize_t), length.out = nvalues)

# E- 
POAL_ydummy_eminus <- as.vector(invlogit(sum_survPOAL$summary["beta[1]","mean"] + (sum_survPOAL$summary["beta[2]","mean"])*xdummy + (sum_survPOAL$summary["beta[3]","mean"])*0 + (sum_survPOAL$summary["beta[4]","mean"])*mean(POAL_data$origin_01)
                    + (sum_survPOAL$summary["beta[5]","mean"])*xdummy*0))
POAL_ydummy_eminus2.5 <- as.vector(invlogit(sum_survPOAL$summary["beta[1]","2.5%"] + (sum_survPOAL$summary["beta[2]","2.5%"])*xdummy + (sum_survPOAL$summary["beta[3]","2.5%"])*0 + (sum_survPOAL$summary["beta[4]","2.5%"])*mean(POAL_data$origin_01)
                                         + (sum_survPOAL$summary["beta[5]","2.5%"])*xdummy*0))
POAL_ydummy_eminus97.5 <- as.vector(invlogit(sum_survPOAL$summary["beta[1]","97.5%"] + (sum_survPOAL$summary["beta[2]","97.5%"])*xdummy + (sum_survPOAL$summary["beta[3]","97.5%"])*0 + (sum_survPOAL$summary["beta[4]","97.5%"])*mean(POAL_data$origin_01)
                                         + (sum_survPOAL$summary["beta[5]","97.5%"])*xdummy*0))

POAL_ydummy_eminus_semin <- as.vector(invlogit((sum_survPOAL$summary["beta[1]","mean"]-sum_survPOAL$summary["beta[1]","se_mean"]) + (sum_survPOAL$summary["beta[2]","mean"]-sum_survPOAL$summary["beta[2]","se_mean"])*xdummy + (sum_survPOAL$summary["beta[3]","mean"]-sum_survPOAL$summary["beta[3]","se_mean"])*0 + (sum_survPOAL$summary["beta[4]","mean"]-sum_survPOAL$summary["beta[4]","se_mean"])*mean(POAL_data$origin_01)
                                              + (sum_survPOAL$summary["beta[5]","mean"]-sum_survPOAL$summary["beta[5]","se_mean"])*xdummy*0))
POAL_ydummy_eminus_semax <- as.vector(invlogit((sum_survPOAL$summary["beta[1]","mean"]+sum_survPOAL$summary["beta[1]","se_mean"]) + (sum_survPOAL$summary["beta[2]","mean"]+sum_survPOAL$summary["beta[2]","se_mean"])*xdummy + (sum_survPOAL$summary["beta[3]","mean"]+sum_survPOAL$summary["beta[3]","se_mean"])*0 + (sum_survPOAL$summary["beta[4]","mean"]+sum_survPOAL$summary["beta[4]","se_mean"])*mean(POAL_data$origin_01)
                                              + (sum_survPOAL$summary["beta[5]","mean"]+sum_survPOAL$summary["beta[5]","se_mean"])*xdummy*0))


POAL_ydummy_eminus_y1 <- as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*0 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01)
                    + mean(post_survPOAL$'beta[5]')*xdummy*0
                    + mean(post_survPOAL$'tau_year[1,1]')))
POAL_ydummy_eminus_y2 <- as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*0 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01)
                    + mean(post_survPOAL$'beta[5]')*xdummy*0
                    + mean(post_survPOAL$'tau_year[1,2]')))
POAL_ydummy_eminus_y3 <- as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + + mean(post_survPOAL$'beta[3]')*0  + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01)
                    + mean(post_survPOAL$'beta[5]')*xdummy*0
                    + mean(post_survPOAL$'tau_year[1,3]')))
POAL_ydummy_eminus_y4 <- as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*0 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01)
                    + mean(post_survPOAL$'beta[5]')*xdummy*0
                    + mean(post_survPOAL$'tau_year[1,4]'))) 
POAL_ydummy_eminus_y5 <- as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*0 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01)
                    + mean(post_survPOAL$'beta[5]')*xdummy*0
                    + mean(post_survPOAL$'tau_year[1,5]')))
POAL_ydummy_eminus_y6 <- as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*0 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01)
                    + mean(post_survPOAL$'beta[5]')*xdummy*0
                    + mean(post_survPOAL$'tau_year[1,6]')))
POAL_ydummy_eminus_y7 <- as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*0 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01)
                    + mean(post_survPOAL$'beta[5]')*xdummy*0
                    + mean(post_survPOAL$'tau_year[1,7]')))
POAL_ydummy_eminus_y8 <- as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*0 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01)
                    + mean(post_survPOAL$'beta[5]')*xdummy*0
                    + mean(post_survPOAL$'tau_year[1,8]')))
POAL_ydummy_eminus_y9 <- as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*0 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01)
                    + mean(post_survPOAL$'beta[5]')*xdummy*0
                    + mean(post_survPOAL$'tau_year[1,9]')))
POAL_ydummy_eminus_y10 <- as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*0 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01)
                    + mean(post_survPOAL$'beta[5]')*xdummy*0
                    + mean(post_survPOAL$'tau_year[1,10]')))
POAL_ydummy_eminus_y11 <- as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*0 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01)
                    + mean(post_survPOAL$'beta[5]')*xdummy*0
                    + mean(post_survPOAL$'tau_year[1,11]')))
POALfitsminus0 <- as.data.frame(cbind(xdummy, POAL_ydummy_eminus, POAL_ydummy_eminus2.5, POAL_ydummy_eminus97.5))
POALfitsminus1 <- as.data.frame(cbind(xdummy,POAL_ydummy_eminus_y1,POAL_ydummy_eminus_y2, POAL_ydummy_eminus_y3, POAL_ydummy_eminus_y4, POAL_ydummy_eminus_y5, POAL_ydummy_eminus_y6, POAL_ydummy_eminus_y7, POAL_ydummy_eminus_y8, POAL_ydummy_eminus_y9, POAL_ydummy_eminus_y10, POAL_ydummy_eminus_y11))
POALfitsminus <- melt(POALfitsminus1, id = "xdummy",
                      measure.vars = c("POAL_ydummy_eminus_y1","POAL_ydummy_eminus_y2", 
                                       "POAL_ydummy_eminus_y3", "POAL_ydummy_eminus_y4", 
                                       "POAL_ydummy_eminus_y5", "POAL_ydummy_eminus_y6", 
                                       "POAL_ydummy_eminus_y7", "POAL_ydummy_eminus_y8", 
                                       "POAL_ydummy_eminus_y9", "POAL_ydummy_eminus_y10",
                                       "POAL_ydummy_eminus_y11")) %>% 
                 mutate(Year = recode(variable, "POAL_ydummy_eminus_y1" = '1',"POAL_ydummy_eminus_y2" = "2", 
                                       "POAL_ydummy_eminus_y3" = "3", "POAL_ydummy_eminus_y4" = "4", 
                                       "POAL_ydummy_eminus_y5" = "5", "POAL_ydummy_eminus_y6" = "6", 
                                       "POAL_ydummy_eminus_y7" = "7", "POAL_ydummy_eminus_y8" = "8", 
                                       "POAL_ydummy_eminus_y9" = "9", "POAL_ydummy_eminus_y10" = "10",
                                       "POAL_ydummy_eminus_y11" = "11") )



# E+
POAL_ydummy_eplus <- as.vector(invlogit(sum_survPOAL$summary["beta[1]","mean"] + (sum_survPOAL$summary["beta[2]","mean"])*xdummy + (sum_survPOAL$summary["beta[3]","mean"])*1 + (sum_survPOAL$summary["beta[4]","mean"])*mean(POAL_data$origin_01)
                                         + (sum_survPOAL$summary["beta[5]","mean"])*xdummy*1))
POAL_ydummy_eplus2.5 <- as.vector(invlogit(sum_survPOAL$summary["beta[1]","2.5%"] + (sum_survPOAL$summary["beta[2]","2.5%"])*xdummy + (sum_survPOAL$summary["beta[3]","2.5%"])*1 + (sum_survPOAL$summary["beta[4]","2.5%"])*mean(POAL_data$origin_01)
                                            + (sum_survPOAL$summary["beta[5]","2.5%"])*xdummy*1))
POAL_ydummy_eplus97.5 <- as.vector(invlogit(sum_survPOAL$summary["beta[1]","97.5%"] + (sum_survPOAL$summary["beta[2]","97.5%"])*xdummy + (sum_survPOAL$summary["beta[3]","97.5%"])*1 + (sum_survPOAL$summary["beta[4]","97.5%"])*mean(POAL_data$origin_01)
                                             + (sum_survPOAL$summary["beta[5]","97.5%"])*xdummy*1))

POAL_ydummy_eplus_semin <- as.vector(invlogit((sum_survPOAL$summary["beta[1]","mean"]-sum_survPOAL$summary["beta[1]","se_mean"]) + (sum_survPOAL$summary["beta[2]","mean"]-sum_survPOAL$summary["beta[2]","se_mean"])*xdummy + (sum_survPOAL$summary["beta[3]","mean"]-sum_survPOAL$summary["beta[3]","se_mean"])*1 + (sum_survPOAL$summary["beta[4]","mean"]-sum_survPOAL$summary["beta[4]","se_mean"])*mean(POAL_data$origin_01)
                                           + (sum_survPOAL$summary["beta[5]","mean"]-sum_survPOAL$summary["beta[5]","se_mean"])*xdummy*1))
POAL_ydummy_eplus_semax <- as.vector(invlogit((sum_survPOAL$summary["beta[1]","mean"]+sum_survPOAL$summary["beta[1]","se_mean"]) + (sum_survPOAL$summary["beta[2]","mean"]+sum_survPOAL$summary["beta[2]","se_mean"])*xdummy + (sum_survPOAL$summary["beta[3]","mean"]+sum_survPOAL$summary["beta[3]","se_mean"])*1 + (sum_survPOAL$summary["beta[4]","mean"]+sum_survPOAL$summary["beta[4]","se_mean"])*mean(POAL_data$origin_01)
                                              + (sum_survPOAL$summary["beta[5]","mean"]+sum_survPOAL$summary["beta[5]","se_mean"])*xdummy*1))

POAL_ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*1 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                     + mean(post_survPOAL$'beta[5]')*xdummy*1
                     + mean(post_survPOAL$'tau_year[2,1]')))
POAL_ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*1 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                     + mean(post_survPOAL$'beta[5]')*xdummy*1
                     + mean(post_survPOAL$'tau_year[2,2')))
POAL_ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*1 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                     + mean(post_survPOAL$'beta[5]')*xdummy*1
                     + mean(post_survPOAL$'tau_year[2,3]')))
POAL_ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*1 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                    + mean(post_survPOAL$'beta[5]')*xdummy*1
                    + mean(post_survPOAL$'tau_year[2,4]')))
POAL_ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*1 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                    + mean(post_survPOAL$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                    + mean(post_survPOAL$'tau_year[2,5]')))
POAL_ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*1 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                    + mean(post_survPOAL$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                    + mean(post_survPOAL$'tau_year[2,6]')))
POAL_ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*1 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                    + mean(post_survPOAL$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                    + mean(post_survPOAL$'tau_year[2,7]')))
POAL_ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*1 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                    + mean(post_survPOAL$'beta[5]')*xdummy*1
                    + mean(post_survPOAL$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
POAL_ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*1 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                    + mean(post_survPOAL$'beta[5]')*xdummy*1
                    + mean(post_survPOAL$'tau_year[2,9]')))
POAL_ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*1 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                    + mean(post_survPOAL$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                    + mean(post_survPOAL$'tau_year[2,10]')))
POAL_ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_survPOAL$'beta[1]') + mean(post_survPOAL$'beta[2]')*xdummy + mean(post_survPOAL$'beta[3]')*1 + mean(post_survPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                    + mean(post_survPOAL$'beta[5]')*xdummy*1
                    + mean(post_survPOAL$'tau_year[2,11]')))
POAL_fitsplus0 <- as.data.frame(cbind(xdummy, POAL_ydummy_eplus, POAL_ydummy_eplus2.5, POAL_ydummy_eplus97.5))  
POALfitsplus1 <- as.data.frame(cbind(xdummy,POAL_ydummy_eplus_y1,POAL_ydummy_eplus_y2, POAL_ydummy_eplus_y3, POAL_ydummy_eplus_y4, POAL_ydummy_eplus_y5, POAL_ydummy_eplus_y6, POAL_ydummy_eplus_y7, POAL_ydummy_eplus_y8, POAL_ydummy_eplus_y9, POAL_ydummy_eplus_y10, POAL_ydummy_eplus_y11))
POALfitsplus <- melt(POALfitsplus1, id = "xdummy",
                                  measure.vars = c("POAL_ydummy_eplus_y1","POAL_ydummy_eplus_y2", 
                                                   "POAL_ydummy_eplus_y3", "POAL_ydummy_eplus_y4", 
                                                   "POAL_ydummy_eplus_y5", "POAL_ydummy_eplus_y6", 
                                                   "POAL_ydummy_eplus_y7", "POAL_ydummy_eplus_y8", 
                                                   "POAL_ydummy_eplus_y9", "POAL_ydummy_eplus_y10",
                                                   "POAL_ydummy_eplus_y11")) %>% 
                  mutate(Year = recode(variable,  "POAL_ydummy_eplus_y1" = '1',"POAL_ydummy_eplus_y2" = "2", 
                                                  "POAL_ydummy_eplus_y3" = "3", "POAL_ydummy_eplus_y4" = "4", 
                                                  "POAL_ydummy_eplus_y5" = "5", "POAL_ydummy_eplus_y6" = "6", 
                                                  "POAL_ydummy_eplus_y7" = "7", "POAL_ydummy_eplus_y8" = "8", 
                                                  "POAL_ydummy_eplus_y9" = "9", "POAL_ydummy_eplus_y10" = "10",
                                                  "POAL_ydummy_eplus_y11" = "11") )
# E+ by year
POALEplusbyyear <- ggplot(data = POALfitsplus) +
geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
geom_point(data = POALsurv_bin1, aes(x = mean_size, y = mean_surv, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
ggtitle("E+ ") + xlab("log(size_t)") + ylab("Prob. of Survival") +
scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = POALfitsplus) +
geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
ggtitle("E+") + xlab("log(size_t)") + ylab("Prob. of Survival")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
  scale_color_manual(values=yearcolors)

# E- by year
POALEminusbyyear <- ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = POALsurv_bin0, aes(x = mean_size, y = mean_surv, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("Prob. of Survival")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()
  
ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("POAL E- Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
  scale_color_manual(values=yearcolors)
  

# effect on mean


POALfits <- as.data.frame(cbind(xdummy, POAL_ydummy_eplus, POAL_ydummy_eplus_semin, POAL_ydummy_eplus_semax, POAL_ydummy_eplus2.5, POAL_ydummy_eplus97.5, POAL_ydummy_eminus, POAL_ydummy_eminus_semin, POAL_ydummy_eminus_semax, POAL_ydummy_eminus2.5, POAL_ydummy_eminus97.5))

# without point
ggplot(data = POALfits) +
  geom_ribbon(aes(x = xdummy, ymin = POAL_ydummy_eplus2.5, ymax = POAL_ydummy_eplus97.5), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = POAL_ydummy_eminus2.5, ymax = POAL_ydummy_eminus97.5), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus), color = "#ff7f00") +
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus), color = "#6a3d9a")+
  ggtitle("POAL Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +theme_classic()

ggplot(data = POALfits) +
  geom_ribbon(aes(x = xdummy, ymin = POAL_ydummy_eplus_semin, ymax = POAL_ydummy_eplus_semax), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = POAL_ydummy_eminus_semin, ymax = POAL_ydummy_eminus_semax), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus), color = "#ff7f00") +
  ggtitle("POAL Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +theme_classic()

# with points
POALmean <- ggplot(data = POALfits) +
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = POALsurv_binmean, aes(x = mean_size, y = mean_surv, color = endo, lwd = samplesize)) +
  ggtitle("Mean Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")


ggplot(data = POALfits) +
  geom_ribbon(aes(x = xdummy, ymin = POAL_ydummy_eplus_semin, ymax = POAL_ydummy_eplus_semax), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = POAL_ydummy_eminus_semin, ymax = POAL_ydummy_eminus_semax), fill = "#ff7f00", alpha = .2)+
   geom_line(aes(x = xdummy, y = POAL_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = POALsurv_binmean, aes(x = mean_size, y = mean_surv, color = endo, lwd = samplesize)) +
  ggtitle("POAL Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +
  scale_color_manual(values = colors2) +theme_classic() + guides(lwd = "none")

POALsurv4panel <- grid.arrange(POALmean, POAL_s,POALEplusbyyear,POALEminusbyyear, ncol= 4)
titlePOALsurv4panel <- annotate_figure(POALsurv4panel, top = "Poa alsodes")

titlePOALsurv4panel


# split up plus and minus with datapoints and add mean line
POAL_minus <- ggplot(data = POALfitsminus1)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = POALfitsminus0, aes(x = xdummy, y = POAL_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = POALsurv_bin0, aes(x = mean_size, y = mean_surv), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


POAL_plus <- ggplot(data = POALfitsplus1)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = POAL_fitsplus0, aes(x = xdummy, y = POAL_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = POALsurv_bin1, aes(x = mean_size, y = mean_surv), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

POALsurv <- grid.arrange(POAL_plus, POAL_minus, ncol= 2)
titlePOALsurv <- annotate_figure(POALsurv, top = "POAL", left = "Probability of Survival")

titlePOALsurv



# with data points
ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .7) +
  geom_point(data = POALsurv_bin0, aes(x = mean_size, y = mean_surv, color = Year)) +
  labs(title = "POAL E- Survival Probability", x = "log(size_t)", y = "Prob. of Survival") +
  scale_color_manual(values=yearcolors)
# without datapoints
ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .7) +
  labs(title = "POAL E- Survival Probability", x = "log(size_t)", y = "Prob. of Survival") +                                                                                                                                                                                                                                                                                                                                                                                                                                     scale_color_manual(values=yearcolors)







######## POSY model fits
# actual data points for probability of survival
POSYsurv_bin0 <- POSY_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  filter(endo == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())

POSYsurv_bin1 <- POSY_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  filter(endo == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())
POSYsurv_bin <- POSY_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())


POSYsurv_binmean <- POSY_data %>% 
  select(surv_t1, logsize_t, endo) %>%
  mutate(endo = recode(endo, "0" = "E-", "1" = "E+", minus = "E-", plus = "E+")) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(POSY_data$logsize_t)
xdummy <- seq(min(POSY_data$logsize_t), max(POSY_data$logsize_t), length.out = nvalues)

# E- 
POSY_ydummy_eminus <- as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*0 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                         + mean(post_survPOSY$'beta[5]')*xdummy*0))
POSY_ydummy_eminus2.5 <- as.vector(invlogit(sum_survPOSY$summary["beta[1]","2.5%"] + (sum_survPOSY$summary["beta[2]","2.5%"])*xdummy + (sum_survPOSY$summary["beta[3]","2.5%"])*0 + (sum_survPOSY$summary["beta[4]","2.5%"])*mean(POSY_data$origin_01)
                                            + (sum_survPOSY$summary["beta[5]","2.5%"])*xdummy*0))
POSY_ydummy_eminus97.5 <- as.vector(invlogit(sum_survPOSY$summary["beta[1]","97.5%"] + (sum_survPOSY$summary["beta[2]","97.5%"])*xdummy + (sum_survPOSY$summary["beta[3]","97.5%"])*0 + (sum_survPOSY$summary["beta[4]","97.5%"])*mean(POSY_data$origin_01)
                                             + (sum_survPOSY$summary["beta[5]","97.5%"])*xdummy*0))

POSY_ydummy_eminus_semin <- as.vector(invlogit((sum_survPOSY$summary["beta[1]","mean"]-sum_survPOSY$summary["beta[1]","se_mean"]) + (sum_survPOSY$summary["beta[2]","mean"]-sum_survPOSY$summary["beta[2]","se_mean"])*xdummy + (sum_survPOSY$summary["beta[3]","mean"]-sum_survPOSY$summary["beta[3]","se_mean"])*0 + (sum_survPOSY$summary["beta[4]","mean"]-sum_survPOSY$summary["beta[4]","se_mean"])*mean(POSY_data$origin_01)
                                               + (sum_survPOSY$summary["beta[5]","mean"]-sum_survPOSY$summary["beta[5]","se_mean"])*xdummy*0))
POSY_ydummy_eminus_semax <- as.vector(invlogit((sum_survPOSY$summary["beta[1]","mean"]+sum_survPOSY$summary["beta[1]","se_mean"]) + (sum_survPOSY$summary["beta[2]","mean"]+sum_survPOSY$summary["beta[2]","se_mean"])*xdummy + (sum_survPOSY$summary["beta[3]","mean"]+sum_survPOSY$summary["beta[3]","se_mean"])*0 + (sum_survPOSY$summary["beta[4]","mean"]+sum_survPOSY$summary["beta[4]","se_mean"])*mean(POSY_data$origin_01)
                                               + (sum_survPOSY$summary["beta[5]","mean"]+sum_survPOSY$summary["beta[5]","se_mean"])*xdummy*0))


POSY_ydummy_eminus_y1 <- as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*0 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_survPOSY$'beta[5]')*xdummy*0
                                            + mean(post_survPOSY$'tau_year[1,1]')))
POSY_ydummy_eminus_y2 <- as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*0 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_survPOSY$'beta[5]')*xdummy*0
                                            + mean(post_survPOSY$'tau_year[1,2]')))
POSY_ydummy_eminus_y3 <- as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + + mean(post_survPOSY$'beta[3]')*0  + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_survPOSY$'beta[5]')*xdummy*0
                                            + mean(post_survPOSY$'tau_year[1,3]')))
POSY_ydummy_eminus_y4 <- as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*0 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_survPOSY$'beta[5]')*xdummy*0
                                            + mean(post_survPOSY$'tau_year[1,4]'))) 
POSY_ydummy_eminus_y5 <- as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*0 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_survPOSY$'beta[5]')*xdummy*0
                                            + mean(post_survPOSY$'tau_year[1,5]')))
POSY_ydummy_eminus_y6 <- as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*0 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_survPOSY$'beta[5]')*xdummy*0
                                            + mean(post_survPOSY$'tau_year[1,6]')))
POSY_ydummy_eminus_y7 <- as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*0 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_survPOSY$'beta[5]')*xdummy*0
                                            + mean(post_survPOSY$'tau_year[1,7]')))
POSY_ydummy_eminus_y8 <- as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*0 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_survPOSY$'beta[5]')*xdummy*0
                                            + mean(post_survPOSY$'tau_year[1,8]')))
POSY_ydummy_eminus_y9 <- as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*0 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_survPOSY$'beta[5]')*xdummy*0
                                            + mean(post_survPOSY$'tau_year[1,9]')))
POSY_ydummy_eminus_y10 <- as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*0 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                             + mean(post_survPOSY$'beta[5]')*xdummy*0
                                             + mean(post_survPOSY$'tau_year[1,10]')))
POSY_ydummy_eminus_y11 <- as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*0 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                             + mean(post_survPOSY$'beta[5]')*xdummy*0
                                             + mean(post_survPOSY$'tau_year[1,11]')))
POSYfitsminus0 <- as.data.frame(cbind(xdummy, POSY_ydummy_eminus))
POSYfitsminus1 <- as.data.frame(cbind(xdummy,POSY_ydummy_eminus_y1,POSY_ydummy_eminus_y2, POSY_ydummy_eminus_y3, POSY_ydummy_eminus_y4, POSY_ydummy_eminus_y5, POSY_ydummy_eminus_y6, POSY_ydummy_eminus_y7, POSY_ydummy_eminus_y8, POSY_ydummy_eminus_y9, POSY_ydummy_eminus_y10, POSY_ydummy_eminus_y11))
POSYfitsminus <- melt(POSYfitsminus1, id = "xdummy",
                      measure.vars = c("POSY_ydummy_eminus_y1","POSY_ydummy_eminus_y2", 
                                       "POSY_ydummy_eminus_y3", "POSY_ydummy_eminus_y4", 
                                       "POSY_ydummy_eminus_y5", "POSY_ydummy_eminus_y6", 
                                       "POSY_ydummy_eminus_y7", "POSY_ydummy_eminus_y8", 
                                       "POSY_ydummy_eminus_y9", "POSY_ydummy_eminus_y10",
                                       "POSY_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "POSY_ydummy_eminus_y1" = '1',"POSY_ydummy_eminus_y2" = "2", 
                       "POSY_ydummy_eminus_y3" = "3", "POSY_ydummy_eminus_y4" = "4", 
                       "POSY_ydummy_eminus_y5" = "5", "POSY_ydummy_eminus_y6" = "6", 
                       "POSY_ydummy_eminus_y7" = "7", "POSY_ydummy_eminus_y8" = "8", 
                       "POSY_ydummy_eminus_y9" = "9", "POSY_ydummy_eminus_y10" = "10",
                       "POSY_ydummy_eminus_y11" = "11") )



# E+
POSY_ydummy_eplus <- as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*1 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                        + mean(post_survPOSY$'beta[5]')*xdummy*1))
POSY_ydummy_eplus2.5 <- as.vector(invlogit(sum_survPOSY$summary["beta[1]","2.5%"] + (sum_survPOSY$summary["beta[2]","2.5%"])*xdummy + (sum_survPOSY$summary["beta[3]","2.5%"])*1 + (sum_survPOSY$summary["beta[4]","2.5%"])*mean(POSY_data$origin_01)
                                            + (sum_survPOSY$summary["beta[5]","2.5%"])*xdummy*1))
POSY_ydummy_eplus97.5 <- as.vector(invlogit(sum_survPOSY$summary["beta[1]","97.5%"] + (sum_survPOSY$summary["beta[2]","97.5%"])*xdummy + (sum_survPOSY$summary["beta[3]","97.5%"])*1 + (sum_survPOSY$summary["beta[4]","97.5%"])*mean(POSY_data$origin_01)
                                             + (sum_survPOSY$summary["beta[5]","97.5%"])*xdummy*1))

POSY_ydummy_eplus_semin <- as.vector(invlogit((sum_survPOSY$summary["beta[1]","mean"]-sum_survPOSY$summary["beta[1]","se_mean"]) + (sum_survPOSY$summary["beta[2]","mean"]-sum_survPOSY$summary["beta[2]","se_mean"])*xdummy + (sum_survPOSY$summary["beta[3]","mean"]-sum_survPOSY$summary["beta[3]","se_mean"])*1 + (sum_survPOSY$summary["beta[4]","mean"]-sum_survPOSY$summary["beta[4]","se_mean"])*mean(POSY_data$origin_01)
                                               + (sum_survPOSY$summary["beta[5]","mean"]-sum_survPOSY$summary["beta[5]","se_mean"])*xdummy*1))
POSY_ydummy_eplus_semax <- as.vector(invlogit((sum_survPOSY$summary["beta[1]","mean"]+sum_survPOSY$summary["beta[1]","se_mean"]) + (sum_survPOSY$summary["beta[2]","mean"]+sum_survPOSY$summary["beta[2]","se_mean"])*xdummy + (sum_survPOSY$summary["beta[3]","mean"]+sum_survPOSY$summary["beta[3]","se_mean"])*1 + (sum_survPOSY$summary["beta[4]","mean"]+sum_survPOSY$summary["beta[4]","se_mean"])*mean(POSY_data$origin_01)
                                               + (sum_survPOSY$summary["beta[5]","mean"]+sum_survPOSY$summary["beta[5]","se_mean"])*xdummy*1))


POSY_ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*1 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_survPOSY$'beta[5]')*xdummy*1
                                            + mean(post_survPOSY$'tau_year[2,1]')))
POSY_ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*1 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_survPOSY$'beta[5]')*xdummy*1
                                            + mean(post_survPOSY$'tau_year[2,2')))
POSY_ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*1 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_survPOSY$'beta[5]')*xdummy*1
                                            + mean(post_survPOSY$'tau_year[2,3]')))
POSY_ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*1 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_survPOSY$'beta[5]')*xdummy*1
                                            + mean(post_survPOSY$'tau_year[2,4]')))
POSY_ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*1 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_survPOSY$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_survPOSY$'tau_year[2,5]')))
POSY_ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*1 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_survPOSY$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_survPOSY$'tau_year[2,6]')))
POSY_ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*1 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_survPOSY$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_survPOSY$'tau_year[2,7]')))
POSY_ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*1 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_survPOSY$'beta[5]')*xdummy*1
                                            + mean(post_survPOSY$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
POSY_ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*1 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_survPOSY$'beta[5]')*xdummy*1
                                            + mean(post_survPOSY$'tau_year[2,9]')))
POSY_ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*1 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                             + mean(post_survPOSY$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_survPOSY$'tau_year[2,10]')))
POSY_ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_survPOSY$'beta[1]') + mean(post_survPOSY$'beta[2]')*xdummy + mean(post_survPOSY$'beta[3]')*1 + mean(post_survPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                             + mean(post_survPOSY$'beta[5]')*xdummy*1
                                             + mean(post_survPOSY$'tau_year[2,11]')))
POSY_fitsplus0 <- as.data.frame(cbind(xdummy, POSY_ydummy_eplus))
POSYfitsplus1 <- as.data.frame(cbind(xdummy,POSY_ydummy_eplus_y1,POSY_ydummy_eplus_y2, POSY_ydummy_eplus_y3, POSY_ydummy_eplus_y4, POSY_ydummy_eplus_y5, POSY_ydummy_eplus_y6, POSY_ydummy_eplus_y7, POSY_ydummy_eplus_y8, POSY_ydummy_eplus_y9, POSY_ydummy_eplus_y10, POSY_ydummy_eplus_y11))
POSYfitsplus <- melt(POSYfitsplus1, id = "xdummy",
                     measure.vars = c("POSY_ydummy_eplus_y1","POSY_ydummy_eplus_y2", 
                                      "POSY_ydummy_eplus_y3", "POSY_ydummy_eplus_y4", 
                                      "POSY_ydummy_eplus_y5", "POSY_ydummy_eplus_y6", 
                                      "POSY_ydummy_eplus_y7", "POSY_ydummy_eplus_y8", 
                                      "POSY_ydummy_eplus_y9", "POSY_ydummy_eplus_y10",
                                      "POSY_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "POSY_ydummy_eplus_y1" = '1',"POSY_ydummy_eplus_y2" = "2", 
                       "POSY_ydummy_eplus_y3" = "3", "POSY_ydummy_eplus_y4" = "4", 
                       "POSY_ydummy_eplus_y5" = "5", "POSY_ydummy_eplus_y6" = "6", 
                       "POSY_ydummy_eplus_y7" = "7", "POSY_ydummy_eplus_y8" = "8", 
                       "POSY_ydummy_eplus_y9" = "9", "POSY_ydummy_eplus_y10" = "10",
                       "POSY_ydummy_eplus_y11" = "11") )

# E+ by year
POSYEplusbyyear <- ggplot(data = POSYfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = POSYsurv_bin1, aes(x = mean_size, y = mean_surv, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("Prob. of Survival") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = POSYfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("Prob. of Survival")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
POSYEminusbyyear <- ggplot(data = POSYfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = POSYsurv_bin0, aes(x = mean_size, y = mean_surv, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("Prob. of Survival")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = POSYfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("POSY E- Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)



# effect on mean

POSYfits <- as.data.frame(cbind(xdummy, POSY_ydummy_eplus, POSY_ydummy_eplus_semin, POSY_ydummy_eplus_semax, POSY_ydummy_eplus2.5, POSY_ydummy_eplus97.5, POSY_ydummy_eminus, POSY_ydummy_eminus_semin, POSY_ydummy_eminus_semax, POSY_ydummy_eminus2.5, POSY_ydummy_eminus97.5))


# without point
ggplot(data = POSYfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("POSY Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = POSYfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = POSYsurv_binmean, aes(x = mean_size, y = mean_surv, color = endo), lwd = 3) +
  ggtitle("POSY Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

POSYmean <- ggplot(data = POSYfits) +
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus), color = "#ff7f00")+ 
  geom_point(data = POSYsurv_binmean, aes(x = mean_size, y = mean_surv, color = endo, lwd = samplesize))+
  ggtitle("Mean Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")

POSYsurv4panel <- grid.arrange(POSYmean, POSY_s,POSYEplusbyyear,POSYEminusbyyear, ncol= 4)
titlePOSYsurv4panel <- annotate_figure(POSYsurv4panel, top = "Poa sylvestris")

titlePOSYsurv4panel

# split up plus and minus with datapoints and add mean line
POSY_minus <- ggplot(data = POSYfitsminus1)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = POSYfitsminus0, aes(x = xdummy, y = POSY_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = POSYsurv_bin0, aes(x = mean_size, y = mean_surv), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


POSY_plus <- ggplot(data = POSYfitsplus1)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = POSY_fitsplus0, aes(x = xdummy, y = POSY_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = POSYsurv_bin1, aes(x = mean_size, y = mean_surv), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

POSYsurv <- grid.arrange(POSY_plus, POSY_minus, ncol= 2)
titlePOSYsurv <- annotate_figure(POSYsurv, top = "POSY", left = "Probability of Survival")

titlePOSYsurv







######## LOAR model fits
# actual data points for probability of survival
LOARsurv_bin0 <- LOAR_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  filter(endo == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())


LOARsurv_bin1 <- LOAR_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  filter(logsize_t>=0) %>% 
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  filter(endo == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())
LOARsurv_bin <- LOAR_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  filter(logsize_t>=0) %>% 
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())


LOARsurv_binmean <- LOAR_data %>% 
  select(surv_t1, logsize_t, endo) %>%
  mutate(endo = recode(endo, "0" = "E-", "1" = "E+", minus = "E-", plus = "E+")) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(LOAR_data$logsize_t)
xdummy <- seq(min(LOAR_data$logsize_t), max(LOAR_data$logsize_t), length.out = nvalues)

# E- 
LOAR_ydummy_eminus <- as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*0 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                         + mean(post_survLOAR$'beta[5]')*xdummy*0))
LOAR_ydummy_eminus2.5 <- as.vector(invlogit(sum_survLOAR$summary["beta[1]","2.5%"] + (sum_survLOAR$summary["beta[2]","2.5%"])*xdummy + (sum_survLOAR$summary["beta[3]","2.5%"])*0 + (sum_survLOAR$summary["beta[4]","2.5%"])*mean(LOAR_data$origin_01)
                                            + (sum_survLOAR$summary["beta[5]","2.5%"])*xdummy*0))
LOAR_ydummy_eminus97.5 <- as.vector(invlogit(sum_survLOAR$summary["beta[1]","97.5%"] + (sum_survLOAR$summary["beta[2]","97.5%"])*xdummy + (sum_survLOAR$summary["beta[3]","97.5%"])*0 + (sum_survLOAR$summary["beta[4]","97.5%"])*mean(LOAR_data$origin_01)
                                             + (sum_survLOAR$summary["beta[5]","97.5%"])*xdummy*0))

LOAR_ydummy_eminus_semin <- as.vector(invlogit((sum_survLOAR$summary["beta[1]","mean"]-sum_survLOAR$summary["beta[1]","se_mean"]) + (sum_survLOAR$summary["beta[2]","mean"]-sum_survLOAR$summary["beta[2]","se_mean"])*xdummy + (sum_survLOAR$summary["beta[3]","mean"]-sum_survLOAR$summary["beta[3]","se_mean"])*0 + (sum_survLOAR$summary["beta[4]","mean"]-sum_survLOAR$summary["beta[4]","se_mean"])*mean(LOAR_data$origin_01)
                                               + (sum_survLOAR$summary["beta[5]","mean"]-sum_survLOAR$summary["beta[5]","se_mean"])*xdummy*0))
LOAR_ydummy_eminus_semax <- as.vector(invlogit((sum_survLOAR$summary["beta[1]","mean"]+sum_survLOAR$summary["beta[1]","se_mean"]) + (sum_survLOAR$summary["beta[2]","mean"]+sum_survLOAR$summary["beta[2]","se_mean"])*xdummy + (sum_survLOAR$summary["beta[3]","mean"]+sum_survLOAR$summary["beta[3]","se_mean"])*0 + (sum_survLOAR$summary["beta[4]","mean"]+sum_survLOAR$summary["beta[4]","se_mean"])*mean(LOAR_data$origin_01)
                                               + (sum_survLOAR$summary["beta[5]","mean"]+sum_survLOAR$summary["beta[5]","se_mean"])*xdummy*0))


LOAR_ydummy_eminus_y1 <- as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*0 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_survLOAR$'beta[5]')*xdummy*0
                                            + mean(post_survLOAR$'tau_year[1,1]')))
LOAR_ydummy_eminus_y2 <- as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*0 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_survLOAR$'beta[5]')*xdummy*0
                                            + mean(post_survLOAR$'tau_year[1,2]')))
LOAR_ydummy_eminus_y3 <- as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + + mean(post_survLOAR$'beta[3]')*0  + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_survLOAR$'beta[5]')*xdummy*0
                                            + mean(post_survLOAR$'tau_year[1,3]')))
LOAR_ydummy_eminus_y4 <- as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*0 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_survLOAR$'beta[5]')*xdummy*0
                                            + mean(post_survLOAR$'tau_year[1,4]'))) 
LOAR_ydummy_eminus_y5 <- as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*0 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_survLOAR$'beta[5]')*xdummy*0
                                            + mean(post_survLOAR$'tau_year[1,5]')))
LOAR_ydummy_eminus_y6 <- as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*0 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_survLOAR$'beta[5]')*xdummy*0
                                            + mean(post_survLOAR$'tau_year[1,6]')))
LOAR_ydummy_eminus_y7 <- as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*0 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_survLOAR$'beta[5]')*xdummy*0
                                            + mean(post_survLOAR$'tau_year[1,7]')))
LOAR_ydummy_eminus_y8 <- as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*0 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_survLOAR$'beta[5]')*xdummy*0
                                            + mean(post_survLOAR$'tau_year[1,8]')))
LOAR_ydummy_eminus_y9 <- as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*0 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_survLOAR$'beta[5]')*xdummy*0
                                            + mean(post_survLOAR$'tau_year[1,9]')))
LOAR_ydummy_eminus_y10 <- as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*0 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                             + mean(post_survLOAR$'beta[5]')*xdummy*0
                                             + mean(post_survLOAR$'tau_year[1,10]')))
LOAR_ydummy_eminus_y11 <- as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*0 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                             + mean(post_survLOAR$'beta[5]')*xdummy*0
                                             + mean(post_survLOAR$'tau_year[1,11]')))
LOARfitsminus0 <- as.data.frame(cbind(xdummy, LOAR_ydummy_eminus))
LOARfitsminus1 <- as.data.frame(cbind(xdummy,LOAR_ydummy_eminus_y1,LOAR_ydummy_eminus_y2, LOAR_ydummy_eminus_y3, LOAR_ydummy_eminus_y4, LOAR_ydummy_eminus_y5, LOAR_ydummy_eminus_y6, LOAR_ydummy_eminus_y7, LOAR_ydummy_eminus_y8, LOAR_ydummy_eminus_y9, LOAR_ydummy_eminus_y10, LOAR_ydummy_eminus_y11))
LOARfitsminus <- melt(LOARfitsminus1, id = "xdummy",
                      measure.vars = c("LOAR_ydummy_eminus_y1","LOAR_ydummy_eminus_y2", 
                                       "LOAR_ydummy_eminus_y3", "LOAR_ydummy_eminus_y4", 
                                       "LOAR_ydummy_eminus_y5", "LOAR_ydummy_eminus_y6", 
                                       "LOAR_ydummy_eminus_y7", "LOAR_ydummy_eminus_y8", 
                                       "LOAR_ydummy_eminus_y9", "LOAR_ydummy_eminus_y10",
                                       "LOAR_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "LOAR_ydummy_eminus_y1" = '1',"LOAR_ydummy_eminus_y2" = "2", 
                       "LOAR_ydummy_eminus_y3" = "3", "LOAR_ydummy_eminus_y4" = "4", 
                       "LOAR_ydummy_eminus_y5" = "5", "LOAR_ydummy_eminus_y6" = "6", 
                       "LOAR_ydummy_eminus_y7" = "7", "LOAR_ydummy_eminus_y8" = "8", 
                       "LOAR_ydummy_eminus_y9" = "9", "LOAR_ydummy_eminus_y10" = "10",
                       "LOAR_ydummy_eminus_y11" = "11") )



# E+
LOAR_ydummy_eplus <- as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*1 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                        + mean(post_survLOAR$'beta[5]')*xdummy*1))
LOAR_ydummy_eplus2.5 <- as.vector(invlogit(sum_survLOAR$summary["beta[1]","2.5%"] + (sum_survLOAR$summary["beta[2]","2.5%"])*xdummy + (sum_survLOAR$summary["beta[3]","2.5%"])*1 + (sum_survLOAR$summary["beta[4]","2.5%"])*mean(LOAR_data$origin_01)
                                            + (sum_survLOAR$summary["beta[5]","2.5%"])*xdummy*1))
LOAR_ydummy_eplus97.5 <- as.vector(invlogit(sum_survLOAR$summary["beta[1]","97.5%"] + (sum_survLOAR$summary["beta[2]","97.5%"])*xdummy + (sum_survLOAR$summary["beta[3]","97.5%"])*1 + (sum_survLOAR$summary["beta[4]","97.5%"])*mean(LOAR_data$origin_01)
                                             + (sum_survLOAR$summary["beta[5]","97.5%"])*xdummy*1))

LOAR_ydummy_eplus_semin <- as.vector(invlogit((sum_survLOAR$summary["beta[1]","mean"]-sum_survLOAR$summary["beta[1]","se_mean"]) + (sum_survLOAR$summary["beta[2]","mean"]-sum_survLOAR$summary["beta[2]","se_mean"])*xdummy + (sum_survLOAR$summary["beta[3]","mean"]-sum_survLOAR$summary["beta[3]","se_mean"])*1 + (sum_survLOAR$summary["beta[4]","mean"]-sum_survLOAR$summary["beta[4]","se_mean"])*mean(LOAR_data$origin_01)
                                               + (sum_survLOAR$summary["beta[5]","mean"]-sum_survLOAR$summary["beta[5]","se_mean"])*xdummy*1))
LOAR_ydummy_eplus_semax <- as.vector(invlogit((sum_survLOAR$summary["beta[1]","mean"]+sum_survLOAR$summary["beta[1]","se_mean"]) + (sum_survLOAR$summary["beta[2]","mean"]+sum_survLOAR$summary["beta[2]","se_mean"])*xdummy + (sum_survLOAR$summary["beta[3]","mean"]+sum_survLOAR$summary["beta[3]","se_mean"])*1 + (sum_survLOAR$summary["beta[4]","mean"]+sum_survLOAR$summary["beta[4]","se_mean"])*mean(LOAR_data$origin_01)
                                               + (sum_survLOAR$summary["beta[5]","mean"]+sum_survLOAR$summary["beta[5]","se_mean"])*xdummy*1))


LOAR_ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*1 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_survLOAR$'beta[5]')*xdummy*1
                                            + mean(post_survLOAR$'tau_year[2,1]')))
LOAR_ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*1 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_survLOAR$'beta[5]')*xdummy*1
                                            + mean(post_survLOAR$'tau_year[2,2')))
LOAR_ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*1 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_survLOAR$'beta[5]')*xdummy*1
                                            + mean(post_survLOAR$'tau_year[2,3]')))
LOAR_ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*1 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_survLOAR$'beta[5]')*xdummy*1
                                            + mean(post_survLOAR$'tau_year[2,4]')))
LOAR_ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*1 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_survLOAR$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_survLOAR$'tau_year[2,5]')))
LOAR_ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*1 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_survLOAR$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_survLOAR$'tau_year[2,6]')))
LOAR_ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*1 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_survLOAR$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_survLOAR$'tau_year[2,7]')))
LOAR_ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*1 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_survLOAR$'beta[5]')*xdummy*1
                                            + mean(post_survLOAR$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
LOAR_ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*1 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_survLOAR$'beta[5]')*xdummy*1
                                            + mean(post_survLOAR$'tau_year[2,9]')))
LOAR_ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*1 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                             + mean(post_survLOAR$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_survLOAR$'tau_year[2,10]')))
LOAR_ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_survLOAR$'beta[1]') + mean(post_survLOAR$'beta[2]')*xdummy + mean(post_survLOAR$'beta[3]')*1 + mean(post_survLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                             + mean(post_survLOAR$'beta[5]')*xdummy*1
                                             + mean(post_survLOAR$'tau_year[2,11]')))
LOAR_fitsplus0 <- as.data.frame(cbind(xdummy, LOAR_ydummy_eplus))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
LOARfitsplus1 <- as.data.frame(cbind(xdummy,LOAR_ydummy_eplus_y1,LOAR_ydummy_eplus_y2, LOAR_ydummy_eplus_y3, LOAR_ydummy_eplus_y4, LOAR_ydummy_eplus_y5, LOAR_ydummy_eplus_y6, LOAR_ydummy_eplus_y7, LOAR_ydummy_eplus_y8, LOAR_ydummy_eplus_y9, LOAR_ydummy_eplus_y10, LOAR_ydummy_eplus_y11))
LOARfitsplus <- melt(LOARfitsplus1, id = "xdummy",
                     measure.vars = c("LOAR_ydummy_eplus_y1","LOAR_ydummy_eplus_y2", 
                                      "LOAR_ydummy_eplus_y3", "LOAR_ydummy_eplus_y4", 
                                      "LOAR_ydummy_eplus_y5", "LOAR_ydummy_eplus_y6", 
                                      "LOAR_ydummy_eplus_y7", "LOAR_ydummy_eplus_y8", 
                                      "LOAR_ydummy_eplus_y9", "LOAR_ydummy_eplus_y10",
                                      "LOAR_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "LOAR_ydummy_eplus_y1" = '1',"LOAR_ydummy_eplus_y2" = "2", 
                       "LOAR_ydummy_eplus_y3" = "3", "LOAR_ydummy_eplus_y4" = "4", 
                       "LOAR_ydummy_eplus_y5" = "5", "LOAR_ydummy_eplus_y6" = "6", 
                       "LOAR_ydummy_eplus_y7" = "7", "LOAR_ydummy_eplus_y8" = "8", 
                       "LOAR_ydummy_eplus_y9" = "9", "LOAR_ydummy_eplus_y10" = "10",
                       "LOAR_ydummy_eplus_y11" = "11") )


# E+ by year
LOAREplusbyyear <- ggplot(data = LOARfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = LOARsurv_bin1, aes(x = mean_size, y = mean_surv, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("Prob. of Survival") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = LOARfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("Prob. of Survival")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
LOAREminusbyyear <- ggplot(data = LOARfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = LOARsurv_bin0, aes(x = mean_size, y = mean_surv, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("Prob. of Survival")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = LOARfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("LOAR E- Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)


# effect on mean


LOARfits <- as.data.frame(cbind(xdummy, LOAR_ydummy_eplus, LOAR_ydummy_eplus_semin, LOAR_ydummy_eplus_semax, LOAR_ydummy_eplus2.5, LOAR_ydummy_eplus97.5, LOAR_ydummy_eminus, LOAR_ydummy_eminus_semin, LOAR_ydummy_eminus_semax, LOAR_ydummy_eminus2.5, LOAR_ydummy_eminus97.5))

# without point
ggplot(data = LOARfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("LOAR Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = LOARfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = LOARsurv_binmean, aes(x = mean_size, y = mean_surv, color = endo), lwd = 3) +
  ggtitle("LOAR Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

LOARmean <- ggplot(data = LOARfits) +
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = LOARsurv_binmean, aes(x = mean_size, y = mean_surv, color = endo, lwd = samplesize)) +
  ggtitle("Mean Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")

LOARsurv4panel <- grid.arrange(LOARmean, LOAR_s,LOAREplusbyyear,LOAREminusbyyear, ncol= 4)
titleLOARsurv4panel <- annotate_figure(LOARsurv4panel, top = "Lolium arundinacea")

titleLOARsurv4panel

# split up plus and minus with datapoints and add mean line
LOAR_minus <- ggplot(data = LOARfitsminus1)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = LOARfitsminus0, aes(x = xdummy, y = LOAR_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = LOARsurv_bin0, aes(x = mean_size, y = mean_surv), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


LOAR_plus <- ggplot(data = LOARfitsplus1)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = LOAR_fitsplus0, aes(x = xdummy, y = LOAR_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = LOARsurv_bin1, aes(x = mean_size, y = mean_surv), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

LOARsurv <- grid.arrange(LOAR_plus, LOAR_minus, ncol= 2)
titleLOARsurv <- annotate_figure(LOARsurv, top = "LOAR", left = "Probability of Survival")

titleLOARsurv







######## FESU model fits
# actual data points for probability of survival
FESUsurv_bin0 <- FESU_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  filter(endo == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())

FESUsurv_bin1 <- FESU_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  filter(endo == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())
FESUsurv_bin <- FESU_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())


FESUsurv_binmean <- FESU_data %>% 
  select(surv_t1, logsize_t, endo) %>%
  mutate(endo = recode(endo, "0" = "E-", "1" = "E+", minus = "E-", plus = "E+")) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(FESU_data$logsize_t)
xdummy <- seq(min(FESU_data$logsize_t), max(FESU_data$logsize_t), length.out = nvalues)

# E- 
FESU_ydummy_eminus <- as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*0 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01)
                                         + mean(post_survFESU$'beta[5]')*xdummy*0))
FESU_ydummy_eminus2.5 <- as.vector(invlogit(sum_survFESU$summary["beta[1]","2.5%"] + (sum_survFESU$summary["beta[2]","2.5%"])*xdummy + (sum_survFESU$summary["beta[3]","2.5%"])*0 + (sum_survFESU$summary["beta[4]","2.5%"])*mean(FESU_data$origin_01)
                                            + (sum_survFESU$summary["beta[5]","2.5%"])*xdummy*0))
FESU_ydummy_eminus97.5 <- as.vector(invlogit(sum_survFESU$summary["beta[1]","97.5%"] + (sum_survFESU$summary["beta[2]","97.5%"])*xdummy + (sum_survFESU$summary["beta[3]","97.5%"])*0 + (sum_survFESU$summary["beta[4]","97.5%"])*mean(FESU_data$origin_01)
                                             + (sum_survFESU$summary["beta[5]","97.5%"])*xdummy*0))

FESU_ydummy_eminus_semin <- as.vector(invlogit((sum_survFESU$summary["beta[1]","mean"]-sum_survFESU$summary["beta[1]","se_mean"]) + (sum_survFESU$summary["beta[2]","mean"]-sum_survFESU$summary["beta[2]","se_mean"])*xdummy + (sum_survFESU$summary["beta[3]","mean"]-sum_survFESU$summary["beta[3]","se_mean"])*0 + (sum_survFESU$summary["beta[4]","mean"]-sum_survFESU$summary["beta[4]","se_mean"])*mean(FESU_data$origin_01)
                                               + (sum_survFESU$summary["beta[5]","mean"]-sum_survFESU$summary["beta[5]","se_mean"])*xdummy*0))
FESU_ydummy_eminus_semax <- as.vector(invlogit((sum_survFESU$summary["beta[1]","mean"]+sum_survFESU$summary["beta[1]","se_mean"]) + (sum_survFESU$summary["beta[2]","mean"]+sum_survFESU$summary["beta[2]","se_mean"])*xdummy + (sum_survFESU$summary["beta[3]","mean"]+sum_survFESU$summary["beta[3]","se_mean"])*0 + (sum_survFESU$summary["beta[4]","mean"]+sum_survFESU$summary["beta[4]","se_mean"])*mean(FESU_data$origin_01)
                                               + (sum_survFESU$summary["beta[5]","mean"]+sum_survFESU$summary["beta[5]","se_mean"])*xdummy*0))


FESU_ydummy_eminus_y1 <- as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*0 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_survFESU$'beta[5]')*xdummy*0
                                            + mean(post_survFESU$'tau_year[1,1]')))
FESU_ydummy_eminus_y2 <- as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*0 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_survFESU$'beta[5]')*xdummy*0
                                            + mean(post_survFESU$'tau_year[1,2]')))
FESU_ydummy_eminus_y3 <- as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + + mean(post_survFESU$'beta[3]')*0  + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_survFESU$'beta[5]')*xdummy*0
                                            + mean(post_survFESU$'tau_year[1,3]')))
FESU_ydummy_eminus_y4 <- as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*0 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_survFESU$'beta[5]')*xdummy*0
                                            + mean(post_survFESU$'tau_year[1,4]'))) 
FESU_ydummy_eminus_y5 <- as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*0 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_survFESU$'beta[5]')*xdummy*0
                                            + mean(post_survFESU$'tau_year[1,5]')))
FESU_ydummy_eminus_y6 <- as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*0 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_survFESU$'beta[5]')*xdummy*0
                                            + mean(post_survFESU$'tau_year[1,6]')))
FESU_ydummy_eminus_y7 <- as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*0 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_survFESU$'beta[5]')*xdummy*0
                                            + mean(post_survFESU$'tau_year[1,7]')))
FESU_ydummy_eminus_y8 <- as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*0 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_survFESU$'beta[5]')*xdummy*0
                                            + mean(post_survFESU$'tau_year[1,8]')))
FESU_ydummy_eminus_y9 <- as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*0 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_survFESU$'beta[5]')*xdummy*0
                                            + mean(post_survFESU$'tau_year[1,9]')))
FESU_ydummy_eminus_y10 <- as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*0 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01)
                                             + mean(post_survFESU$'beta[5]')*xdummy*0
                                             + mean(post_survFESU$'tau_year[1,10]')))
FESU_ydummy_eminus_y11 <- as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*0 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01)
                                             + mean(post_survFESU$'beta[5]')*xdummy*0
                                             + mean(post_survFESU$'tau_year[1,11]')))
FESUfitsminus0 <- as.data.frame(cbind(xdummy, FESU_ydummy_eminus))
FESUfitsminus1 <- as.data.frame(cbind(xdummy,FESU_ydummy_eminus_y1,FESU_ydummy_eminus_y2, FESU_ydummy_eminus_y3, FESU_ydummy_eminus_y4, FESU_ydummy_eminus_y5, FESU_ydummy_eminus_y6, FESU_ydummy_eminus_y7, FESU_ydummy_eminus_y8, FESU_ydummy_eminus_y9, FESU_ydummy_eminus_y10, FESU_ydummy_eminus_y11))
FESUfitsminus <- melt(FESUfitsminus1, id = "xdummy",
                      measure.vars = c("FESU_ydummy_eminus_y1","FESU_ydummy_eminus_y2", 
                                       "FESU_ydummy_eminus_y3", "FESU_ydummy_eminus_y4", 
                                       "FESU_ydummy_eminus_y5", "FESU_ydummy_eminus_y6", 
                                       "FESU_ydummy_eminus_y7", "FESU_ydummy_eminus_y8", 
                                       "FESU_ydummy_eminus_y9", "FESU_ydummy_eminus_y10",
                                       "FESU_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "FESU_ydummy_eminus_y1" = '1',"FESU_ydummy_eminus_y2" = "2", 
                       "FESU_ydummy_eminus_y3" = "3", "FESU_ydummy_eminus_y4" = "4", 
                       "FESU_ydummy_eminus_y5" = "5", "FESU_ydummy_eminus_y6" = "6", 
                       "FESU_ydummy_eminus_y7" = "7", "FESU_ydummy_eminus_y8" = "8", 
                       "FESU_ydummy_eminus_y9" = "9", "FESU_ydummy_eminus_y10" = "10",
                       "FESU_ydummy_eminus_y11" = "11") )



# E+
FESU_ydummy_eplus <- as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*1 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01)
                                        + mean(post_survFESU$'beta[5]')*xdummy*1))
FESU_ydummy_eplus2.5 <- as.vector(invlogit(sum_survFESU$summary["beta[1]","2.5%"] + (sum_survFESU$summary["beta[2]","2.5%"])*xdummy + (sum_survFESU$summary["beta[3]","2.5%"])*1 + (sum_survFESU$summary["beta[4]","2.5%"])*mean(FESU_data$origin_01)
                                            + (sum_survFESU$summary["beta[5]","2.5%"])*xdummy*1))
FESU_ydummy_eplus97.5 <- as.vector(invlogit(sum_survFESU$summary["beta[1]","97.5%"] + (sum_survFESU$summary["beta[2]","97.5%"])*xdummy + (sum_survFESU$summary["beta[3]","97.5%"])*1 + (sum_survFESU$summary["beta[4]","97.5%"])*mean(FESU_data$origin_01)
                                             + (sum_survFESU$summary["beta[5]","97.5%"])*xdummy*1))

FESU_ydummy_eplus_semin <- as.vector(invlogit((sum_survFESU$summary["beta[1]","mean"]-sum_survFESU$summary["beta[1]","se_mean"]) + (sum_survFESU$summary["beta[2]","mean"]-sum_survFESU$summary["beta[2]","se_mean"])*xdummy + (sum_survFESU$summary["beta[3]","mean"]-sum_survFESU$summary["beta[3]","se_mean"])*1 + (sum_survFESU$summary["beta[4]","mean"]-sum_survFESU$summary["beta[4]","se_mean"])*mean(FESU_data$origin_01)
                                               + (sum_survFESU$summary["beta[5]","mean"]-sum_survFESU$summary["beta[5]","se_mean"])*xdummy*1))
FESU_ydummy_eplus_semax <- as.vector(invlogit((sum_survFESU$summary["beta[1]","mean"]+sum_survFESU$summary["beta[1]","se_mean"]) + (sum_survFESU$summary["beta[2]","mean"]+sum_survFESU$summary["beta[2]","se_mean"])*xdummy + (sum_survFESU$summary["beta[3]","mean"]+sum_survFESU$summary["beta[3]","se_mean"])*1 + (sum_survFESU$summary["beta[4]","mean"]+sum_survFESU$summary["beta[4]","se_mean"])*mean(FESU_data$origin_01)
                                               + (sum_survFESU$summary["beta[5]","mean"]+sum_survFESU$summary["beta[5]","se_mean"])*xdummy*1))


FESU_ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*1 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_survFESU$'beta[5]')*xdummy*1
                                            + mean(post_survFESU$'tau_year[2,1]')))
FESU_ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*1 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_survFESU$'beta[5]')*xdummy*1
                                            + mean(post_survFESU$'tau_year[2,2')))
FESU_ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*1 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_survFESU$'beta[5]')*xdummy*1
                                            + mean(post_survFESU$'tau_year[2,3]')))
FESU_ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*1 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_survFESU$'beta[5]')*xdummy*1
                                            + mean(post_survFESU$'tau_year[2,4]')))
FESU_ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*1 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_survFESU$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_survFESU$'tau_year[2,5]')))
FESU_ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*1 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_survFESU$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_survFESU$'tau_year[2,6]')))
FESU_ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*1 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_survFESU$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_survFESU$'tau_year[2,7]')))
FESU_ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*1 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_survFESU$'beta[5]')*xdummy*1
                                            + mean(post_survFESU$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
FESU_ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*1 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_survFESU$'beta[5]')*xdummy*1
                                            + mean(post_survFESU$'tau_year[2,9]')))
FESU_ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*1 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                             + mean(post_survFESU$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_survFESU$'tau_year[2,10]')))
FESU_ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_survFESU$'beta[1]') + mean(post_survFESU$'beta[2]')*xdummy + mean(post_survFESU$'beta[3]')*1 + mean(post_survFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                             + mean(post_survFESU$'beta[5]')*xdummy*1
                                             + mean(post_survFESU$'tau_year[2,11]')))
FESU_fitsplus0 <- as.data.frame(cbind(xdummy, FESU_ydummy_eplus))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
FESUfitsplus1 <- as.data.frame(cbind(xdummy,FESU_ydummy_eplus_y1,FESU_ydummy_eplus_y2, FESU_ydummy_eplus_y3, FESU_ydummy_eplus_y4, FESU_ydummy_eplus_y5, FESU_ydummy_eplus_y6, FESU_ydummy_eplus_y7, FESU_ydummy_eplus_y8, FESU_ydummy_eplus_y9, FESU_ydummy_eplus_y10, FESU_ydummy_eplus_y11))
FESUfitsplus <- melt(FESUfitsplus1, id = "xdummy",
                     measure.vars = c("FESU_ydummy_eplus_y1","FESU_ydummy_eplus_y2", 
                                      "FESU_ydummy_eplus_y3", "FESU_ydummy_eplus_y4", 
                                      "FESU_ydummy_eplus_y5", "FESU_ydummy_eplus_y6", 
                                      "FESU_ydummy_eplus_y7", "FESU_ydummy_eplus_y8", 
                                      "FESU_ydummy_eplus_y9", "FESU_ydummy_eplus_y10",
                                      "FESU_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "FESU_ydummy_eplus_y1" = '1',"FESU_ydummy_eplus_y2" = "2", 
                       "FESU_ydummy_eplus_y3" = "3", "FESU_ydummy_eplus_y4" = "4", 
                       "FESU_ydummy_eplus_y5" = "5", "FESU_ydummy_eplus_y6" = "6", 
                       "FESU_ydummy_eplus_y7" = "7", "FESU_ydummy_eplus_y8" = "8", 
                       "FESU_ydummy_eplus_y9" = "9", "FESU_ydummy_eplus_y10" = "10",
                       "FESU_ydummy_eplus_y11" = "11") )




# E+ by year
FESUEplusbyyear <- ggplot(data = FESUfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = FESUsurv_bin1, aes(x = mean_size, y = mean_surv, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("Prob. of Survival") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = FESUfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("Prob. of Survival")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
FESUEminusbyyear <- ggplot(data = FESUfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = FESUsurv_bin0, aes(x = mean_size, y = mean_surv, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("Prob. of Survival")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = FESUfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("FESU E- Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)



# effect on mean


FESUfits <- as.data.frame(cbind(xdummy, FESU_ydummy_eplus, FESU_ydummy_eplus_semin, FESU_ydummy_eplus_semax, FESU_ydummy_eplus2.5, FESU_ydummy_eplus97.5, FESU_ydummy_eminus, FESU_ydummy_eminus_semin, FESU_ydummy_eminus_semax, FESU_ydummy_eminus2.5, FESU_ydummy_eminus97.5))

# without point
ggplot(data = FESUfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("FESU Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = FESUfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = FESUsurv_binmean, aes(x = mean_size, y = mean_surv, color = endo), lwd = 3) +
  ggtitle("FESU Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

FESUmean <- ggplot(data = FESUfits) +
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = FESUsurv_binmean, aes(x = mean_size, y = mean_surv, color = endo, lwd = samplesize)) +
  ggtitle("Mean Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")


FESUsurv4panel <- grid.arrange(FESUmean, FESU_s,FESUEplusbyyear,FESUEminusbyyear, ncol= 4)
titleFESUsurv4panel <- annotate_figure(FESUsurv4panel, top = "Festuca subverticillata")

titleFESUsurv4panel

# split up plus and minus with datapoints and add mean line
FESU_minus <- ggplot(data = FESUfitsminus1)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = FESUfitsminus0, aes(x = xdummy, y = FESU_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = FESUsurv_bin0, aes(x = mean_size, y = mean_surv), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


FESU_plus <- ggplot(data = FESUfitsplus1)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = FESU_fitsplus0, aes(x = xdummy, y = FESU_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = FESUsurv_bin1, aes(x = mean_size, y = mean_surv), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

FESUsurv <- grid.arrange(FESU_plus, FESU_minus, ncol= 2)
titleFESUsurv <- annotate_figure(FESUsurv, top = "FESU", left = "Probability of Survival")

titleFESUsurv






######## ELVI model fits
# actual data points for probability of survival
ELVIsurv_bin0 <- ELVI_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  filter(endo == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())

ELVIsurv_bin1 <- ELVI_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  filter(endo == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())
ELVIsurv_bin <- ELVI_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())


ELVIsurv_binmean <- ELVI_data %>% 
  select(surv_t1, logsize_t, endo) %>%
  mutate(endo = recode(endo, "0" = "E-", "1" = "E+", minus = "E-", plus = "E+")) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(ELVI_data$logsize_t)
xdummy <- seq(min(ELVI_data$logsize_t), max(ELVI_data$logsize_t), length.out = nvalues)

# E- 
ELVI_ydummy_eminus <- as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*0 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                         + mean(post_survELVI$'beta[5]')*xdummy*0))
ELVI_ydummy_eminus2.5 <- as.vector(invlogit(sum_survELVI$summary["beta[1]","2.5%"] + (sum_survELVI$summary["beta[2]","2.5%"])*xdummy + (sum_survELVI$summary["beta[3]","2.5%"])*0 + (sum_survELVI$summary["beta[4]","2.5%"])*mean(ELVI_data$origin_01)
                                            + (sum_survELVI$summary["beta[5]","2.5%"])*xdummy*0))
ELVI_ydummy_eminus97.5 <- as.vector(invlogit(sum_survELVI$summary["beta[1]","97.5%"] + (sum_survELVI$summary["beta[2]","97.5%"])*xdummy + (sum_survELVI$summary["beta[3]","97.5%"])*0 + (sum_survELVI$summary["beta[4]","97.5%"])*mean(ELVI_data$origin_01)
                                             + (sum_survELVI$summary["beta[5]","97.5%"])*xdummy*0))

ELVI_ydummy_eminus_semin <- as.vector(invlogit((sum_survELVI$summary["beta[1]","mean"]-sum_survELVI$summary["beta[1]","se_mean"]) + (sum_survELVI$summary["beta[2]","mean"]-sum_survELVI$summary["beta[2]","se_mean"])*xdummy + (sum_survELVI$summary["beta[3]","mean"]-sum_survELVI$summary["beta[3]","se_mean"])*0 + (sum_survELVI$summary["beta[4]","mean"]-sum_survELVI$summary["beta[4]","se_mean"])*mean(ELVI_data$origin_01)
                                               + (sum_survELVI$summary["beta[5]","mean"]-sum_survELVI$summary["beta[5]","se_mean"])*xdummy*0))
ELVI_ydummy_eminus_semax <- as.vector(invlogit((sum_survELVI$summary["beta[1]","mean"]+sum_survELVI$summary["beta[1]","se_mean"]) + (sum_survELVI$summary["beta[2]","mean"]+sum_survELVI$summary["beta[2]","se_mean"])*xdummy + (sum_survELVI$summary["beta[3]","mean"]+sum_survELVI$summary["beta[3]","se_mean"])*0 + (sum_survELVI$summary["beta[4]","mean"]+sum_survELVI$summary["beta[4]","se_mean"])*mean(ELVI_data$origin_01)
                                               + (sum_survELVI$summary["beta[5]","mean"]+sum_survELVI$summary["beta[5]","se_mean"])*xdummy*0))


ELVI_ydummy_eminus_y1 <- as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*0 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_survELVI$'beta[5]')*xdummy*0
                                            + mean(post_survELVI$'tau_year[1,1]')))
ELVI_ydummy_eminus_y2 <- as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*0 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_survELVI$'beta[5]')*xdummy*0
                                            + mean(post_survELVI$'tau_year[1,2]')))
ELVI_ydummy_eminus_y3 <- as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + + mean(post_survELVI$'beta[3]')*0  + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_survELVI$'beta[5]')*xdummy*0
                                            + mean(post_survELVI$'tau_year[1,3]')))
ELVI_ydummy_eminus_y4 <- as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*0 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_survELVI$'beta[5]')*xdummy*0
                                            + mean(post_survELVI$'tau_year[1,4]'))) 
ELVI_ydummy_eminus_y5 <- as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*0 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_survELVI$'beta[5]')*xdummy*0
                                            + mean(post_survELVI$'tau_year[1,5]')))
ELVI_ydummy_eminus_y6 <- as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*0 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_survELVI$'beta[5]')*xdummy*0
                                            + mean(post_survELVI$'tau_year[1,6]')))
ELVI_ydummy_eminus_y7 <- as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*0 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_survELVI$'beta[5]')*xdummy*0
                                            + mean(post_survELVI$'tau_year[1,7]')))
ELVI_ydummy_eminus_y8 <- as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*0 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_survELVI$'beta[5]')*xdummy*0
                                            + mean(post_survELVI$'tau_year[1,8]')))
ELVI_ydummy_eminus_y9 <- as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*0 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_survELVI$'beta[5]')*xdummy*0
                                            + mean(post_survELVI$'tau_year[1,9]')))
ELVI_ydummy_eminus_y10 <- as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*0 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                             + mean(post_survELVI$'beta[5]')*xdummy*0
                                             + mean(post_survELVI$'tau_year[1,10]')))
ELVI_ydummy_eminus_y11 <- as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*0 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                             + mean(post_survELVI$'beta[5]')*xdummy*0
                                             + mean(post_survELVI$'tau_year[1,11]')))
ELVIfitsminus0 <- as.data.frame(cbind(xdummy, ELVI_ydummy_eminus))
ELVIfitsminus1 <- as.data.frame(cbind(xdummy,ELVI_ydummy_eminus_y1,ELVI_ydummy_eminus_y2, ELVI_ydummy_eminus_y3, ELVI_ydummy_eminus_y4, ELVI_ydummy_eminus_y5, ELVI_ydummy_eminus_y6, ELVI_ydummy_eminus_y7, ELVI_ydummy_eminus_y8, ELVI_ydummy_eminus_y9, ELVI_ydummy_eminus_y10, ELVI_ydummy_eminus_y11))
ELVIfitsminus <- melt(ELVIfitsminus1, id = "xdummy",
                      measure.vars = c("ELVI_ydummy_eminus_y1","ELVI_ydummy_eminus_y2", 
                                       "ELVI_ydummy_eminus_y3", "ELVI_ydummy_eminus_y4", 
                                       "ELVI_ydummy_eminus_y5", "ELVI_ydummy_eminus_y6", 
                                       "ELVI_ydummy_eminus_y7", "ELVI_ydummy_eminus_y8", 
                                       "ELVI_ydummy_eminus_y9", "ELVI_ydummy_eminus_y10",
                                       "ELVI_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "ELVI_ydummy_eminus_y1" = '1',"ELVI_ydummy_eminus_y2" = "2", 
                       "ELVI_ydummy_eminus_y3" = "3", "ELVI_ydummy_eminus_y4" = "4", 
                       "ELVI_ydummy_eminus_y5" = "5", "ELVI_ydummy_eminus_y6" = "6", 
                       "ELVI_ydummy_eminus_y7" = "7", "ELVI_ydummy_eminus_y8" = "8", 
                       "ELVI_ydummy_eminus_y9" = "9", "ELVI_ydummy_eminus_y10" = "10",
                       "ELVI_ydummy_eminus_y11" = "11") )



# E+
ELVI_ydummy_eplus <- as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*1 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                        + mean(post_survELVI$'beta[5]')*xdummy*1))
ELVI_ydummy_eplus2.5 <- as.vector(invlogit(sum_survELVI$summary["beta[1]","2.5%"] + (sum_survELVI$summary["beta[2]","2.5%"])*xdummy + (sum_survELVI$summary["beta[3]","2.5%"])*1 + (sum_survELVI$summary["beta[4]","2.5%"])*mean(ELVI_data$origin_01)
                                            + (sum_survELVI$summary["beta[5]","2.5%"])*xdummy*1))
ELVI_ydummy_eplus97.5 <- as.vector(invlogit(sum_survELVI$summary["beta[1]","97.5%"] + (sum_survELVI$summary["beta[2]","97.5%"])*xdummy + (sum_survELVI$summary["beta[3]","97.5%"])*1 + (sum_survELVI$summary["beta[4]","97.5%"])*mean(ELVI_data$origin_01)
                                             + (sum_survELVI$summary["beta[5]","97.5%"])*xdummy*1))

ELVI_ydummy_eplus_semin <- as.vector(invlogit((sum_survELVI$summary["beta[1]","mean"]-sum_survELVI$summary["beta[1]","se_mean"]) + (sum_survELVI$summary["beta[2]","mean"]-sum_survELVI$summary["beta[2]","se_mean"])*xdummy + (sum_survELVI$summary["beta[3]","mean"]-sum_survELVI$summary["beta[3]","se_mean"])*1 + (sum_survELVI$summary["beta[4]","mean"]-sum_survELVI$summary["beta[4]","se_mean"])*mean(ELVI_data$origin_01)
                                               + (sum_survELVI$summary["beta[5]","mean"]-sum_survELVI$summary["beta[5]","se_mean"])*xdummy*1))
ELVI_ydummy_eplus_semax <- as.vector(invlogit((sum_survELVI$summary["beta[1]","mean"]+sum_survELVI$summary["beta[1]","se_mean"]) + (sum_survELVI$summary["beta[2]","mean"]+sum_survELVI$summary["beta[2]","se_mean"])*xdummy + (sum_survELVI$summary["beta[3]","mean"]+sum_survELVI$summary["beta[3]","se_mean"])*1 + (sum_survELVI$summary["beta[4]","mean"]+sum_survELVI$summary["beta[4]","se_mean"])*mean(ELVI_data$origin_01)
                                               + (sum_survELVI$summary["beta[5]","mean"]+sum_survELVI$summary["beta[5]","se_mean"])*xdummy*1))


ELVI_ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*1 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_survELVI$'beta[5]')*xdummy*1
                                            + mean(post_survELVI$'tau_year[2,1]')))
ELVI_ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*1 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_survELVI$'beta[5]')*xdummy*1
                                            + mean(post_survELVI$'tau_year[2,2')))
ELVI_ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*1 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_survELVI$'beta[5]')*xdummy*1
                                            + mean(post_survELVI$'tau_year[2,3]')))
ELVI_ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*1 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_survELVI$'beta[5]')*xdummy*1
                                            + mean(post_survELVI$'tau_year[2,4]')))
ELVI_ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*1 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_survELVI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_survELVI$'tau_year[2,5]')))
ELVI_ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*1 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_survELVI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_survELVI$'tau_year[2,6]')))
ELVI_ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*1 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_survELVI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_survELVI$'tau_year[2,7]')))
ELVI_ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*1 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_survELVI$'beta[5]')*xdummy*1
                                            + mean(post_survELVI$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
ELVI_ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*1 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_survELVI$'beta[5]')*xdummy*1
                                            + mean(post_survELVI$'tau_year[2,9]')))
ELVI_ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*1 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                             + mean(post_survELVI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_survELVI$'tau_year[2,10]')))
ELVI_ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_survELVI$'beta[1]') + mean(post_survELVI$'beta[2]')*xdummy + mean(post_survELVI$'beta[3]')*1 + mean(post_survELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                             + mean(post_survELVI$'beta[5]')*xdummy*1
                                             + mean(post_survELVI$'tau_year[2,11]')))
ELVI_fitsplus0 <- as.data.frame(cbind(xdummy, ELVI_ydummy_eplus))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
ELVIfitsplus1 <- as.data.frame(cbind(xdummy,ELVI_ydummy_eplus_y1,ELVI_ydummy_eplus_y2, ELVI_ydummy_eplus_y3, ELVI_ydummy_eplus_y4, ELVI_ydummy_eplus_y5, ELVI_ydummy_eplus_y6, ELVI_ydummy_eplus_y7, ELVI_ydummy_eplus_y8, ELVI_ydummy_eplus_y9, ELVI_ydummy_eplus_y10, ELVI_ydummy_eplus_y11))
ELVIfitsplus <- melt(ELVIfitsplus1, id = "xdummy",
                     measure.vars = c("ELVI_ydummy_eplus_y1","ELVI_ydummy_eplus_y2", 
                                      "ELVI_ydummy_eplus_y3", "ELVI_ydummy_eplus_y4", 
                                      "ELVI_ydummy_eplus_y5", "ELVI_ydummy_eplus_y6", 
                                      "ELVI_ydummy_eplus_y7", "ELVI_ydummy_eplus_y8", 
                                      "ELVI_ydummy_eplus_y9", "ELVI_ydummy_eplus_y10",
                                      "ELVI_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "ELVI_ydummy_eplus_y1" = '1',"ELVI_ydummy_eplus_y2" = "2", 
                       "ELVI_ydummy_eplus_y3" = "3", "ELVI_ydummy_eplus_y4" = "4", 
                       "ELVI_ydummy_eplus_y5" = "5", "ELVI_ydummy_eplus_y6" = "6", 
                       "ELVI_ydummy_eplus_y7" = "7", "ELVI_ydummy_eplus_y8" = "8", 
                       "ELVI_ydummy_eplus_y9" = "9", "ELVI_ydummy_eplus_y10" = "10",
                       "ELVI_ydummy_eplus_y11" = "11") )


# E+ by year
ELVIEplusbyyear <- ggplot(data = ELVIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = ELVIsurv_bin1, aes(x = mean_size, y = mean_surv, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("Prob. of Survival") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = ELVIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("Prob. of Survival")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
ELVIEminusbyyear <- ggplot(data = ELVIfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = ELVIsurv_bin0, aes(x = mean_size, y = mean_surv, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("Prob. of Survival")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = ELVIfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("ELVI E- Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)




# effect on mean


ELVIfits <- as.data.frame(cbind(xdummy, ELVI_ydummy_eplus, ELVI_ydummy_eplus_semin, ELVI_ydummy_eplus_semax, ELVI_ydummy_eplus2.5, ELVI_ydummy_eplus97.5, ELVI_ydummy_eminus, ELVI_ydummy_eminus_semin, ELVI_ydummy_eminus_semax, ELVI_ydummy_eminus2.5, ELVI_ydummy_eminus97.5))

# without point
ggplot(data = ELVIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("ELVI Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = ELVIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = ELVIsurv_binmean, aes(x = mean_size, y = mean_surv, color = endo), lwd = 3) +
  ggtitle("ELVI Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

ELVImean <- ggplot(data = ELVIfits) +
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = ELVIsurv_binmean, aes(x = mean_size, y = mean_surv, color = endo, lwd = samplesize)) +
  ggtitle("Mean Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")


ELVIsurv4panel <- grid.arrange(ELVImean, ELVI_s,ELVIEplusbyyear,ELVIEminusbyyear, ncol= 4)
titleELVIsurv4panel <- annotate_figure(ELVIsurv4panel, top = "Elymus virginicus")

titleELVIsurv4panel

# split up plus and minus with datapoints and add mean line
ELVI_minus <- ggplot(data = ELVIfitsminus1)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = ELVIfitsminus0, aes(x = xdummy, y = ELVI_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = ELVIsurv_bin0, aes(x = mean_size, y = mean_surv), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


ELVI_plus <- ggplot(data = ELVIfitsplus1)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = ELVI_fitsplus0, aes(x = xdummy, y = ELVI_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = ELVIsurv_bin1, aes(x = mean_size, y = mean_surv), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

ELVIsurv <- grid.arrange(ELVI_plus, ELVI_minus, ncol= 2)
titleELVIsurv <- annotate_figure(ELVIsurv, top = "ELVI", left = "Probability of Survival")

titleELVIsurv






######## ELRI model fits
# actual data points for probability of survival
ELRIsurv_bin0 <- ELRI_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  filter(endo == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())

ELRIsurv_bin1 <- ELRI_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  filter(endo == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())
ELRIsurv_bin <- ELRI_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())


ELRIsurv_binmean <- ELRI_data %>% 
  select(surv_t1, logsize_t, endo) %>%
  mutate(endo = recode(endo, "0" = "E-", "1" = "E+", minus = "E-", plus = "E+")) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(ELRI_data$logsize_t)
xdummy <- seq(min(ELRI_data$logsize_t), max(ELRI_data$logsize_t), length.out = nvalues)

# E- 
ELRI_ydummy_eminus <- as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*0 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                         + mean(post_survELRI$'beta[5]')*xdummy*0))
ELRI_ydummy_eminus2.5 <- as.vector(invlogit(sum_survELRI$summary["beta[1]","2.5%"] + (sum_survELRI$summary["beta[2]","2.5%"])*xdummy + (sum_survELRI$summary["beta[3]","2.5%"])*0 + (sum_survELRI$summary["beta[4]","2.5%"])*mean(ELRI_data$origin_01)
                                            + (sum_survELRI$summary["beta[5]","2.5%"])*xdummy*0))
ELRI_ydummy_eminus97.5 <- as.vector(invlogit(sum_survELRI$summary["beta[1]","97.5%"] + (sum_survELRI$summary["beta[2]","97.5%"])*xdummy + (sum_survELRI$summary["beta[3]","97.5%"])*0 + (sum_survELRI$summary["beta[4]","97.5%"])*mean(ELRI_data$origin_01)
                                             + (sum_survELRI$summary["beta[5]","97.5%"])*xdummy*0))

ELRI_ydummy_eminus_semin <- as.vector(invlogit((sum_survELRI$summary["beta[1]","mean"]-sum_survELRI$summary["beta[1]","se_mean"]) + (sum_survELRI$summary["beta[2]","mean"]-sum_survELRI$summary["beta[2]","se_mean"])*xdummy + (sum_survELRI$summary["beta[3]","mean"]-sum_survELRI$summary["beta[3]","se_mean"])*0 + (sum_survELRI$summary["beta[4]","mean"]-sum_survELRI$summary["beta[4]","se_mean"])*mean(ELRI_data$origin_01)
                                               + (sum_survELRI$summary["beta[5]","mean"]-sum_survELRI$summary["beta[5]","se_mean"])*xdummy*0))
ELRI_ydummy_eminus_semax <- as.vector(invlogit((sum_survELRI$summary["beta[1]","mean"]+sum_survELRI$summary["beta[1]","se_mean"]) + (sum_survELRI$summary["beta[2]","mean"]+sum_survELRI$summary["beta[2]","se_mean"])*xdummy + (sum_survELRI$summary["beta[3]","mean"]+sum_survELRI$summary["beta[3]","se_mean"])*0 + (sum_survELRI$summary["beta[4]","mean"]+sum_survELRI$summary["beta[4]","se_mean"])*mean(ELRI_data$origin_01)
                                               + (sum_survELRI$summary["beta[5]","mean"]+sum_survELRI$summary["beta[5]","se_mean"])*xdummy*0))


ELRI_ydummy_eminus_y1 <- as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*0 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_survELRI$'beta[5]')*xdummy*0
                                            + mean(post_survELRI$'tau_year[1,1]')))
ELRI_ydummy_eminus_y2 <- as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*0 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_survELRI$'beta[5]')*xdummy*0
                                            + mean(post_survELRI$'tau_year[1,2]')))
ELRI_ydummy_eminus_y3 <- as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + + mean(post_survELRI$'beta[3]')*0  + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_survELRI$'beta[5]')*xdummy*0
                                            + mean(post_survELRI$'tau_year[1,3]')))
ELRI_ydummy_eminus_y4 <- as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*0 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_survELRI$'beta[5]')*xdummy*0
                                            + mean(post_survELRI$'tau_year[1,4]'))) 
ELRI_ydummy_eminus_y5 <- as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*0 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_survELRI$'beta[5]')*xdummy*0
                                            + mean(post_survELRI$'tau_year[1,5]')))
ELRI_ydummy_eminus_y6 <- as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*0 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_survELRI$'beta[5]')*xdummy*0
                                            + mean(post_survELRI$'tau_year[1,6]')))
ELRI_ydummy_eminus_y7 <- as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*0 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_survELRI$'beta[5]')*xdummy*0
                                            + mean(post_survELRI$'tau_year[1,7]')))
ELRI_ydummy_eminus_y8 <- as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*0 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_survELRI$'beta[5]')*xdummy*0
                                            + mean(post_survELRI$'tau_year[1,8]')))
ELRI_ydummy_eminus_y9 <- as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*0 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_survELRI$'beta[5]')*xdummy*0
                                            + mean(post_survELRI$'tau_year[1,9]')))
ELRI_ydummy_eminus_y10 <- as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*0 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                             + mean(post_survELRI$'beta[5]')*xdummy*0
                                             + mean(post_survELRI$'tau_year[1,10]')))
ELRI_ydummy_eminus_y11 <- as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*0 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                             + mean(post_survELRI$'beta[5]')*xdummy*0
                                             + mean(post_survELRI$'tau_year[1,11]')))
ELRIfitsminus0 <- as.data.frame(cbind(xdummy, ELRI_ydummy_eminus))
ELRIfitsminus1 <- as.data.frame(cbind(xdummy,ELRI_ydummy_eminus_y1,ELRI_ydummy_eminus_y2, ELRI_ydummy_eminus_y3, ELRI_ydummy_eminus_y4, ELRI_ydummy_eminus_y5, ELRI_ydummy_eminus_y6, ELRI_ydummy_eminus_y7, ELRI_ydummy_eminus_y8, ELRI_ydummy_eminus_y9, ELRI_ydummy_eminus_y10, ELRI_ydummy_eminus_y11))
ELRIfitsminus <- melt(ELRIfitsminus1, id = "xdummy",
                      measure.vars = c("ELRI_ydummy_eminus_y1","ELRI_ydummy_eminus_y2", 
                                       "ELRI_ydummy_eminus_y3", "ELRI_ydummy_eminus_y4", 
                                       "ELRI_ydummy_eminus_y5", "ELRI_ydummy_eminus_y6", 
                                       "ELRI_ydummy_eminus_y7", "ELRI_ydummy_eminus_y8", 
                                       "ELRI_ydummy_eminus_y9", "ELRI_ydummy_eminus_y10",
                                       "ELRI_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "ELRI_ydummy_eminus_y1" = '1',"ELRI_ydummy_eminus_y2" = "2", 
                       "ELRI_ydummy_eminus_y3" = "3", "ELRI_ydummy_eminus_y4" = "4", 
                       "ELRI_ydummy_eminus_y5" = "5", "ELRI_ydummy_eminus_y6" = "6", 
                       "ELRI_ydummy_eminus_y7" = "7", "ELRI_ydummy_eminus_y8" = "8", 
                       "ELRI_ydummy_eminus_y9" = "9", "ELRI_ydummy_eminus_y10" = "10",
                       "ELRI_ydummy_eminus_y11" = "11") )



# E+
ELRI_ydummy_eplus <- as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*1 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                        + mean(post_survELRI$'beta[5]')*xdummy*1))
ELRI_ydummy_eplus2.5 <- as.vector(invlogit(sum_survELRI$summary["beta[1]","2.5%"] + (sum_survELRI$summary["beta[2]","2.5%"])*xdummy + (sum_survELRI$summary["beta[3]","2.5%"])*1 + (sum_survELRI$summary["beta[4]","2.5%"])*mean(ELRI_data$origin_01)
                                            + (sum_survELRI$summary["beta[5]","2.5%"])*xdummy*1))
ELRI_ydummy_eplus97.5 <- as.vector(invlogit(sum_survELRI$summary["beta[1]","97.5%"] + (sum_survELRI$summary["beta[2]","97.5%"])*xdummy + (sum_survELRI$summary["beta[3]","97.5%"])*1 + (sum_survELRI$summary["beta[4]","97.5%"])*mean(ELRI_data$origin_01)
                                             + (sum_survELRI$summary["beta[5]","97.5%"])*xdummy*1))

ELRI_ydummy_eplus_semin <- as.vector(invlogit((sum_survELRI$summary["beta[1]","mean"]-sum_survELRI$summary["beta[1]","se_mean"]) + (sum_survELRI$summary["beta[2]","mean"]-sum_survELRI$summary["beta[2]","se_mean"])*xdummy + (sum_survELRI$summary["beta[3]","mean"]-sum_survELRI$summary["beta[3]","se_mean"])*1 + (sum_survELRI$summary["beta[4]","mean"]-sum_survELRI$summary["beta[4]","se_mean"])*mean(ELRI_data$origin_01)
                                               + (sum_survELRI$summary["beta[5]","mean"]-sum_survELRI$summary["beta[5]","se_mean"])*xdummy*1))
ELRI_ydummy_eplus_semax <- as.vector(invlogit((sum_survELRI$summary["beta[1]","mean"]+sum_survELRI$summary["beta[1]","se_mean"]) + (sum_survELRI$summary["beta[2]","mean"]+sum_survELRI$summary["beta[2]","se_mean"])*xdummy + (sum_survELRI$summary["beta[3]","mean"]+sum_survELRI$summary["beta[3]","se_mean"])*1 + (sum_survELRI$summary["beta[4]","mean"]+sum_survELRI$summary["beta[4]","se_mean"])*mean(ELRI_data$origin_01)
                                               + (sum_survELRI$summary["beta[5]","mean"]+sum_survELRI$summary["beta[5]","se_mean"])*xdummy*1))


ELRI_ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*1 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_survELRI$'beta[5]')*xdummy*1
                                            + mean(post_survELRI$'tau_year[2,1]')))
ELRI_ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*1 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_survELRI$'beta[5]')*xdummy*1
                                            + mean(post_survELRI$'tau_year[2,2')))
ELRI_ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*1 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_survELRI$'beta[5]')*xdummy*1
                                            + mean(post_survELRI$'tau_year[2,3]')))
ELRI_ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*1 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_survELRI$'beta[5]')*xdummy*1
                                            + mean(post_survELRI$'tau_year[2,4]')))
ELRI_ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*1 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_survELRI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_survELRI$'tau_year[2,5]')))
ELRI_ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*1 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_survELRI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_survELRI$'tau_year[2,6]')))
ELRI_ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*1 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_survELRI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_survELRI$'tau_year[2,7]')))
ELRI_ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*1 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_survELRI$'beta[5]')*xdummy*1
                                            + mean(post_survELRI$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
ELRI_ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*1 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_survELRI$'beta[5]')*xdummy*1
                                            + mean(post_survELRI$'tau_year[2,9]')))
ELRI_ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*1 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                             + mean(post_survELRI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_survELRI$'tau_year[2,10]')))
ELRI_ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_survELRI$'beta[1]') + mean(post_survELRI$'beta[2]')*xdummy + mean(post_survELRI$'beta[3]')*1 + mean(post_survELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                             + mean(post_survELRI$'beta[5]')*xdummy*1
                                             + mean(post_survELRI$'tau_year[2,11]')))
ELRI_fitsplus0 <- as.data.frame(cbind(xdummy, ELRI_ydummy_eplus))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
ELRIfitsplus1 <- as.data.frame(cbind(xdummy,ELRI_ydummy_eplus_y1,ELRI_ydummy_eplus_y2, ELRI_ydummy_eplus_y3, ELRI_ydummy_eplus_y4, ELRI_ydummy_eplus_y5, ELRI_ydummy_eplus_y6, ELRI_ydummy_eplus_y7, ELRI_ydummy_eplus_y8, ELRI_ydummy_eplus_y9, ELRI_ydummy_eplus_y10, ELRI_ydummy_eplus_y11))
ELRIfitsplus <- melt(ELRIfitsplus1, id = "xdummy",
                     measure.vars = c("ELRI_ydummy_eplus_y1","ELRI_ydummy_eplus_y2", 
                                      "ELRI_ydummy_eplus_y3", "ELRI_ydummy_eplus_y4", 
                                      "ELRI_ydummy_eplus_y5", "ELRI_ydummy_eplus_y6", 
                                      "ELRI_ydummy_eplus_y7", "ELRI_ydummy_eplus_y8", 
                                      "ELRI_ydummy_eplus_y9", "ELRI_ydummy_eplus_y10",
                                      "ELRI_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "ELRI_ydummy_eplus_y1" = '1',"ELRI_ydummy_eplus_y2" = "2", 
                       "ELRI_ydummy_eplus_y3" = "3", "ELRI_ydummy_eplus_y4" = "4", 
                       "ELRI_ydummy_eplus_y5" = "5", "ELRI_ydummy_eplus_y6" = "6", 
                       "ELRI_ydummy_eplus_y7" = "7", "ELRI_ydummy_eplus_y8" = "8", 
                       "ELRI_ydummy_eplus_y9" = "9", "ELRI_ydummy_eplus_y10" = "10",
                       "ELRI_ydummy_eplus_y11" = "11") )

# E+ by year
ELRIEplusbyyear <- ggplot(data = ELRIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = ELRIsurv_bin1, aes(x = mean_size, y = mean_surv, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("Prob. of Survival") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = ELRIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("Prob. of Survival")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
ELRIEminusbyyear <- ggplot(data = ELRIfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = ELRIsurv_bin0, aes(x = mean_size, y = mean_surv, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("Prob. of Survival")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = ELRIfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("ELRI E- Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)




# effect on mean


ELRIfits <- as.data.frame(cbind(xdummy, ELRI_ydummy_eplus, ELRI_ydummy_eplus_semin, ELRI_ydummy_eplus_semax, ELRI_ydummy_eplus2.5, ELRI_ydummy_eplus97.5, ELRI_ydummy_eminus, ELRI_ydummy_eminus_semin, ELRI_ydummy_eminus_semax, ELRI_ydummy_eminus2.5, ELRI_ydummy_eminus97.5))


# without point
ggplot(data = ELRIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("ELRI Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = ELRIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = ELRIsurv_binmean, aes(x = mean_size, y = mean_surv, color = endo), lwd = 3) +
  ggtitle("ELRI Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

ELRImean <- ggplot(data = ELRIfits) +
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = ELRIsurv_binmean, aes(x = mean_size, y = mean_surv, color = endo, lwd = samplesize)) +
  ggtitle("Mean Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")

ELRIsurv4panel <- grid.arrange(ELRImean, ELRI_s,ELRIEplusbyyear,ELRIEminusbyyear, ncol= 4)
titleELRIsurv4panel <- annotate_figure(ELRIsurv4panel, top = "Elymus riparius")

titleELRIsurv4panel

# split up plus and minus with datapoints and add mean line
ELRI_minus <- ggplot(data = ELRIfitsminus1)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = ELRIfitsminus0, aes(x = xdummy, y = ELRI_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = ELRIsurv_bin0, aes(x = mean_size, y = mean_surv), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


ELRI_plus <- ggplot(data = ELRIfitsplus1)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = ELRI_fitsplus0, aes(x = xdummy, y = ELRI_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = ELRIsurv_bin1, aes(x = mean_size, y = mean_surv), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

ELRIsurv <- grid.arrange(ELRI_plus, ELRI_minus, ncol= 2)
titleELRIsurv <- annotate_figure(ELRIsurv, top = "ELRI", left = "Probability of Survival")

titleELRIsurv





########### AGPE model fits
# actual data points for probability of survival
AGPEsurv_bin0 <- AGPE_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  filter(endo == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())

AGPEsurv_bin1 <- AGPE_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  filter(endo == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())
AGPEsurv_bin <- AGPE_data %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, "0" = 0, "1" = 1, minus = 0, plus = 1)) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())


AGPEsurv_binmean <- AGPE_data %>% 
  select(surv_t1, logsize_t, endo) %>%
  mutate(endo = recode(endo, "0" = "E-", "1" = "E+", minus = "E-", plus = "E+")) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(AGPE_data$logsize_t)
xdummy <- seq(min(AGPE_data$logsize_t), max(AGPE_data$logsize_t), length.out = nvalues)

# E- 
AGPE_ydummy_eminus <- as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*0 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                         + mean(post_survAGPE$'beta[5]')*xdummy*0))
AGPE_ydummy_eminus2.5 <- as.vector(invlogit(sum_survAGPE$summary["beta[1]","2.5%"] + (sum_survAGPE$summary["beta[2]","2.5%"])*xdummy + (sum_survAGPE$summary["beta[3]","2.5%"])*0 + (sum_survAGPE$summary["beta[4]","2.5%"])*mean(AGPE_data$origin_01)
                                            + (sum_survAGPE$summary["beta[5]","2.5%"])*xdummy*0))
AGPE_ydummy_eminus97.5 <- as.vector(invlogit(sum_survAGPE$summary["beta[1]","97.5%"] + (sum_survAGPE$summary["beta[2]","97.5%"])*xdummy + (sum_survAGPE$summary["beta[3]","97.5%"])*0 + (sum_survAGPE$summary["beta[4]","97.5%"])*mean(AGPE_data$origin_01)
                                             + (sum_survAGPE$summary["beta[5]","97.5%"])*xdummy*0))

AGPE_ydummy_eminus_semin <- as.vector(invlogit((sum_survAGPE$summary["beta[1]","mean"]-sum_survAGPE$summary["beta[1]","se_mean"]) + (sum_survAGPE$summary["beta[2]","mean"]-sum_survAGPE$summary["beta[2]","se_mean"])*xdummy + (sum_survAGPE$summary["beta[3]","mean"]-sum_survAGPE$summary["beta[3]","se_mean"])*0 + (sum_survAGPE$summary["beta[4]","mean"]-sum_survAGPE$summary["beta[4]","se_mean"])*mean(AGPE_data$origin_01)
                                               + (sum_survAGPE$summary["beta[5]","mean"]-sum_survAGPE$summary["beta[5]","se_mean"])*xdummy*0))
AGPE_ydummy_eminus_semax <- as.vector(invlogit((sum_survAGPE$summary["beta[1]","mean"]+sum_survAGPE$summary["beta[1]","se_mean"]) + (sum_survAGPE$summary["beta[2]","mean"]+sum_survAGPE$summary["beta[2]","se_mean"])*xdummy + (sum_survAGPE$summary["beta[3]","mean"]+sum_survAGPE$summary["beta[3]","se_mean"])*0 + (sum_survAGPE$summary["beta[4]","mean"]+sum_survAGPE$summary["beta[4]","se_mean"])*mean(AGPE_data$origin_01)
                                               + (sum_survAGPE$summary["beta[5]","mean"]+sum_survAGPE$summary["beta[5]","se_mean"])*xdummy*0))


AGPE_ydummy_eminus_y1 <- as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*0 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_survAGPE$'beta[5]')*xdummy*0
                                            + mean(post_survAGPE$'tau_year[1,1]')))
AGPE_ydummy_eminus_y2 <- as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*0 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_survAGPE$'beta[5]')*xdummy*0
                                            + mean(post_survAGPE$'tau_year[1,2]')))
AGPE_ydummy_eminus_y3 <- as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + + mean(post_survAGPE$'beta[3]')*0  + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_survAGPE$'beta[5]')*xdummy*0
                                            + mean(post_survAGPE$'tau_year[1,3]')))
AGPE_ydummy_eminus_y4 <- as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*0 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_survAGPE$'beta[5]')*xdummy*0
                                            + mean(post_survAGPE$'tau_year[1,4]'))) 
AGPE_ydummy_eminus_y5 <- as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*0 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_survAGPE$'beta[5]')*xdummy*0
                                            + mean(post_survAGPE$'tau_year[1,5]')))
AGPE_ydummy_eminus_y6 <- as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*0 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_survAGPE$'beta[5]')*xdummy*0
                                            + mean(post_survAGPE$'tau_year[1,6]')))
AGPE_ydummy_eminus_y7 <- as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*0 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_survAGPE$'beta[5]')*xdummy*0
                                            + mean(post_survAGPE$'tau_year[1,7]')))
AGPE_ydummy_eminus_y8 <- as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*0 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_survAGPE$'beta[5]')*xdummy*0
                                            + mean(post_survAGPE$'tau_year[1,8]')))
AGPE_ydummy_eminus_y9 <- as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*0 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_survAGPE$'beta[5]')*xdummy*0
                                            + mean(post_survAGPE$'tau_year[1,9]')))
AGPE_ydummy_eminus_y10 <- as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*0 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                             + mean(post_survAGPE$'beta[5]')*xdummy*0
                                             + mean(post_survAGPE$'tau_year[1,10]')))
AGPE_ydummy_eminus_y11 <- as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*0 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                             + mean(post_survAGPE$'beta[5]')*xdummy*0
                                             + mean(post_survAGPE$'tau_year[1,11]')))
AGPEfitsminus0 <- as.data.frame(cbind(xdummy, AGPE_ydummy_eminus))
AGPEfitsminus1 <- as.data.frame(cbind(xdummy,AGPE_ydummy_eminus_y1,AGPE_ydummy_eminus_y2, AGPE_ydummy_eminus_y3, AGPE_ydummy_eminus_y4, AGPE_ydummy_eminus_y5, AGPE_ydummy_eminus_y6, AGPE_ydummy_eminus_y7, AGPE_ydummy_eminus_y8, AGPE_ydummy_eminus_y9, AGPE_ydummy_eminus_y10, AGPE_ydummy_eminus_y11))
AGPEfitsminus <- melt(AGPEfitsminus1, id = "xdummy",
                      measure.vars = c("AGPE_ydummy_eminus_y1","AGPE_ydummy_eminus_y2", 
                                       "AGPE_ydummy_eminus_y3", "AGPE_ydummy_eminus_y4", 
                                       "AGPE_ydummy_eminus_y5", "AGPE_ydummy_eminus_y6", 
                                       "AGPE_ydummy_eminus_y7", "AGPE_ydummy_eminus_y8", 
                                       "AGPE_ydummy_eminus_y9", "AGPE_ydummy_eminus_y10",
                                       "AGPE_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "AGPE_ydummy_eminus_y1" = '1',"AGPE_ydummy_eminus_y2" = "2", 
                       "AGPE_ydummy_eminus_y3" = "3", "AGPE_ydummy_eminus_y4" = "4", 
                       "AGPE_ydummy_eminus_y5" = "5", "AGPE_ydummy_eminus_y6" = "6", 
                       "AGPE_ydummy_eminus_y7" = "7", "AGPE_ydummy_eminus_y8" = "8", 
                       "AGPE_ydummy_eminus_y9" = "9", "AGPE_ydummy_eminus_y10" = "10",
                       "AGPE_ydummy_eminus_y11" = "11") )



# E+
AGPE_ydummy_eplus <- as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*1 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                        + mean(post_survAGPE$'beta[5]')*xdummy*1))
AGPE_ydummy_eplus2.5 <- as.vector(invlogit(sum_survAGPE$summary["beta[1]","2.5%"] + (sum_survAGPE$summary["beta[2]","2.5%"])*xdummy + (sum_survAGPE$summary["beta[3]","2.5%"])*1 + (sum_survAGPE$summary["beta[4]","2.5%"])*mean(AGPE_data$origin_01)
                                            + (sum_survAGPE$summary["beta[5]","2.5%"])*xdummy*1))
AGPE_ydummy_eplus97.5 <- as.vector(invlogit(sum_survAGPE$summary["beta[1]","97.5%"] + (sum_survAGPE$summary["beta[2]","97.5%"])*xdummy + (sum_survAGPE$summary["beta[3]","97.5%"])*1 + (sum_survAGPE$summary["beta[4]","97.5%"])*mean(AGPE_data$origin_01)
                                             + (sum_survAGPE$summary["beta[5]","97.5%"])*xdummy*1))

AGPE_ydummy_eplus_semin <- as.vector(invlogit((sum_survAGPE$summary["beta[1]","mean"]-sum_survAGPE$summary["beta[1]","se_mean"]) + (sum_survAGPE$summary["beta[2]","mean"]-sum_survAGPE$summary["beta[2]","se_mean"])*xdummy + (sum_survAGPE$summary["beta[3]","mean"]-sum_survAGPE$summary["beta[3]","se_mean"])*1 + (sum_survAGPE$summary["beta[4]","mean"]-sum_survAGPE$summary["beta[4]","se_mean"])*mean(AGPE_data$origin_01)
                                               + (sum_survAGPE$summary["beta[5]","mean"]-sum_survAGPE$summary["beta[5]","se_mean"])*xdummy*1))
AGPE_ydummy_eplus_semax <- as.vector(invlogit((sum_survAGPE$summary["beta[1]","mean"]+sum_survAGPE$summary["beta[1]","se_mean"]) + (sum_survAGPE$summary["beta[2]","mean"]+sum_survAGPE$summary["beta[2]","se_mean"])*xdummy + (sum_survAGPE$summary["beta[3]","mean"]+sum_survAGPE$summary["beta[3]","se_mean"])*1 + (sum_survAGPE$summary["beta[4]","mean"]+sum_survAGPE$summary["beta[4]","se_mean"])*mean(AGPE_data$origin_01)
                                               + (sum_survAGPE$summary["beta[5]","mean"]+sum_survAGPE$summary["beta[5]","se_mean"])*xdummy*1))


AGPE_ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*1 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_survAGPE$'beta[5]')*xdummy*1
                                            + mean(post_survAGPE$'tau_year[2,1]')))
AGPE_ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*1 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_survAGPE$'beta[5]')*xdummy*1
                                            + mean(post_survAGPE$'tau_year[2,2')))
AGPE_ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*1 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_survAGPE$'beta[5]')*xdummy*1
                                            + mean(post_survAGPE$'tau_year[2,3]')))
AGPE_ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*1 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_survAGPE$'beta[5]')*xdummy*1
                                            + mean(post_survAGPE$'tau_year[2,4]')))
AGPE_ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*1 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_survAGPE$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_survAGPE$'tau_year[2,5]')))
AGPE_ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*1 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_survAGPE$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_survAGPE$'tau_year[2,6]')))
AGPE_ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*1 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_survAGPE$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_survAGPE$'tau_year[2,7]')))
AGPE_ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*1 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_survAGPE$'beta[5]')*xdummy*1
                                            + mean(post_survAGPE$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
AGPE_ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*1 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_survAGPE$'beta[5]')*xdummy*1
                                            + mean(post_survAGPE$'tau_year[2,9]')))
AGPE_ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*1 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                             + mean(post_survAGPE$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_survAGPE$'tau_year[2,10]')))
AGPE_ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_survAGPE$'beta[1]') + mean(post_survAGPE$'beta[2]')*xdummy + mean(post_survAGPE$'beta[3]')*1 + mean(post_survAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                             + mean(post_survAGPE$'beta[5]')*xdummy*1
                                             + mean(post_survAGPE$'tau_year[2,11]')))
AGPE_fitsplus0 <- as.data.frame(cbind(xdummy, AGPE_ydummy_eplus))     
AGPEfitsplus1 <- as.data.frame(cbind(xdummy,AGPE_ydummy_eplus_y1,AGPE_ydummy_eplus_y2, AGPE_ydummy_eplus_y3, AGPE_ydummy_eplus_y4, AGPE_ydummy_eplus_y5, AGPE_ydummy_eplus_y6, AGPE_ydummy_eplus_y7, AGPE_ydummy_eplus_y8, AGPE_ydummy_eplus_y9, AGPE_ydummy_eplus_y10, AGPE_ydummy_eplus_y11))
AGPEfitsplus <- melt(AGPEfitsplus1, id = "xdummy",
                     measure.vars = c("AGPE_ydummy_eplus_y1","AGPE_ydummy_eplus_y2", 
                                      "AGPE_ydummy_eplus_y3", "AGPE_ydummy_eplus_y4", 
                                      "AGPE_ydummy_eplus_y5", "AGPE_ydummy_eplus_y6", 
                                      "AGPE_ydummy_eplus_y7", "AGPE_ydummy_eplus_y8", 
                                      "AGPE_ydummy_eplus_y9", "AGPE_ydummy_eplus_y10",
                                      "AGPE_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "AGPE_ydummy_eplus_y1" = '1',"AGPE_ydummy_eplus_y2" = "2", 
                       "AGPE_ydummy_eplus_y3" = "3", "AGPE_ydummy_eplus_y4" = "4", 
                       "AGPE_ydummy_eplus_y5" = "5", "AGPE_ydummy_eplus_y6" = "6", 
                       "AGPE_ydummy_eplus_y7" = "7", "AGPE_ydummy_eplus_y8" = "8", 
                       "AGPE_ydummy_eplus_y9" = "9", "AGPE_ydummy_eplus_y10" = "10",
                       "AGPE_ydummy_eplus_y11" = "11") )

# E+ by year
AGPEEplusbyyear <- ggplot(data = AGPEfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = AGPEsurv_bin1, aes(x = mean_size, y = mean_surv, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("Prob. of Survival") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = AGPEfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("Prob. of Survival")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
AGPEEminusbyyear <- ggplot(data = AGPEfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = AGPEsurv_bin0, aes(x = mean_size, y = mean_surv, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("Prob. of Survival")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = AGPEfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("AGPE E- Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)




# effect on mean


AGPEfits <- as.data.frame(cbind(xdummy, AGPE_ydummy_eplus, AGPE_ydummy_eplus_semin, AGPE_ydummy_eplus_semax, AGPE_ydummy_eplus2.5, AGPE_ydummy_eplus97.5, AGPE_ydummy_eminus, AGPE_ydummy_eminus_semin, AGPE_ydummy_eminus_semax, AGPE_ydummy_eminus2.5, AGPE_ydummy_eminus97.5))


# without point
ggplot(data = AGPEfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("AGPE Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = AGPEfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = AGPEsurv_binmean, aes(x = mean_size, y = mean_surv, color = endo), lwd = 3) +
  ggtitle("AGPE Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

AGPEmean <- ggplot(data = AGPEfits) +
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = AGPEsurv_binmean, aes(x = mean_size, y = mean_surv, color = endo, lwd = samplesize)) +
  ggtitle("Mean Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")


AGPEsurv4panel <- grid.arrange(AGPEmean, AGPE_s,AGPEEplusbyyear,AGPEEminusbyyear, ncol= 4)
titleAGPEsurv4panel <- annotate_figure(AGPEsurv4panel, top = "Agrostis perennans")

titleAGPEsurv4panel

# split up plus and minus with datapoints and add mean line
AGPE_minus <- ggplot(data = AGPEfitsminus1)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = AGPEfitsminus0, aes(x = xdummy, y = AGPE_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = AGPEsurv_bin0, aes(x = mean_size, y = mean_surv), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


AGPE_plus <- ggplot(data = AGPEfitsplus1)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = AGPE_fitsplus0, aes(x = xdummy, y = AGPE_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = AGPEsurv_bin1, aes(x = mean_size, y = mean_surv), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

AGPEsurv <- grid.arrange(AGPE_plus, AGPE_minus, ncol= 2)
titleAGPEsurv <- annotate_figure(AGPEsurv, top = "AGPE", left = "Probability of Survival")

titleAGPEsurv




survfits <- grid.arrange(titleAGPEsurv,titleELRIsurv,titleELVIsurv,titleFESUsurv,titlePOALsurv,titlePOSYsurv, nrow = 3)
survtitle <- annotate_figure(survfits, top = "Survival Probability")
survtitle

surv4panel <- grid.arrange(titleAGPEsurv4panel,titleELRIsurv4panel,titleELVIsurv4panel,titleFESUsurv4panel,titleLOARsurv4panel,titlePOALsurv4panel,titlePOSYsurv4panel, nrow = 7)
surv4paneltitle <- annotate_figure(surv4panel, top = "Survival Probability")

surv4paneltitle






######### plots of Growth model fits#####

######## POAL model fits
# actual data points for size_t1
POALgrow_bin0 <- POAL_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  mutate(endo_01 = recode(endo_01, minus = 0, plus = 1)) %>% 
  filter(endo_01 == 0)


POALgrow_bin1 <- POAL_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) 

POALgrow_bin0mean <- POAL_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 8)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())



POALgrow_bin1mean <- POAL_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 8)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())






POALgrow_binmean <- POAL_data %>% 
  select(size_t1, logsize_t, endo_01) %>%
  mutate(endo_01 = as.character(endo_01)) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(POAL_data$logsize_t)
xdummy <- seq(min(POAL_data$logsize_t), max(POAL_data$logsize_t), length.out = nvalues)

# E- 
POAL_ydummy_eminus <- as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*0 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                         + mean(post_growPOAL$'beta[5]')*xdummy*0))
POAL_ydummy_eminus_y1 <- as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*0 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_growPOAL$'beta[5]')*xdummy*0
                                            + mean(post_growPOAL$'tau_year[1,1]')))
POAL_ydummy_eminus_y2 <- as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*0 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_growPOAL$'beta[5]')*xdummy*0
                                            + mean(post_growPOAL$'tau_year[1,2]')))
POAL_ydummy_eminus_y3 <- as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + + mean(post_growPOAL$'beta[3]')*0  + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_growPOAL$'beta[5]')*xdummy*0
                                            + mean(post_growPOAL$'tau_year[1,3]')))
POAL_ydummy_eminus_y4 <- as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*0 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_growPOAL$'beta[5]')*xdummy*0
                                            + mean(post_growPOAL$'tau_year[1,4]'))) 
POAL_ydummy_eminus_y5 <- as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*0 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_growPOAL$'beta[5]')*xdummy*0
                                            + mean(post_growPOAL$'tau_year[1,5]')))
POAL_ydummy_eminus_y6 <- as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*0 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_growPOAL$'beta[5]')*xdummy*0
                                            + mean(post_growPOAL$'tau_year[1,6]')))
POAL_ydummy_eminus_y7 <- as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*0 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_growPOAL$'beta[5]')*xdummy*0
                                            + mean(post_growPOAL$'tau_year[1,7]')))
POAL_ydummy_eminus_y8 <- as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*0 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_growPOAL$'beta[5]')*xdummy*0
                                            + mean(post_growPOAL$'tau_year[1,8]')))
POAL_ydummy_eminus_y9 <- as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*0 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_growPOAL$'beta[5]')*xdummy*0
                                            + mean(post_growPOAL$'tau_year[1,9]')))
POAL_ydummy_eminus_y10 <- as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*0 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                             + mean(post_growPOAL$'beta[5]')*xdummy*0
                                             + mean(post_growPOAL$'tau_year[1,10]')))
POAL_ydummy_eminus_y11 <- as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*0 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                             + mean(post_growPOAL$'beta[5]')*xdummy*0
                                             + mean(post_growPOAL$'tau_year[1,11]')))
POALfitsminus0 <- as.data.frame(cbind(xdummy, POAL_ydummy_eminus))
POALfitsminus1 <- as.data.frame(cbind(xdummy,POAL_ydummy_eminus_y1,POAL_ydummy_eminus_y2, POAL_ydummy_eminus_y3, POAL_ydummy_eminus_y4, POAL_ydummy_eminus_y5, POAL_ydummy_eminus_y6, POAL_ydummy_eminus_y7, POAL_ydummy_eminus_y8, POAL_ydummy_eminus_y9, POAL_ydummy_eminus_y10, POAL_ydummy_eminus_y11))
POALfitsminus <- melt(POALfitsminus1, id = "xdummy",
                      measure.vars = c("POAL_ydummy_eminus_y1","POAL_ydummy_eminus_y2", 
                                       "POAL_ydummy_eminus_y3", "POAL_ydummy_eminus_y4", 
                                       "POAL_ydummy_eminus_y5", "POAL_ydummy_eminus_y6", 
                                       "POAL_ydummy_eminus_y7", "POAL_ydummy_eminus_y8", 
                                       "POAL_ydummy_eminus_y9", "POAL_ydummy_eminus_y10",
                                       "POAL_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "POAL_ydummy_eminus_y1" = '1',"POAL_ydummy_eminus_y2" = "2", 
                       "POAL_ydummy_eminus_y3" = "3", "POAL_ydummy_eminus_y4" = "4", 
                       "POAL_ydummy_eminus_y5" = "5", "POAL_ydummy_eminus_y6" = "6", 
                       "POAL_ydummy_eminus_y7" = "7", "POAL_ydummy_eminus_y8" = "8", 
                       "POAL_ydummy_eminus_y9" = "9", "POAL_ydummy_eminus_y10" = "10",
                       "POAL_ydummy_eminus_y11" = "11") )



# E+
POAL_ydummy_eplus <- as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*1 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                        + mean(post_growPOAL$'beta[5]')*xdummy*1))
POAL_ydummy_eplus_y1 <-  as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*1 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_growPOAL$'beta[5]')*xdummy*1
                                            + mean(post_growPOAL$'tau_year[2,1]')))
POAL_ydummy_eplus_y2 <-  as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*1 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_growPOAL$'beta[5]')*xdummy*1
                                            + mean(post_growPOAL$'tau_year[2,2')))
POAL_ydummy_eplus_y3 <-  as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*1 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_growPOAL$'beta[5]')*xdummy*1
                                            + mean(post_growPOAL$'tau_year[2,3]')))
POAL_ydummy_eplus_y4 <-  as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*1 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_growPOAL$'beta[5]')*xdummy*1
                                            + mean(post_growPOAL$'tau_year[2,4]')))
POAL_ydummy_eplus_y5 <-  as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*1 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_growPOAL$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_growPOAL$'tau_year[2,5]')))
POAL_ydummy_eplus_y6 <-  as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*1 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_growPOAL$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_growPOAL$'tau_year[2,6]')))
POAL_ydummy_eplus_y7 <-  as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*1 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_growPOAL$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_growPOAL$'tau_year[2,7]')))
POAL_ydummy_eplus_y8 <-  as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*1 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_growPOAL$'beta[5]')*xdummy*1
                                            + mean(post_growPOAL$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
POAL_ydummy_eplus_y9 <-  as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*1 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_growPOAL$'beta[5]')*xdummy*1
                                            + mean(post_growPOAL$'tau_year[2,9]')))
POAL_ydummy_eplus_y10 <-  as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*1 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                             + mean(post_growPOAL$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_growPOAL$'tau_year[2,10]')))
POAL_ydummy_eplus_y11 <-  as.vector(exp(mean(post_growPOAL$'beta[1]') + mean(post_growPOAL$'beta[2]')*xdummy + mean(post_growPOAL$'beta[3]')*1 + mean(post_growPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                             + mean(post_growPOAL$'beta[5]')*xdummy*1
                                             + mean(post_growPOAL$'tau_year[2,11]')))
POAL_fitsplus0 <- as.data.frame(cbind(xdummy, POAL_ydummy_eplus))
POALfitsplus1 <- as.data.frame(cbind(xdummy,POAL_ydummy_eplus_y1,POAL_ydummy_eplus_y2, POAL_ydummy_eplus_y3, POAL_ydummy_eplus_y4, POAL_ydummy_eplus_y5, POAL_ydummy_eplus_y6, POAL_ydummy_eplus_y7, POAL_ydummy_eplus_y8, POAL_ydummy_eplus_y9, POAL_ydummy_eplus_y10, POAL_ydummy_eplus_y11))
POALfitsplus <- melt(POALfitsplus1, id = "xdummy",
                     measure.vars = c("POAL_ydummy_eplus_y1","POAL_ydummy_eplus_y2", 
                                      "POAL_ydummy_eplus_y3", "POAL_ydummy_eplus_y4", 
                                      "POAL_ydummy_eplus_y5", "POAL_ydummy_eplus_y6", 
                                      "POAL_ydummy_eplus_y7", "POAL_ydummy_eplus_y8", 
                                      "POAL_ydummy_eplus_y9", "POAL_ydummy_eplus_y10",
                                      "POAL_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "POAL_ydummy_eplus_y1" = '1',"POAL_ydummy_eplus_y2" = "2", 
                       "POAL_ydummy_eplus_y3" = "3", "POAL_ydummy_eplus_y4" = "4", 
                       "POAL_ydummy_eplus_y5" = "5", "POAL_ydummy_eplus_y6" = "6", 
                       "POAL_ydummy_eplus_y7" = "7", "POAL_ydummy_eplus_y8" = "8", 
                       "POAL_ydummy_eplus_y9" = "9", "POAL_ydummy_eplus_y10" = "10",
                       "POAL_ydummy_eplus_y11" = "11") )

ggplot(data = POALfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8)+
  geom_point(data = POALgrow_bin1, aes(x = logsize_t, y = size_t1), lwd = 2) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
  ggtitle("POAL E+ Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth")  
scale_color_manual(values=yearcolors)

ggplot(data = POALfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("POAL E+ Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# effect on mean


POALfits <- as.data.frame(cbind(xdummy, POAL_ydummy_eminus, POAL_ydummy_eplus))
POALfits <- melt(POALfits, id = "xdummy", measure.vars = c("POAL_ydummy_eplus", "POAL_ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "POAL_ydummy_eminus" = "E-", "POAL_ydummy_eplus" = "E+"))

# without point
ggplot(data = POALfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("POAL Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = POALfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = POALgrow_binmean, aes(x = mean_size, y = mean_grow, color = endo_01), lwd = 3) +
  ggtitle("POAL Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth") +   
  scale_color_manual(values=colors2)


# split up plus and minus with datapoints and add mean line
POAL_minus <- ggplot(data = POALfitsminus1)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = POALfitsminus0, aes(x = xdummy, y = POAL_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = POALgrow_bin0mean, aes(x = mean_size, y = mean_grow), color ="#ff7f00") +
  theme_classic()+ ylim(0,80)+
  labs(title = "E-", x = "log(size_t)", y = "")


POAL_plus <- ggplot(data = POALfitsplus1)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = POAL_fitsplus0, aes(x = xdummy, y = POAL_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = POALgrow_bin1mean, aes(x = mean_size, y = mean_grow), color ="#6a3d9a") +
  theme_classic()+ ylim(0,80)+
  labs(title = "E+", x = "log(size_t)", y = "")

# E+ by year
POALEplusbyyear <- ggplot(data = POALfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = POALgrow_bin1mean, aes(x = mean_size, y = mean_grow, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("size_t1") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = POALfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("size_t1")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
POALEminusbyyear <- ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = POALgrow_bin0mean, aes(x = mean_size, y = mean_grow, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("size_t1")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("POAL E- Growth Probability") + xlab("log(size_t)") + ylab("size_t1")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)


# effect on mean


POALfits <- as.data.frame(cbind(xdummy, POAL_ydummy_eplus,POAL_ydummy_eminus))

# without point
ggplot(data = POALfits) +
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus), color = "#ff7f00") +
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus), color = "#6a3d9a")+
  ggtitle("POAL Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +theme_classic()

ggplot(data = POALfits) +
  geom_ribbon(aes(x = xdummy, ymin = POAL_ydummy_eplus_semin, ymax = POAL_ydummy_eplus_semax), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = POAL_ydummy_eminus_semin, ymax = POAL_ydummy_eminus_semax), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus), color = "#ff7f00") +
  ggtitle("POAL Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +theme_classic()

# with points
POALmean <- ggplot(data = POALfits) +
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = POALgrow_binmean, aes(x = mean_size, y = mean_grow, color = endo_01, lwd = samplesize)) +
  ggtitle("Mean Growth") + xlab("log(size_t)") + ylab("size_t1") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")




POALgrow <- grid.arrange(POAL_plus, POAL_minus, ncol= 2)


titlePOALgrow <- annotate_figure(POALgrow, top = "POAL", left = "size_t1")

titlePOALgrow



# with data points
ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .7) +
  geom_point(data = POALgrow_bin0, aes(x = mean_size, y = mean_grow, color = Year)) +
  labs(title = "POAL E- Growth Probability", x = "log(size_t)", y = "size_t1") +
  scale_color_manual(values=yearcolors)
# without datapoints
ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .7) +
  labs(title = "POAL E- Growth Probability", x = "log(size_t)", y = "size_t1") +                                                                                                                                                                                                                                                                                                                                                                                                                                     scale_color_manual(values=yearcolors)






POALgrow4panel <- grid.arrange(POALmean, POAL_g,POALEplusbyyear,POALEminusbyyear, ncol= 4)
titlePOALgrow4panel <- annotate_figure(POALgrow4panel, top = "Poa alsodes")

titlePOALgrow4panel


######## POSY model fits
# actual data points for size_t1
POSYgrow_bin0 <- POSY_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) +
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())

POSYgrow_bin1 <- POSY_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  mutate(endo_01 = recode(endo_01, minus = 0, plus = 1)) %>% 
  filter(endo_01 == 1) 

POALgrow_bin1 <- POAL_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) 

POSYgrow_bin0mean <- POSY_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 8)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())



POSYgrow_bin1mean <- POSY_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 8)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())






POSYgrow_bin <- POSY_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())


POSYgrow_binmean <- POSY_data %>% 
  select(size_t1, logsize_t, endo_01) %>%
  mutate(endo_01 = as.character(endo_01)) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())



POSYgrow_binmean_0 <- POSY_data %>% 
  select(size_t1, logsize_t, endo_01) %>%
  mutate(size_bin = cut_interval(logsize_t, n = 30)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n()) %>% 
  filter(endo_01 == "E-")



POSYgrow_binmean_1 <- POSY_data %>% 
  select(size_t1, logsize_t, endo_01) %>%
  mutate(endo_01 = recode(endo_01, minus = 'E-', plus = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 30)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n()) %>% 
  filter(endo_01 == "E+")



# fitted model line
nvalues <- length(POSY_data$logsize_t)
xdummy <- seq(min(POSY_data$logsize_t), max(POSY_data$logsize_t), length.out = nvalues)

# E- 
POSY_ydummy_eminus <- as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*0 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                         + mean(post_growPOSY$'beta[5]')*xdummy*0))
POSY_ydummy_eminus_y1 <- as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*0 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_growPOSY$'beta[5]')*xdummy*0
                                            + mean(post_growPOSY$'tau_year[1,1]')))
POSY_ydummy_eminus_y2 <- as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*0 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_growPOSY$'beta[5]')*xdummy*0
                                            + mean(post_growPOSY$'tau_year[1,2]')))
POSY_ydummy_eminus_y3 <- as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + + mean(post_growPOSY$'beta[3]')*0  + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_growPOSY$'beta[5]')*xdummy*0
                                            + mean(post_growPOSY$'tau_year[1,3]')))
POSY_ydummy_eminus_y4 <- as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*0 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_growPOSY$'beta[5]')*xdummy*0
                                            + mean(post_growPOSY$'tau_year[1,4]'))) 
POSY_ydummy_eminus_y5 <- as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*0 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_growPOSY$'beta[5]')*xdummy*0
                                            + mean(post_growPOSY$'tau_year[1,5]')))
POSY_ydummy_eminus_y6 <- as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*0 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_growPOSY$'beta[5]')*xdummy*0
                                            + mean(post_growPOSY$'tau_year[1,6]')))
POSY_ydummy_eminus_y7 <- as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*0 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_growPOSY$'beta[5]')*xdummy*0
                                            + mean(post_growPOSY$'tau_year[1,7]')))
POSY_ydummy_eminus_y8 <- as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*0 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_growPOSY$'beta[5]')*xdummy*0
                                            + mean(post_growPOSY$'tau_year[1,8]')))
POSY_ydummy_eminus_y9 <- as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*0 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_growPOSY$'beta[5]')*xdummy*0
                                            + mean(post_growPOSY$'tau_year[1,9]')))
POSY_ydummy_eminus_y10 <- as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*0 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                             + mean(post_growPOSY$'beta[5]')*xdummy*0
                                             + mean(post_growPOSY$'tau_year[1,10]')))
POSY_ydummy_eminus_y11 <- as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*0 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                             + mean(post_growPOSY$'beta[5]')*xdummy*0
                                             + mean(post_growPOSY$'tau_year[1,11]')))
POSYfitsminus0 <- as.data.frame(cbind(xdummy, POSY_ydummy_eminus))
POSYfitsminus1 <- as.data.frame(cbind(xdummy,POSY_ydummy_eminus_y1,POSY_ydummy_eminus_y2, POSY_ydummy_eminus_y3, POSY_ydummy_eminus_y4, POSY_ydummy_eminus_y5, POSY_ydummy_eminus_y6, POSY_ydummy_eminus_y7, POSY_ydummy_eminus_y8, POSY_ydummy_eminus_y9, POSY_ydummy_eminus_y10, POSY_ydummy_eminus_y11))
POSYfitsminus <- melt(POSYfitsminus1, id = "xdummy",
                      measure.vars = c("POSY_ydummy_eminus_y1","POSY_ydummy_eminus_y2", 
                                       "POSY_ydummy_eminus_y3", "POSY_ydummy_eminus_y4", 
                                       "POSY_ydummy_eminus_y5", "POSY_ydummy_eminus_y6", 
                                       "POSY_ydummy_eminus_y7", "POSY_ydummy_eminus_y8", 
                                       "POSY_ydummy_eminus_y9", "POSY_ydummy_eminus_y10",
                                       "POSY_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "POSY_ydummy_eminus_y1" = '1',"POSY_ydummy_eminus_y2" = "2", 
                       "POSY_ydummy_eminus_y3" = "3", "POSY_ydummy_eminus_y4" = "4", 
                       "POSY_ydummy_eminus_y5" = "5", "POSY_ydummy_eminus_y6" = "6", 
                       "POSY_ydummy_eminus_y7" = "7", "POSY_ydummy_eminus_y8" = "8", 
                       "POSY_ydummy_eminus_y9" = "9", "POSY_ydummy_eminus_y10" = "10",
                       "POSY_ydummy_eminus_y11" = "11") )



# E+
POSY_ydummy_eplus <- as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*1 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                        + mean(post_growPOSY$'beta[5]')*xdummy*1))
POSY_ydummy_eplus_y1 <-  as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*1 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_growPOSY$'beta[5]')*xdummy*1
                                            + mean(post_growPOSY$'tau_year[2,1]')))
POSY_ydummy_eplus_y2 <-  as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*1 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_growPOSY$'beta[5]')*xdummy*1
                                            + mean(post_growPOSY$'tau_year[2,2')))
POSY_ydummy_eplus_y3 <-  as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*1 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_growPOSY$'beta[5]')*xdummy*1
                                            + mean(post_growPOSY$'tau_year[2,3]')))
POSY_ydummy_eplus_y4 <-  as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*1 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_growPOSY$'beta[5]')*xdummy*1
                                            + mean(post_growPOSY$'tau_year[2,4]')))
POSY_ydummy_eplus_y5 <-  as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*1 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_growPOSY$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_growPOSY$'tau_year[2,5]')))
POSY_ydummy_eplus_y6 <-  as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*1 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_growPOSY$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_growPOSY$'tau_year[2,6]')))
POSY_ydummy_eplus_y7 <-  as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*1 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_growPOSY$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_growPOSY$'tau_year[2,7]')))
POSY_ydummy_eplus_y8 <-  as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*1 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_growPOSY$'beta[5]')*xdummy*1
                                            + mean(post_growPOSY$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
POSY_ydummy_eplus_y9 <-  as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*1 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_growPOSY$'beta[5]')*xdummy*1
                                            + mean(post_growPOSY$'tau_year[2,9]')))
POSY_ydummy_eplus_y10 <-  as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*1 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                             + mean(post_growPOSY$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_growPOSY$'tau_year[2,10]')))
POSY_ydummy_eplus_y11 <-  as.vector(exp(mean(post_growPOSY$'beta[1]') + mean(post_growPOSY$'beta[2]')*xdummy + mean(post_growPOSY$'beta[3]')*1 + mean(post_growPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                             + mean(post_growPOSY$'beta[5]')*xdummy*1
                                             + mean(post_growPOSY$'tau_year[2,11]')))
POSY_fitsplus0 <- as.data.frame(cbind(xdummy, POSY_ydummy_eplus))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
POSYfitsplus1 <- as.data.frame(cbind(xdummy,POSY_ydummy_eplus_y1,POSY_ydummy_eplus_y2, POSY_ydummy_eplus_y3, POSY_ydummy_eplus_y4, POSY_ydummy_eplus_y5, POSY_ydummy_eplus_y6, POSY_ydummy_eplus_y7, POSY_ydummy_eplus_y8, POSY_ydummy_eplus_y9, POSY_ydummy_eplus_y10, POSY_ydummy_eplus_y11))
POSYfitsplus <- melt(POSYfitsplus1, id = "xdummy",
                     measure.vars = c("POSY_ydummy_eplus_y1","POSY_ydummy_eplus_y2", 
                                      "POSY_ydummy_eplus_y3", "POSY_ydummy_eplus_y4", 
                                      "POSY_ydummy_eplus_y5", "POSY_ydummy_eplus_y6", 
                                      "POSY_ydummy_eplus_y7", "POSY_ydummy_eplus_y8", 
                                      "POSY_ydummy_eplus_y9", "POSY_ydummy_eplus_y10",
                                      "POSY_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "POSY_ydummy_eplus_y1" = '1',"POSY_ydummy_eplus_y2" = "2", 
                       "POSY_ydummy_eplus_y3" = "3", "POSY_ydummy_eplus_y4" = "4", 
                       "POSY_ydummy_eplus_y5" = "5", "POSY_ydummy_eplus_y6" = "6", 
                       "POSY_ydummy_eplus_y7" = "7", "POSY_ydummy_eplus_y8" = "8", 
                       "POSY_ydummy_eplus_y9" = "9", "POSY_ydummy_eplus_y10" = "10",
                       "POSY_ydummy_eplus_y11" = "11") )

ggplot(data = POSYfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = POSYgrow_bin1, aes(x = mean_size, y = mean_grow, color = Year), lwd = 2) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("POSY E+ Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth")  
scale_color_manual(values=yearcolors)

ggplot(data = POSYfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("POSY E+ Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# effect on mean


POSYfits <- as.data.frame(cbind(xdummy, POSY_ydummy_eminus, POSY_ydummy_eplus))
POSYfits <- melt(POSYfits, id = "xdummy", measure.vars = c("POSY_ydummy_eplus", "POSY_ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "POSY_ydummy_eminus" = "E-", "POSY_ydummy_eplus" = "E+"))

# without point
ggplot(data = POSYfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("POSY Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = POSYfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = POSYgrow_binmean, aes(x = mean_size, y = mean_grow, color = endo_01), lwd = 3) +
  ggtitle("POSY Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth") +   
  scale_color_manual(values=colors2)


# split up plus and minus with datapoints and add mean line
POSY_minus <- ggplot(data = POSYfitsminus1)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = POSYfitsminus0, aes(x = xdummy, y = POSY_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = POSYgrow_bin0mean, aes(x = mean_size, y = mean_grow), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "") + ylim(0,80)


POSY_plus <- ggplot(data = POSYfitsplus1)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = POSY_fitsplus0, aes(x = xdummy, y = POSY_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = POSYgrow_bin1mean, aes(x = mean_size, y = mean_grow), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")+ ylim(0,80)

POSYgrow <- grid.arrange(POSY_plus, POSY_minus, ncol= 2)
titlePOSYgrow <- annotate_figure(POSYgrow, top = "POSY", left = "size_t1")

titlePOSYgrow


# E+ by year
POSYEplusbyyear <- ggplot(data = POSYfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = POSYgrow_bin1mean, aes(x = mean_size, y = mean_grow, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("size_t1") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = POSYfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("size_t1")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
POSYEminusbyyear <- ggplot(data = POSYfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = POSYgrow_bin0mean, aes(x = mean_size, y = mean_grow, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("size_t1")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = POSYfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("POSY E- Growth Probability") + xlab("log(size_t)") + ylab("size_t1")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)



# effect on mean

POSYfits <- as.data.frame(cbind(xdummy, POSY_ydummy_eplus,POSY_ydummy_eminus))


# without point
ggplot(data = POSYfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("POSY Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = POSYfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = POSYsurv_binmean, aes(x = mean_size, y = mean_grow, color = endo_01), lwd = 3) +
  ggtitle("POSY Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +   
  scale_color_manual(values=colors2)

POSYmean <- ggplot(data = POSYfits) +
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus), color = "#ff7f00")+ 
  geom_point(data = POSYgrow_binmean, aes(x = mean_size, y = mean_grow, color = endo_01, lwd = samplesize))+
  ggtitle("Mean Growth") + xlab("log(size_t)") + ylab("size_t1") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")



POSYgrow4panel <- grid.arrange(POSYmean, POSY_g,POSYEplusbyyear,POSYEminusbyyear, ncol= 4)
titlePOSYgrow4panel <- annotate_figure(POSYgrow4panel, top = "Poa sylvestris")

titlePOSYgrow4panel



######## LOAR model fits
# actual data points for size_t1
LOARgrow_bin0 <- LOAR_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>%  
  mutate(size_bin = cut_interval(logsize_t, n = 15))
group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())

LOARgrow_bin1 <- LOAR_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(logsize_t>=0) %>% 
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())


LOARgrow_bin0mean <- LOAR_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 8)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())



LOARgrow_bin1mean <- LOAR_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 8)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())



LOARgrow_bin <- LOAR_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(logsize_t>=0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())


LOARgrow_binmean <- LOAR_data %>% 
  select(size_t1, logsize_t, endo_01) %>%
  mutate(endo_01 = as.character(endo_01)) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(LOAR_data$logsize_t)
xdummy <- seq(min(LOAR_data$logsize_t), max(LOAR_data$logsize_t), length.out = nvalues)

# E- 
LOAR_ydummy_eminus <- as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*0 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                         + mean(post_growLOAR$'beta[5]')*xdummy*0))
LOAR_ydummy_eminus_y1 <- as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*0 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_growLOAR$'beta[5]')*xdummy*0
                                            + mean(post_growLOAR$'tau_year[1,1]')))
LOAR_ydummy_eminus_y2 <- as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*0 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_growLOAR$'beta[5]')*xdummy*0
                                            + mean(post_growLOAR$'tau_year[1,2]')))
LOAR_ydummy_eminus_y3 <- as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + + mean(post_growLOAR$'beta[3]')*0  + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_growLOAR$'beta[5]')*xdummy*0
                                            + mean(post_growLOAR$'tau_year[1,3]')))
LOAR_ydummy_eminus_y4 <- as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*0 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_growLOAR$'beta[5]')*xdummy*0
                                            + mean(post_growLOAR$'tau_year[1,4]'))) 
LOAR_ydummy_eminus_y5 <- as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*0 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_growLOAR$'beta[5]')*xdummy*0
                                            + mean(post_growLOAR$'tau_year[1,5]')))
LOAR_ydummy_eminus_y6 <- as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*0 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_growLOAR$'beta[5]')*xdummy*0
                                            + mean(post_growLOAR$'tau_year[1,6]')))
LOAR_ydummy_eminus_y7 <- as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*0 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_growLOAR$'beta[5]')*xdummy*0
                                            + mean(post_growLOAR$'tau_year[1,7]')))
LOAR_ydummy_eminus_y8 <- as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*0 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_growLOAR$'beta[5]')*xdummy*0
                                            + mean(post_growLOAR$'tau_year[1,8]')))
LOAR_ydummy_eminus_y9 <- as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*0 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_growLOAR$'beta[5]')*xdummy*0
                                            + mean(post_growLOAR$'tau_year[1,9]')))
LOAR_ydummy_eminus_y10 <- as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*0 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                             + mean(post_growLOAR$'beta[5]')*xdummy*0
                                             + mean(post_growLOAR$'tau_year[1,10]')))
LOAR_ydummy_eminus_y11 <- as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*0 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                             + mean(post_growLOAR$'beta[5]')*xdummy*0
                                             + mean(post_growLOAR$'tau_year[1,11]')))
LOARfitsminus0 <- as.data.frame(cbind(xdummy, LOAR_ydummy_eminus))
LOARfitsminus1 <- as.data.frame(cbind(xdummy,LOAR_ydummy_eminus_y1,LOAR_ydummy_eminus_y2, LOAR_ydummy_eminus_y3, LOAR_ydummy_eminus_y4, LOAR_ydummy_eminus_y5, LOAR_ydummy_eminus_y6, LOAR_ydummy_eminus_y7, LOAR_ydummy_eminus_y8, LOAR_ydummy_eminus_y9, LOAR_ydummy_eminus_y10, LOAR_ydummy_eminus_y11))
LOARfitsminus <- melt(LOARfitsminus1, id = "xdummy",
                      measure.vars = c("LOAR_ydummy_eminus_y1","LOAR_ydummy_eminus_y2", 
                                       "LOAR_ydummy_eminus_y3", "LOAR_ydummy_eminus_y4", 
                                       "LOAR_ydummy_eminus_y5", "LOAR_ydummy_eminus_y6", 
                                       "LOAR_ydummy_eminus_y7", "LOAR_ydummy_eminus_y8", 
                                       "LOAR_ydummy_eminus_y9", "LOAR_ydummy_eminus_y10",
                                       "LOAR_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "LOAR_ydummy_eminus_y1" = '1',"LOAR_ydummy_eminus_y2" = "2", 
                       "LOAR_ydummy_eminus_y3" = "3", "LOAR_ydummy_eminus_y4" = "4", 
                       "LOAR_ydummy_eminus_y5" = "5", "LOAR_ydummy_eminus_y6" = "6", 
                       "LOAR_ydummy_eminus_y7" = "7", "LOAR_ydummy_eminus_y8" = "8", 
                       "LOAR_ydummy_eminus_y9" = "9", "LOAR_ydummy_eminus_y10" = "10",
                       "LOAR_ydummy_eminus_y11" = "11") )



# E+
LOAR_ydummy_eplus <- as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*1 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                        + mean(post_growLOAR$'beta[5]')*xdummy*1))
LOAR_ydummy_eplus_y1 <-  as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*1 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_growLOAR$'beta[5]')*xdummy*1
                                            + mean(post_growLOAR$'tau_year[2,1]')))
LOAR_ydummy_eplus_y2 <-  as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*1 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_growLOAR$'beta[5]')*xdummy*1
                                            + mean(post_growLOAR$'tau_year[2,2')))
LOAR_ydummy_eplus_y3 <-  as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*1 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_growLOAR$'beta[5]')*xdummy*1
                                            + mean(post_growLOAR$'tau_year[2,3]')))
LOAR_ydummy_eplus_y4 <-  as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*1 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_growLOAR$'beta[5]')*xdummy*1
                                            + mean(post_growLOAR$'tau_year[2,4]')))
LOAR_ydummy_eplus_y5 <-  as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*1 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_growLOAR$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_growLOAR$'tau_year[2,5]')))
LOAR_ydummy_eplus_y6 <-  as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*1 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_growLOAR$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_growLOAR$'tau_year[2,6]')))
LOAR_ydummy_eplus_y7 <-  as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*1 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_growLOAR$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_growLOAR$'tau_year[2,7]')))
LOAR_ydummy_eplus_y8 <-  as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*1 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_growLOAR$'beta[5]')*xdummy*1
                                            + mean(post_growLOAR$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
LOAR_ydummy_eplus_y9 <-  as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*1 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_growLOAR$'beta[5]')*xdummy*1
                                            + mean(post_growLOAR$'tau_year[2,9]')))
LOAR_ydummy_eplus_y10 <-  as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*1 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                             + mean(post_growLOAR$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_growLOAR$'tau_year[2,10]')))
LOAR_ydummy_eplus_y11 <-  as.vector(exp(mean(post_growLOAR$'beta[1]') + mean(post_growLOAR$'beta[2]')*xdummy + mean(post_growLOAR$'beta[3]')*1 + mean(post_growLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                             + mean(post_growLOAR$'beta[5]')*xdummy*1
                                             + mean(post_growLOAR$'tau_year[2,11]')))
LOAR_fitsplus0 <- as.data.frame(cbind(xdummy, LOAR_ydummy_eplus))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
LOARfitsplus1 <- as.data.frame(cbind(xdummy,LOAR_ydummy_eplus_y1,LOAR_ydummy_eplus_y2, LOAR_ydummy_eplus_y3, LOAR_ydummy_eplus_y4, LOAR_ydummy_eplus_y5, LOAR_ydummy_eplus_y6, LOAR_ydummy_eplus_y7, LOAR_ydummy_eplus_y8, LOAR_ydummy_eplus_y9, LOAR_ydummy_eplus_y10, LOAR_ydummy_eplus_y11))
LOARfitsplus <- melt(LOARfitsplus1, id = "xdummy",
                     measure.vars = c("LOAR_ydummy_eplus_y1","LOAR_ydummy_eplus_y2", 
                                      "LOAR_ydummy_eplus_y3", "LOAR_ydummy_eplus_y4", 
                                      "LOAR_ydummy_eplus_y5", "LOAR_ydummy_eplus_y6", 
                                      "LOAR_ydummy_eplus_y7", "LOAR_ydummy_eplus_y8", 
                                      "LOAR_ydummy_eplus_y9", "LOAR_ydummy_eplus_y10",
                                      "LOAR_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "LOAR_ydummy_eplus_y1" = '1',"LOAR_ydummy_eplus_y2" = "2", 
                       "LOAR_ydummy_eplus_y3" = "3", "LOAR_ydummy_eplus_y4" = "4", 
                       "LOAR_ydummy_eplus_y5" = "5", "LOAR_ydummy_eplus_y6" = "6", 
                       "LOAR_ydummy_eplus_y7" = "7", "LOAR_ydummy_eplus_y8" = "8", 
                       "LOAR_ydummy_eplus_y9" = "9", "LOAR_ydummy_eplus_y10" = "10",
                       "LOAR_ydummy_eplus_y11" = "11") )

ggplot(data = LOARfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = LOARgrow_bin1, aes(x = mean_size, y = mean_grow, color = Year), lwd = 2) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("LOAR E+ Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth")  
scale_color_manual(values=yearcolors)

ggplot(data = LOARfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("LOAR E+ Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# effect on mean


LOARfits <- as.data.frame(cbind(xdummy, LOAR_ydummy_eminus, LOAR_ydummy_eplus))
LOARfits <- melt(LOARfits, id = "xdummy", measure.vars = c("LOAR_ydummy_eplus", "LOAR_ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "LOAR_ydummy_eminus" = "E-", "LOAR_ydummy_eplus" = "E+"))

# without point
ggplot(data = LOARfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("LOAR Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = LOARfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = LOARgrow_binmean, aes(x = mean_size, y = mean_grow, color = endo_01), lwd = 3) +
  ggtitle("LOAR Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth") +   
  scale_color_manual(values=colors2)


# split up plus and minus with datapoints and add mean line
LOAR_minus <- ggplot(data = LOARfitsminus1)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = LOARfitsminus0, aes(x = xdummy, y = LOAR_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = LOARgrow_bin0mean, aes(x = mean_size, y = mean_grow), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


LOAR_plus <- ggplot(data = LOARfitsplus1)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = LOAR_fitsplus0, aes(x = xdummy, y = LOAR_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = LOARgrow_bin1mean, aes(x = mean_size, y = mean_grow), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

LOARgrow <- grid.arrange(LOAR_plus, LOAR_minus, ncol= 2)
titleLOARgrow <- annotate_figure(LOARgrow, top = "LOAR", left = "size_t1")

titleLOARgrow


# E+ by year
LOAREplusbyyear <- ggplot(data = LOARfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = LOARgrow_bin1mean, aes(x = mean_size, y = mean_grow, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("size_t1") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = LOARfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("size_t1")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
LOAREminusbyyear <- ggplot(data = LOARfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = LOARgrow_bin0mean, aes(x = mean_size, y = mean_grow, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("size_t1")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = LOARfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("LOAR E- Growth Probability") + xlab("log(size_t)") + ylab("size_t1")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)



# effect on mean

LOARfits <- as.data.frame(cbind(xdummy, LOAR_ydummy_eplus, LOAR_ydummy_eminus))


# without point
ggplot(data = LOARfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("LOAR Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = LOARfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = LOARsurv_binmean, aes(x = mean_size, y = mean_grow, color = endo), lwd = 3) +
  ggtitle("LOAR Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +   
  scale_color_manual(values=colors2)

LOARmean <- ggplot(data = LOARfits) +
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus), color = "#ff7f00")+ 
  geom_point(data = LOARgrow_binmean, aes(x = mean_size, y = mean_grow, color = endo_01, lwd = samplesize))+
  ggtitle("Mean Growth") + xlab("log(size_t)") + ylab("size_t1") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")




LOARgrow4panel <- grid.arrange(LOARmean, LOAR_g,LOAREplusbyyear,LOAREminusbyyear, ncol= 4)
titleLOARgrow4panel <- annotate_figure(LOARgrow4panel, top = "Lolium arundinacea")

titleLOARgrow4panel





######## FESU model fits
# actual data points for size_t1
FESUgrow_bin0 <- FESU_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())

FESUgrow_bin1 <- FESU_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())

FESUgrow_bin0mean <- FESU_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 8)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())



FESUgrow_bin1mean <- FESU_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 8)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())



FESUgrow_bin <- FESU_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())


FESUgrow_binmean <- FESU_data %>% 
  select(size_t1, logsize_t, endo_01) %>%
  mutate(endo_01 = as.character(endo_01)) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(FESU_data$logsize_t)
xdummy <- seq(min(FESU_data$logsize_t), max(FESU_data$logsize_t), length.out = nvalues)

# E- 
FESU_ydummy_eminus <- as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*0 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01)
                                         + mean(post_growFESU$'beta[5]')*xdummy*0))
FESU_ydummy_eminus_y1 <- as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*0 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_growFESU$'beta[5]')*xdummy*0
                                            + mean(post_growFESU$'tau_year[1,1]')))
FESU_ydummy_eminus_y2 <- as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*0 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_growFESU$'beta[5]')*xdummy*0
                                            + mean(post_growFESU$'tau_year[1,2]')))
FESU_ydummy_eminus_y3 <- as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + + mean(post_growFESU$'beta[3]')*0  + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_growFESU$'beta[5]')*xdummy*0
                                            + mean(post_growFESU$'tau_year[1,3]')))
FESU_ydummy_eminus_y4 <- as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*0 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_growFESU$'beta[5]')*xdummy*0
                                            + mean(post_growFESU$'tau_year[1,4]'))) 
FESU_ydummy_eminus_y5 <- as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*0 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_growFESU$'beta[5]')*xdummy*0
                                            + mean(post_growFESU$'tau_year[1,5]')))
FESU_ydummy_eminus_y6 <- as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*0 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_growFESU$'beta[5]')*xdummy*0
                                            + mean(post_growFESU$'tau_year[1,6]')))
FESU_ydummy_eminus_y7 <- as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*0 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_growFESU$'beta[5]')*xdummy*0
                                            + mean(post_growFESU$'tau_year[1,7]')))
FESU_ydummy_eminus_y8 <- as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*0 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_growFESU$'beta[5]')*xdummy*0
                                            + mean(post_growFESU$'tau_year[1,8]')))
FESU_ydummy_eminus_y9 <- as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*0 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_growFESU$'beta[5]')*xdummy*0
                                            + mean(post_growFESU$'tau_year[1,9]')))
FESU_ydummy_eminus_y10 <- as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*0 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01)
                                             + mean(post_growFESU$'beta[5]')*xdummy*0
                                             + mean(post_growFESU$'tau_year[1,10]')))
FESU_ydummy_eminus_y11 <- as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*0 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01)
                                             + mean(post_growFESU$'beta[5]')*xdummy*0
                                             + mean(post_growFESU$'tau_year[1,11]')))
FESUfitsminus0 <- as.data.frame(cbind(xdummy, FESU_ydummy_eminus))
FESUfitsminus1 <- as.data.frame(cbind(xdummy,FESU_ydummy_eminus_y1,FESU_ydummy_eminus_y2, FESU_ydummy_eminus_y3, FESU_ydummy_eminus_y4, FESU_ydummy_eminus_y5, FESU_ydummy_eminus_y6, FESU_ydummy_eminus_y7, FESU_ydummy_eminus_y8, FESU_ydummy_eminus_y9, FESU_ydummy_eminus_y10, FESU_ydummy_eminus_y11))
FESUfitsminus <- melt(FESUfitsminus1, id = "xdummy",
                      measure.vars = c("FESU_ydummy_eminus_y1","FESU_ydummy_eminus_y2", 
                                       "FESU_ydummy_eminus_y3", "FESU_ydummy_eminus_y4", 
                                       "FESU_ydummy_eminus_y5", "FESU_ydummy_eminus_y6", 
                                       "FESU_ydummy_eminus_y7", "FESU_ydummy_eminus_y8", 
                                       "FESU_ydummy_eminus_y9", "FESU_ydummy_eminus_y10",
                                       "FESU_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "FESU_ydummy_eminus_y1" = '1',"FESU_ydummy_eminus_y2" = "2", 
                       "FESU_ydummy_eminus_y3" = "3", "FESU_ydummy_eminus_y4" = "4", 
                       "FESU_ydummy_eminus_y5" = "5", "FESU_ydummy_eminus_y6" = "6", 
                       "FESU_ydummy_eminus_y7" = "7", "FESU_ydummy_eminus_y8" = "8", 
                       "FESU_ydummy_eminus_y9" = "9", "FESU_ydummy_eminus_y10" = "10",
                       "FESU_ydummy_eminus_y11" = "11") )



# E+
FESU_ydummy_eplus <- as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*1 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01)
                                        + mean(post_growFESU$'beta[5]')*xdummy*1))
FESU_ydummy_eplus_y1 <-  as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*1 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_growFESU$'beta[5]')*xdummy*1
                                            + mean(post_growFESU$'tau_year[2,1]')))
FESU_ydummy_eplus_y2 <-  as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*1 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_growFESU$'beta[5]')*xdummy*1
                                            + mean(post_growFESU$'tau_year[2,2')))
FESU_ydummy_eplus_y3 <-  as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*1 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_growFESU$'beta[5]')*xdummy*1
                                            + mean(post_growFESU$'tau_year[2,3]')))
FESU_ydummy_eplus_y4 <-  as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*1 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_growFESU$'beta[5]')*xdummy*1
                                            + mean(post_growFESU$'tau_year[2,4]')))
FESU_ydummy_eplus_y5 <-  as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*1 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_growFESU$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_growFESU$'tau_year[2,5]')))
FESU_ydummy_eplus_y6 <-  as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*1 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_growFESU$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_growFESU$'tau_year[2,6]')))
FESU_ydummy_eplus_y7 <-  as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*1 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_growFESU$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_growFESU$'tau_year[2,7]')))
FESU_ydummy_eplus_y8 <-  as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*1 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_growFESU$'beta[5]')*xdummy*1
                                            + mean(post_growFESU$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
FESU_ydummy_eplus_y9 <-  as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*1 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_growFESU$'beta[5]')*xdummy*1
                                            + mean(post_growFESU$'tau_year[2,9]')))
FESU_ydummy_eplus_y10 <-  as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*1 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                             + mean(post_growFESU$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_growFESU$'tau_year[2,10]')))
FESU_ydummy_eplus_y11 <-  as.vector(exp(mean(post_growFESU$'beta[1]') + mean(post_growFESU$'beta[2]')*xdummy + mean(post_growFESU$'beta[3]')*1 + mean(post_growFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                             + mean(post_growFESU$'beta[5]')*xdummy*1
                                             + mean(post_growFESU$'tau_year[2,11]')))
FESU_fitsplus0 <- as.data.frame(cbind(xdummy, FESU_ydummy_eplus))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
FESUfitsplus1 <- as.data.frame(cbind(xdummy,FESU_ydummy_eplus_y1,FESU_ydummy_eplus_y2, FESU_ydummy_eplus_y3, FESU_ydummy_eplus_y4, FESU_ydummy_eplus_y5, FESU_ydummy_eplus_y6, FESU_ydummy_eplus_y7, FESU_ydummy_eplus_y8, FESU_ydummy_eplus_y9, FESU_ydummy_eplus_y10, FESU_ydummy_eplus_y11))
FESUfitsplus <- melt(FESUfitsplus1, id = "xdummy",
                     measure.vars = c("FESU_ydummy_eplus_y1","FESU_ydummy_eplus_y2", 
                                      "FESU_ydummy_eplus_y3", "FESU_ydummy_eplus_y4", 
                                      "FESU_ydummy_eplus_y5", "FESU_ydummy_eplus_y6", 
                                      "FESU_ydummy_eplus_y7", "FESU_ydummy_eplus_y8", 
                                      "FESU_ydummy_eplus_y9", "FESU_ydummy_eplus_y10",
                                      "FESU_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "FESU_ydummy_eplus_y1" = '1',"FESU_ydummy_eplus_y2" = "2", 
                       "FESU_ydummy_eplus_y3" = "3", "FESU_ydummy_eplus_y4" = "4", 
                       "FESU_ydummy_eplus_y5" = "5", "FESU_ydummy_eplus_y6" = "6", 
                       "FESU_ydummy_eplus_y7" = "7", "FESU_ydummy_eplus_y8" = "8", 
                       "FESU_ydummy_eplus_y9" = "9", "FESU_ydummy_eplus_y10" = "10",
                       "FESU_ydummy_eplus_y11" = "11") )

ggplot(data = FESUfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = FESUgrow_bin1, aes(x = mean_size, y = mean_grow, color = Year), lwd = 2) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("FESU E+ Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth")  
scale_color_manual(values=yearcolors)

ggplot(data = FESUfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("FESU E+ Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# effect on mean


FESUfits <- as.data.frame(cbind(xdummy, FESU_ydummy_eminus, FESU_ydummy_eplus))
FESUfits <- melt(FESUfits, id = "xdummy", measure.vars = c("FESU_ydummy_eplus", "FESU_ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "FESU_ydummy_eminus" = "E-", "FESU_ydummy_eplus" = "E+"))

# without point
ggplot(data = FESUfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("FESU Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = FESUfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = FESUgrow_binmean, aes(x = mean_size, y = mean_grow, color = endo_01), lwd = 3) +
  ggtitle("FESU Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth") +   
  scale_color_manual(values=colors2)


# split up plus and minus with datapoints and add mean line
FESU_minus <- ggplot(data = FESUfitsminus1)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = FESUfitsminus0, aes(x = xdummy, y = FESU_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = FESUgrow_bin0mean, aes(x = mean_size, y = mean_grow), color ="#ff7f00") +
  theme_classic()+ ylim(0,40)+
  labs(title = "E-", x = "log(size_t)", y = "")


FESU_plus <- ggplot(data = FESUfitsplus1)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = FESU_fitsplus0, aes(x = xdummy, y = FESU_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = FESUgrow_bin1mean, aes(x = mean_size, y = mean_grow), color ="#6a3d9a") +
  theme_classic()+ ylim(0,40)+
  labs(title = "E+", x = "log(size_t)", y = "")

FESUgrow <- grid.arrange(FESU_plus, FESU_minus, ncol= 2)
titleFESUgrow <- annotate_figure(FESUgrow, top = "FESU", left = "size_t1")

titleFESUgrow



# E+ by year
FESUEplusbyyear <- ggplot(data = FESUfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = FESUgrow_bin1mean, aes(x = mean_size, y = mean_grow, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("size_t1") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = FESUfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("size_t1")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
FESUEminusbyyear <- ggplot(data = FESUfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = FESUgrow_bin0mean, aes(x = mean_size, y = mean_grow, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("size_t1")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = FESUfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("FESU E- Growth Probability") + xlab("log(size_t)") + ylab("size_t1")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)



# effect on mean

FESUfits <- as.data.frame(cbind(xdummy, FESU_ydummy_eplus,FESU_ydummy_eminus))


# without point
ggplot(data = FESUfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("FESU Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = FESUfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = FESUsurv_binmean, aes(x = mean_size, y = mean_grow, color = endo), lwd = 3) +
  ggtitle("FESU Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +   
  scale_color_manual(values=colors2)

FESUmean <- ggplot(data = FESUfits) +
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus), color = "#ff7f00")+ 
  geom_point(data = FESUgrow_binmean, aes(x = mean_size, y = mean_grow, color = endo_01, lwd = samplesize))+
  ggtitle("Mean Growth") + xlab("log(size_t)") + ylab("size_t1") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")




FESUgrow4panel <- grid.arrange(FESUmean, FESU_g,FESUEplusbyyear,FESUEminusbyyear, ncol= 4)
titleFESUgrow4panel <- annotate_figure(FESUgrow4panel, top = "Festuca subverticillata")

titleFESUgrow4panel




######## ELVI model fits
# actual data points for size_t1
ELVIgrow_bin0 <- ELVI_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())

ELVIgrow_bin1 <- ELVI_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())

ELVIgrow_bin0mean <- ELVI_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 8)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())



ELVIgrow_bin1mean <- ELVI_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 8)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())



ELVIgrow_bin <- ELVI_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())


ELVIgrow_binmean <- ELVI_data %>% 
  select(size_t1, logsize_t, endo_01) %>%
  mutate(endo_01 = as.character(endo_01)) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(ELVI_data$logsize_t)
xdummy <- seq(min(ELVI_data$logsize_t), max(ELVI_data$logsize_t), length.out = nvalues)

# E- 
ELVI_ydummy_eminus <- as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*0 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                         + mean(post_growELVI$'beta[5]')*xdummy*0))
ELVI_ydummy_eminus_y1 <- as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*0 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_growELVI$'beta[5]')*xdummy*0
                                            + mean(post_growELVI$'tau_year[1,1]')))
ELVI_ydummy_eminus_y2 <- as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*0 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_growELVI$'beta[5]')*xdummy*0
                                            + mean(post_growELVI$'tau_year[1,2]')))
ELVI_ydummy_eminus_y3 <- as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + + mean(post_growELVI$'beta[3]')*0  + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_growELVI$'beta[5]')*xdummy*0
                                            + mean(post_growELVI$'tau_year[1,3]')))
ELVI_ydummy_eminus_y4 <- as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*0 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_growELVI$'beta[5]')*xdummy*0
                                            + mean(post_growELVI$'tau_year[1,4]'))) 
ELVI_ydummy_eminus_y5 <- as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*0 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_growELVI$'beta[5]')*xdummy*0
                                            + mean(post_growELVI$'tau_year[1,5]')))
ELVI_ydummy_eminus_y6 <- as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*0 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_growELVI$'beta[5]')*xdummy*0
                                            + mean(post_growELVI$'tau_year[1,6]')))
ELVI_ydummy_eminus_y7 <- as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*0 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_growELVI$'beta[5]')*xdummy*0
                                            + mean(post_growELVI$'tau_year[1,7]')))
ELVI_ydummy_eminus_y8 <- as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*0 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_growELVI$'beta[5]')*xdummy*0
                                            + mean(post_growELVI$'tau_year[1,8]')))
ELVI_ydummy_eminus_y9 <- as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*0 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_growELVI$'beta[5]')*xdummy*0
                                            + mean(post_growELVI$'tau_year[1,9]')))
ELVI_ydummy_eminus_y10 <- as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*0 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                             + mean(post_growELVI$'beta[5]')*xdummy*0
                                             + mean(post_growELVI$'tau_year[1,10]')))
ELVI_ydummy_eminus_y11 <- as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*0 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                             + mean(post_growELVI$'beta[5]')*xdummy*0
                                             + mean(post_growELVI$'tau_year[1,11]')))
ELVIfitsminus0 <- as.data.frame(cbind(xdummy, ELVI_ydummy_eminus))
ELVIfitsminus1 <- as.data.frame(cbind(xdummy,ELVI_ydummy_eminus_y1,ELVI_ydummy_eminus_y2, ELVI_ydummy_eminus_y3, ELVI_ydummy_eminus_y4, ELVI_ydummy_eminus_y5, ELVI_ydummy_eminus_y6, ELVI_ydummy_eminus_y7, ELVI_ydummy_eminus_y8, ELVI_ydummy_eminus_y9, ELVI_ydummy_eminus_y10, ELVI_ydummy_eminus_y11))
ELVIfitsminus <- melt(ELVIfitsminus1, id = "xdummy",
                      measure.vars = c("ELVI_ydummy_eminus_y1","ELVI_ydummy_eminus_y2", 
                                       "ELVI_ydummy_eminus_y3", "ELVI_ydummy_eminus_y4", 
                                       "ELVI_ydummy_eminus_y5", "ELVI_ydummy_eminus_y6", 
                                       "ELVI_ydummy_eminus_y7", "ELVI_ydummy_eminus_y8", 
                                       "ELVI_ydummy_eminus_y9", "ELVI_ydummy_eminus_y10",
                                       "ELVI_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "ELVI_ydummy_eminus_y1" = '1',"ELVI_ydummy_eminus_y2" = "2", 
                       "ELVI_ydummy_eminus_y3" = "3", "ELVI_ydummy_eminus_y4" = "4", 
                       "ELVI_ydummy_eminus_y5" = "5", "ELVI_ydummy_eminus_y6" = "6", 
                       "ELVI_ydummy_eminus_y7" = "7", "ELVI_ydummy_eminus_y8" = "8", 
                       "ELVI_ydummy_eminus_y9" = "9", "ELVI_ydummy_eminus_y10" = "10",
                       "ELVI_ydummy_eminus_y11" = "11") )



# E+
ELVI_ydummy_eplus <- as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*1 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                        + mean(post_growELVI$'beta[5]')*xdummy*1))
ELVI_ydummy_eplus_y1 <-  as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*1 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_growELVI$'beta[5]')*xdummy*1
                                            + mean(post_growELVI$'tau_year[2,1]')))
ELVI_ydummy_eplus_y2 <-  as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*1 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_growELVI$'beta[5]')*xdummy*1
                                            + mean(post_growELVI$'tau_year[2,2')))
ELVI_ydummy_eplus_y3 <-  as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*1 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_growELVI$'beta[5]')*xdummy*1
                                            + mean(post_growELVI$'tau_year[2,3]')))
ELVI_ydummy_eplus_y4 <-  as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*1 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_growELVI$'beta[5]')*xdummy*1
                                            + mean(post_growELVI$'tau_year[2,4]')))
ELVI_ydummy_eplus_y5 <-  as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*1 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_growELVI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_growELVI$'tau_year[2,5]')))
ELVI_ydummy_eplus_y6 <-  as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*1 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_growELVI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_growELVI$'tau_year[2,6]')))
ELVI_ydummy_eplus_y7 <-  as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*1 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_growELVI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_growELVI$'tau_year[2,7]')))
ELVI_ydummy_eplus_y8 <-  as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*1 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_growELVI$'beta[5]')*xdummy*1
                                            + mean(post_growELVI$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
ELVI_ydummy_eplus_y9 <-  as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*1 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_growELVI$'beta[5]')*xdummy*1
                                            + mean(post_growELVI$'tau_year[2,9]')))
ELVI_ydummy_eplus_y10 <-  as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*1 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                             + mean(post_growELVI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_growELVI$'tau_year[2,10]')))
ELVI_ydummy_eplus_y11 <-  as.vector(exp(mean(post_growELVI$'beta[1]') + mean(post_growELVI$'beta[2]')*xdummy + mean(post_growELVI$'beta[3]')*1 + mean(post_growELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                             + mean(post_growELVI$'beta[5]')*xdummy*1
                                             + mean(post_growELVI$'tau_year[2,11]')))
ELVI_fitsplus0 <- as.data.frame(cbind(xdummy, ELVI_ydummy_eplus))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
ELVIfitsplus1 <- as.data.frame(cbind(xdummy,ELVI_ydummy_eplus_y1,ELVI_ydummy_eplus_y2, ELVI_ydummy_eplus_y3, ELVI_ydummy_eplus_y4, ELVI_ydummy_eplus_y5, ELVI_ydummy_eplus_y6, ELVI_ydummy_eplus_y7, ELVI_ydummy_eplus_y8, ELVI_ydummy_eplus_y9, ELVI_ydummy_eplus_y10, ELVI_ydummy_eplus_y11))
ELVIfitsplus <- melt(ELVIfitsplus1, id = "xdummy",
                     measure.vars = c("ELVI_ydummy_eplus_y1","ELVI_ydummy_eplus_y2", 
                                      "ELVI_ydummy_eplus_y3", "ELVI_ydummy_eplus_y4", 
                                      "ELVI_ydummy_eplus_y5", "ELVI_ydummy_eplus_y6", 
                                      "ELVI_ydummy_eplus_y7", "ELVI_ydummy_eplus_y8", 
                                      "ELVI_ydummy_eplus_y9", "ELVI_ydummy_eplus_y10",
                                      "ELVI_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "ELVI_ydummy_eplus_y1" = '1',"ELVI_ydummy_eplus_y2" = "2", 
                       "ELVI_ydummy_eplus_y3" = "3", "ELVI_ydummy_eplus_y4" = "4", 
                       "ELVI_ydummy_eplus_y5" = "5", "ELVI_ydummy_eplus_y6" = "6", 
                       "ELVI_ydummy_eplus_y7" = "7", "ELVI_ydummy_eplus_y8" = "8", 
                       "ELVI_ydummy_eplus_y9" = "9", "ELVI_ydummy_eplus_y10" = "10",
                       "ELVI_ydummy_eplus_y11" = "11") )

ggplot(data = ELVIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = ELVIgrow_bin1, aes(x = mean_size, y = mean_grow, color = Year), lwd = 2) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("ELVI E+ Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth")  
scale_color_manual(values=yearcolors)

ggplot(data = ELVIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("ELVI E+ Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# effect on mean


ELVIfits <- as.data.frame(cbind(xdummy, ELVI_ydummy_eminus, ELVI_ydummy_eplus))
ELVIfits <- melt(ELVIfits, id = "xdummy", measure.vars = c("ELVI_ydummy_eplus", "ELVI_ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ELVI_ydummy_eminus" = "E-", "ELVI_ydummy_eplus" = "E+"))

# without point
ggplot(data = ELVIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("ELVI Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = ELVIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = ELVIgrow_binmean, aes(x = mean_size, y = mean_grow, color = endo_01), lwd = 3) +
  ggtitle("ELVI Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth") +   
  scale_color_manual(values=colors2)


# split up plus and minus with datapoints and add mean line
ELVI_minus <- ggplot(data = ELVIfitsminus1)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = ELVIfitsminus0, aes(x = xdummy, y = ELVI_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = ELVIgrow_bin0mean, aes(x = mean_size, y = mean_grow), color ="#ff7f00") +
  theme_classic()+ ylim(0,30)+
  labs(title = "E-", x = "log(size_t)", y = "")


ELVI_plus <- ggplot(data = ELVIfitsplus1)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = ELVI_fitsplus0, aes(x = xdummy, y = ELVI_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = ELVIgrow_bin1mean, aes(x = mean_size, y = mean_grow), color ="#6a3d9a") +
  theme_classic()+ylim(0,30)+
  labs(title = "E+", x = "log(size_t)", y = "")

ELVIgrow <- grid.arrange(ELVI_plus, ELVI_minus, ncol= 2)
titleELVIgrow <- annotate_figure(ELVIgrow, top = "ELVI", left = "size_t1")

titleELVIgrow




# E+ by year
ELVIEplusbyyear <- ggplot(data = ELVIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = ELVIgrow_bin1mean, aes(x = mean_size, y = mean_grow, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("size_t1") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = ELVIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("size_t1")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
ELVIEminusbyyear <- ggplot(data = ELVIfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = ELVIgrow_bin0mean, aes(x = mean_size, y = mean_grow, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("size_t1")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = ELVIfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("ELVI E- Growth Probability") + xlab("log(size_t)") + ylab("size_t1")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)



# effect on mean

ELVIfits <- as.data.frame(cbind(xdummy, ELVI_ydummy_eplus,ELVI_ydummy_eminus))


# without point
ggplot(data = ELVIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("ELVI Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = ELVIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = ELVIsurv_binmean, aes(x = mean_size, y = mean_grow, color = endo), lwd = 3) +
  ggtitle("ELVI Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +   
  scale_color_manual(values=colors2)

ELVImean <- ggplot(data = ELVIfits) +
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus), color = "#ff7f00")+ 
  geom_point(data = ELVIgrow_binmean, aes(x = mean_size, y = mean_grow, color = endo_01, lwd = samplesize))+
  ggtitle("Mean Growth") + xlab("log(size_t)") + ylab("size_t1") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")




ELVIgrow4panel <- grid.arrange(ELVImean, ELVI_g,ELVIEplusbyyear,ELVIEminusbyyear, ncol= 4)
titleELVIgrow4panel <- annotate_figure(ELVIgrow4panel, top = "Elymus virginicus")

titleELVIgrow4panel



######## ELRI model fits
# actual data points for size_t1
ELRIgrow_bin0 <- ELRI_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())

ELRIgrow_bin1 <- ELRI_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())

ELRIgrow_bin0mean <- ELRI_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 8)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())



ELRIgrow_bin1mean <- ELRI_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 8)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())



ELRIgrow_bin <- ELRI_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())


ELRIgrow_binmean <- ELRI_data %>% 
  select(size_t1, logsize_t, endo_01) %>%
  mutate(endo_01 =as.character(endo_01)) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(ELRI_data$logsize_t)
xdummy <- seq(min(ELRI_data$logsize_t), max(ELRI_data$logsize_t), length.out = nvalues)

# E- 
ELRI_ydummy_eminus <- as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*0 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                         + mean(post_growELRI$'beta[5]')*xdummy*0))
ELRI_ydummy_eminus_y1 <- as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*0 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_growELRI$'beta[5]')*xdummy*0
                                            + mean(post_growELRI$'tau_year[1,1]')))
ELRI_ydummy_eminus_y2 <- as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*0 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_growELRI$'beta[5]')*xdummy*0
                                            + mean(post_growELRI$'tau_year[1,2]')))
ELRI_ydummy_eminus_y3 <- as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + + mean(post_growELRI$'beta[3]')*0  + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_growELRI$'beta[5]')*xdummy*0
                                            + mean(post_growELRI$'tau_year[1,3]')))
ELRI_ydummy_eminus_y4 <- as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*0 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_growELRI$'beta[5]')*xdummy*0
                                            + mean(post_growELRI$'tau_year[1,4]'))) 
ELRI_ydummy_eminus_y5 <- as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*0 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_growELRI$'beta[5]')*xdummy*0
                                            + mean(post_growELRI$'tau_year[1,5]')))
ELRI_ydummy_eminus_y6 <- as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*0 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_growELRI$'beta[5]')*xdummy*0
                                            + mean(post_growELRI$'tau_year[1,6]')))
ELRI_ydummy_eminus_y7 <- as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*0 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_growELRI$'beta[5]')*xdummy*0
                                            + mean(post_growELRI$'tau_year[1,7]')))
ELRI_ydummy_eminus_y8 <- as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*0 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_growELRI$'beta[5]')*xdummy*0
                                            + mean(post_growELRI$'tau_year[1,8]')))
ELRI_ydummy_eminus_y9 <- as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*0 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_growELRI$'beta[5]')*xdummy*0
                                            + mean(post_growELRI$'tau_year[1,9]')))
ELRI_ydummy_eminus_y10 <- as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*0 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                             + mean(post_growELRI$'beta[5]')*xdummy*0
                                             + mean(post_growELRI$'tau_year[1,10]')))
ELRI_ydummy_eminus_y11 <- as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*0 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                             + mean(post_growELRI$'beta[5]')*xdummy*0
                                             + mean(post_growELRI$'tau_year[1,11]')))
ELRIfitsminus0 <- as.data.frame(cbind(xdummy, ELRI_ydummy_eminus))
ELRIfitsminus1 <- as.data.frame(cbind(xdummy,ELRI_ydummy_eminus_y1,ELRI_ydummy_eminus_y2, ELRI_ydummy_eminus_y3, ELRI_ydummy_eminus_y4, ELRI_ydummy_eminus_y5, ELRI_ydummy_eminus_y6, ELRI_ydummy_eminus_y7, ELRI_ydummy_eminus_y8, ELRI_ydummy_eminus_y9, ELRI_ydummy_eminus_y10, ELRI_ydummy_eminus_y11))
ELRIfitsminus <- melt(ELRIfitsminus1, id = "xdummy",
                      measure.vars = c("ELRI_ydummy_eminus_y1","ELRI_ydummy_eminus_y2", 
                                       "ELRI_ydummy_eminus_y3", "ELRI_ydummy_eminus_y4", 
                                       "ELRI_ydummy_eminus_y5", "ELRI_ydummy_eminus_y6", 
                                       "ELRI_ydummy_eminus_y7", "ELRI_ydummy_eminus_y8", 
                                       "ELRI_ydummy_eminus_y9", "ELRI_ydummy_eminus_y10",
                                       "ELRI_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "ELRI_ydummy_eminus_y1" = '1',"ELRI_ydummy_eminus_y2" = "2", 
                       "ELRI_ydummy_eminus_y3" = "3", "ELRI_ydummy_eminus_y4" = "4", 
                       "ELRI_ydummy_eminus_y5" = "5", "ELRI_ydummy_eminus_y6" = "6", 
                       "ELRI_ydummy_eminus_y7" = "7", "ELRI_ydummy_eminus_y8" = "8", 
                       "ELRI_ydummy_eminus_y9" = "9", "ELRI_ydummy_eminus_y10" = "10",
                       "ELRI_ydummy_eminus_y11" = "11") )



# E+
ELRI_ydummy_eplus <- as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*1 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                        + mean(post_growELRI$'beta[5]')*xdummy*1))
ELRI_ydummy_eplus_y1 <-  as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*1 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_growELRI$'beta[5]')*xdummy*1
                                            + mean(post_growELRI$'tau_year[2,1]')))
ELRI_ydummy_eplus_y2 <-  as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*1 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_growELRI$'beta[5]')*xdummy*1
                                            + mean(post_growELRI$'tau_year[2,2')))
ELRI_ydummy_eplus_y3 <-  as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*1 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_growELRI$'beta[5]')*xdummy*1
                                            + mean(post_growELRI$'tau_year[2,3]')))
ELRI_ydummy_eplus_y4 <-  as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*1 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_growELRI$'beta[5]')*xdummy*1
                                            + mean(post_growELRI$'tau_year[2,4]')))
ELRI_ydummy_eplus_y5 <-  as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*1 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_growELRI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_growELRI$'tau_year[2,5]')))
ELRI_ydummy_eplus_y6 <-  as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*1 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_growELRI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_growELRI$'tau_year[2,6]')))
ELRI_ydummy_eplus_y7 <-  as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*1 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_growELRI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_growELRI$'tau_year[2,7]')))
ELRI_ydummy_eplus_y8 <-  as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*1 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_growELRI$'beta[5]')*xdummy*1
                                            + mean(post_growELRI$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
ELRI_ydummy_eplus_y9 <-  as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*1 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_growELRI$'beta[5]')*xdummy*1
                                            + mean(post_growELRI$'tau_year[2,9]')))
ELRI_ydummy_eplus_y10 <-  as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*1 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                             + mean(post_growELRI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_growELRI$'tau_year[2,10]')))
ELRI_ydummy_eplus_y11 <-  as.vector(exp(mean(post_growELRI$'beta[1]') + mean(post_growELRI$'beta[2]')*xdummy + mean(post_growELRI$'beta[3]')*1 + mean(post_growELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                             + mean(post_growELRI$'beta[5]')*xdummy*1
                                             + mean(post_growELRI$'tau_year[2,11]')))
ELRI_fitsplus0 <- as.data.frame(cbind(xdummy, ELRI_ydummy_eplus))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
ELRIfitsplus1 <- as.data.frame(cbind(xdummy,ELRI_ydummy_eplus_y1,ELRI_ydummy_eplus_y2, ELRI_ydummy_eplus_y3, ELRI_ydummy_eplus_y4, ELRI_ydummy_eplus_y5, ELRI_ydummy_eplus_y6, ELRI_ydummy_eplus_y7, ELRI_ydummy_eplus_y8, ELRI_ydummy_eplus_y9, ELRI_ydummy_eplus_y10, ELRI_ydummy_eplus_y11))
ELRIfitsplus <- melt(ELRIfitsplus1, id = "xdummy",
                     measure.vars = c("ELRI_ydummy_eplus_y1","ELRI_ydummy_eplus_y2", 
                                      "ELRI_ydummy_eplus_y3", "ELRI_ydummy_eplus_y4", 
                                      "ELRI_ydummy_eplus_y5", "ELRI_ydummy_eplus_y6", 
                                      "ELRI_ydummy_eplus_y7", "ELRI_ydummy_eplus_y8", 
                                      "ELRI_ydummy_eplus_y9", "ELRI_ydummy_eplus_y10",
                                      "ELRI_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "ELRI_ydummy_eplus_y1" = '1',"ELRI_ydummy_eplus_y2" = "2", 
                       "ELRI_ydummy_eplus_y3" = "3", "ELRI_ydummy_eplus_y4" = "4", 
                       "ELRI_ydummy_eplus_y5" = "5", "ELRI_ydummy_eplus_y6" = "6", 
                       "ELRI_ydummy_eplus_y7" = "7", "ELRI_ydummy_eplus_y8" = "8", 
                       "ELRI_ydummy_eplus_y9" = "9", "ELRI_ydummy_eplus_y10" = "10",
                       "ELRI_ydummy_eplus_y11" = "11") )

ggplot(data = ELRIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = ELRIgrow_bin1, aes(x = mean_size, y = mean_grow, color = Year), lwd = 2) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("ELRI E+ Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth")  
scale_color_manual(values=yearcolors)

ggplot(data = ELRIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("ELRI E+ Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# effect on mean


ELRIfits <- as.data.frame(cbind(xdummy, ELRI_ydummy_eminus, ELRI_ydummy_eplus))
ELRIfits <- melt(ELRIfits, id = "xdummy", measure.vars = c("ELRI_ydummy_eplus", "ELRI_ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ELRI_ydummy_eminus" = "E-", "ELRI_ydummy_eplus" = "E+"))

# without point
ggplot(data = ELRIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("ELRI Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = ELRIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = ELRIgrow_binmean, aes(x = mean_size, y = mean_grow, color = endo_01), lwd = 3) +
  ggtitle("ELRI Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth") +   
  scale_color_manual(values=colors2)


# split up plus and minus with datapoints and add mean line
ELRI_minus <- ggplot(data = ELRIfitsminus1)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = ELRIfitsminus0, aes(x = xdummy, y = ELRI_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = ELRIgrow_bin0mean, aes(x = mean_size, y = mean_grow), color ="#ff7f00") +
  theme_classic()+ylim(0,40)+
  labs(title = "E-", x = "log(size_t)", y = "")


ELRI_plus <- ggplot(data = ELRIfitsplus1)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = ELRI_fitsplus0, aes(x = xdummy, y = ELRI_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = ELRIgrow_bin1mean, aes(x = mean_size, y = mean_grow), color ="#6a3d9a") +
  theme_classic()+ylim(0,40)+
  labs(title = "E+", x = "log(size_t)", y = "")

ELRIgrow <- grid.arrange(ELRI_plus, ELRI_minus, ncol= 2)
titleELRIgrow <- annotate_figure(ELRIgrow, top = "ELRI", left = "size_t1")

titleELRIgrow




# E+ by year
ELRIEplusbyyear <- ggplot(data = ELRIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = ELRIgrow_bin1mean, aes(x = mean_size, y = mean_grow, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("size_t1") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = ELRIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("size_t1")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
ELRIEminusbyyear <- ggplot(data = ELRIfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = ELRIgrow_bin0mean, aes(x = mean_size, y = mean_grow, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("size_t1")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = ELRIfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("ELRI E- Growth Probability") + xlab("log(size_t)") + ylab("size_t1")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)



# effect on mean

ELRIfits <- as.data.frame(cbind(xdummy, ELRI_ydummy_eplus, ELRI_ydummy_eminus))


# without point
ggplot(data = ELRIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("ELRI Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = ELRIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = ELRIsurv_binmean, aes(x = mean_size, y = mean_grow, color = endo), lwd = 3) +
  ggtitle("ELRI Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +   
  scale_color_manual(values=colors2)

ELRImean <- ggplot(data = ELRIfits) +
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus), color = "#ff7f00")+ 
  geom_point(data = ELRIgrow_binmean, aes(x = mean_size, y = mean_grow, color = endo_01, lwd = samplesize))+
  ggtitle("Mean Growth") + xlab("log(size_t)") + ylab("size_t1") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")





ELRIgrow4panel <- grid.arrange(ELRImean, ELRI_g,ELRIEplusbyyear,ELRIEminusbyyear, ncol= 4)
titleELRIgrow4panel <- annotate_figure(ELRIgrow4panel, top = "Elymus riparius")

titleELRIgrow4panel


########### AGPE model fits
# actual data points for size_t1
AGPEgrow_bin0 <- AGPE_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())

AGPEgrow_bin1 <- AGPE_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())

AGPEgrow_bin0mean <- AGPE_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 8)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())



AGPEgrow_bin1mean <- AGPE_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 8)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())



AGPEgrow_bin <- AGPE_data %>% 
  select(size_t1, logsize_t, year_t, endo_01) %>%
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())


AGPEgrow_binmean <- AGPE_data %>% 
  select(size_t1, logsize_t, endo_01) %>%
  mutate(endo_01 = as.character(endo_01)) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_grow = mean(size_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(AGPE_data$logsize_t)
xdummy <- seq(min(AGPE_data$logsize_t), max(AGPE_data$logsize_t), length.out = nvalues)

# E- 
AGPE_ydummy_eminus <- as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*0 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                         + mean(post_growAGPE$'beta[5]')*xdummy*0))
AGPE_ydummy_eminus_y1 <- as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*0 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_growAGPE$'beta[5]')*xdummy*0
                                            + mean(post_growAGPE$'tau_year[1,1]')))
AGPE_ydummy_eminus_y2 <- as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*0 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_growAGPE$'beta[5]')*xdummy*0
                                            + mean(post_growAGPE$'tau_year[1,2]')))
AGPE_ydummy_eminus_y3 <- as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + + mean(post_growAGPE$'beta[3]')*0  + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_growAGPE$'beta[5]')*xdummy*0
                                            + mean(post_growAGPE$'tau_year[1,3]')))
AGPE_ydummy_eminus_y4 <- as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*0 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_growAGPE$'beta[5]')*xdummy*0
                                            + mean(post_growAGPE$'tau_year[1,4]'))) 
AGPE_ydummy_eminus_y5 <- as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*0 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_growAGPE$'beta[5]')*xdummy*0
                                            + mean(post_growAGPE$'tau_year[1,5]')))
AGPE_ydummy_eminus_y6 <- as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*0 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_growAGPE$'beta[5]')*xdummy*0
                                            + mean(post_growAGPE$'tau_year[1,6]')))
AGPE_ydummy_eminus_y7 <- as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*0 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_growAGPE$'beta[5]')*xdummy*0
                                            + mean(post_growAGPE$'tau_year[1,7]')))
AGPE_ydummy_eminus_y8 <- as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*0 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_growAGPE$'beta[5]')*xdummy*0
                                            + mean(post_growAGPE$'tau_year[1,8]')))
AGPE_ydummy_eminus_y9 <- as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*0 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_growAGPE$'beta[5]')*xdummy*0
                                            + mean(post_growAGPE$'tau_year[1,9]')))
AGPE_ydummy_eminus_y10 <- as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*0 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                             + mean(post_growAGPE$'beta[5]')*xdummy*0
                                             + mean(post_growAGPE$'tau_year[1,10]')))
AGPE_ydummy_eminus_y11 <- as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*0 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                             + mean(post_growAGPE$'beta[5]')*xdummy*0
                                             + mean(post_growAGPE$'tau_year[1,11]')))
AGPEfitsminus0 <- as.data.frame(cbind(xdummy, AGPE_ydummy_eminus))
AGPEfitsminus1 <- as.data.frame(cbind(xdummy,AGPE_ydummy_eminus_y1,AGPE_ydummy_eminus_y2, AGPE_ydummy_eminus_y3, AGPE_ydummy_eminus_y4, AGPE_ydummy_eminus_y5, AGPE_ydummy_eminus_y6, AGPE_ydummy_eminus_y7, AGPE_ydummy_eminus_y8, AGPE_ydummy_eminus_y9, AGPE_ydummy_eminus_y10, AGPE_ydummy_eminus_y11))
AGPEfitsminus <- melt(AGPEfitsminus1, id = "xdummy",
                      measure.vars = c("AGPE_ydummy_eminus_y1","AGPE_ydummy_eminus_y2", 
                                       "AGPE_ydummy_eminus_y3", "AGPE_ydummy_eminus_y4", 
                                       "AGPE_ydummy_eminus_y5", "AGPE_ydummy_eminus_y6", 
                                       "AGPE_ydummy_eminus_y7", "AGPE_ydummy_eminus_y8", 
                                       "AGPE_ydummy_eminus_y9", "AGPE_ydummy_eminus_y10",
                                       "AGPE_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "AGPE_ydummy_eminus_y1" = '1',"AGPE_ydummy_eminus_y2" = "2", 
                       "AGPE_ydummy_eminus_y3" = "3", "AGPE_ydummy_eminus_y4" = "4", 
                       "AGPE_ydummy_eminus_y5" = "5", "AGPE_ydummy_eminus_y6" = "6", 
                       "AGPE_ydummy_eminus_y7" = "7", "AGPE_ydummy_eminus_y8" = "8", 
                       "AGPE_ydummy_eminus_y9" = "9", "AGPE_ydummy_eminus_y10" = "10",
                       "AGPE_ydummy_eminus_y11" = "11") )



# E+
AGPE_ydummy_eplus <- as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*1 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                        + mean(post_growAGPE$'beta[5]')*xdummy*1))
AGPE_ydummy_eplus_y1 <-  as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*1 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_growAGPE$'beta[5]')*xdummy*1
                                            + mean(post_growAGPE$'tau_year[2,1]')))
AGPE_ydummy_eplus_y2 <-  as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*1 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_growAGPE$'beta[5]')*xdummy*1
                                            + mean(post_growAGPE$'tau_year[2,2')))
AGPE_ydummy_eplus_y3 <-  as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*1 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_growAGPE$'beta[5]')*xdummy*1
                                            + mean(post_growAGPE$'tau_year[2,3]')))
AGPE_ydummy_eplus_y4 <-  as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*1 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_growAGPE$'beta[5]')*xdummy*1
                                            + mean(post_growAGPE$'tau_year[2,4]')))
AGPE_ydummy_eplus_y5 <-  as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*1 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_growAGPE$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_growAGPE$'tau_year[2,5]')))
AGPE_ydummy_eplus_y6 <-  as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*1 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_growAGPE$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_growAGPE$'tau_year[2,6]')))
AGPE_ydummy_eplus_y7 <-  as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*1 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_growAGPE$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_growAGPE$'tau_year[2,7]')))
AGPE_ydummy_eplus_y8 <-  as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*1 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_growAGPE$'beta[5]')*xdummy*1
                                            + mean(post_growAGPE$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
AGPE_ydummy_eplus_y9 <-  as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*1 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_growAGPE$'beta[5]')*xdummy*1
                                            + mean(post_growAGPE$'tau_year[2,9]')))
AGPE_ydummy_eplus_y10 <-  as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*1 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                             + mean(post_growAGPE$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_growAGPE$'tau_year[2,10]')))
AGPE_ydummy_eplus_y11 <-  as.vector(exp(mean(post_growAGPE$'beta[1]') + mean(post_growAGPE$'beta[2]')*xdummy + mean(post_growAGPE$'beta[3]')*1 + mean(post_growAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                             + mean(post_growAGPE$'beta[5]')*xdummy*1
                                             + mean(post_growAGPE$'tau_year[2,11]')))
AGPE_fitsplus0 <- as.data.frame(cbind(xdummy, AGPE_ydummy_eplus))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
AGPEfitsplus1 <- as.data.frame(cbind(xdummy,AGPE_ydummy_eplus_y1,AGPE_ydummy_eplus_y2, AGPE_ydummy_eplus_y3, AGPE_ydummy_eplus_y4, AGPE_ydummy_eplus_y5, AGPE_ydummy_eplus_y6, AGPE_ydummy_eplus_y7, AGPE_ydummy_eplus_y8, AGPE_ydummy_eplus_y9, AGPE_ydummy_eplus_y10, AGPE_ydummy_eplus_y11))
AGPEfitsplus <- melt(AGPEfitsplus1, id = "xdummy",
                     measure.vars = c("AGPE_ydummy_eplus_y1","AGPE_ydummy_eplus_y2", 
                                      "AGPE_ydummy_eplus_y3", "AGPE_ydummy_eplus_y4", 
                                      "AGPE_ydummy_eplus_y5", "AGPE_ydummy_eplus_y6", 
                                      "AGPE_ydummy_eplus_y7", "AGPE_ydummy_eplus_y8", 
                                      "AGPE_ydummy_eplus_y9", "AGPE_ydummy_eplus_y10",
                                      "AGPE_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "AGPE_ydummy_eplus_y1" = '1',"AGPE_ydummy_eplus_y2" = "2", 
                       "AGPE_ydummy_eplus_y3" = "3", "AGPE_ydummy_eplus_y4" = "4", 
                       "AGPE_ydummy_eplus_y5" = "5", "AGPE_ydummy_eplus_y6" = "6", 
                       "AGPE_ydummy_eplus_y7" = "7", "AGPE_ydummy_eplus_y8" = "8", 
                       "AGPE_ydummy_eplus_y9" = "9", "AGPE_ydummy_eplus_y10" = "10",
                       "AGPE_ydummy_eplus_y11" = "11") )

ggplot(data = AGPEfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = AGPEgrow_bin1, aes(x = mean_size, y = mean_grow, color = Year), lwd = 2) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("AGPE E+ Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth")  
scale_color_manual(values=yearcolors)

ggplot(data = AGPEfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("AGPE E+ Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# effect on mean


AGPEfits <- as.data.frame(cbind(xdummy, AGPE_ydummy_eminus, AGPE_ydummy_eplus))
AGPEfits <- melt(AGPEfits, id = "xdummy", measure.vars = c("AGPE_ydummy_eplus", "AGPE_ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "AGPE_ydummy_eminus" = "E-", "AGPE_ydummy_eplus" = "E+"))

# without point
ggplot(data = AGPEfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("AGPE Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = AGPEfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = AGPEgrow_binmean, aes(x = mean_size, y = mean_grow, color = endo_01), lwd = 3) +
  ggtitle("AGPE Growth Probability") + xlab("log(size_t)") + ylab("Prob. of Growth") +   
  scale_color_manual(values=colors2)


# split up plus and minus with datapoints and add mean line
AGPE_minus <- ggplot(data = AGPEfitsminus1)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = AGPEfitsminus0, aes(x = xdummy, y = AGPE_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = AGPEgrow_bin0mean, aes(x = mean_size, y = mean_grow), color ="#ff7f00") +
  theme_classic()+ ylim(0,50)+
  labs(title = "E-", x = "log(size_t)", y = "")


AGPE_plus <- ggplot(data = AGPEfitsplus1)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = AGPE_fitsplus0, aes(x = xdummy, y = AGPE_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = AGPEgrow_bin1mean, aes(x = mean_size, y = mean_grow), color ="#6a3d9a") +
  theme_classic()+ ylim(0,50)+
  labs(title = "E+", x = "log(size_t)", y = "")

AGPEgrow <- grid.arrange(AGPE_plus, AGPE_minus, ncol= 2)
titleAGPEgrow <- annotate_figure(AGPEgrow, top = "AGPE", left = "size_t1")

titleAGPEgrow




# E+ by year
AGPEEplusbyyear <- ggplot(data = AGPEfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = AGPEgrow_bin1mean, aes(x = mean_size, y = mean_grow, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("size_t1") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = AGPEfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("size_t1")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
AGPEEminusbyyear <- ggplot(data = AGPEfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = AGPEgrow_bin0mean, aes(x = mean_size, y = mean_grow, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("size_t1")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = AGPEfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("AGPE E- Growth Probability") + xlab("log(size_t)") + ylab("size_t1")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)



# effect on mean

AGPEfits <- as.data.frame(cbind(xdummy, AGPE_ydummy_eplus, AGPE_ydummy_eminus))


# without point
ggplot(data = AGPEfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("AGPE Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = AGPEfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = AGPEsurv_binmean, aes(x = mean_size, y = mean_grow, color = endo), lwd = 3) +
  ggtitle("AGPE Growth Probability") + xlab("log(size_t)") + ylab("size_t1") +   
  scale_color_manual(values=colors2)

AGPEmean <- ggplot(data = AGPEfits) +
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus), color = "#ff7f00")+ 
  geom_point(data = AGPEgrow_binmean, aes(x = mean_size, y = mean_grow, color = endo_01, lwd = samplesize))+
  ggtitle("Mean Growth") + xlab("log(size_t)") + ylab("size_t1") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")





AGPEgrow4panel <- grid.arrange(AGPEmean, AGPE_g,AGPEEplusbyyear,AGPEEminusbyyear, ncol= 4)
titleAGPEgrow4panel <- annotate_figure(AGPEgrow4panel, top = "Agrostis perennans")

titleAGPEgrow4panel




growfits <- grid.arrange(titleAGPEgrow,titleELRIgrow,titleELVIgrow,titleFESUgrow,titlePOALgrow,titlePOSYgrow, nrow = 3)
growtitle <- annotate_figure(growfits, top = "Growth Probability")
growtitle

grow4panel <- grid.arrange(titleAGPEgrow4panel,titleELRIgrow4panel,titleELVIgrow4panel,titleFESUgrow4panel,titleLOARgrow4panel,titlePOALgrow4panel,titlePOSYgrow4panel, nrow = 7)
grow4paneltitle <- annotate_figure(grow4panel, top = "Growth Probability")

grow4paneltitle






##### Flowering Model fits #####
######## POAL model fits
# actual data points for probability of flowering status
POALflw_bin0 <- POAL_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())

POALflw_bin1 <- POAL_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())
POALflw_bin <- POAL_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())


POALflw_binmean <- POAL_data %>% 
  select(flw_stat_t1, logsize_t, endo_01) %>%
  mutate(endo_01 = recode(endo_01, "0" = "E-", "1" = "E+", minus = 'E-', plus = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(POAL_data$logsize_t)
xdummy <- seq(min(POAL_data$logsize_t), max(POAL_data$logsize_t), length.out = nvalues)

# E- 
POAL_ydummy_eminus <- as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*0 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                         + mean(post_flwPOAL$'beta[5]')*xdummy*0))
POAL_ydummy_eminus_y1 <- as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*0 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*0
                                            + mean(post_flwPOAL$'tau_year[1,1]')))
POAL_ydummy_eminus_y2 <- as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*0 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*0
                                            + mean(post_flwPOAL$'tau_year[1,2]')))
POAL_ydummy_eminus_y3 <- as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + + mean(post_flwPOAL$'beta[3]')*0  + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*0
                                            + mean(post_flwPOAL$'tau_year[1,3]')))
POAL_ydummy_eminus_y4 <- as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*0 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*0
                                            + mean(post_flwPOAL$'tau_year[1,4]'))) 
POAL_ydummy_eminus_y5 <- as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*0 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*0
                                            + mean(post_flwPOAL$'tau_year[1,5]')))
POAL_ydummy_eminus_y6 <- as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*0 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*0
                                            + mean(post_flwPOAL$'tau_year[1,6]')))
POAL_ydummy_eminus_y7 <- as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*0 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*0
                                            + mean(post_flwPOAL$'tau_year[1,7]')))
POAL_ydummy_eminus_y8 <- as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*0 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*0
                                            + mean(post_flwPOAL$'tau_year[1,8]')))
POAL_ydummy_eminus_y9 <- as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*0 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*0
                                            + mean(post_flwPOAL$'tau_year[1,9]')))
POAL_ydummy_eminus_y10 <- as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*0 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                             + mean(post_flwPOAL$'beta[5]')*xdummy*0
                                             + mean(post_flwPOAL$'tau_year[1,10]')))
POAL_ydummy_eminus_y11 <- as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*0 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                             + mean(post_flwPOAL$'beta[5]')*xdummy*0
                                             + mean(post_flwPOAL$'tau_year[1,11]')))
POALfitsminus0 <- as.data.frame(cbind(xdummy, POAL_ydummy_eminus))
POALfitsminus1 <- as.data.frame(cbind(xdummy,POAL_ydummy_eminus_y1,POAL_ydummy_eminus_y2, POAL_ydummy_eminus_y3, POAL_ydummy_eminus_y4, POAL_ydummy_eminus_y5, POAL_ydummy_eminus_y6, POAL_ydummy_eminus_y7, POAL_ydummy_eminus_y8, POAL_ydummy_eminus_y9, POAL_ydummy_eminus_y10, POAL_ydummy_eminus_y11))
POALfitsminus <- melt(POALfitsminus1, id = "xdummy",
                      measure.vars = c("POAL_ydummy_eminus_y1","POAL_ydummy_eminus_y2", 
                                       "POAL_ydummy_eminus_y3", "POAL_ydummy_eminus_y4", 
                                       "POAL_ydummy_eminus_y5", "POAL_ydummy_eminus_y6", 
                                       "POAL_ydummy_eminus_y7", "POAL_ydummy_eminus_y8", 
                                       "POAL_ydummy_eminus_y9", "POAL_ydummy_eminus_y10",
                                       "POAL_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "POAL_ydummy_eminus_y1" = '1',"POAL_ydummy_eminus_y2" = "2", 
                       "POAL_ydummy_eminus_y3" = "3", "POAL_ydummy_eminus_y4" = "4", 
                       "POAL_ydummy_eminus_y5" = "5", "POAL_ydummy_eminus_y6" = "6", 
                       "POAL_ydummy_eminus_y7" = "7", "POAL_ydummy_eminus_y8" = "8", 
                       "POAL_ydummy_eminus_y9" = "9", "POAL_ydummy_eminus_y10" = "10",
                       "POAL_ydummy_eminus_y11" = "11") )



# E+
POAL_ydummy_eplus <- as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*1 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01)
                                        + mean(post_flwPOAL$'beta[5]')*xdummy*1))
POAL_ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*1 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*1
                                            + mean(post_flwPOAL$'tau_year[2,1]')))
POAL_ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*1 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*1
                                            + mean(post_flwPOAL$'tau_year[2,2')))
POAL_ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*1 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*1
                                            + mean(post_flwPOAL$'tau_year[2,3]')))
POAL_ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*1 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*1
                                            + mean(post_flwPOAL$'tau_year[2,4]')))
POAL_ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*1 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_flwPOAL$'tau_year[2,5]')))
POAL_ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*1 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_flwPOAL$'tau_year[2,6]')))
POAL_ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*1 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_flwPOAL$'tau_year[2,7]')))
POAL_ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*1 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*1
                                            + mean(post_flwPOAL$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
POAL_ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*1 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                            + mean(post_flwPOAL$'beta[5]')*xdummy*1
                                            + mean(post_flwPOAL$'tau_year[2,9]')))
POAL_ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*1 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                             + mean(post_flwPOAL$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_flwPOAL$'tau_year[2,10]')))
POAL_ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_flwPOAL$'beta[1]') + mean(post_flwPOAL$'beta[2]')*xdummy + mean(post_flwPOAL$'beta[3]')*1 + mean(post_flwPOAL$'beta[4]')*mean(POAL_data$origin_01) 
                                             + mean(post_flwPOAL$'beta[5]')*xdummy*1
                                             + mean(post_flwPOAL$'tau_year[2,11]')))
POAL_fitsplus0 <- as.data.frame(cbind(xdummy, POAL_ydummy_eplus))  
POALfitsplus1 <- as.data.frame(cbind(xdummy,POAL_ydummy_eplus_y1,POAL_ydummy_eplus_y2, POAL_ydummy_eplus_y3, POAL_ydummy_eplus_y4, POAL_ydummy_eplus_y5, POAL_ydummy_eplus_y6, POAL_ydummy_eplus_y7, POAL_ydummy_eplus_y8, POAL_ydummy_eplus_y9, POAL_ydummy_eplus_y10, POAL_ydummy_eplus_y11))
POALfitsplus <- melt(POALfitsplus1, id = "xdummy",
                     measure.vars = c("POAL_ydummy_eplus_y1","POAL_ydummy_eplus_y2", 
                                      "POAL_ydummy_eplus_y3", "POAL_ydummy_eplus_y4", 
                                      "POAL_ydummy_eplus_y5", "POAL_ydummy_eplus_y6", 
                                      "POAL_ydummy_eplus_y7", "POAL_ydummy_eplus_y8", 
                                      "POAL_ydummy_eplus_y9", "POAL_ydummy_eplus_y10",
                                      "POAL_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "POAL_ydummy_eplus_y1" = '1',"POAL_ydummy_eplus_y2" = "2", 
                       "POAL_ydummy_eplus_y3" = "3", "POAL_ydummy_eplus_y4" = "4", 
                       "POAL_ydummy_eplus_y5" = "5", "POAL_ydummy_eplus_y6" = "6", 
                       "POAL_ydummy_eplus_y7" = "7", "POAL_ydummy_eplus_y8" = "8", 
                       "POAL_ydummy_eplus_y9" = "9", "POAL_ydummy_eplus_y10" = "10",
                       "POAL_ydummy_eplus_y11" = "11") )

ggplot(data = POALfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = POALflw_bin1, aes(x = mean_size, y = mean_flw, color = Year), lwd = 2) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("POAL E+ flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ")  
scale_color_manual(values=yearcolors)

ggplot(data = POALfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("POAL E+ flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# effect on mean


POALfits <- as.data.frame(cbind(xdummy, POAL_ydummy_eminus_woyear, POAL_ydummy_eplus_woyear))
POALfits <- melt(POALfits, id = "xdummy", measure.vars = c("POAL_ydummy_eplus_woyear", "POAL_ydummy_eminus_woyear")) %>% 
  mutate(Endo = recode(variable, "POAL_ydummy_eminus_woyear" = "E-", "POAL_ydummy_eplus_woyear" = "E+"))

# without point
ggplot(data = POALfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("POAL flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = POALfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = POALflw_binmean, aes(x = mean_size, y = mean_flw, color = endo_01), lwd = 3) +
  ggtitle("POAL flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ") +   
  scale_color_manual(values=colors2)


# split up plus and minus with datapoints and add mean line
POAL_minus <- ggplot(data = POALfitsminus1)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = POALfitsminus0, aes(x = xdummy, y = POAL_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = POALflw_bin0, aes(x = mean_size, y = mean_flw), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


POAL_plus <- ggplot(data = POALfitsplus1)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = POAL_fitsplus0, aes(x = xdummy, y = POAL_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = POALflw_bin1, aes(x = mean_size, y = mean_flw), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

POALflw <- grid.arrange(POAL_plus, POAL_minus, ncol= 2)
titlePOALflw <- annotate_figure(POALflw, top = "POAL", left = "Probability of flowering ")

titlePOALflw

# E+ by year
POALEplusbyyear <- ggplot(data = POALfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = POALflw_bin1, aes(x = mean_size, y = mean_flw, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("Prob. of Flowering") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = POALfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("Prob. of Flowering")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
POALEminusbyyear <- ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = POALflw_bin0, aes(x = mean_size, y = mean_flw, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("Prob. of Flowering")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("POAL E- Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)


# effect on mean


POALfits <- as.data.frame(cbind(xdummy, POAL_ydummy_eplus, POAL_ydummy_eminus))

# without point
ggplot(data = POALfits) +
  geom_ribbon(aes(x = xdummy, ymin = POAL_ydummy_eplus2.5, ymax = POAL_ydummy_eplus97.5), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = POAL_ydummy_eminus2.5, ymax = POAL_ydummy_eminus97.5), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus), color = "#ff7f00") +
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus), color = "#6a3d9a")+
  ggtitle("POAL Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +theme_classic()

ggplot(data = POALfits) +
  geom_ribbon(aes(x = xdummy, ymin = POAL_ydummy_eplus_semin, ymax = POAL_ydummy_eplus_semax), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = POAL_ydummy_eminus_semin, ymax = POAL_ydummy_eminus_semax), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus), color = "#ff7f00") +
  ggtitle("POAL Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +theme_classic()

# with points
POALmean <- ggplot(data = POALfits) +
  geom_line(aes(x = xdummy, y = POAL_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = POAL_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = POALflw_binmean, aes(x = mean_size, y = mean_flw, color = endo_01, lwd = samplesize)) +
  ggtitle("Mean Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")



POALflw4panel <- grid.arrange(POALmean, POAL_f,POALEplusbyyear,POALEminusbyyear, ncol= 4)
titlePOALflw4panel <- annotate_figure(POALflw4panel, top = "Poa alsodes")

titlePOALflw4panel

# with data points
ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .7) +
  geom_point(data = POALflw_bin0, aes(x = mean_size, y = mean_flw, color = Year)) +
  labs(title = "POAL E- flowering  Probability", x = "log(size_t)", y = "Prob. of flowering ") +
  scale_color_manual(values=yearcolors)
# without datapoints
ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .7) +
  labs(title = "POAL E- flowering  Probability", x = "log(size_t)", y = "Prob. of flowering ") +                                                                                                                                                                                                                                                                                                                                                                                                                                     scale_color_manual(values=yearcolors)







######## POSY model fits
# actual data points for probability of flowering 
POSYflw_bin0 <- POSY_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())

POSYflw_bin1 <- POSY_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())
POSYflw_bin <- POSY_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())


POSYflw_binmean <- POSY_data %>% 
  select(flw_stat_t1, logsize_t, endo_01) %>%
  mutate(endo_01 = recode(endo_01, "0" = "E-", "1" = "E+", minus = 'E-', plus = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(POSY_data$logsize_t)
xdummy <- seq(min(POSY_data$logsize_t), max(POSY_data$logsize_t), length.out = nvalues)

# E- 
POSY_ydummy_eminus <- as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*0 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                         + mean(post_flwPOSY$'beta[5]')*xdummy*0))
POSY_ydummy_eminus_y1 <- as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*0 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*0
                                            + mean(post_flwPOSY$'tau_year[1,1]')))
POSY_ydummy_eminus_y2 <- as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*0 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*0
                                            + mean(post_flwPOSY$'tau_year[1,2]')))
POSY_ydummy_eminus_y3 <- as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + + mean(post_flwPOSY$'beta[3]')*0  + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*0
                                            + mean(post_flwPOSY$'tau_year[1,3]')))
POSY_ydummy_eminus_y4 <- as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*0 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*0
                                            + mean(post_flwPOSY$'tau_year[1,4]'))) 
POSY_ydummy_eminus_y5 <- as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*0 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*0
                                            + mean(post_flwPOSY$'tau_year[1,5]')))
POSY_ydummy_eminus_y6 <- as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*0 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*0
                                            + mean(post_flwPOSY$'tau_year[1,6]')))
POSY_ydummy_eminus_y7 <- as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*0 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*0
                                            + mean(post_flwPOSY$'tau_year[1,7]')))
POSY_ydummy_eminus_y8 <- as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*0 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*0
                                            + mean(post_flwPOSY$'tau_year[1,8]')))
POSY_ydummy_eminus_y9 <- as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*0 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*0
                                            + mean(post_flwPOSY$'tau_year[1,9]')))
POSY_ydummy_eminus_y10 <- as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*0 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                             + mean(post_flwPOSY$'beta[5]')*xdummy*0
                                             + mean(post_flwPOSY$'tau_year[1,10]')))
POSY_ydummy_eminus_y11 <- as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*0 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                             + mean(post_flwPOSY$'beta[5]')*xdummy*0
                                             + mean(post_flwPOSY$'tau_year[1,11]')))
POSYfitsminus0 <- as.data.frame(cbind(xdummy, POSY_ydummy_eminus))
POSYfitsminus1 <- as.data.frame(cbind(xdummy,POSY_ydummy_eminus_y1,POSY_ydummy_eminus_y2, POSY_ydummy_eminus_y3, POSY_ydummy_eminus_y4, POSY_ydummy_eminus_y5, POSY_ydummy_eminus_y6, POSY_ydummy_eminus_y7, POSY_ydummy_eminus_y8, POSY_ydummy_eminus_y9, POSY_ydummy_eminus_y10, POSY_ydummy_eminus_y11))
POSYfitsminus <- melt(POSYfitsminus1, id = "xdummy",
                      measure.vars = c("POSY_ydummy_eminus_y1","POSY_ydummy_eminus_y2", 
                                       "POSY_ydummy_eminus_y3", "POSY_ydummy_eminus_y4", 
                                       "POSY_ydummy_eminus_y5", "POSY_ydummy_eminus_y6", 
                                       "POSY_ydummy_eminus_y7", "POSY_ydummy_eminus_y8", 
                                       "POSY_ydummy_eminus_y9", "POSY_ydummy_eminus_y10",
                                       "POSY_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "POSY_ydummy_eminus_y1" = '1',"POSY_ydummy_eminus_y2" = "2", 
                       "POSY_ydummy_eminus_y3" = "3", "POSY_ydummy_eminus_y4" = "4", 
                       "POSY_ydummy_eminus_y5" = "5", "POSY_ydummy_eminus_y6" = "6", 
                       "POSY_ydummy_eminus_y7" = "7", "POSY_ydummy_eminus_y8" = "8", 
                       "POSY_ydummy_eminus_y9" = "9", "POSY_ydummy_eminus_y10" = "10",
                       "POSY_ydummy_eminus_y11" = "11") )



# E+
POSY_ydummy_eplus <- as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*1 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01)
                                        + mean(post_flwPOSY$'beta[5]')*xdummy*1))
POSY_ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*1 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*1
                                            + mean(post_flwPOSY$'tau_year[2,1]')))
POSY_ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*1 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*1
                                            + mean(post_flwPOSY$'tau_year[2,2')))
POSY_ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*1 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*1
                                            + mean(post_flwPOSY$'tau_year[2,3]')))
POSY_ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*1 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*1
                                            + mean(post_flwPOSY$'tau_year[2,4]')))
POSY_ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*1 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_flwPOSY$'tau_year[2,5]')))
POSY_ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*1 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_flwPOSY$'tau_year[2,6]')))
POSY_ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*1 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_flwPOSY$'tau_year[2,7]')))
POSY_ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*1 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*1
                                            + mean(post_flwPOSY$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
POSY_ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*1 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                            + mean(post_flwPOSY$'beta[5]')*xdummy*1
                                            + mean(post_flwPOSY$'tau_year[2,9]')))
POSY_ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*1 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                             + mean(post_flwPOSY$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_flwPOSY$'tau_year[2,10]')))
POSY_ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_flwPOSY$'beta[1]') + mean(post_flwPOSY$'beta[2]')*xdummy + mean(post_flwPOSY$'beta[3]')*1 + mean(post_flwPOSY$'beta[4]')*mean(POSY_data$origin_01) 
                                             + mean(post_flwPOSY$'beta[5]')*xdummy*1
                                             + mean(post_flwPOSY$'tau_year[2,11]')))
POSY_fitsplus0 <- as.data.frame(cbind(xdummy, POSY_ydummy_eplus))
POSYfitsplus1 <- as.data.frame(cbind(xdummy,POSY_ydummy_eplus_y1,POSY_ydummy_eplus_y2, POSY_ydummy_eplus_y3, POSY_ydummy_eplus_y4, POSY_ydummy_eplus_y5, POSY_ydummy_eplus_y6, POSY_ydummy_eplus_y7, POSY_ydummy_eplus_y8, POSY_ydummy_eplus_y9, POSY_ydummy_eplus_y10, POSY_ydummy_eplus_y11))
POSYfitsplus <- melt(POSYfitsplus1, id = "xdummy",
                     measure.vars = c("POSY_ydummy_eplus_y1","POSY_ydummy_eplus_y2", 
                                      "POSY_ydummy_eplus_y3", "POSY_ydummy_eplus_y4", 
                                      "POSY_ydummy_eplus_y5", "POSY_ydummy_eplus_y6", 
                                      "POSY_ydummy_eplus_y7", "POSY_ydummy_eplus_y8", 
                                      "POSY_ydummy_eplus_y9", "POSY_ydummy_eplus_y10",
                                      "POSY_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "POSY_ydummy_eplus_y1" = '1',"POSY_ydummy_eplus_y2" = "2", 
                       "POSY_ydummy_eplus_y3" = "3", "POSY_ydummy_eplus_y4" = "4", 
                       "POSY_ydummy_eplus_y5" = "5", "POSY_ydummy_eplus_y6" = "6", 
                       "POSY_ydummy_eplus_y7" = "7", "POSY_ydummy_eplus_y8" = "8", 
                       "POSY_ydummy_eplus_y9" = "9", "POSY_ydummy_eplus_y10" = "10",
                       "POSY_ydummy_eplus_y11" = "11") )

ggplot(data = POSYfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = POSYflw_bin1, aes(x = mean_size, y = mean_flw, color = Year), lwd = 2) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("POSY E+ flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ")  
scale_color_manual(values=yearcolors)

ggplot(data = POSYfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("POSY E+ flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# effect on mean


POSYfits <- as.data.frame(cbind(xdummy, POSY_ydummy_eminus, POSY_ydummy_eplus))
POSYfits <- melt(POSYfits, id = "xdummy", measure.vars = c("POSY_ydummy_eplus", "POSY_ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "POSY_ydummy_eminus" = "E-", "POSY_ydummy_eplus" = "E+"))

# without point
ggplot(data = POSYfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("POSY flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = POSYfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = POSYflw_binmean, aes(x = mean_size, y = mean_flw, color = endo_01), lwd = 3) +
  ggtitle("POSY flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ") +   
  scale_color_manual(values=colors2)


# E+ by year
POSYEplusbyyear <- ggplot(data = POSYfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = POSYflw_bin1, aes(x = mean_size, y = mean_flw, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("Prob. of Flowering") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = POSYfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("Prob. of Flowering")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
POSYEminusbyyear <- ggplot(data = POSYfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = POSYflw_bin0, aes(x = mean_size, y = mean_flw, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("Prob. of Flowering")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = POSYfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("POSY E- Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)


# effect on mean


POSYfits <- as.data.frame(cbind(xdummy, POSY_ydummy_eplus,  POSY_ydummy_eminus))

# without point
ggplot(data = POSYfits) +
  geom_ribbon(aes(x = xdummy, ymin = POSY_ydummy_eplus2.5, ymax = POSY_ydummy_eplus97.5), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = POSY_ydummy_eminus2.5, ymax = POSY_ydummy_eminus97.5), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus), color = "#ff7f00") +
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus), color = "#6a3d9a")+
  ggtitle("POSY Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +theme_classic()

ggplot(data = POSYfits) +
  geom_ribbon(aes(x = xdummy, ymin = POSY_ydummy_eplus_semin, ymax = POSY_ydummy_eplus_semax), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = POSY_ydummy_eminus_semin, ymax = POSY_ydummy_eminus_semax), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus), color = "#ff7f00") +
  ggtitle("POSY Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +theme_classic()

# with points
POSYmean <- ggplot(data = POSYfits) +
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = POSYflw_binmean, aes(x = mean_size, y = mean_flw, color = endo_01, lwd = samplesize)) +
  ggtitle("Mean Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")



POSYflw4panel <- grid.arrange(POSYmean, POSY_f,POSYEplusbyyear,POSYEminusbyyear, ncol= 4)
titlePOSYflw4panel <- annotate_figure(POSYflw4panel, top = "Poa sylvestris")

titlePOSYflw4panel


# split up plus and minus with datapoints and add mean line
POSY_minus <- ggplot(data = POSYfitsminus1)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = POSYfitsminus0, aes(x = xdummy, y = POSY_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = POSYflw_bin0, aes(x = mean_size, y = mean_flw), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


POSY_plus <- ggplot(data = POSYfitsplus1)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = POSY_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = POSY_fitsplus0, aes(x = xdummy, y = POSY_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = POSYflw_bin1, aes(x = mean_size, y = mean_flw), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

POSYflw <- grid.arrange(POSY_plus, POSY_minus, ncol= 2)
titlePOSYflw <- annotate_figure(POSYflw, top = "POSY", left = "Probability of flowering ")

titlePOSYflw







######## LOAR model fits
# actual data points for probability of flowering 
LOARflw_bin0 <- LOAR_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(!is.na(logsize_t)) %>% 
  filter(logsize_t>0) %>% 
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>%  
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())

LOARflw_bin1 <- LOAR_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(logsize_t>=0) %>% 
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())
LOARflw_bin <- LOAR_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(logsize_t>=0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())


LOARflw_binmean <- LOAR_data %>% 
  select(flw_stat_t1, logsize_t, endo_01) %>%
  mutate(endo_01 = recode(endo_01, "0" = "E-", "1" = "E+", minus = 'E-', plus = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(LOAR_data$logsize_t)
xdummy <- seq(min(LOAR_data$logsize_t), max(LOAR_data$logsize_t), length.out = nvalues)

# E- 
LOAR_ydummy_eminus <- as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*0 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                         + mean(post_flwLOAR$'beta[5]')*xdummy*0))
LOAR_ydummy_eminus_y1 <- as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*0 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*0
                                            + mean(post_flwLOAR$'tau_year[1,1]')))
LOAR_ydummy_eminus_y2 <- as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*0 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*0
                                            + mean(post_flwLOAR$'tau_year[1,2]')))
LOAR_ydummy_eminus_y3 <- as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + + mean(post_flwLOAR$'beta[3]')*0  + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*0
                                            + mean(post_flwLOAR$'tau_year[1,3]')))
LOAR_ydummy_eminus_y4 <- as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*0 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*0
                                            + mean(post_flwLOAR$'tau_year[1,4]'))) 
LOAR_ydummy_eminus_y5 <- as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*0 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*0
                                            + mean(post_flwLOAR$'tau_year[1,5]')))
LOAR_ydummy_eminus_y6 <- as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*0 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*0
                                            + mean(post_flwLOAR$'tau_year[1,6]')))
LOAR_ydummy_eminus_y7 <- as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*0 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*0
                                            + mean(post_flwLOAR$'tau_year[1,7]')))
LOAR_ydummy_eminus_y8 <- as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*0 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*0
                                            + mean(post_flwLOAR$'tau_year[1,8]')))
LOAR_ydummy_eminus_y9 <- as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*0 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*0
                                            + mean(post_flwLOAR$'tau_year[1,9]')))
LOAR_ydummy_eminus_y10 <- as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*0 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                             + mean(post_flwLOAR$'beta[5]')*xdummy*0
                                             + mean(post_flwLOAR$'tau_year[1,10]')))
LOAR_ydummy_eminus_y11 <- as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*0 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                             + mean(post_flwLOAR$'beta[5]')*xdummy*0
                                             + mean(post_flwLOAR$'tau_year[1,11]')))
LOARfitsminus0 <- as.data.frame(cbind(xdummy, LOAR_ydummy_eminus))
LOARfitsminus1 <- as.data.frame(cbind(xdummy,LOAR_ydummy_eminus_y1,LOAR_ydummy_eminus_y2, LOAR_ydummy_eminus_y3, LOAR_ydummy_eminus_y4, LOAR_ydummy_eminus_y5, LOAR_ydummy_eminus_y6, LOAR_ydummy_eminus_y7, LOAR_ydummy_eminus_y8, LOAR_ydummy_eminus_y9, LOAR_ydummy_eminus_y10, LOAR_ydummy_eminus_y11))
LOARfitsminus <- melt(LOARfitsminus1, id = "xdummy",
                      measure.vars = c("LOAR_ydummy_eminus_y1","LOAR_ydummy_eminus_y2", 
                                       "LOAR_ydummy_eminus_y3", "LOAR_ydummy_eminus_y4", 
                                       "LOAR_ydummy_eminus_y5", "LOAR_ydummy_eminus_y6", 
                                       "LOAR_ydummy_eminus_y7", "LOAR_ydummy_eminus_y8", 
                                       "LOAR_ydummy_eminus_y9", "LOAR_ydummy_eminus_y10",
                                       "LOAR_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "LOAR_ydummy_eminus_y1" = '1',"LOAR_ydummy_eminus_y2" = "2", 
                       "LOAR_ydummy_eminus_y3" = "3", "LOAR_ydummy_eminus_y4" = "4", 
                       "LOAR_ydummy_eminus_y5" = "5", "LOAR_ydummy_eminus_y6" = "6", 
                       "LOAR_ydummy_eminus_y7" = "7", "LOAR_ydummy_eminus_y8" = "8", 
                       "LOAR_ydummy_eminus_y9" = "9", "LOAR_ydummy_eminus_y10" = "10",
                       "LOAR_ydummy_eminus_y11" = "11") )



# E+
LOAR_ydummy_eplus <- as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*1 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01)
                                        + mean(post_flwLOAR$'beta[5]')*xdummy*1))
LOAR_ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*1 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*1
                                            + mean(post_flwLOAR$'tau_year[2,1]')))
LOAR_ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*1 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*1
                                            + mean(post_flwLOAR$'tau_year[2,2')))
LOAR_ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*1 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*1
                                            + mean(post_flwLOAR$'tau_year[2,3]')))
LOAR_ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*1 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*1
                                            + mean(post_flwLOAR$'tau_year[2,4]')))
LOAR_ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*1 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_flwLOAR$'tau_year[2,5]')))
LOAR_ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*1 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_flwLOAR$'tau_year[2,6]')))
LOAR_ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*1 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_flwLOAR$'tau_year[2,7]')))
LOAR_ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*1 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*1
                                            + mean(post_flwLOAR$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
LOAR_ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*1 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                            + mean(post_flwLOAR$'beta[5]')*xdummy*1
                                            + mean(post_flwLOAR$'tau_year[2,9]')))
LOAR_ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*1 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                             + mean(post_flwLOAR$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_flwLOAR$'tau_year[2,10]')))
LOAR_ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_flwLOAR$'beta[1]') + mean(post_flwLOAR$'beta[2]')*xdummy + mean(post_flwLOAR$'beta[3]')*1 + mean(post_flwLOAR$'beta[4]')*mean(LOAR_data$origin_01) 
                                             + mean(post_flwLOAR$'beta[5]')*xdummy*1
                                             + mean(post_flwLOAR$'tau_year[2,11]')))
LOAR_fitsplus0 <- as.data.frame(cbind(xdummy, LOAR_ydummy_eplus))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
LOARfitsplus1 <- as.data.frame(cbind(xdummy,LOAR_ydummy_eplus_y1,LOAR_ydummy_eplus_y2, LOAR_ydummy_eplus_y3, LOAR_ydummy_eplus_y4, LOAR_ydummy_eplus_y5, LOAR_ydummy_eplus_y6, LOAR_ydummy_eplus_y7, LOAR_ydummy_eplus_y8, LOAR_ydummy_eplus_y9, LOAR_ydummy_eplus_y10, LOAR_ydummy_eplus_y11))
LOARfitsplus <- melt(LOARfitsplus1, id = "xdummy",
                     measure.vars = c("LOAR_ydummy_eplus_y1","LOAR_ydummy_eplus_y2", 
                                      "LOAR_ydummy_eplus_y3", "LOAR_ydummy_eplus_y4", 
                                      "LOAR_ydummy_eplus_y5", "LOAR_ydummy_eplus_y6", 
                                      "LOAR_ydummy_eplus_y7", "LOAR_ydummy_eplus_y8", 
                                      "LOAR_ydummy_eplus_y9", "LOAR_ydummy_eplus_y10",
                                      "LOAR_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "LOAR_ydummy_eplus_y1" = '1',"LOAR_ydummy_eplus_y2" = "2", 
                       "LOAR_ydummy_eplus_y3" = "3", "LOAR_ydummy_eplus_y4" = "4", 
                       "LOAR_ydummy_eplus_y5" = "5", "LOAR_ydummy_eplus_y6" = "6", 
                       "LOAR_ydummy_eplus_y7" = "7", "LOAR_ydummy_eplus_y8" = "8", 
                       "LOAR_ydummy_eplus_y9" = "9", "LOAR_ydummy_eplus_y10" = "10",
                       "LOAR_ydummy_eplus_y11" = "11") )

ggplot(data = LOARfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = LOARflw_bin1, aes(x = mean_size, y = mean_flw, color = Year), lwd = 2) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("LOAR E+ flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ")  
scale_color_manual(values=yearcolors)

ggplot(data = LOARfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("LOAR E+ flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# effect on mean


LOARfits <- as.data.frame(cbind(xdummy, LOAR_ydummy_eminus, LOAR_ydummy_eplus))
LOARfits <- melt(LOARfits, id = "xdummy", measure.vars = c("LOAR_ydummy_eplus", "LOAR_ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "LOAR_ydummy_eminus" = "E-", "LOAR_ydummy_eplus" = "E+"))

# without point
ggplot(data = LOARfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("LOAR flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = LOARfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = LOARflw_binmean, aes(x = mean_size, y = mean_flw, color = endo_01), lwd = 3) +
  ggtitle("LOAR flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ") +   
  scale_color_manual(values=colors2)

# E+ by year
LOAREplusbyyear <- ggplot(data = LOARfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = LOARflw_bin1, aes(x = mean_size, y = mean_flw, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("Prob. of Flowering") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = LOARfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("Prob. of Flowering")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
LOAREminusbyyear <- ggplot(data = LOARfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = LOARflw_bin0, aes(x = mean_size, y = mean_flw, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("Prob. of Flowering")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = LOARfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("LOAR E- Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)


# effect on mean


LOARfits <- as.data.frame(cbind(xdummy, LOAR_ydummy_eplus, LOAR_ydummy_eminus))

# without point
ggplot(data = LOARfits) +
  geom_ribbon(aes(x = xdummy, ymin = LOAR_ydummy_eplus2.5, ymax = LOAR_ydummy_eplus97.5), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = LOAR_ydummy_eminus2.5, ymax = LOAR_ydummy_eminus97.5), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus), color = "#ff7f00") +
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus), color = "#6a3d9a")+
  ggtitle("LOAR Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +theme_classic()

ggplot(data = LOARfits) +
  geom_ribbon(aes(x = xdummy, ymin = LOAR_ydummy_eplus_semin, ymax = LOAR_ydummy_eplus_semax), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = LOAR_ydummy_eminus_semin, ymax = LOAR_ydummy_eminus_semax), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus), color = "#ff7f00") +
  ggtitle("LOAR Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +theme_classic()

# with points
LOARmean <- ggplot(data = LOARfits) +
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = LOARflw_binmean, aes(x = mean_size, y = mean_flw, color = endo_01, lwd = samplesize)) +
  ggtitle("Mean Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")



LOARflw4panel <- grid.arrange(LOARmean, LOAR_f,LOAREplusbyyear,LOAREminusbyyear, ncol= 4)
titleLOARflw4panel <- annotate_figure(LOARflw4panel, top = "Lolium arundinaceae")

titleLOARflw4panel

# split up plus and minus with datapoints and add mean line
LOAR_minus <- ggplot(data = LOARfitsminus1)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = LOARfitsminus0, aes(x = xdummy, y = LOAR_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = LOARflw_bin0, aes(x = mean_size, y = mean_flw), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


LOAR_plus <- ggplot(data = LOARfitsplus1)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = LOAR_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = LOAR_fitsplus0, aes(x = xdummy, y = LOAR_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = LOARflw_bin1, aes(x = mean_size, y = mean_flw), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

LOARflw <- grid.arrange(LOAR_plus, LOAR_minus, ncol= 2)
titleLOARflw <- annotate_figure(LOARflw, top = "LOAR", left = "Probability of flowering ")

titleLOARflw







######## FESU model fits
# actual data points for probability of flowering 
FESUflw_bin0 <- FESU_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())

FESUflw_bin1 <- FESU_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())
FESUflw_bin <- FESU_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())


FESUflw_binmean <- FESU_data %>% 
  select(flw_stat_t1, logsize_t, endo_01) %>%
  mutate(endo_01 = recode(endo_01, "0" = "E-", "1" = "E+", minus = 'E-', plus = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(FESU_data$logsize_t)
xdummy <- seq(min(FESU_data$logsize_t), max(FESU_data$logsize_t), length.out = nvalues)

# E- 
FESU_ydummy_eminus <- as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*0 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01)
                                         + mean(post_flwFESU$'beta[5]')*xdummy*0))
FESU_ydummy_eminus_y1 <- as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*0 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_flwFESU$'beta[5]')*xdummy*0
                                            + mean(post_flwFESU$'tau_year[1,1]')))
FESU_ydummy_eminus_y2 <- as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*0 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_flwFESU$'beta[5]')*xdummy*0
                                            + mean(post_flwFESU$'tau_year[1,2]')))
FESU_ydummy_eminus_y3 <- as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + + mean(post_flwFESU$'beta[3]')*0  + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_flwFESU$'beta[5]')*xdummy*0
                                            + mean(post_flwFESU$'tau_year[1,3]')))
FESU_ydummy_eminus_y4 <- as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*0 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_flwFESU$'beta[5]')*xdummy*0
                                            + mean(post_flwFESU$'tau_year[1,4]'))) 
FESU_ydummy_eminus_y5 <- as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*0 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_flwFESU$'beta[5]')*xdummy*0
                                            + mean(post_flwFESU$'tau_year[1,5]')))
FESU_ydummy_eminus_y6 <- as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*0 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_flwFESU$'beta[5]')*xdummy*0
                                            + mean(post_flwFESU$'tau_year[1,6]')))
FESU_ydummy_eminus_y7 <- as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*0 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_flwFESU$'beta[5]')*xdummy*0
                                            + mean(post_flwFESU$'tau_year[1,7]')))
FESU_ydummy_eminus_y8 <- as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*0 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_flwFESU$'beta[5]')*xdummy*0
                                            + mean(post_flwFESU$'tau_year[1,8]')))
FESU_ydummy_eminus_y9 <- as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*0 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01)
                                            + mean(post_flwFESU$'beta[5]')*xdummy*0
                                            + mean(post_flwFESU$'tau_year[1,9]')))
FESU_ydummy_eminus_y10 <- as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*0 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01)
                                             + mean(post_flwFESU$'beta[5]')*xdummy*0
                                             + mean(post_flwFESU$'tau_year[1,10]')))
# FESU_ydummy_eminus_y11 <- as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*0 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01)
                                             + mean(post_flwFESU$'beta[5]')*xdummy*0
                                             + mean(post_flwFESU$'tau_year[1,11]')))
FESUfitsminus0 <- as.data.frame(cbind(xdummy, FESU_ydummy_eminus))
FESUfitsminus1 <- as.data.frame(cbind(xdummy,FESU_ydummy_eminus_y1,FESU_ydummy_eminus_y2, FESU_ydummy_eminus_y3, FESU_ydummy_eminus_y4, FESU_ydummy_eminus_y5, FESU_ydummy_eminus_y6, FESU_ydummy_eminus_y7, FESU_ydummy_eminus_y8, FESU_ydummy_eminus_y9, FESU_ydummy_eminus_y10))
FESUfitsminus <- melt(FESUfitsminus1, id = "xdummy",
                      measure.vars = c("FESU_ydummy_eminus_y1","FESU_ydummy_eminus_y2", 
                                       "FESU_ydummy_eminus_y3", "FESU_ydummy_eminus_y4", 
                                       "FESU_ydummy_eminus_y5", "FESU_ydummy_eminus_y6", 
                                       "FESU_ydummy_eminus_y7", "FESU_ydummy_eminus_y8", 
                                       "FESU_ydummy_eminus_y9", "FESU_ydummy_eminus_y10"
                                       )) %>% 
  mutate(Year = recode(variable, "FESU_ydummy_eminus_y1" = '1',"FESU_ydummy_eminus_y2" = "2", 
                       "FESU_ydummy_eminus_y3" = "3", "FESU_ydummy_eminus_y4" = "4", 
                       "FESU_ydummy_eminus_y5" = "5", "FESU_ydummy_eminus_y6" = "6", 
                       "FESU_ydummy_eminus_y7" = "7", "FESU_ydummy_eminus_y8" = "8", 
                       "FESU_ydummy_eminus_y9" = "9", "FESU_ydummy_eminus_y10" = "10",
                       "FESU_ydummy_eminus_y11" = "11") )



# E+
FESU_ydummy_eplus <- as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*1 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01)
                                        + mean(post_flwFESU$'beta[5]')*xdummy*1))
FESU_ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*1 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_flwFESU$'beta[5]')*xdummy*1
                                            + mean(post_flwFESU$'tau_year[2,1]')))
FESU_ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*1 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_flwFESU$'beta[5]')*xdummy*1
                                            + mean(post_flwFESU$'tau_year[2,2')))
FESU_ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*1 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_flwFESU$'beta[5]')*xdummy*1
                                            + mean(post_flwFESU$'tau_year[2,3]')))
FESU_ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*1 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_flwFESU$'beta[5]')*xdummy*1
                                            + mean(post_flwFESU$'tau_year[2,4]')))
FESU_ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*1 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_flwFESU$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_flwFESU$'tau_year[2,5]')))
FESU_ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*1 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_flwFESU$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_flwFESU$'tau_year[2,6]')))
FESU_ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*1 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_flwFESU$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_flwFESU$'tau_year[2,7]')))
FESU_ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*1 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_flwFESU$'beta[5]')*xdummy*1
                                            + mean(post_flwFESU$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
FESU_ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*1 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                            + mean(post_flwFESU$'beta[5]')*xdummy*1
                                            + mean(post_flwFESU$'tau_year[2,9]')))
FESU_ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*1 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                             + mean(post_flwFESU$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_flwFESU$'tau_year[2,10]')))
# FESU_ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_flwFESU$'beta[1]') + mean(post_flwFESU$'beta[2]')*xdummy + mean(post_flwFESU$'beta[3]')*1 + mean(post_flwFESU$'beta[4]')*mean(FESU_data$origin_01) 
                                             + mean(post_flwFESU$'beta[5]')*xdummy*1
                                             + mean(post_flwFESU$'tau_year[2,11]')))
FESU_fitsplus0 <- as.data.frame(cbind(xdummy, FESU_ydummy_eplus))
FESUfitsplus1 <- as.data.frame(cbind(xdummy,FESU_ydummy_eplus_y1,FESU_ydummy_eplus_y2, FESU_ydummy_eplus_y3, FESU_ydummy_eplus_y4, FESU_ydummy_eplus_y5, FESU_ydummy_eplus_y6, FESU_ydummy_eplus_y7, FESU_ydummy_eplus_y8, FESU_ydummy_eplus_y9, FESU_ydummy_eplus_y10))
FESUfitsplus <- melt(FESUfitsplus1, id = "xdummy",
                     measure.vars = c("FESU_ydummy_eplus_y1","FESU_ydummy_eplus_y2", 
                                      "FESU_ydummy_eplus_y3", "FESU_ydummy_eplus_y4", 
                                      "FESU_ydummy_eplus_y5", "FESU_ydummy_eplus_y6", 
                                      "FESU_ydummy_eplus_y7", "FESU_ydummy_eplus_y8", 
                                      "FESU_ydummy_eplus_y9", "FESU_ydummy_eplus_y10"
                                      )) %>% 
  mutate(Year = recode(variable,  "FESU_ydummy_eplus_y1" = '1',"FESU_ydummy_eplus_y2" = "2", 
                       "FESU_ydummy_eplus_y3" = "3", "FESU_ydummy_eplus_y4" = "4", 
                       "FESU_ydummy_eplus_y5" = "5", "FESU_ydummy_eplus_y6" = "6", 
                       "FESU_ydummy_eplus_y7" = "7", "FESU_ydummy_eplus_y8" = "8", 
                       "FESU_ydummy_eplus_y9" = "9", "FESU_ydummy_eplus_y10" = "10",
                       "FESU_ydummy_eplus_y11" = "11") )

ggplot(data = FESUfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = FESUflw_bin1, aes(x = mean_size, y = mean_flw, color = Year), lwd = 2) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("FESU E+ flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ")  
scale_color_manual(values=yearcolors)

ggplot(data = FESUfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("FESU E+ flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# effect on mean


FESUfits <- as.data.frame(cbind(xdummy, FESU_ydummy_eminus, FESU_ydummy_eplus))
FESUfits <- melt(FESUfits, id = "xdummy", measure.vars = c("FESU_ydummy_eplus", "FESU_ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "FESU_ydummy_eminus" = "E-", "FESU_ydummy_eplus" = "E+"))

# without point
ggplot(data = FESUfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("FESU flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = FESUfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = FESUflw_binmean, aes(x = mean_size, y = mean_flw, color = endo_01), lwd = 3) +
  ggtitle("FESU flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ") +   
  scale_color_manual(values=colors2)

# E+ by year
FESUEplusbyyear <- ggplot(data = FESUfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = FESUflw_bin1, aes(x = mean_size, y = mean_flw, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("Prob. of Flowering") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = FESUfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("Prob. of Flowering")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
FESUEminusbyyear <- ggplot(data = FESUfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = FESUflw_bin0, aes(x = mean_size, y = mean_flw, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("Prob. of Flowering")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = FESUfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("FESU E- Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)


# effect on mean


FESUfits <- as.data.frame(cbind(xdummy, FESU_ydummy_eplus, FESU_ydummy_eminus))

# without point
ggplot(data = FESUfits) +
  geom_ribbon(aes(x = xdummy, ymin = FESU_ydummy_eplus2.5, ymax = FESU_ydummy_eplus97.5), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = FESU_ydummy_eminus2.5, ymax = FESU_ydummy_eminus97.5), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus), color = "#ff7f00") +
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus), color = "#6a3d9a")+
  ggtitle("FESU Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +theme_classic()

ggplot(data = FESUfits) +
  geom_ribbon(aes(x = xdummy, ymin = FESU_ydummy_eplus_semin, ymax = FESU_ydummy_eplus_semax), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = FESU_ydummy_eminus_semin, ymax = FESU_ydummy_eminus_semax), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus), color = "#ff7f00") +
  ggtitle("FESU Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +theme_classic()

# with points
FESUmean <- ggplot(data = FESUfits) +
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = FESUflw_binmean, aes(x = mean_size, y = mean_flw, color = endo_01, lwd = samplesize)) +
  ggtitle("Mean Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")


FESUflw4panel <- grid.arrange(FESUmean, FESU_f,FESUEplusbyyear,FESUEminusbyyear, ncol= 4)
titleFESUflw4panel <- annotate_figure(FESUflw4panel, top = "Festuca subverticillata")

titleFESUflw4panel

# split up plus and minus with datapoints and add mean line
FESU_minus <- ggplot(data = FESUfitsminus1)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = FESUfitsminus0, aes(x = xdummy, y = FESU_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = FESUflw_bin0, aes(x = mean_size, y = mean_flw), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


FESU_plus <- ggplot(data = FESUfitsplus1)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = FESU_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = FESU_fitsplus0, aes(x = xdummy, y = FESU_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = FESUflw_bin1, aes(x = mean_size, y = mean_flw), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

FESUflw <- grid.arrange(FESU_plus, FESU_minus, ncol= 2)
titleFESUflw <- annotate_figure(FESUflw, top = "FESU", left = "Probability of flowering ")

titleFESUflw






######## ELVI model fits
# actual data points for probability of flowering 
ELVIflw_bin0 <- ELVI_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())

ELVIflw_bin1 <- ELVI_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())
ELVIflw_bin <- ELVI_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())


ELVIflw_binmean <- ELVI_data %>% 
  select(flw_stat_t1, logsize_t, endo_01) %>%
  mutate(endo_01 = recode(endo_01, "0" = "E-", "1" = "E+", minus = 'E-', plus = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(ELVI_data$logsize_t)
xdummy <- seq(min(ELVI_data$logsize_t), max(ELVI_data$logsize_t), length.out = nvalues)

# E- 
ELVI_ydummy_eminus <- as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*0 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                         + mean(post_flwELVI$'beta[5]')*xdummy*0))
ELVI_ydummy_eminus_y1 <- as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*0 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_flwELVI$'beta[5]')*xdummy*0
                                            + mean(post_flwELVI$'tau_year[1,1]')))
ELVI_ydummy_eminus_y2 <- as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*0 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_flwELVI$'beta[5]')*xdummy*0
                                            + mean(post_flwELVI$'tau_year[1,2]')))
ELVI_ydummy_eminus_y3 <- as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + + mean(post_flwELVI$'beta[3]')*0  + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_flwELVI$'beta[5]')*xdummy*0
                                            + mean(post_flwELVI$'tau_year[1,3]')))
ELVI_ydummy_eminus_y4 <- as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*0 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_flwELVI$'beta[5]')*xdummy*0
                                            + mean(post_flwELVI$'tau_year[1,4]'))) 
ELVI_ydummy_eminus_y5 <- as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*0 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_flwELVI$'beta[5]')*xdummy*0
                                            + mean(post_flwELVI$'tau_year[1,5]')))
ELVI_ydummy_eminus_y6 <- as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*0 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_flwELVI$'beta[5]')*xdummy*0
                                            + mean(post_flwELVI$'tau_year[1,6]')))
ELVI_ydummy_eminus_y7 <- as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*0 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_flwELVI$'beta[5]')*xdummy*0
                                            + mean(post_flwELVI$'tau_year[1,7]')))
ELVI_ydummy_eminus_y8 <- as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*0 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_flwELVI$'beta[5]')*xdummy*0
                                            + mean(post_flwELVI$'tau_year[1,8]')))
ELVI_ydummy_eminus_y9 <- as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*0 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                            + mean(post_flwELVI$'beta[5]')*xdummy*0
                                            + mean(post_flwELVI$'tau_year[1,9]')))
ELVI_ydummy_eminus_y10 <- as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*0 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                             + mean(post_flwELVI$'beta[5]')*xdummy*0
                                             + mean(post_flwELVI$'tau_year[1,10]')))
ELVI_ydummy_eminus_y11 <- as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*0 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                             + mean(post_flwELVI$'beta[5]')*xdummy*0
                                             + mean(post_flwELVI$'tau_year[1,11]')))
ELVIfitsminus0 <- as.data.frame(cbind(xdummy, ELVI_ydummy_eminus))
ELVIfitsminus1 <- as.data.frame(cbind(xdummy,ELVI_ydummy_eminus_y1,ELVI_ydummy_eminus_y2, ELVI_ydummy_eminus_y3, ELVI_ydummy_eminus_y4, ELVI_ydummy_eminus_y5, ELVI_ydummy_eminus_y6, ELVI_ydummy_eminus_y7, ELVI_ydummy_eminus_y8, ELVI_ydummy_eminus_y9, ELVI_ydummy_eminus_y10, ELVI_ydummy_eminus_y11))
ELVIfitsminus <- melt(ELVIfitsminus1, id = "xdummy",
                      measure.vars = c("ELVI_ydummy_eminus_y1","ELVI_ydummy_eminus_y2", 
                                       "ELVI_ydummy_eminus_y3", "ELVI_ydummy_eminus_y4", 
                                       "ELVI_ydummy_eminus_y5", "ELVI_ydummy_eminus_y6", 
                                       "ELVI_ydummy_eminus_y7", "ELVI_ydummy_eminus_y8", 
                                       "ELVI_ydummy_eminus_y9", "ELVI_ydummy_eminus_y10",
                                       "ELVI_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "ELVI_ydummy_eminus_y1" = '1',"ELVI_ydummy_eminus_y2" = "2", 
                       "ELVI_ydummy_eminus_y3" = "3", "ELVI_ydummy_eminus_y4" = "4", 
                       "ELVI_ydummy_eminus_y5" = "5", "ELVI_ydummy_eminus_y6" = "6", 
                       "ELVI_ydummy_eminus_y7" = "7", "ELVI_ydummy_eminus_y8" = "8", 
                       "ELVI_ydummy_eminus_y9" = "9", "ELVI_ydummy_eminus_y10" = "10",
                       "ELVI_ydummy_eminus_y11" = "11") )



# E+
ELVI_ydummy_eplus <- as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*1 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01)
                                        + mean(post_flwELVI$'beta[5]')*xdummy*1))
ELVI_ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*1 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_flwELVI$'beta[5]')*xdummy*1
                                            + mean(post_flwELVI$'tau_year[2,1]')))
ELVI_ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*1 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_flwELVI$'beta[5]')*xdummy*1
                                            + mean(post_flwELVI$'tau_year[2,2')))
ELVI_ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*1 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_flwELVI$'beta[5]')*xdummy*1
                                            + mean(post_flwELVI$'tau_year[2,3]')))
ELVI_ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*1 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_flwELVI$'beta[5]')*xdummy*1
                                            + mean(post_flwELVI$'tau_year[2,4]')))
ELVI_ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*1 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_flwELVI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_flwELVI$'tau_year[2,5]')))
ELVI_ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*1 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_flwELVI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_flwELVI$'tau_year[2,6]')))
ELVI_ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*1 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_flwELVI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_flwELVI$'tau_year[2,7]')))
ELVI_ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*1 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_flwELVI$'beta[5]')*xdummy*1
                                            + mean(post_flwELVI$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
ELVI_ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*1 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                            + mean(post_flwELVI$'beta[5]')*xdummy*1
                                            + mean(post_flwELVI$'tau_year[2,9]')))
ELVI_ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*1 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                             + mean(post_flwELVI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_flwELVI$'tau_year[2,10]')))
ELVI_ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_flwELVI$'beta[1]') + mean(post_flwELVI$'beta[2]')*xdummy + mean(post_flwELVI$'beta[3]')*1 + mean(post_flwELVI$'beta[4]')*mean(ELVI_data$origin_01) 
                                             + mean(post_flwELVI$'beta[5]')*xdummy*1
                                             + mean(post_flwELVI$'tau_year[2,11]')))
ELVI_fitsplus0 <- as.data.frame(cbind(xdummy, ELVI_ydummy_eplus))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
ELVIfitsplus1 <- as.data.frame(cbind(xdummy,ELVI_ydummy_eplus_y1,ELVI_ydummy_eplus_y2, ELVI_ydummy_eplus_y3, ELVI_ydummy_eplus_y4, ELVI_ydummy_eplus_y5, ELVI_ydummy_eplus_y6, ELVI_ydummy_eplus_y7, ELVI_ydummy_eplus_y8, ELVI_ydummy_eplus_y9, ELVI_ydummy_eplus_y10, ELVI_ydummy_eplus_y11))
ELVIfitsplus <- melt(ELVIfitsplus1, id = "xdummy",
                     measure.vars = c("ELVI_ydummy_eplus_y1","ELVI_ydummy_eplus_y2", 
                                      "ELVI_ydummy_eplus_y3", "ELVI_ydummy_eplus_y4", 
                                      "ELVI_ydummy_eplus_y5", "ELVI_ydummy_eplus_y6", 
                                      "ELVI_ydummy_eplus_y7", "ELVI_ydummy_eplus_y8", 
                                      "ELVI_ydummy_eplus_y9", "ELVI_ydummy_eplus_y10",
                                      "ELVI_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "ELVI_ydummy_eplus_y1" = '1',"ELVI_ydummy_eplus_y2" = "2", 
                       "ELVI_ydummy_eplus_y3" = "3", "ELVI_ydummy_eplus_y4" = "4", 
                       "ELVI_ydummy_eplus_y5" = "5", "ELVI_ydummy_eplus_y6" = "6", 
                       "ELVI_ydummy_eplus_y7" = "7", "ELVI_ydummy_eplus_y8" = "8", 
                       "ELVI_ydummy_eplus_y9" = "9", "ELVI_ydummy_eplus_y10" = "10",
                       "ELVI_ydummy_eplus_y11" = "11") )

ggplot(data = ELVIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = ELVIflw_bin1, aes(x = mean_size, y = mean_flw, color = Year), lwd = 2) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("ELVI E+ flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ")  
scale_color_manual(values=yearcolors)

ggplot(data = ELVIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("ELVI E+ flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# effect on mean


ELVIfits <- as.data.frame(cbind(xdummy, ELVI_ydummy_eminus, ELVI_ydummy_eplus))
ELVIfits <- melt(ELVIfits, id = "xdummy", measure.vars = c("ELVI_ydummy_eplus", "ELVI_ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ELVI_ydummy_eminus" = "E-", "ELVI_ydummy_eplus" = "E+"))

# without point
ggplot(data = ELVIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("ELVI flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = ELVIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = ELVIflw_binmean, aes(x = mean_size, y = mean_flw, color = endo_01), lwd = 3) +
  ggtitle("ELVI flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ") +   
  scale_color_manual(values=colors2)

# E+ by year
ELVIEplusbyyear <- ggplot(data = ELVIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = ELVIflw_bin1, aes(x = mean_size, y = mean_flw, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("Prob. of Flowering") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = ELVIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("Prob. of Flowering")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
ELVIEminusbyyear <- ggplot(data = ELVIfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = ELVIflw_bin0, aes(x = mean_size, y = mean_flw, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("Prob. of Flowering")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = ELVIfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("ELVI E- Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)


# effect on mean


ELVIfits <- as.data.frame(cbind(xdummy, ELVI_ydummy_eplus, ELVI_ydummy_eminus))

# without point
ggplot(data = ELVIfits) +
  geom_ribbon(aes(x = xdummy, ymin = ELVI_ydummy_eplus2.5, ymax = ELVI_ydummy_eplus97.5), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = ELVI_ydummy_eminus2.5, ymax = ELVI_ydummy_eminus97.5), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus), color = "#ff7f00") +
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus), color = "#6a3d9a")+
  ggtitle("ELVI Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +theme_classic()

ggplot(data = ELVIfits) +
  geom_ribbon(aes(x = xdummy, ymin = ELVI_ydummy_eplus_semin, ymax = ELVI_ydummy_eplus_semax), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = ELVI_ydummy_eminus_semin, ymax = ELVI_ydummy_eminus_semax), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus), color = "#ff7f00") +
  ggtitle("ELVI Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +theme_classic()

# with points
ELVImean <- ggplot(data = ELVIfits) +
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = ELVIflw_binmean, aes(x = mean_size, y = mean_flw, color = endo_01, lwd = samplesize)) +
  ggtitle("Mean Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")


ELVIflw4panel <- grid.arrange(ELVImean, ELVI_f,ELVIEplusbyyear,ELVIEminusbyyear, ncol= 4)
titleELVIflw4panel <- annotate_figure(ELVIflw4panel, top = "Elymus virginicus")

titleELVIflw4panel

# split up plus and minus with datapoints and add mean line
ELVI_minus <- ggplot(data = ELVIfitsminus1)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = ELVIfitsminus0, aes(x = xdummy, y = ELVI_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = ELVIflw_bin0, aes(x = mean_size, y = mean_flw), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


ELVI_plus <- ggplot(data = ELVIfitsplus1)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELVI_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = ELVI_fitsplus0, aes(x = xdummy, y = ELVI_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = ELVIflw_bin1, aes(x = mean_size, y = mean_flw), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

ELVIflw <- grid.arrange(ELVI_plus, ELVI_minus, ncol= 2)
titleELVIflw <- annotate_figure(ELVIflw, top = "ELVI", left = "Probability of flowering ")

titleELVIflw






######## ELRI model fits
# actual data points for probability of flowering 
ELRIflw_bin0 <- ELRI_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())

ELRIflw_bin1 <- ELRI_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())
ELRIflw_bin <- ELRI_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())


ELRIflw_binmean <- ELRI_data %>% 
  select(flw_stat_t1, logsize_t, endo_01) %>%
  mutate(endo_01 = recode(endo_01, "0" = "E-", "1" = "E+", minus = 'E-', plus = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(ELRI_data$logsize_t)
xdummy <- seq(min(ELRI_data$logsize_t), max(ELRI_data$logsize_t), length.out = nvalues)

# E- 
ELRI_ydummy_eminus <- as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*0 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                         + mean(post_flwELRI$'beta[5]')*xdummy*0))
ELRI_ydummy_eminus_y1 <- as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*0 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_flwELRI$'beta[5]')*xdummy*0
                                            + mean(post_flwELRI$'tau_year[1,1]')))
ELRI_ydummy_eminus_y2 <- as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*0 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_flwELRI$'beta[5]')*xdummy*0
                                            + mean(post_flwELRI$'tau_year[1,2]')))
ELRI_ydummy_eminus_y3 <- as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + + mean(post_flwELRI$'beta[3]')*0  + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_flwELRI$'beta[5]')*xdummy*0
                                            + mean(post_flwELRI$'tau_year[1,3]')))
ELRI_ydummy_eminus_y4 <- as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*0 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_flwELRI$'beta[5]')*xdummy*0
                                            + mean(post_flwELRI$'tau_year[1,4]'))) 
ELRI_ydummy_eminus_y5 <- as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*0 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_flwELRI$'beta[5]')*xdummy*0
                                            + mean(post_flwELRI$'tau_year[1,5]')))
ELRI_ydummy_eminus_y6 <- as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*0 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_flwELRI$'beta[5]')*xdummy*0
                                            + mean(post_flwELRI$'tau_year[1,6]')))
ELRI_ydummy_eminus_y7 <- as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*0 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_flwELRI$'beta[5]')*xdummy*0
                                            + mean(post_flwELRI$'tau_year[1,7]')))
ELRI_ydummy_eminus_y8 <- as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*0 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_flwELRI$'beta[5]')*xdummy*0
                                            + mean(post_flwELRI$'tau_year[1,8]')))
ELRI_ydummy_eminus_y9 <- as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*0 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                            + mean(post_flwELRI$'beta[5]')*xdummy*0
                                            + mean(post_flwELRI$'tau_year[1,9]')))
ELRI_ydummy_eminus_y10 <- as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*0 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                             + mean(post_flwELRI$'beta[5]')*xdummy*0
                                             + mean(post_flwELRI$'tau_year[1,10]')))
ELRI_ydummy_eminus_y11 <- as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*0 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                             + mean(post_flwELRI$'beta[5]')*xdummy*0
                                             + mean(post_flwELRI$'tau_year[1,11]')))
ELRIfitsminus0 <- as.data.frame(cbind(xdummy, ELRI_ydummy_eminus))
ELRIfitsminus1 <- as.data.frame(cbind(xdummy,ELRI_ydummy_eminus_y1,ELRI_ydummy_eminus_y2, ELRI_ydummy_eminus_y3, ELRI_ydummy_eminus_y4, ELRI_ydummy_eminus_y5, ELRI_ydummy_eminus_y6, ELRI_ydummy_eminus_y7, ELRI_ydummy_eminus_y8, ELRI_ydummy_eminus_y9, ELRI_ydummy_eminus_y10, ELRI_ydummy_eminus_y11))
ELRIfitsminus <- melt(ELRIfitsminus1, id = "xdummy",
                      measure.vars = c("ELRI_ydummy_eminus_y1","ELRI_ydummy_eminus_y2", 
                                       "ELRI_ydummy_eminus_y3", "ELRI_ydummy_eminus_y4", 
                                       "ELRI_ydummy_eminus_y5", "ELRI_ydummy_eminus_y6", 
                                       "ELRI_ydummy_eminus_y7", "ELRI_ydummy_eminus_y8", 
                                       "ELRI_ydummy_eminus_y9", "ELRI_ydummy_eminus_y10",
                                       "ELRI_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "ELRI_ydummy_eminus_y1" = '1',"ELRI_ydummy_eminus_y2" = "2", 
                       "ELRI_ydummy_eminus_y3" = "3", "ELRI_ydummy_eminus_y4" = "4", 
                       "ELRI_ydummy_eminus_y5" = "5", "ELRI_ydummy_eminus_y6" = "6", 
                       "ELRI_ydummy_eminus_y7" = "7", "ELRI_ydummy_eminus_y8" = "8", 
                       "ELRI_ydummy_eminus_y9" = "9", "ELRI_ydummy_eminus_y10" = "10",
                       "ELRI_ydummy_eminus_y11" = "11") )



# E+
ELRI_ydummy_eplus <- as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*1 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01)
                                        + mean(post_flwELRI$'beta[5]')*xdummy*1))
ELRI_ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*1 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_flwELRI$'beta[5]')*xdummy*1
                                            + mean(post_flwELRI$'tau_year[2,1]')))
ELRI_ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*1 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_flwELRI$'beta[5]')*xdummy*1
                                            + mean(post_flwELRI$'tau_year[2,2')))
ELRI_ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*1 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_flwELRI$'beta[5]')*xdummy*1
                                            + mean(post_flwELRI$'tau_year[2,3]')))
ELRI_ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*1 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_flwELRI$'beta[5]')*xdummy*1
                                            + mean(post_flwELRI$'tau_year[2,4]')))
ELRI_ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*1 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_flwELRI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_flwELRI$'tau_year[2,5]')))
ELRI_ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*1 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_flwELRI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_flwELRI$'tau_year[2,6]')))
ELRI_ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*1 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_flwELRI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_flwELRI$'tau_year[2,7]')))
ELRI_ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*1 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_flwELRI$'beta[5]')*xdummy*1
                                            + mean(post_flwELRI$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
ELRI_ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*1 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                            + mean(post_flwELRI$'beta[5]')*xdummy*1
                                            + mean(post_flwELRI$'tau_year[2,9]')))
ELRI_ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*1 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                             + mean(post_flwELRI$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_flwELRI$'tau_year[2,10]')))
ELRI_ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_flwELRI$'beta[1]') + mean(post_flwELRI$'beta[2]')*xdummy + mean(post_flwELRI$'beta[3]')*1 + mean(post_flwELRI$'beta[4]')*mean(ELRI_data$origin_01) 
                                             + mean(post_flwELRI$'beta[5]')*xdummy*1
                                             + mean(post_flwELRI$'tau_year[2,11]')))
ELRI_fitsplus0 <- as.data.frame(cbind(xdummy, ELRI_ydummy_eplus))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
ELRIfitsplus1 <- as.data.frame(cbind(xdummy,ELRI_ydummy_eplus_y1,ELRI_ydummy_eplus_y2, ELRI_ydummy_eplus_y3, ELRI_ydummy_eplus_y4, ELRI_ydummy_eplus_y5, ELRI_ydummy_eplus_y6, ELRI_ydummy_eplus_y7, ELRI_ydummy_eplus_y8, ELRI_ydummy_eplus_y9, ELRI_ydummy_eplus_y10, ELRI_ydummy_eplus_y11))
ELRIfitsplus <- melt(ELRIfitsplus1, id = "xdummy",
                     measure.vars = c("ELRI_ydummy_eplus_y1","ELRI_ydummy_eplus_y2", 
                                      "ELRI_ydummy_eplus_y3", "ELRI_ydummy_eplus_y4", 
                                      "ELRI_ydummy_eplus_y5", "ELRI_ydummy_eplus_y6", 
                                      "ELRI_ydummy_eplus_y7", "ELRI_ydummy_eplus_y8", 
                                      "ELRI_ydummy_eplus_y9", "ELRI_ydummy_eplus_y10",
                                      "ELRI_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "ELRI_ydummy_eplus_y1" = '1',"ELRI_ydummy_eplus_y2" = "2", 
                       "ELRI_ydummy_eplus_y3" = "3", "ELRI_ydummy_eplus_y4" = "4", 
                       "ELRI_ydummy_eplus_y5" = "5", "ELRI_ydummy_eplus_y6" = "6", 
                       "ELRI_ydummy_eplus_y7" = "7", "ELRI_ydummy_eplus_y8" = "8", 
                       "ELRI_ydummy_eplus_y9" = "9", "ELRI_ydummy_eplus_y10" = "10",
                       "ELRI_ydummy_eplus_y11" = "11") )

ggplot(data = ELRIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = ELRIflw_bin1, aes(x = mean_size, y = mean_flw, color = Year), lwd = 2) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("ELRI E+ flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ")  
scale_color_manual(values=yearcolors)

ggplot(data = ELRIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("ELRI E+ flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# effect on mean


ELRIfits <- as.data.frame(cbind(xdummy, ELRI_ydummy_eminus, ELRI_ydummy_eplus))
ELRIfits <- melt(ELRIfits, id = "xdummy", measure.vars = c("ELRI_ydummy_eplus", "ELRI_ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ELRI_ydummy_eminus" = "E-", "ELRI_ydummy_eplus" = "E+"))

# without point
ggplot(data = ELRIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("ELRI flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = ELRIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = ELRIflw_binmean, aes(x = mean_size, y = mean_flw, color = endo_01), lwd = 3) +
  ggtitle("ELRI flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ") +   
  scale_color_manual(values=colors2)

# E+ by year
ELRIEplusbyyear <- ggplot(data = ELRIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = ELRIflw_bin1, aes(x = mean_size, y = mean_flw, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("Prob. of Flowering") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = ELRIfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("Prob. of Flowering")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
ELRIEminusbyyear <- ggplot(data = ELRIfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = ELRIflw_bin0, aes(x = mean_size, y = mean_flw, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("Prob. of Flowering")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = ELRIfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("ELRI E- Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)


# effect on mean


ELRIfits <- as.data.frame(cbind(xdummy, ELRI_ydummy_eplus, ELRI_ydummy_eminus))

# without point
ggplot(data = ELRIfits) +
  geom_ribbon(aes(x = xdummy, ymin = ELRI_ydummy_eplus2.5, ymax = ELRI_ydummy_eplus97.5), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = ELRI_ydummy_eminus2.5, ymax = ELRI_ydummy_eminus97.5), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus), color = "#ff7f00") +
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus), color = "#6a3d9a")+
  ggtitle("ELRI Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +theme_classic()

ggplot(data = ELRIfits) +
  geom_ribbon(aes(x = xdummy, ymin = ELRI_ydummy_eplus_semin, ymax = ELRI_ydummy_eplus_semax), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = ELRI_ydummy_eminus_semin, ymax = ELRI_ydummy_eminus_semax), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus), color = "#ff7f00") +
  ggtitle("ELRI Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +theme_classic()

# with points
ELRImean <- ggplot(data = ELRIfits) +
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = ELRIflw_binmean, aes(x = mean_size, y = mean_flw, color = endo_01, lwd = samplesize)) +
  ggtitle("Mean Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")



ELRIflw4panel <- grid.arrange(ELRImean, ELRI_f,ELRIEplusbyyear,ELRIEminusbyyear, ncol= 4)
titleELRIflw4panel <- annotate_figure(ELRIflw4panel, top = "Elymus riparius")

titleELRIflw4panel
# split up plus and minus with datapoints and add mean line
ELRI_minus <- ggplot(data = ELRIfitsminus1)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = ELRIfitsminus0, aes(x = xdummy, y = ELRI_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = ELRIflw_bin0, aes(x = mean_size, y = mean_flw), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


ELRI_plus <- ggplot(data = ELRIfitsplus1)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = ELRI_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = ELRI_fitsplus0, aes(x = xdummy, y = ELRI_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = ELRIflw_bin1, aes(x = mean_size, y = mean_flw), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

ELRIflw <- grid.arrange(ELRI_plus, ELRI_minus, ncol= 2)
titleELRIflw <- annotate_figure(ELRIflw, top = "ELRI", left = "Probability of flowering ")

titleELRIflw





########### AGPE model fits
# actual data points for probability of flowering 
AGPEflw_bin0 <- AGPE_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())

AGPEflw_bin1 <- AGPE_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  filter(endo_01 == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())
AGPEflw_bin <- AGPE_data %>% 
  select(flw_stat_t1, logsize_t, year_t, endo_01) %>%
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())


AGPEflw_binmean <- AGPE_data %>% 
  select(flw_stat_t1, logsize_t, endo_01) %>%
  mutate(endo_01 = recode(endo_01, "0" = "E-", "1" = "E+", minus = 'E-', plus = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo_01) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_stat_t1,na.rm=T),
            samplesize = n())
# fitted model line
nvalues <- length(AGPE_data$logsize_t)
xdummy <- seq(min(AGPE_data$logsize_t), max(AGPE_data$logsize_t), length.out = nvalues)

# E- 
AGPE_ydummy_eminus <- as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*0 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                         + mean(post_flwAGPE$'beta[5]')*xdummy*0))
AGPE_ydummy_eminus_y1 <- as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*0 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*0
                                            + mean(post_flwAGPE$'tau_year[1,1]')))
AGPE_ydummy_eminus_y2 <- as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*0 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*0
                                            + mean(post_flwAGPE$'tau_year[1,2]')))
AGPE_ydummy_eminus_y3 <- as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + + mean(post_flwAGPE$'beta[3]')*0  + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*0
                                            + mean(post_flwAGPE$'tau_year[1,3]')))
AGPE_ydummy_eminus_y4 <- as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*0 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*0
                                            + mean(post_flwAGPE$'tau_year[1,4]'))) 
AGPE_ydummy_eminus_y5 <- as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*0 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*0
                                            + mean(post_flwAGPE$'tau_year[1,5]')))
AGPE_ydummy_eminus_y6 <- as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*0 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*0
                                            + mean(post_flwAGPE$'tau_year[1,6]')))
AGPE_ydummy_eminus_y7 <- as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*0 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*0
                                            + mean(post_flwAGPE$'tau_year[1,7]')))
AGPE_ydummy_eminus_y8 <- as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*0 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*0
                                            + mean(post_flwAGPE$'tau_year[1,8]')))
AGPE_ydummy_eminus_y9 <- as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*0 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*0
                                            + mean(post_flwAGPE$'tau_year[1,9]')))
AGPE_ydummy_eminus_y10 <- as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*0 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                             + mean(post_flwAGPE$'beta[5]')*xdummy*0
                                             + mean(post_flwAGPE$'tau_year[1,10]')))
AGPE_ydummy_eminus_y11 <- as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*0 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                             + mean(post_flwAGPE$'beta[5]')*xdummy*0
                                             + mean(post_flwAGPE$'tau_year[1,11]')))
AGPEfitsminus0 <- as.data.frame(cbind(xdummy, AGPE_ydummy_eminus))
AGPEfitsminus1 <- as.data.frame(cbind(xdummy,AGPE_ydummy_eminus_y1,AGPE_ydummy_eminus_y2, AGPE_ydummy_eminus_y3, AGPE_ydummy_eminus_y4, AGPE_ydummy_eminus_y5, AGPE_ydummy_eminus_y6, AGPE_ydummy_eminus_y7, AGPE_ydummy_eminus_y8, AGPE_ydummy_eminus_y9, AGPE_ydummy_eminus_y10, AGPE_ydummy_eminus_y11))
AGPEfitsminus <- melt(AGPEfitsminus1, id = "xdummy",
                      measure.vars = c("AGPE_ydummy_eminus_y1","AGPE_ydummy_eminus_y2", 
                                       "AGPE_ydummy_eminus_y3", "AGPE_ydummy_eminus_y4", 
                                       "AGPE_ydummy_eminus_y5", "AGPE_ydummy_eminus_y6", 
                                       "AGPE_ydummy_eminus_y7", "AGPE_ydummy_eminus_y8", 
                                       "AGPE_ydummy_eminus_y9", "AGPE_ydummy_eminus_y10",
                                       "AGPE_ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "AGPE_ydummy_eminus_y1" = '1',"AGPE_ydummy_eminus_y2" = "2", 
                       "AGPE_ydummy_eminus_y3" = "3", "AGPE_ydummy_eminus_y4" = "4", 
                       "AGPE_ydummy_eminus_y5" = "5", "AGPE_ydummy_eminus_y6" = "6", 
                       "AGPE_ydummy_eminus_y7" = "7", "AGPE_ydummy_eminus_y8" = "8", 
                       "AGPE_ydummy_eminus_y9" = "9", "AGPE_ydummy_eminus_y10" = "10",
                       "AGPE_ydummy_eminus_y11" = "11") )



# E+
AGPE_ydummy_eplus <- as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*1 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01)
                                        + mean(post_flwAGPE$'beta[5]')*xdummy*1))
AGPE_ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*1 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*1
                                            + mean(post_flwAGPE$'tau_year[2,1]')))
AGPE_ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*1 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*1
                                            + mean(post_flwAGPE$'tau_year[2,2')))
AGPE_ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*1 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*1
                                            + mean(post_flwAGPE$'tau_year[2,3]')))
AGPE_ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*1 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*1
                                            + mean(post_flwAGPE$'tau_year[2,4]')))
AGPE_ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*1 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
                                            + mean(post_flwAGPE$'tau_year[2,5]')))
AGPE_ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*1 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
                                            + mean(post_flwAGPE$'tau_year[2,6]')))
AGPE_ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*1 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                                            + mean(post_flwAGPE$'tau_year[2,7]')))
AGPE_ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*1 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*1
                                            + mean(post_flwAGPE$'tau_year[2,8]')))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
AGPE_ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*1 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                            + mean(post_flwAGPE$'beta[5]')*xdummy*1
                                            + mean(post_flwAGPE$'tau_year[2,9]')))
AGPE_ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*1 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                             + mean(post_flwAGPE$'beta[5]')*xdummy*1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
                                             + mean(post_flwAGPE$'tau_year[2,10]')))
AGPE_ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_flwAGPE$'beta[1]') + mean(post_flwAGPE$'beta[2]')*xdummy + mean(post_flwAGPE$'beta[3]')*1 + mean(post_flwAGPE$'beta[4]')*mean(AGPE_data$origin_01) 
                                             + mean(post_flwAGPE$'beta[5]')*xdummy*1
                                             + mean(post_flwAGPE$'tau_year[2,11]')))
AGPE_fitsplus0 <- as.data.frame(cbind(xdummy, AGPE_ydummy_eplus))     
AGPEfitsplus1 <- as.data.frame(cbind(xdummy,AGPE_ydummy_eplus_y1,AGPE_ydummy_eplus_y2, AGPE_ydummy_eplus_y3, AGPE_ydummy_eplus_y4, AGPE_ydummy_eplus_y5, AGPE_ydummy_eplus_y6, AGPE_ydummy_eplus_y7, AGPE_ydummy_eplus_y8, AGPE_ydummy_eplus_y9, AGPE_ydummy_eplus_y10, AGPE_ydummy_eplus_y11))
AGPEfitsplus <- melt(AGPEfitsplus1, id = "xdummy",
                     measure.vars = c("AGPE_ydummy_eplus_y1","AGPE_ydummy_eplus_y2", 
                                      "AGPE_ydummy_eplus_y3", "AGPE_ydummy_eplus_y4", 
                                      "AGPE_ydummy_eplus_y5", "AGPE_ydummy_eplus_y6", 
                                      "AGPE_ydummy_eplus_y7", "AGPE_ydummy_eplus_y8", 
                                      "AGPE_ydummy_eplus_y9", "AGPE_ydummy_eplus_y10",
                                      "AGPE_ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable,  "AGPE_ydummy_eplus_y1" = '1',"AGPE_ydummy_eplus_y2" = "2", 
                       "AGPE_ydummy_eplus_y3" = "3", "AGPE_ydummy_eplus_y4" = "4", 
                       "AGPE_ydummy_eplus_y5" = "5", "AGPE_ydummy_eplus_y6" = "6", 
                       "AGPE_ydummy_eplus_y7" = "7", "AGPE_ydummy_eplus_y8" = "8", 
                       "AGPE_ydummy_eplus_y9" = "9", "AGPE_ydummy_eplus_y10" = "10",
                       "AGPE_ydummy_eplus_y11" = "11") )

ggplot(data = AGPEfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = AGPEflw_bin1, aes(x = mean_size, y = mean_flw, color = Year), lwd = 2) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("AGPE E+ flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ")  
scale_color_manual(values=yearcolors)

ggplot(data = AGPEfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("AGPE E+ flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# effect on mean


AGPEfits <- as.data.frame(cbind(xdummy, AGPE_ydummy_eminus, AGPE_ydummy_eplus))
AGPEfits <- melt(AGPEfits, id = "xdummy", measure.vars = c("AGPE_ydummy_eplus", "AGPE_ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "AGPE_ydummy_eminus" = "E-", "AGPE_ydummy_eplus" = "E+"))

# without point
ggplot(data = AGPEfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("AGPE flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = AGPEfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = AGPEflw_binmean, aes(x = mean_size, y = mean_flw, color = endo_01), lwd = 3) +
  ggtitle("AGPE flowering  Probability") + xlab("log(size_t)") + ylab("Prob. of flowering ") +   
  scale_color_manual(values=colors2)

# E+ by year
AGPEEplusbyyear <- ggplot(data = AGPEfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .8) +
  geom_point(data = AGPEflw_bin1, aes(x = mean_size, y = mean_flw, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E+ ") + xlab("log(size_t)") + ylab("Prob. of Flowering") +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none") +theme_classic()

ggplot(data = AGPEfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("E+") + xlab("log(size_t)") + ylab("Prob. of Flowering")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)

# E- by year
AGPEEminusbyyear <- ggplot(data = AGPEfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = AGPEflw_bin0, aes(x = mean_size, y = mean_flw, color = Year, lwd = samplesize)) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
  ggtitle("E-") + xlab("log(size_t)") + ylab("Prob. of Flowering")  +
  scale_color_brewer(palette = "Paired") + guides(lwd = "none")+theme_classic()

ggplot(data = AGPEfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("AGPE E- Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering")                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
scale_color_manual(values=yearcolors)


# effect on mean


AGPEfits <- as.data.frame(cbind(xdummy, AGPE_ydummy_eplus, AGPE_ydummy_eminus))

# without point
ggplot(data = AGPEfits) +
  geom_ribbon(aes(x = xdummy, ymin = AGPE_ydummy_eplus2.5, ymax = AGPE_ydummy_eplus97.5), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = AGPE_ydummy_eminus2.5, ymax = AGPE_ydummy_eminus97.5), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus), color = "#ff7f00") +
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus), color = "#6a3d9a")+
  ggtitle("AGPE Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +theme_classic()

ggplot(data = AGPEfits) +
  geom_ribbon(aes(x = xdummy, ymin = AGPE_ydummy_eplus_semin, ymax = AGPE_ydummy_eplus_semax), fill = "#6a3d9a", alpha = .2) +
  geom_ribbon(aes(x = xdummy, ymin = AGPE_ydummy_eminus_semin, ymax = AGPE_ydummy_eminus_semax), fill = "#ff7f00", alpha = .2)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus), color = "#ff7f00") +
  ggtitle("AGPE Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +theme_classic()

# with points
AGPEmean <- ggplot(data = AGPEfits) +
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus), color = "#6a3d9a")+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus), color = "#ff7f00") +
  geom_point(data = AGPEflw_binmean, aes(x = mean_size, y = mean_flw, color = endo_01, lwd = samplesize)) +
  ggtitle("Mean Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +
  scale_color_manual(values = colors2) + theme_classic() + guides(lwd = "none")



AGPEflw4panel <- grid.arrange(AGPEmean, AGPE_f,AGPEEplusbyyear,AGPEEminusbyyear, ncol= 4)
titleAGPEflw4panel <- annotate_figure(AGPEflw4panel, top = "Agrostis perennans")

titleAGPEflw4panel

# split up plus and minus with datapoints and add mean line
AGPE_minus <- ggplot(data = AGPEfitsminus1)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y1), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y2), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y3), col = "#ff7f00", alpha = .5) +
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y4), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y5), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y6), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y7), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y8), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y9), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y10), col = "#ff7f00", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eminus_y11), col = "#ff7f00", alpha = .5)+
  geom_line(data = AGPEfitsminus0, aes(x = xdummy, y = AGPE_ydummy_eminus), col = "#ff7f00", lwd = 1) +
  geom_point(data = AGPEflw_bin0, aes(x = mean_size, y = mean_flw), color ="#ff7f00") +
  theme_classic()+
  labs(title = "E-", x = "log(size_t)", y = "")


AGPE_plus <- ggplot(data = AGPEfitsplus1)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y1), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y2), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y3), col = "#6a3d9a", alpha = .5) +
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y4), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y5), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y6), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y7), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y8), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y9), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y10), col = "#6a3d9a", alpha = .5)+
  geom_line(aes(x = xdummy, y = AGPE_ydummy_eplus_y11), col = "#6a3d9a", alpha = .5)+
  geom_line(data = AGPE_fitsplus0, aes(x = xdummy, y = AGPE_ydummy_eplus), col = "#6a3d9a", lwd = 1) +
  geom_point(data = AGPEflw_bin1, aes(x = mean_size, y = mean_flw), color ="#6a3d9a") +
  theme_classic()+
  labs(title = "E+", x = "log(size_t)", y = "")

AGPEflw <- grid.arrange(AGPE_plus, AGPE_minus, ncol= 2)
titleAGPEflw <- annotate_figure(AGPEflw, top = "AGPE", left = "Probability of flowering ")

titleAGPEflw




flwfits <- grid.arrange(titleAGPEflw,titleELRIflw,titleELVIflw,titleFESUflw,titlePOALflw,titlePOSYflw, nrow = 3)
flwtitle <- annotate_figure(flwfits, top = "Flowering  Probability")
flwtitle

flw4panel <- grid.arrange(titleAGPEflw4panel,titleELRIflw4panel,titleELVIflw4panel,titleFESUflw4panel,titleLOARflw4panel,titlePOALflw4panel,titlePOSYflw4panel, nrow = 7)
flw4paneltitle <- annotate_figure(flw4panel, top = "Flowering Probability")

flw4paneltitle

