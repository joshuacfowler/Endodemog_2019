## Title: Grass endophyte vital rate variance with a bayesian framework
## Purpose: Load R environment for endodemog_figure.R including packages and vital rate model outputs 
## for Survival, Growth and Fertility models 
## Authors: Josh and Tom

## Packages for making graphs and handling model outputs
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
library(gmodels)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

#  source in raw output from models
# filepaths
joshpath <- "/Users/joshuacfowler/Dropbox/EndodemogData/"
tompath <- "C:/Users/tm9/Dropbox/EndodemogData/"

## Read in the survival model output for all species
survPOAL <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_surv_POAL.rds")
survPOSY <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_surv_POSY.rds")
survLOAR <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_surv_LOAR.rds")
survELVI <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_surv_ELVI.rds")
survELRI <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_surv_ELRI.rds")
survFESU <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_surv_FESU.rds")
survAGPE <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_surv_AGPE.rds")


## Read in the flowering model output for all species
flwPOAL <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_flw_POAL.rds")
flwPOSY <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_flw_POSY.rds")
flwLOAR <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_flw_LOAR.rds")
flwELVI <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_flw_ELVI.rds")
flwELRI <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_flw_ELRI.rds")
flwFESU <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_flw_FESU.rds")
flwAGPE <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_flw_AGPE.rds")


## Read in the growth model output for all species
growPOAL <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_grow_POAL.rds")
growPOSY <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_grow_POSY.rds")
growLOAR <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_grow_LOAR.rds")
growELVI <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_grow_ELVI.rds")
growELRI <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_grow_ELRI.rds")
growFESU <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_grow_FESU.rds")
growAGPE <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_grow_AGPE.rds")

## Read in the fertility model output for all species
fertPOAL <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_fert_POAL.rds")
fertPOSY <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_fert_POSY.rds")
fertLOAR <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_fert_LOAR.rds")
fertELVI <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_fert_ELVI.rds")
fertELRI <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_fert_ELRI.rds")
fertFESU <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_fert_FESU.rds")
fertAGPE <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_fert_AGPE.rds")

## Read in the spikelet model output for all species
spikePOAL <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_spike_POAL.rds")
spikePOSY <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_spike_POSY.rds")
spikeLOAR <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_spike_LOAR.rds")
spikeELVI <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_spike_ELVI.rds")
spikeELRI <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_spike_ELRI.rds")
spikeFESU <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_spike_FESU.rds")
spikeAGPE <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_spike_AGPE.rds")


## Read in the seed mean model output for all species
seedPOAL <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_POAL.rds")
seedPOSY <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_POSY.rds")
seedLOAR <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_LOAR.rds")
seedELVI <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_ELVI.rds")
seedELRI <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_ELRI.rds")
seedFESU <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_FESU.rds")
seedAGPE <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_AGPE.rds")


## Read in the seed to seedling model output for all species
s_to_sPOAL <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_POAL.rds")
s_to_sPOSY <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_POSY.rds")
s_to_sLOAR <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_LOAR.rds")
s_to_sELVI <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_ELVI.rds")
s_to_sELRI <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_ELRI.rds")
s_to_sFESU <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_FESU.rds")
s_to_sAGPE <- read_rds(path = "~/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_AGPE.rds")


# save posteriors and within dataframes

post_survPOAL <- rstan::extract(survPOAL)
post_survPOSY <- rstan::extract(survPOSY)
post_survLOAR <- rstan::extract(survLOAR)
post_survELVI <- rstan::extract(survELVI)
post_survELRI <- rstan::extract(survELRI)
post_survFESU <- rstan::extract(survFESU)
post_survAGPE <- rstan::extract(survAGPE)

post_flwPOAL <- rstan::extract(flwPOAL)
post_flwPOSY <- rstan::extract(flwPOSY)
post_flwLOAR <- rstan::extract(flwLOAR)
post_flwELVI <- rstan::extract(flwELVI)
post_flwELRI <- rstan::extract(flwELRI)
post_flwFESU <- rstan::extract(flwFESU)
post_flwAGPE <- rstan::extract(flwAGPE)

post_growPOAL <- rstan::extract(growPOAL)
post_growPOSY <- rstan::extract(growPOSY)
post_growLOAR <- rstan::extract(growLOAR)
post_growELVI <- rstan::extract(growELVI)
post_growELRI <- rstan::extract(growELRI)
post_growFESU <- rstan::extract(growFESU)
post_growAGPE <- rstan::extract(growAGPE)

post_fertPOAL <- rstan::extract(fertPOAL)
post_fertPOSY <- rstan::extract(fertPOSY)
post_fertLOAR <- rstan::extract(fertLOAR)
post_fertELVI <- rstan::extract(fertELVI)
post_fertELRI <- rstan::extract(fertELRI)
post_fertFESU <- rstan::extract(fertFESU)
post_fertAGPE <- rstan::extract(fertAGPE)

post_spikePOAL <- rstan::extract(spikePOAL)
post_spikePOSY <- rstan::extract(spikePOSY)
post_spikeLOAR <- rstan::extract(spikeLOAR)
post_spikeELVI <- rstan::extract(spikeELVI)
post_spikeELRI <- rstan::extract(spikeELRI)
post_spikeFESU <- rstan::extract(spikeFESU)
post_spikeAGPE <- rstan::extract(spikeAGPE)

post_seedPOAL <- rstan::extract(seedPOAL)
post_seedPOSY <- rstan::extract(seedPOSY)
post_seedLOAR <- rstan::extract(seedLOAR)
post_seedELVI <- rstan::extract(seedELVI)
post_seedELRI <- rstan::extract(seedELRI)
post_seedFESU <- rstan::extract(seedFESU)
post_seedAGPE <- rstan::extract(seedAGPE)


post_s_to_sPOAL <- rstan::extract(s_to_sPOAL)
post_s_to_sPOSY <- rstan::extract(s_to_sPOSY)
post_s_to_sLOAR <- rstan::extract(s_to_sLOAR)
post_s_to_sELVI <- rstan::extract(s_to_sELVI)
post_s_to_sELRI <- rstan::extract(s_to_sELRI)
post_s_to_sFESU <- rstan::extract(s_to_sFESU)
post_s_to_sAGPE <- rstan::extract(s_to_sAGPE)





