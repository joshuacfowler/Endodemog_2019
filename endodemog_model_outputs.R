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

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

#  source in raw output from models
## Read in the survival model output for all species
survPOAL <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_POAL.rds")
survPOSY <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_POSY.rds")
survLOAR <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_LOAR.rds")
survELVI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_ELVI.rds")
survELRI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_ELRI.rds")
survFESU <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_FESU.rds")
survAGPE <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_AGPE.rds")


## Read in the flowering model output for all species
flwPOAL <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_POAL.rds")
flwPOSY <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_POSY.rds")
flwLOAR <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_LOAR.rds")
flwELVI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_ELVI.rds")
flwELRI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_ELRI.rds")
flwFESU <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_FESU.rds")
flwAGPE <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_AGPE.rds")


## Read in the growth model output for all species
growPOAL <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_grow_POAL.rds")
growPOSY <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_grow_POSY.rds")
growLOAR <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_grow_LOAR.rds")
growELVI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_grow_ELVI.rds")
growELRI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_grow_ELRI.rds")
growFESU <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_grow_FESU.rds")
growAGPE <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_grow_AGPE.rds")

## Read in the fertility model output for all species
fertPOAL <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_fert_POAL.rds")
fertPOSY <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_fert_POSY.rds")
fertLOAR <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_fert_LOAR.rds")
fertELVI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_fert_ELVI.rds")
fertELRI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_fert_ELRI.rds")
fertFESU <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_fert_FESU.rds")
fertAGPE <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_fert_AGPE.rds")

## Read in the spikelet model output for all species
spikePOAL <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_POAL.rds")
spikePOSY <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_POSY.rds")
spikeLOAR <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_LOAR.rds")
spikeELVI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_ELVI.rds")
spikeELRI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_ELRI.rds")
spikeFESU <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_FESU.rds")
spikeAGPE <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_AGPE.rds")


## Read in the seed mean model output for all species
seedPOAL <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_POAL.rds")
seedPOSY <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_POSY.rds")
seedLOAR <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_LOAR.rds")
# seedELVI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_ELVI.rds")
# seedELRI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_ELRI.rds")
seedFESU <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_FESU.rds")
seedAGPE <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_AGPE.rds")


## Read in the seed to seedling model output for all species
s_to_sPOAL <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_POAL.rds")
s_to_sPOSY <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_POSY.rds")
s_to_sLOAR <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_LOAR.rds")
# s_to_sELVI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_ELVI.rds")
# s_to_sELRI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_ELRI.rds")
s_to_sFESU <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_FESU.rds")
s_to_sAGPE <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_AGPE.rds")


# save posteriors and within dataframes
params = c("beta", "tau_year", "sigma_e", "tau_plot")

post_survPOAL <- as.data.frame(survPOAL, pars = params)
post_survPOSY <- as.data.frame(survPOSY, pars = params)
post_survLOAR <- as.data.frame(survLOAR, pars = params)
post_survELVI <- as.data.frame(survELVI, pars = params)
post_survELRI <- as.data.frame(survELRI, pars = params)
post_survFESU <- as.data.frame(survFESU, pars = params)
post_survAGPE <- as.data.frame(survAGPE, pars = params)

post_flwPOAL <- as.data.frame(flwPOAL, pars = params)
post_flwPOSY <- as.data.frame(flwPOSY, pars = params)
post_flwLOAR <- as.data.frame(flwLOAR, pars = params)
post_flwELVI <- as.data.frame(flwELVI, pars = params)
post_flwELRI <- as.data.frame(flwELRI, pars = params)
post_flwFESU <- as.data.frame(flwFESU, pars = params)
post_flwAGPE <- as.data.frame(flwAGPE, pars = params)

post_growPOAL <- as.data.frame(growPOAL, pars = params)
post_growPOSY <- as.data.frame(growPOSY, pars = params)
post_growLOAR <- as.data.frame(growLOAR, pars = params)
post_growELVI <- as.data.frame(growELVI, pars = params)
post_growELRI <- as.data.frame(growELRI, pars = params)
post_growFESU <- as.data.frame(growFESU, pars = params)
post_growAGPE <- as.data.frame(growAGPE, pars = params)

post_fertPOAL <- as.data.frame(fertPOAL, pars = params)
post_fertPOSY <- as.data.frame(fertPOSY, pars = params)
post_fertLOAR <- as.data.frame(fertLOAR, pars = params)
post_fertELVI <- as.data.frame(fertELVI, pars = params)
post_fertELRI <- as.data.frame(fertELRI, pars = params)
post_fertFESU <- as.data.frame(fertFESU, pars = params)
post_fertAGPE <- as.data.frame(fertAGPE, pars = params)

post_spikePOAL <- as.data.frame(spikePOAL, pars = params)
post_spikePOSY <- as.data.frame(spikePOSY, pars = params)
post_spikeLOAR <- as.data.frame(spikeLOAR, pars = params)
post_spikeELVI <- as.data.frame(spikeELVI, pars = params)
post_spikeELRI <- as.data.frame(spikeELRI, pars = params)
post_spikeFESU <- as.data.frame(spikeFESU, pars = params)
post_spikeAGPE <- as.data.frame(spikeAGPE, pars = params)

post_seedPOAL <- as.data.frame(seedPOAL)
post_seedPOSY <- as.data.frame(seedPOSY)
post_seedLOAR <- as.data.frame(seedLOAR)
# post_seedELVI <- as.data.frame(seedELVI)
# post_seedELRI <- as.data.frame(seedELRI)
post_seedFESU <- as.data.frame(seedFESU)
post_seedAGPE <- as.data.frame(seedAGPE)


post_s_to_sPOAL <- as.data.frame(s_to_sPOAL)
post_s_to_sPOSY <- as.data.frame(s_to_sPOSY)
post_s_to_sLOAR <- as.data.frame(s_to_sLOAR)
# post_s_to_sELVI <- as.data.frame(s_to_sELVI)
# post_s_to_sELRI <- as.data.frame(s_to_sELRI)
post_s_to_sFESU <- as.data.frame(s_to_sFESU)
post_s_to_sAGPE <- as.data.frame(s_to_sAGPE)





