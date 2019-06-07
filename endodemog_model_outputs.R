## Title: Grass endophyte vital rate variance with a bayesian framework
## Purpose: Load R environment for endodemog_figure.R including packages and vital rate model outputs 
## for Survival, Growth and Fertility models 
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
survLOAR <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_LOAR_withplot.rds")
survELVI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_ELVI_withplot.rds")
survELRI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_ELRI_withplot.rds")
survFESU <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_FESU_withplot.rds")
survAGPE <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_AGPE_withplot.rds")


## Read in the flowering model output for all species
flwPOAL <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_POAL_withplot.rds")
flwPOSY <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_POSY_withplot.rds")
flwLOAR <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_LOAR_withplot.rds")
flwELVI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_ELVI_withplot.rds")
flwELRI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_ELRI_withplot.rds")
flwFESU <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_FESU_withplot.rds")
flwAGPE <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_AGPE_withplot.rds")


## Read in the growth model output for all species
growPOAL <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_grow_POAL_withplot.rds")
growPOSY <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_grow_POSY_withplot.rds")
growLOAR <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_grow_LOAR_withplot.rds")
growELVI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_grow_ELVI_withplot.rds")
growELRI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_grow_ELRI_withplot.rds")
growFESU <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_grow_FESU_withplot.rds")
growAGPE <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_grow_AGPE_withplot.rds")

## Read in the fertility model output for all species
fertPOAL <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_fert_POAL_withplot.rds")
fertPOSY <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_fert_POSY_withplot.rds")
fertLOAR <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_fert_LOAR_withplot.rds")
fertELVI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_fert_ELVI_withplot.rds")
fertELRI <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_fert_ELRI_withplot.rds")
fertFESU <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_fert_FESU_withplot.rds")
fertAGPE <- read_rds(path = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_fert_AGPE_withplot.rds")


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

