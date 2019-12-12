## Title: Grass endophyte population model with a bayesian framework
## Purpose: Assemble matrix model from vital rate estimates
## Authors: Joshua and Tom
## Updated: 11/19/2019
#############################################################
library(tidyverse)
library(rstan)
library(StanHeaders)
library(popbio)
library(bbmle)
library(SPEI)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

#############################################################################################
####### Data manipulation to prepare data as lists for Stan models------------------
#############################################################################################
# filepaths
joshpath <- "/Users/joshuacfowler/Dropbox/EndodemogData/"
tompath <- "C:/Users/tm9/Dropbox/EndodemogData/"
# Growth data lists are generated in the endodemog_data_processing.R file
# within the section titled "Preparing datalists for Growth Kernel"
source("endodemog_data_processing.R")
# endo_demog <- read.csv("~/Dropbox/EndodemogData/fulldataplusmetadata/endo_demog_long.csv")
# loar <- endo_demog %>% 
#   filter(species == "LOAR") %>% 
#   mutate(log_size_t = log(size_t),
#          log_size_t1 = log(size_t1),
#          flow_t = as.integer(seed_t > 0))
# 
# # DATA ISSUES:
# # I noticed that there are about 10 plants with a size of 0 tillers. They are:
# filter(loar,size_t==0)
# # I also noticed that there are a few instances where survival was > 1:
# filter(loar,surv_t1 > 1)
# # We'll need to cut these for the analysis.
# loar <- loar %>% filter(size_t > 0,
#                         surv_t1 <= 1)

# 
# #  Read in the model outputs ----------------------------------------------
source("endodemog_model_outputs.R")

# survLOAR <- read_rds(path = paste(joshpath, "Model_Runs/endodemog_surv_LOAR.rds", sep = ""))
# 
# growLOAR <- read_rds(path = paste(joshpath, "Model_Runs/endodemog_grow_LOAR.rds", sep = ""))
# flwLOAR <- read_rds(path =  paste(joshpath, "Model_Runs/endodemog_flw_LOAR.rds", sep = ""))
# 
# fertLOAR <- read_rds(path =  paste(joshpath, "Model_Runs/endodemog_fert_LOAR.rds", sep = ""))
# 
# spikeLOAR <- read_rds(path = paste(joshpath, "Model_Runs/endodemog_spike_LOAR.rds", sep = ""))
# seedLOAR <- read_rds(path = paste(joshpath, "Model_Runs/endodemog_seed_mean_LOAR.rds", sep = ""))
# s_to_sLOAR <- read_rds(path = paste(joshpath, "Model_Runs/endodemog_s_to_s_LOAR.rds", sep = ""))

#############################################################################################
# Assembling matrix model -------------------------------------------------
#############################################################################################

# collect all the parameters into one vector

# This is the linear predictor for surv, grow, flw, fert, spike
# mu[n] = beta[1] + beta[2]*logsize_t[n] + beta[3]*endo_01[n] +beta[4]*origin_01[n]
# + tau_year[endo_index[n],year_t[n]]
# + tau_plot[plot[n]];
#

# creates vector for no mean endo effect and no endo variance effect for recruits
params_0m_0v_rec <- function(surv, grow,flw, fert, spike, s_to_s, seed, data){
  params <- c()
  params[1:11] <- mean(surv$beta[,1]) + mean(surv$beta[,4]) + colMeans(surv$tau_year[,1,1:11]); names(params)[1:11] <- paste0("surv_intercept_y", 1:11)
  params[12] <- mean(surv$beta[,2]); names(params)[12] <- "surv_slope"
  params[13:23] <-  mean(grow$beta[,1]) + mean(grow$beta[,4]) + colMeans(grow$tau_year[,1,1:11]); names(params)[13:23] <- paste0("grow_intercept_y", 1:11)
  params[24] <- mean(grow$beta[,2]); names(params)[24] <- "grow_slope"
  params[25:35] <- mean(flw$beta[,1]) + mean(flw$beta[,4]) + colMeans(flw$tau_year[,1,1:11]); names(params)[25:35] <- paste0("flw_intercept_y", 1:11)
  params[36] <- mean(flw$beta[,2]); names(params)[36] <- "flw_slope"
  params[37:47] <- mean(fert$beta[,1]) + mean(fert$beta[,4]) + colMeans(fert$tau_year[,1,1:11]); names(params)[37:47] <- paste0("fert_intercept_y", 1:11)
  params[48] <- mean(fert$beta[,2]); names(params)[48] <- "fert_slope"
  params[49:59] <- mean(spike$beta[,1]) + mean(spike$beta[,4]) + colMeans(spike$tau_year[,1,1:11]); names(params)[49:59] <- paste0("spike_intercept_y", 1:11)
  params[60] <- mean(spike$beta[,2]); names(params)[60] <- "spike_slope"
  params[61:71] <- mean(s_to_s$beta[,1]) + colMeans(s_to_s$tau_year[,1,1:11]); names(params)[61:71] <- paste0("s_to_s_intercept_y", 1:11)
  params[72] <- mean(seed$mu_seed); names(params)[72] <- "mu_seed"
  params[73] <- mean(grow$phi); names(params)[73] <- "grow_phi"
  params[74] <- mean(fert$phi); names(params)[74] <- "fert_phi"
  params[75] <- mean(spike$phi); names(params)[75] <- "spike_phi"
  params[76] <- 1; names(params)[76]<-"min_size"
  params[77] <- quantile(data$size_t1,0.95,na.rm=T); names(params)[77]<-"max_size"
  return(params)
}

# creates vector for mean endo effect and no endo variance effect for recruits
params_1m_0v_rec <- function(surv, grow,flw, fert, spike, s_to_s, seed, data){
  params <- c()
  params[1:11] <- mean(surv$beta[,1]) + mean(surv$beta[,3]) + mean(surv$beta[,4]) + colMeans(surv$tau_year[,1,1:11]); names(params)[1:11] <- paste0("surv_intercept_y", 1:11)
  params[12] <- mean(surv$beta[,2]); names(params)[12] <- "surv_slope"
  params[13:23] <-  mean(grow$beta[,1]) + mean(grow$beta[,3]) + mean(grow$beta[,4]) + colMeans(grow$tau_year[,1,1:11]); names(params)[13:23] <- paste0("grow_intercept_y", 1:11)
  params[24] <- mean(grow$beta[,2]); names(params)[24] <- "grow_slope"
  params[25:35] <- mean(flw$beta[,1]) + mean(flw$beta[,3]) + mean(flw$beta[,4]) + colMeans(flw$tau_year[,1,1:11]); names(params)[25:35] <- paste0("flw_intercept_y", 1:11)
  params[36] <- mean(flw$beta[,2]); names(params)[36] <- "flw_slope"
  params[37:47] <- mean(fert$beta[,1]) + mean(fert$beta[,3]) + mean(fert$beta[,4]) + colMeans(fert$tau_year[,1,1:11]); names(params)[37:47] <- paste0("fert_intercept_y", 1:11)
  params[48] <- mean(fert$beta[,2]); names(params)[48] <- "fert_slope"
  params[49:59] <- mean(spike$beta[,1]) + mean(spike$beta[,3]) + mean(spike$beta[,4]) + colMeans(spike$tau_year[,1,1:11]); names(params)[49:59] <- paste0("spike_intercept_y", 1:11)
  params[60] <- mean(spike$beta[,2]); names(params)[60] <- "spike_slope"
  params[61:71] <- mean(s_to_s$beta[,1]) + mean(s_to_s$beta[,2]) + colMeans(s_to_s$tau_year[,1,1:11]); names(params)[61:71] <- paste0("s_to_s_intercept_y", 1:11)
  params[72] <- mean(seed$mu_seed); names(params)[72] <- "mu_seed"
  params[73] <- mean(grow$phi); names(params)[73] <- "grow_phi"
  params[74] <- mean(fert$phi); names(params)[74] <- "fert_phi"
  params[75] <- mean(spike$phi); names(params)[75] <- "spike_phi"
  params[76] <- 1; names(params)[76]<-"min_size"
  params[77] <- quantile(data$size_t1,0.95,na.rm=T); names(params)[77]<-"max_size"
  return(params)
}

# creates vector for no mean endo effect and with endo variance effect for recruits
params_0m_1v_rec <- function(surv, grow,flw, fert, spike, s_to_s, seed, data){
  params <- c()
  params[1:11] <- mean(surv$beta[,1]) + mean(surv$beta[,4]) + colMeans(surv$tau_year[,2,1:11]); names(params)[1:11] <- paste0("surv_intercept_y", 1:11)
  params[12] <- mean(surv$beta[,2]); names(params)[12] <- "surv_slope"
  params[13:23] <-  mean(grow$beta[,1]) + mean(grow$beta[,4]) + colMeans(grow$tau_year[,2,1:11]); names(params)[13:23] <- paste0("grow_intercept_y", 1:11)
  params[24] <- mean(grow$beta[,2]); names(params)[24] <- "grow_slope"
  params[25:35] <- mean(flw$beta[,1]) + mean(flw$beta[,4]) + colMeans(flw$tau_year[,2,1:11]); names(params)[25:35] <- paste0("flw_intercept_y", 1:11)
  params[36] <- mean(flw$beta[,2]); names(params)[36] <- "flw_slope"
  params[37:47] <- mean(fert$beta[,1]) + mean(fert$beta[,4]) + colMeans(fert$tau_year[,2,1:11]); names(params)[37:47] <- paste0("fert_intercept_y", 1:11)
  params[48] <- mean(fert$beta[,2]); names(params)[48] <- "fert_slope"
  params[49:59] <- mean(spike$beta[,1]) + mean(spike$beta[,4]) + colMeans(spike$tau_year[,2,1:11]); names(params)[49:59] <- paste0("spike_intercept_y", 1:11)
  params[60] <- mean(spike$beta[,2]); names(params)[60] <- "spike_slope"
  params[61:71] <- mean(s_to_s$beta[,1]) + colMeans(s_to_s$tau_year[,2,1:11]); names(params)[61:71] <- paste0("s_to_s_intercept_y", 1:11)
  params[72] <- mean(seed$mu_seed); names(params)[72] <- "mu_seed"
  params[73] <- mean(grow$phi); names(params)[73] <- "grow_phi"
  params[74] <- mean(fert$phi); names(params)[74] <- "fert_phi"
  params[75] <- mean(spike$phi); names(params)[75] <- "spike_phi"
  params[76] <- 1; names(params)[76]<-"min_size"
  params[77] <- quantile(data$size_t1,0.95,na.rm=T); names(params)[77]<-"max_size"
  return(params)
}

# creates vector for mean endo effect and with endo variance effect for recruits
params_1m_1v_rec <- function(surv, grow,flw, fert, spike, s_to_s, seed, data){
  params <- c()
  params[1:11] <- mean(surv$beta[,1]) + mean(surv$beta[,3]) + mean(surv$beta[,4]) + colMeans(surv$tau_year[,2,1:11]); names(params)[1:11] <- paste0("surv_intercept_y", 1:11)
  params[12] <- mean(surv$beta[,2]); names(params)[12] <- "surv_slope"
  params[13:23] <-  mean(grow$beta[,1]) + mean(grow$beta[,3]) + mean(grow$beta[,4]) + colMeans(grow$tau_year[,2,1:11]); names(params)[13:23] <- paste0("grow_intercept_y", 1:11)
  params[24] <- mean(grow$beta[,2]); names(params)[24] <- "grow_slope"
  params[25:35] <- mean(flw$beta[,1]) + mean(flw$beta[,3]) + mean(flw$beta[,4]) + colMeans(flw$tau_year[,2,1:11]); names(params)[25:35] <- paste0("flw_intercept_y", 1:11)
  params[36] <- mean(flw$beta[,2]); names(params)[36] <- "flw_slope"
  params[37:47] <- mean(fert$beta[,1]) + mean(fert$beta[,3]) + mean(fert$beta[,4]) + colMeans(fert$tau_year[,2,1:11]); names(params)[37:47] <- paste0("fert_intercept_y", 1:11)
  params[48] <- mean(fert$beta[,2]); names(params)[48] <- "fert_slope"
  params[49:59] <- mean(spike$beta[,1]) + mean(spike$beta[,3]) + mean(spike$beta[,4]) + colMeans(spike$tau_year[,2,1:11]); names(params)[49:59] <- paste0("spike_intercept_y", 1:11)
  params[60] <- mean(spike$beta[,2]); names(params)[60] <- "spike_slope"
  params[61:71] <- mean(s_to_s$beta[,1]) + mean(s_to_s$beta[,2]) + colMeans(s_to_s$tau_year[,2,1:11]); names(params)[61:71] <- paste0("s_to_s_intercept_y", 1:11)
  params[72] <- mean(seed$mu_seed); names(params)[72] <- "mu_seed"
  params[73] <- mean(grow$phi); names(params)[73] <- "grow_phi"
  params[74] <- mean(fert$phi); names(params)[74] <- "fert_phi"
  params[75] <- mean(spike$phi); names(params)[75] <- "spike_phi"
  params[76] <- 1; names(params)[76]<-"min_size"
  params[77] <- quantile(data$size_t1,0.95,na.rm=T); names(params)[77]<-"max_size"
  return(params)
}
# params[28] <- quantile(data$size_t1,0.95,na.rm=T); names(params)[28]<-"max_size"


agpe_0m_0v_rec_params <- params_0m_0v_rec(surv = post_survAGPE,
                                          grow = post_growAGPE, 
                                          flw = post_flwAGPE, 
                                          fert = post_fertAGPE,
                                          spike = post_spikeAGPE,
                                          s_to_s = post_s_to_sAGPE, 
                                          seed = post_seedAGPE, 
                                          data = AGPE_surv_data)
elri_0m_0v_rec_params <- params_0m_0v_rec(surv = post_survELRI,
                                          grow = post_growELRI, 
                                          flw = post_flwELRI, 
                                          fert = post_fertELRI,
                                          spike = post_spikeELRI,
                                          s_to_s = post_s_to_sELRI, 
                                          seed = post_seedELRI, 
                                          data = ELRI_surv_data)
elvi_0m_0v_rec_params <- params_0m_0v_rec(surv = post_survELVI,
                                          grow = post_growELVI, 
                                          flw = post_flwELVI, 
                                          fert = post_fertELVI,
                                          spike = post_spikeELVI,
                                          s_to_s = post_s_to_sELVI, 
                                          seed = post_seedELVI, 
                                          data = ELVI_surv_data)
fesu_0m_0v_rec_params <- params_0m_0v_rec(surv = post_survFESU,
                                          grow = post_growFESU, 
                                          flw = post_flwFESU, 
                                          fert = post_fertFESU,
                                          spike = post_spikeFESU,
                                          s_to_s = post_s_to_sFESU, 
                                          seed = post_seedFESU, 
                                          data = FESU_surv_data)
loar_0m_0v_rec_params <- params_0m_0v_rec(surv = post_survLOAR,
                                          grow = post_growLOAR, 
                                          flw = post_flwLOAR, 
                                          fert = post_fertLOAR,
                                          spike = post_spikeLOAR,
                                          s_to_s = post_s_to_sLOAR, 
                                          seed = post_seedLOAR, 
                                          data = LOAR_surv_data)
poal_0m_0v_rec_params <- params_0m_0v_rec(surv = post_survPOAL,
                                          grow = post_growPOAL, 
                                          flw = post_flwPOAL, 
                                          fert = post_fertPOAL,
                                          spike = post_spikePOAL,
                                          s_to_s = post_s_to_sPOAL, 
                                          seed = post_seedPOAL, 
                                          data = POAL_surv_data)
posy_0m_0v_rec_params <- params_0m_0v_rec(surv = post_survPOSY,
                                          grow = post_growPOSY, 
                                          flw = post_flwPOSY, 
                                          fert = post_fertPOSY,
                                          spike = post_spikePOSY,
                                          s_to_s = post_s_to_sPOSY, 
                                          seed = post_seedPOSY, 
                                          data = POSY_surv_data)





agpe_1m_0v_rec_params <- params_1m_0v_rec(surv = post_survAGPE,
                                          grow = post_growAGPE, 
                                          flw = post_flwAGPE, 
                                          fert = post_fertAGPE,
                                          spike = post_spikeAGPE,
                                          s_to_s = post_s_to_sAGPE, 
                                          seed = post_seedAGPE, 
                                          data = AGPE_surv_data)
elri_1m_0v_rec_params <- params_1m_0v_rec(surv = post_survELRI,
                                          grow = post_growELRI, 
                                          flw = post_flwELRI, 
                                          fert = post_fertELRI,
                                          spike = post_spikeELRI,
                                          s_to_s = post_s_to_sELRI, 
                                          seed = post_seedELRI, 
                                          data = ELRI_surv_data)
elvi_1m_0v_rec_params <- params_1m_0v_rec(surv = post_survELVI,
                                          grow = post_growELVI, 
                                          flw = post_flwELVI, 
                                          fert = post_fertELVI,
                                          spike = post_spikeELVI,
                                          s_to_s = post_s_to_sELVI, 
                                          seed = post_seedELVI, 
                                          data = ELVI_surv_data)
fesu_1m_0v_rec_params <- params_1m_0v_rec(surv = post_survFESU,
                                          grow = post_growFESU, 
                                          flw = post_flwFESU, 
                                          fert = post_fertFESU,
                                          spike = post_spikeFESU,
                                          s_to_s = post_s_to_sFESU, 
                                          seed = post_seedFESU, 
                                          data = FESU_surv_data)
loar_1m_0v_rec_params <- params_1m_0v_rec(surv = post_survLOAR,
                                          grow = post_growLOAR, 
                                          flw = post_flwLOAR, 
                                          fert = post_fertLOAR,
                                          spike = post_spikeLOAR,
                                          s_to_s = post_s_to_sLOAR, 
                                          seed = post_seedLOAR, 
                                          data = LOAR_surv_data)
poal_1m_0v_rec_params <- params_1m_0v_rec(surv = post_survPOAL,
                                          grow = post_growPOAL, 
                                          flw = post_flwPOAL, 
                                          fert = post_fertPOAL,
                                          spike = post_spikePOAL,
                                          s_to_s = post_s_to_sPOAL, 
                                          seed = post_seedPOAL, 
                                          data = POAL_surv_data)
posy_1m_0v_rec_params <- params_1m_0v_rec(surv = post_survPOSY,
                                          grow = post_growPOSY, 
                                          flw = post_flwPOSY, 
                                          fert = post_fertPOSY,
                                          spike = post_spikePOSY,
                                          s_to_s = post_s_to_sPOSY, 
                                          seed = post_seedPOSY, 
                                          data = POSY_surv_data)





agpe_0m_1v_rec_params <- params_0m_1v_rec(surv = post_survAGPE,
                                          grow = post_growAGPE, 
                                          flw = post_flwAGPE, 
                                          fert = post_fertAGPE,
                                          spike = post_spikeAGPE,
                                          s_to_s = post_s_to_sAGPE, 
                                          seed = post_seedAGPE, 
                                          data = AGPE_surv_data)
elri_0m_1v_rec_params <- params_0m_1v_rec(surv = post_survELRI,
                                          grow = post_growELRI, 
                                          flw = post_flwELRI, 
                                          fert = post_fertELRI,
                                          spike = post_spikeELRI,
                                          s_to_s = post_s_to_sELRI, 
                                          seed = post_seedELRI, 
                                          data = ELRI_surv_data)
elvi_0m_1v_rec_params <- params_0m_1v_rec(surv = post_survELVI,
                                          grow = post_growELVI, 
                                          flw = post_flwELVI, 
                                          fert = post_fertELVI,
                                          spike = post_spikeELVI,
                                          s_to_s = post_s_to_sELVI, 
                                          seed = post_seedELVI, 
                                          data = ELVI_surv_data)
fesu_0m_1v_rec_params <- params_0m_1v_rec(surv = post_survFESU,
                                          grow = post_growFESU, 
                                          flw = post_flwFESU, 
                                          fert = post_fertFESU,
                                          spike = post_spikeFESU,
                                          s_to_s = post_s_to_sFESU, 
                                          seed = post_seedFESU, 
                                          data = FESU_surv_data)
loar_0m_1v_rec_params <- params_0m_1v_rec(surv = post_survLOAR,
                                          grow = post_growLOAR, 
                                          flw = post_flwLOAR, 
                                          fert = post_fertLOAR,
                                          spike = post_spikeLOAR,
                                          s_to_s = post_s_to_sLOAR, 
                                          seed = post_seedLOAR, 
                                          data = LOAR_surv_data)
poal_0m_1v_rec_params <- params_0m_1v_rec(surv = post_survPOAL,
                                          grow = post_growPOAL, 
                                          flw = post_flwPOAL, 
                                          fert = post_fertPOAL,
                                          spike = post_spikePOAL,
                                          s_to_s = post_s_to_sPOAL, 
                                          seed = post_seedPOAL, 
                                          data = POAL_surv_data)
posy_0m_1v_rec_params <- params_0m_1v_rec(surv = post_survPOSY,
                                          grow = post_growPOSY, 
                                          flw = post_flwPOSY, 
                                          fert = post_fertPOSY,
                                          spike = post_spikePOSY,
                                          s_to_s = post_s_to_sPOSY, 
                                          seed = post_seedPOSY, 
                                          data = POSY_surv_data)



agpe_1m_1v_rec_params <- params_1m_1v_rec(surv = post_survAGPE,
                                          grow = post_growAGPE, 
                                          flw = post_flwAGPE, 
                                          fert = post_fertAGPE,
                                          spike = post_spikeAGPE,
                                          s_to_s = post_s_to_sAGPE, 
                                          seed = post_seedAGPE, 
                                          data = AGPE_surv_data)
elri_1m_1v_rec_params <- params_1m_1v_rec(surv = post_survELRI,
                                          grow = post_growELRI, 
                                          flw = post_flwELRI, 
                                          fert = post_fertELRI,
                                          spike = post_spikeELRI,
                                          s_to_s = post_s_to_sELRI, 
                                          seed = post_seedELRI, 
                                          data = ELRI_surv_data)
elvi_1m_1v_rec_params <- params_1m_1v_rec(surv = post_survELVI,
                                          grow = post_growELVI, 
                                          flw = post_flwELVI, 
                                          fert = post_fertELVI,
                                          spike = post_spikeELVI,
                                          s_to_s = post_s_to_sELVI, 
                                          seed = post_seedELVI, 
                                          data = ELVI_surv_data)
fesu_1m_1v_rec_params <- params_1m_1v_rec(surv = post_survFESU,
                                          grow = post_growFESU, 
                                          flw = post_flwFESU, 
                                          fert = post_fertFESU,
                                          spike = post_spikeFESU,
                                          s_to_s = post_s_to_sFESU, 
                                          seed = post_seedFESU, 
                                          data = FESU_surv_data)
loar_1m_1v_rec_params <- params_1m_1v_rec(surv = post_survLOAR,
                                          grow = post_growLOAR, 
                                          flw = post_flwLOAR, 
                                          fert = post_fertLOAR,
                                          spike = post_spikeLOAR,
                                          s_to_s = post_s_to_sLOAR, 
                                          seed = post_seedLOAR, 
                                          data = LOAR_surv_data)
poal_1m_1v_rec_params <- params_1m_1v_rec(surv = post_survPOAL,
                                          grow = post_growPOAL, 
                                          flw = post_flwPOAL, 
                                          fert = post_fertPOAL,
                                          spike = post_spikePOAL,
                                          s_to_s = post_s_to_sPOAL, 
                                          seed = post_seedPOAL, 
                                          data = POAL_surv_data)
posy_1m_1v_rec_params <- params_1m_1v_rec(surv = post_survPOSY,
                                          grow = post_growPOSY, 
                                          flw = post_flwPOSY, 
                                          fert = post_fertPOSY,
                                          spike = post_spikePOSY,
                                          s_to_s = post_s_to_sPOSY, 
                                          seed = post_seedPOSY, 
                                          data = POSY_surv_data)

# 
# 
# getparams <- function(surv, grow, flw, fert, spike, seed, s_to_s, data){
# params <- c()
# params[1] <- lapply(rstan::extract(surv, pars = "beta[1]"), FUN = mean); names(params)[1] <- "surv_beta1"
# params[2] <- lapply(rstan::extract(surv, pars = "beta[2]"), FUN = mean); names(params)[2] <- "surv_beta2"
# params[3] <- lapply(rstan::extract(surv, pars = "beta[3]"), FUN = mean); names(params)[3] <- "surv_beta3"
# params[4] <- lapply(rstan::extract(surv, pars = "beta[4]"), FUN = mean); names(params)[4] <- "surv_beta4"
# # collect the growth parameters
# params[5] <- lapply(rstan::extract(grow, pars = "beta[1]"), FUN = mean); names(params)[5] <- "grow_beta1"
# params[6] <- lapply(rstan::extract(grow, pars = "beta[2]"), FUN = mean); names(params)[6] <- "grow_beta2"
# params[7] <- lapply(rstan::extract(grow, pars = "beta[3]"), FUN = mean); names(params)[7] <- "grow_beta3"
# params[8] <- lapply(rstan::extract(grow, pars = "beta[4]"), FUN = mean); names(params)[8] <- "grow_beta4"
# params[9] <- lapply(rstan::extract(grow, pars = "phi"), FUN = mean); names(params)[9] <- "grow_phi"
# # collect the flowering parameters
# params[10] <- lapply(rstan::extract(flw, pars = "beta[1]"), FUN = mean); names(params)[10] <- "flw_beta1"
# params[11] <- lapply(rstan::extract(flw, pars = "beta[2]"), FUN = mean); names(params)[11] <- "flw_beta2"
# params[12] <- lapply(rstan::extract(flw, pars = "beta[3]"), FUN = mean); names(params)[12] <- "flw_beta3"
# params[13] <- lapply(rstan::extract(flw, pars = "beta[4]"), FUN = mean); names(params)[13] <- "flw_beta4"
# # collect the fertility parameters
# params[14] <- lapply(rstan::extract(fert, pars = "beta[1]"), FUN = mean); names(params)[14] <- "fert_beta1"
# params[15] <- lapply(rstan::extract(fert, pars = "beta[2]"), FUN = mean); names(params)[15] <- "fert_beta2"
# params[16] <- lapply(rstan::extract(fert, pars = "beta[3]"), FUN = mean); names(params)[16] <- "fert_beta3"
# params[17] <- lapply(rstan::extract(fert, pars = "beta[4]"), FUN = mean); names(params)[17] <- "fert_beta4"
# params[18] <- lapply(rstan::extract(fert, pars = "phi"), FUN = mean); names(params)[18] <- "fert_phi"
# # collect the spikelets/infl parameters
# params[19] <- lapply(rstan::extract(spike, pars = "beta[1]"), FUN = mean); names(params)[19] <- "spike_beta1"
# params[20] <- lapply(rstan::extract(spike, pars = "beta[2]"), FUN = mean); names(params)[20] <- "spike_beta2"
# params[21] <- lapply(rstan::extract(spike, pars = "beta[3]"), FUN = mean); names(params)[21] <- "spike_beta3"
# params[22] <- lapply(rstan::extract(spike, pars = "beta[4]"), FUN = mean); names(params)[22] <- "spike_beta4"
# params[23] <- lapply(rstan::extract(spike, pars = "phi"), FUN = mean); names(params)[23] <- "spike_phi"
# # collect the seed/spikelet parameters
# params[24] <- lapply(rstan::extract(seed, pars = "mu_seed"), FUN = mean); names(params)[24] <- "mu_seed"
# # collect the recruitment parameters
# params[25] <- lapply(rstan::extract(s_to_s, pars = "beta[1]"), FUN = mean); names(params)[25] <- "s_to_s_beta1"
# params[26] <- lapply(rstan::extract(s_to_s, pars = "beta[2]"), FUN = mean); names(params)[26] <- "s_to_s_beta2"
# # min and max size
# params[27] <- 1; names(params)[27]<-"min_size"
# params[28] <- quantile(data$size_t1,0.95,na.rm=T); names(params)[28]<-"max_size"
# # note that I define max size as the 95TH pctile of observed sized. The very max sizes observed often have very poor 
# # replication, and I find that these few observations (and the corresponding vital rate predictions) can have a strong
# # influence on the results. So this approach is conservative, but you can always experiment with this.
# # collect the endo specific year parameters
# # surv
# params[29] <- lapply(rstan::extract(surv, pars = "tau_year[1,1]"), FUN = mean); names(params)[29] <- "surv_eminus_y1"
# params[30] <- lapply(rstan::extract(surv, pars = "tau_year[1,2]"), FUN = mean); names(params)[30] <- "surv_eminus_y2"
# params[31] <- lapply(rstan::extract(surv, pars = "tau_year[1,3]"), FUN = mean); names(params)[31] <- "surv_eminus_y3"
# params[32] <- lapply(rstan::extract(surv, pars = "tau_year[1,4]"), FUN = mean); names(params)[32] <- "surv_eminus_y4"
# params[33] <- lapply(rstan::extract(surv, pars = "tau_year[1,5]"), FUN = mean); names(params)[33] <- "surv_eminus_y5"
# params[34] <- lapply(rstan::extract(surv, pars = "tau_year[1,6]"), FUN = mean); names(params)[34] <- "surv_eminus_y6"
# params[35] <- lapply(rstan::extract(surv, pars = "tau_year[1,7]"), FUN = mean); names(params)[35] <- "surv_eminus_y7"
# params[36] <- lapply(rstan::extract(surv, pars = "tau_year[1,8]"), FUN = mean); names(params)[36] <- "surv_eminus_y8"
# params[37] <- lapply(rstan::extract(surv, pars = "tau_year[1,9]"), FUN = mean); names(params)[37] <- "surv_eminus_y9"
# params[38] <- lapply(rstan::extract(surv, pars = "tau_year[1,10]"), FUN = mean); names(params)[38] <- "surv_eminus_y10"
# params[39] <- lapply(rstan::extract(surv, pars = "tau_year[1,11]"), FUN = mean); names(params)[39] <- "surv_eminus_y11"
# 
# params[40] <- lapply(rstan::extract(surv, pars = "tau_year[2,1]"), FUN = mean); names(params)[40] <- "surv_eplus_y1"
# params[41] <- lapply(rstan::extract(surv, pars = "tau_year[2,2]"), FUN = mean); names(params)[41] <- "surv_eplus_y2"
# params[42] <- lapply(rstan::extract(surv, pars = "tau_year[2,3]"), FUN = mean); names(params)[42] <- "surv_eplus_y3"
# params[43] <- lapply(rstan::extract(surv, pars = "tau_year[2,4]"), FUN = mean); names(params)[43] <- "surv_eplus_y4"
# params[44] <- lapply(rstan::extract(surv, pars = "tau_year[2,5]"), FUN = mean); names(params)[44] <- "surv_eplus_y5"
# params[45] <- lapply(rstan::extract(surv, pars = "tau_year[2,6]"), FUN = mean); names(params)[45] <- "surv_eplus_y6"
# params[46] <- lapply(rstan::extract(surv, pars = "tau_year[2,7]"), FUN = mean); names(params)[46] <- "surv_eplus_y7"
# params[47] <- lapply(rstan::extract(surv, pars = "tau_year[2,8]"), FUN = mean); names(params)[47] <- "surv_eplus_y8"
# params[48] <- lapply(rstan::extract(surv, pars = "tau_year[2,9]"), FUN = mean); names(params)[48] <- "surv_eplus_y9"
# params[49] <- lapply(rstan::extract(surv, pars = "tau_year[2,10]"), FUN = mean); names(params)[49] <- "surv_eplus_y10"
# params[50] <- lapply(rstan::extract(surv, pars = "tau_year[2,11]"), FUN = mean); names(params)[50] <- "surv_eplus_y11"
# 
# # grow
# params[51] <- lapply(rstan::extract(grow, pars = "tau_year[1,1]"), FUN = mean); names(params)[51] <- "grow_eminus_y1"
# params[52] <- lapply(rstan::extract(grow, pars = "tau_year[1,2]"), FUN = mean); names(params)[52] <- "grow_eminus_y2"
# params[53] <- lapply(rstan::extract(grow, pars = "tau_year[1,3]"), FUN = mean); names(params)[53] <- "grow_eminus_y3"
# params[54] <- lapply(rstan::extract(grow, pars = "tau_year[1,4]"), FUN = mean); names(params)[54] <- "grow_eminus_y4"
# params[55] <- lapply(rstan::extract(grow, pars = "tau_year[1,5]"), FUN = mean); names(params)[55] <- "grow_eminus_y5"
# params[56] <- lapply(rstan::extract(grow, pars = "tau_year[1,6]"), FUN = mean); names(params)[56] <- "grow_eminus_y6"
# params[57] <- lapply(rstan::extract(grow, pars = "tau_year[1,7]"), FUN = mean); names(params)[57] <- "grow_eminus_y7"
# params[58] <- lapply(rstan::extract(grow, pars = "tau_year[1,8]"), FUN = mean); names(params)[58] <- "grow_eminus_y8"
# params[59] <- lapply(rstan::extract(grow, pars = "tau_year[1,9]"), FUN = mean); names(params)[59] <- "grow_eminus_y9"
# params[60] <- lapply(rstan::extract(grow, pars = "tau_year[1,10]"), FUN = mean); names(params)[60] <- "grow_eminus_y10"
# params[61] <- lapply(rstan::extract(grow, pars = "tau_year[1,11]"), FUN = mean); names(params)[61] <- "grow_eminus_y11"
# 
# params[62] <- lapply(rstan::extract(grow, pars = "tau_year[2,1]"), FUN = mean); names(params)[62] <- "grow_eplus_y1"
# params[63] <- lapply(rstan::extract(grow, pars = "tau_year[2,2]"), FUN = mean); names(params)[63] <- "grow_eplus_y2"
# params[64] <- lapply(rstan::extract(grow, pars = "tau_year[2,3]"), FUN = mean); names(params)[64] <- "grow_eplus_y3"
# params[65] <- lapply(rstan::extract(grow, pars = "tau_year[2,4]"), FUN = mean); names(params)[65] <- "grow_eplus_y4"
# params[66] <- lapply(rstan::extract(grow, pars = "tau_year[2,5]"), FUN = mean); names(params)[66] <- "grow_eplus_y5"
# params[67] <- lapply(rstan::extract(grow, pars = "tau_year[2,6]"), FUN = mean); names(params)[67] <- "grow_eplus_y6"
# params[68] <- lapply(rstan::extract(grow, pars = "tau_year[2,7]"), FUN = mean); names(params)[68] <- "grow_eplus_y7"
# params[69] <- lapply(rstan::extract(grow, pars = "tau_year[2,8]"), FUN = mean); names(params)[69] <- "grow_eplus_y8"
# params[70] <- lapply(rstan::extract(grow, pars = "tau_year[2,9]"), FUN = mean); names(params)[70] <- "grow_eplus_y9"
# params[71] <- lapply(rstan::extract(grow, pars = "tau_year[2,10]"), FUN = mean); names(params)[71] <- "grow_eplus_y10"
# params[72] <- lapply(rstan::extract(grow, pars = "tau_year[2,11]"), FUN = mean); names(params)[72] <- "grow_eplus_y11"
# 
# # flw
# params[73] <- lapply(rstan::extract(flw, pars = "tau_year[1,1]"), FUN = mean); names(params)[73] <- "flw_eminus_y1"
# params[74] <- lapply(rstan::extract(flw, pars = "tau_year[1,2]"), FUN = mean); names(params)[74] <- "flw_eminus_y2"
# params[75] <- lapply(rstan::extract(flw, pars = "tau_year[1,3]"), FUN = mean); names(params)[75] <- "flw_eminus_y3"
# params[76] <- lapply(rstan::extract(flw, pars = "tau_year[1,4]"), FUN = mean); names(params)[76] <- "flw_eminus_y4"
# params[77] <- lapply(rstan::extract(flw, pars = "tau_year[1,5]"), FUN = mean); names(params)[77] <- "flw_eminus_y5"
# params[78] <- lapply(rstan::extract(flw, pars = "tau_year[1,6]"), FUN = mean); names(params)[78] <- "flw_eminus_y6"
# params[79] <- lapply(rstan::extract(flw, pars = "tau_year[1,7]"), FUN = mean); names(params)[79] <- "flw_eminus_y7"
# params[80] <- lapply(rstan::extract(flw, pars = "tau_year[1,8]"), FUN = mean); names(params)[80] <- "flw_eminus_y8"
# params[81] <- lapply(rstan::extract(flw, pars = "tau_year[1,9]"), FUN = mean); names(params)[81] <- "flw_eminus_y9"
# params[82] <- lapply(rstan::extract(flw, pars = "tau_year[1,10]"), FUN = mean); names(params)[82] <- "flw_eminus_y10"
# params[83] <- lapply(rstan::extract(flw, pars = "tau_year[1,11]"), FUN = mean); names(params)[83] <- "flw_eminus_y11"
# 
# params[84] <- lapply(rstan::extract(flw, pars = "tau_year[2,1]"), FUN = mean); names(params)[84] <- "flw_eplus_y1"
# params[85] <- lapply(rstan::extract(flw, pars = "tau_year[2,2]"), FUN = mean); names(params)[85] <- "flw_eplus_y2"
# params[86] <- lapply(rstan::extract(flw, pars = "tau_year[2,3]"), FUN = mean); names(params)[86] <- "flw_eplus_y3"
# params[87] <- lapply(rstan::extract(flw, pars = "tau_year[2,4]"), FUN = mean); names(params)[87] <- "flw_eplus_y4"
# params[88] <- lapply(rstan::extract(flw, pars = "tau_year[2,5]"), FUN = mean); names(params)[88] <- "flw_eplus_y5"
# params[89] <- lapply(rstan::extract(flw, pars = "tau_year[2,6]"), FUN = mean); names(params)[89] <- "flw_eplus_y6"
# params[90] <- lapply(rstan::extract(flw, pars = "tau_year[2,7]"), FUN = mean); names(params)[90] <- "flw_eplus_y7"
# params[91] <- lapply(rstan::extract(flw, pars = "tau_year[2,8]"), FUN = mean); names(params)[91] <- "flw_eplus_y8"
# params[92] <- lapply(rstan::extract(flw, pars = "tau_year[2,9]"), FUN = mean); names(params)[92] <- "flw_eplus_y9"
# params[93] <- lapply(rstan::extract(flw, pars = "tau_year[2,10]"), FUN = mean); names(params)[93] <- "flw_eplus_y10"
# params[94] <- lapply(rstan::extract(flw, pars = "tau_year[2,11]"), FUN = mean); names(params)[94] <- "flw_eplus_y11"
# 
# # fert
# params[95] <- lapply(rstan::extract(fert, pars = "tau_year[1,1]"), FUN = mean); names(params)[95] <- "fert_eminus_y1"
# params[96] <- lapply(rstan::extract(fert, pars = "tau_year[1,2]"), FUN = mean); names(params)[96] <- "fert_eminus_y2"
# params[97] <- lapply(rstan::extract(fert, pars = "tau_year[1,3]"), FUN = mean); names(params)[97] <- "fert_eminus_y3"
# params[98] <- lapply(rstan::extract(fert, pars = "tau_year[1,4]"), FUN = mean); names(params)[98] <- "fert_eminus_y4"
# params[99] <- lapply(rstan::extract(fert, pars = "tau_year[1,5]"), FUN = mean); names(params)[99] <- "fert_eminus_y5"
# params[100] <- lapply(rstan::extract(fert, pars = "tau_year[1,6]"), FUN = mean); names(params)[100] <- "fert_eminus_y6"
# params[101] <- lapply(rstan::extract(fert, pars = "tau_year[1,7]"), FUN = mean); names(params)[101] <- "fert_eminus_y7"
# params[102] <- lapply(rstan::extract(fert, pars = "tau_year[1,8]"), FUN = mean); names(params)[102] <- "fert_eminus_y8"
# params[103] <- lapply(rstan::extract(fert, pars = "tau_year[1,9]"), FUN = mean); names(params)[103] <- "fert_eminus_y9"
# params[104] <- lapply(rstan::extract(fert, pars = "tau_year[1,10]"), FUN = mean); names(params)[104] <- "fert_eminus_y10"
# params[105] <- lapply(rstan::extract(fert, pars = "tau_year[1,11]"), FUN = mean); names(params)[105] <- "fert_eminus_y11"
# 
# params[106] <- lapply(rstan::extract(fert, pars = "tau_year[2,1]"), FUN = mean); names(params)[106] <- "fert_eplus_y1"
# params[107] <- lapply(rstan::extract(fert, pars = "tau_year[2,2]"), FUN = mean); names(params)[107] <- "fert_eplus_y2"
# params[108] <- lapply(rstan::extract(fert, pars = "tau_year[2,3]"), FUN = mean); names(params)[108] <- "fert_eplus_y3"
# params[109] <- lapply(rstan::extract(fert, pars = "tau_year[2,4]"), FUN = mean); names(params)[109] <- "fert_eplus_y4"
# params[110] <- lapply(rstan::extract(fert, pars = "tau_year[2,5]"), FUN = mean); names(params)[110] <- "fert_eplus_y5"
# params[111] <- lapply(rstan::extract(fert, pars = "tau_year[2,6]"), FUN = mean); names(params)[111] <- "fert_eplus_y6"
# params[112] <- lapply(rstan::extract(fert, pars = "tau_year[2,7]"), FUN = mean); names(params)[112] <- "fert_eplus_y7"
# params[113] <- lapply(rstan::extract(fert, pars = "tau_year[2,8]"), FUN = mean); names(params)[113] <- "fert_eplus_y8"
# params[114] <- lapply(rstan::extract(fert, pars = "tau_year[2,9]"), FUN = mean); names(params)[114] <- "fert_eplus_y9"
# params[115] <- lapply(rstan::extract(fert, pars = "tau_year[2,10]"), FUN = mean); names(params)[115] <- "fert_eplus_y10"
# params[116] <- lapply(rstan::extract(fert, pars = "tau_year[2,11]"), FUN = mean); names(params)[116] <- "fert_eplus_y11"
# 
# # spike
# params[117] <- lapply(rstan::extract(spike, pars = "tau_year[1,1]"), FUN = mean); names(params)[117] <- "spike_eminus_y1"
# params[118] <- lapply(rstan::extract(spike, pars = "tau_year[1,2]"), FUN = mean); names(params)[118] <- "spike_eminus_y2"
# params[119] <- lapply(rstan::extract(spike, pars = "tau_year[1,3]"), FUN = mean); names(params)[119] <- "spike_eminus_y3"
# params[120] <- lapply(rstan::extract(spike, pars = "tau_year[1,4]"), FUN = mean); names(params)[120] <- "spike_eminus_y4"
# params[121] <- lapply(rstan::extract(spike, pars = "tau_year[1,5]"), FUN = mean); names(params)[121] <- "spike_eminus_y5"
# params[122] <- lapply(rstan::extract(spike, pars = "tau_year[1,6]"), FUN = mean); names(params)[122] <- "spike_eminus_y6"
# params[123] <- lapply(rstan::extract(spike, pars = "tau_year[1,7]"), FUN = mean); names(params)[123] <- "spike_eminus_y7"
# params[124] <- lapply(rstan::extract(spike, pars = "tau_year[1,8]"), FUN = mean); names(params)[124] <- "spike_eminus_y8"
# params[125] <- lapply(rstan::extract(spike, pars = "tau_year[1,9]"), FUN = mean); names(params)[125] <- "spike_eminus_y9"
# params[126] <- lapply(rstan::extract(spike, pars = "tau_year[1,10]"), FUN = mean); names(params)[126] <- "spike_eminus_y10"
# params[127] <- lapply(rstan::extract(spike, pars = "tau_year[1,11]"), FUN = mean); names(params)[127] <- "spike_eminus_y11"
# 
# params[128] <- lapply(rstan::extract(spike, pars = "tau_year[2,1]"), FUN = mean); names(params)[128] <- "spike_eplus_y1"
# params[129] <- lapply(rstan::extract(spike, pars = "tau_year[2,2]"), FUN = mean); names(params)[129] <- "spike_eplus_y2"
# params[130] <- lapply(rstan::extract(spike, pars = "tau_year[2,3]"), FUN = mean); names(params)[130] <- "spike_eplus_y3"
# params[131] <- lapply(rstan::extract(spike, pars = "tau_year[2,4]"), FUN = mean); names(params)[131] <- "spike_eplus_y4"
# params[132] <- lapply(rstan::extract(spike, pars = "tau_year[2,5]"), FUN = mean); names(params)[132] <- "spike_eplus_y5"
# params[133] <- lapply(rstan::extract(spike, pars = "tau_year[2,6]"), FUN = mean); names(params)[133] <- "spike_eplus_y6"
# params[134] <- lapply(rstan::extract(spike, pars = "tau_year[2,7]"), FUN = mean); names(params)[134] <- "spike_eplus_y7"
# params[135] <- lapply(rstan::extract(spike, pars = "tau_year[2,8]"), FUN = mean); names(params)[135] <- "spike_eplus_y8"
# params[136] <- lapply(rstan::extract(spike, pars = "tau_year[2,9]"), FUN = mean); names(params)[136] <- "spike_eplus_y9"
# params[137] <- lapply(rstan::extract(spike, pars = "tau_year[2,10]"), FUN = mean); names(params)[137] <- "spike_eplus_y10"
# params[138] <- lapply(rstan::extract(spike, pars = "tau_year[2,11]"), FUN = mean); names(params)[138] <- "spike_eplus_y11"
# 
# # s_to_s
# params[139] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,1]"), FUN = mean); names(params)[139] <- "s_to_s_eminus_y1"
# params[140] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,2]"), FUN = mean); names(params)[140] <- "s_to_s_eminus_y2"
# params[141] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,3]"), FUN = mean); names(params)[141] <- "s_to_s_eminus_y3"
# params[142] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,4]"), FUN = mean); names(params)[142] <- "s_to_s_eminus_y4"
# params[143] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,5]"), FUN = mean); names(params)[143] <- "s_to_s_eminus_y5"
# params[144] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,6]"), FUN = mean); names(params)[144] <- "s_to_s_eminus_y6"
# params[145] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,7]"), FUN = mean); names(params)[145] <- "s_to_s_eminus_y7"
# params[146] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,8]"), FUN = mean); names(params)[146] <- "s_to_s_eminus_y8"
# params[147] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,9]"), FUN = mean); names(params)[147] <- "s_to_s_eminus_y9"
# params[148] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,10]"), FUN = mean); names(params)[148] <- "s_to_s_eminus_y10"
# params[149] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,11]"), FUN = mean); names(params)[149] <- "s_to_s_eminus_y11"
# 
# params[150] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,1]"), FUN = mean); names(params)[150] <- "s_to_s_eplus_y1"
# params[151] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,2]"), FUN = mean); names(params)[151] <- "s_to_s_eplus_y2"
# params[152] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,3]"), FUN = mean); names(params)[152] <- "s_to_s_eplus_y3"
# params[153] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,4]"), FUN = mean); names(params)[153] <- "s_to_s_eplus_y4"
# params[154] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,5]"), FUN = mean); names(params)[154] <- "s_to_s_eplus_y5"
# params[155] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,6]"), FUN = mean); names(params)[155] <- "s_to_s_eplus_y6"
# params[156] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,7]"), FUN = mean); names(params)[156] <- "s_to_s_eplus_y7"
# params[157] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,8]"), FUN = mean); names(params)[157] <- "s_to_s_eplus_y8"
# params[158] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,9]"), FUN = mean); names(params)[158] <- "s_to_s_eplus_y9"
# params[159] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,10]"), FUN = mean); names(params)[159] <- "s_to_s_eplus_y10"
# params[160] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,11]"), FUN = mean); names(params)[160] <- "s_to_s_eplus_y11"
# 
# params <- unlist(params)
# return(params)
# }
# 
# agpe_params <- getparams(surv = survAGPE, grow = growAGPE, flw = flwAGPE, fert = fertAGPE, spike = spikeAGPE, seed = seedAGPE, s_to_s = s_to_sAGPE, data = AGPE_surv_data)
# elri_params <- getparams(surv = survELRI, grow = growELRI, flw = flwELRI, fert = fertELRI, spike = spikeELRI, seed = seedELRI, s_to_s = s_to_sELRI, data = ELRI_surv_data)
# elvi_params <- getparams(surv = survELVI, grow = growELVI, flw = flwELVI, fert = fertELVI, spike = spikeELVI, seed = seedELVI, s_to_s = s_to_sELVI, data = ELVI_surv_data)
# fesu_params <- getparams(surv = survFESU, grow = growFESU, flw = flwFESU, fert = fertFESU, spike = spikeFESU, seed = seedFESU, s_to_s = s_to_sFESU, data = FESU_surv_data)
# loar_params <- getparams(surv = survLOAR, grow = growLOAR, flw = flwLOAR, fert = fertLOAR, spike = spikeLOAR, seed = seedLOAR, s_to_s = s_to_sLOAR, data = LOAR_surv_data)
# poal_params <- getparams(surv = survPOAL, grow = growPOAL, flw = flwPOAL, fert = fertPOAL, spike = spikePOAL, seed = seedPOAL, s_to_s = s_to_sPOAL, data = POAL_surv_data)
# posy_params <- getparams(surv = survPOSY, grow = growPOSY, flw = flwPOSY, fert = fertPOSY, spike = spikePOSY, seed = seedPOSY, s_to_s = s_to_sPOSY, data = POSY_surv_data)



# define functions that will be used to populate projection matrix
#SURVIVAL AT SIZE X.
# currently this is fitting as if the E- is the intercept
sx<-function(x,params, year){
  invlogit(params[paste0("surv_intercept_y", year)] + params["surv_slope"]*log(x))
}
sx(3, agpe_1m_1v_rec_params, year = 1:11)
#PROBABILITY OF GROWTH FROM SIZE X TO Y
#This function truncates the density asscociation with x==0 and x>x.max
gxy<-function(x,y,params, year){
  grow.mean<-params[paste0("grow_intercept_y", year)] + params["grow_slope"]*log(x)
  pr_grow<-dnbinom(x=y,mu=exp(grow.mean), size = params["grow_phi"])
  truncLower<-dnbinom(x=0,mu=exp(grow.mean),size = params["grow_phi"])
  truncUpper<-sum(dnbinom(x=(params["max_size"]+1):500, mu=exp(grow.mean),size = params["grow_phi"]))
  return(pr_grow/(1-(truncLower+truncUpper)))
}
# gxy(3,4, agpe_1m_1v_rec_params, year = 1)

#SURVIVAL*GROWTH
# This function generates the transition matrix by multiplying growth and survival.
pxy<-function(x,y,params, year){
  sx(x,params,year) * gxy(x,y,params, year)
}
# pxy(x = 1,y = 4, params = agpe_0m_0v_rec_params, year = 1)

#FERTILITY--returns number of seedlings, which we will assume (for now) to be 1-tiller, produced by size X
fx<-function(x,params, year){
  flower <- invlogit(params[paste0("flw_intercept_y", year)] + params["flw_slope"]*log(x))
  fert <- exp(params[paste0("fert_intercept_y", year)] + params["fert_slope"]*log(x))
  spike <- exp(params[paste0("spike_intercept_y", year)] + params["spike_slope"]*log(x))
  seed.mean <- params["mu_seed"]
  recruitment <- invlogit(params[paste0("s_to_s_intercept_y", year)])
  seedlings <- flower * fert * spike * seed.mean * recruitment
  return(seedlings)
}
# fx(3, agpe_1m_1v_rec_params, year = 7)


## note from Tom: we should think about adding a reproductive delay
# finally, here is the function that takes the parameter vector and assembles the matrix model from all of the pieces
bigmatrix<-function(params, year){   
  matdim<-params["max_size"]## matrix dimension
  y <- 1:params["max_size"]## size (tiller number) associated with each class
  # Fertility matrix
  Fmat <- matrix(0,matdim,matdim)
  # all seedlings get dumped into top row (1-tiller)
  Fmat[1,]<-fx(x = y, params,year)
  
  # Growth/survival transition matrix
  Tmat <-matrix(0,matdim,matdim)
  # Filling the transition matrix 
  Tmat<-t(outer(y,y,pxy,params, year))
  
  # Put it all together
  # sum the Tmat & Fmat to get the whole matrix
  MPMmat<-Tmat + Fmat
  return(list(MPMmat = MPMmat, Tmat = Tmat, Fmat = Fmat))
}
lambda(bigmatrix(loar_0m_1v_rec_params, year = 2)$MPMmat)




## population growth rate (eigenalaysis of the projection matrix)
# matrix <- bigmatrix(loar_params)
# calculating the matrices
agpe_mat <- bigmatrix(agpe_params)
elri_mat <- bigmatrix(elri_params)
elvi_mat <- bigmatrix(elvi_params)
fesu_mat <- bigmatrix(fesu_params)
loar_mat <- bigmatrix(loar_params)
poal_mat <- bigmatrix(poal_params)
posy_mat <- bigmatrix(posy_params)

# saving the yearly growth rates to vectors 
yearly <- function(params, nyear = 11){
  yearly_eminus <- c()
  yearly_eplus <- c()
  yearly_eminus_rec <- c()
  yearly_eplus_rec <- c()
  year_t1 <- c(2008:2018)
  # eminus original
  for(i in 1:nyear){
    # eminus original
    yearly_eminus[i] <- lambda(bigmatrix(params = params, endo = FALSE, recruit = FALSE, year = i)$MPMmat)
    # eplus original
    yearly_eplus[i] <- lambda(bigmatrix(params = params, endo = TRUE, recruit = FALSE, year = i)$MPMmat)
    # eminus recruit
    yearly_eminus_rec[i] <- lambda(bigmatrix(params = params, endo = FALSE, recruit = TRUE, year = i)$MPMmat)
    # eplus recruit
    yearly_eplus_rec[i] <- lambda(bigmatrix(params = params, endo = TRUE, recruit = TRUE, year = i)$MPMmat)
    
  }

  lambdas <- as_tibble(cbind(yearly_eminus, yearly_eplus, yearly_eminus_rec, yearly_eplus_rec, year_t1))

  return(lambdas)
}
yearly(fesu_params)


# some starter histograms
# eminus = orange, eplus = purple
ggplot(data = yearly(agpe_params))+
  geom_histogram(aes(yearly_eminus_rec), bins = 8, fill = "#ff7f00", alpha = .6) +
  geom_histogram(aes(yearly_eplus_rec), bins = 8, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") +labs(title = "AGPE") + theme_classic(base_size = 20)

ggplot(data = yearly(elri_params))+
  geom_histogram(aes(yearly_eminus_rec), bins = 8, fill = "#ff7f00", alpha = .6) +
  geom_histogram(aes(yearly_eplus_rec), bins = 8, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") +labs(title = "ELRI") + theme_classic(base_size = 20)

ggplot(data = yearly(elvi_params))+
  geom_histogram(aes(yearly_eminus_rec), bins = 8, fill = "#ff7f00", alpha = .6) +
  geom_histogram(aes(yearly_eplus_rec), bins = 8, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") +labs(title = "ELVI") + theme_classic(base_size = 20)

ggplot(data = yearly(fesu_params))+
  geom_histogram(aes(yearly_eminus_rec), bins = 8, fill = "#ff7f00", alpha = .6) +
  geom_histogram(aes(yearly_eplus_rec), bins = 8, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") +labs(title = "FESU") + theme_classic(base_size = 20)

ggplot(data = yearly(loar_params))+
  geom_histogram(aes(yearly_eminus_rec), bins = 8, fill = "#ff7f00", alpha = .6) +
  geom_histogram(aes(yearly_eplus_rec), bins = 8, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") +labs(title = "LOAR") + theme_classic(base_size = 20)

ggplot(data = yearly(poal_params))+
  geom_histogram(aes(yearly_eminus_rec), bins = 8, fill = "#ff7f00", alpha = .6) +
  geom_histogram(aes(yearly_eplus_rec), bins = 8, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") +labs(title = "POAL") + theme_classic(base_size = 20)

ggplot(data = yearly(posy_params))+
  geom_histogram(aes(yearly_eminus_rec), bins = 8, fill = "#ff7f00", alpha = .6) +
  geom_histogram(aes(yearly_eplus_rec), bins = 8, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") +labs(title = "POSY") + theme_classic(base_size = 20) 

# option of a density plot
ggplot(data = yearly(posy_params))+
  geom_density(aes(yearly_eminus_rec), bw = .1, fill = "#ff7f00", alpha = .6) +
  geom_density(aes(yearly_eplus_rec), bw = .1, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") + labs(title = "POSY") +theme_classic(base_size = 20)

# some starter time series
ggplot(data = yearly(agpe_params))+
  geom_point(aes(x = year_t1, y = yearly_eminus_rec),size = 4, col = "#ff7f00") +
  geom_point(aes(x = year_t1, y = yearly_eplus_rec),size = 4, col = "#6a3d9a") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(y = "population growth rate", x = "year_t1") + ggtitle("AGPE") + theme_classic(base_size = 20) + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))
ggplot(data = yearly(elri_params))+
  geom_point(aes(x = year_t1, y = yearly_eminus_rec),size = 4, col = "#ff7f00") +
  geom_point(aes(x = year_t1, y = yearly_eplus_rec),size = 4, col = "#6a3d9a") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(y = "population growth rate", x = "year_t1") + ggtitle("ELRI") + theme_classic(base_size = 20) + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))
ggplot(data = yearly(elvi_params))+
  geom_point(aes(x = year_t1, y = yearly_eminus_rec),size = 4, col = "#ff7f00") +
  geom_point(aes(x = year_t1, y = yearly_eplus_rec),size = 4, col = "#6a3d9a") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(y = "population growth rate", x = "year_t1") + ggtitle("ELVI") + theme_classic(base_size = 20) + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))

ggplot(data = yearly(fesu_params))+
  geom_point(aes(x = year_t1, y = yearly_eminus_rec),size = 4, col = "#ff7f00") +
  geom_point(aes(x = year_t1, y = yearly_eplus_rec),size = 4, col = "#6a3d9a") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(y = "population growth rate", x = "year_t1") + ggtitle("FESU") + theme_classic(base_size = 20) + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))

ggplot(data = yearly(loar_params))+
  geom_point(aes(x = year_t1, y = yearly_eminus_rec),size = 4, col = "#ff7f00") +
  geom_point(aes(x = year_t1, y = yearly_eplus_rec),size = 4, col = "#6a3d9a") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(y = "population growth rate", x = "year_t1") + ggtitle("LOAR") + theme_classic(base_size = 20) + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))

ggplot(data = yearly(poal_params))+
  geom_point(aes(x = year_t1, y = yearly_eminus_rec),size = 4, col = "#ff7f00") +
  geom_point(aes(x = year_t1, y = yearly_eplus_rec),size = 4, col = "#6a3d9a") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(y = "population growth rate", x = "year_t1") + ggtitle("POAL") + theme_classic(base_size = 20) + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))

ggplot(data = yearly(posy_params))+
  geom_point(aes(x = year_t1, y = yearly_eminus_rec),size = 4, col = "#ff7f00") +
  geom_point(aes(x = year_t1, y = yearly_eplus_rec),size = 4, col = "#6a3d9a") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(y = "population growth rate", x = "year_t1") + ggtitle("POSY") + theme_classic(base_size = 20) + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))



# I'm gonna read in env. data here
# data downloaded from PRISM, daily ppt, tmin, tmean, tmax for GPS point 39.235900000000,-86.218100000000 from Jan 1, 2006 to Jan 1, 2019
climate <- read_csv(file = "~/Dropbox/EndodemogData/PRISMClimateData_BrownCo.csv") %>% 
  mutate(year = year(Date), month = month(Date), day = day(Date)) %>% 
  rename(ppt = `ppt (mm)`, tmean = `tmean (degrees C)`) %>% 
  mutate(site_lat = 39.235900000000, site_long = -86.218100000000)
  

AGPE_climate <- climate %>% 
  mutate(census_month = 9, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>% 
  filter(climate_year != 2006) %>% 
  group_by(climate_year) %>% 
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean))
AGPE <- yearly(agpe_params) %>% 
  mutate(species = "AGPE") %>% 
  merge(AGPE_climate, by.x = c("year_t1"), by.y = c("climate_year"))
ELRI_climate <- climate %>%
  mutate(census_month = 7, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>%
  filter(climate_year != 2006) %>%
  group_by(climate_year) %>%
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean))
ELRI <- yearly(elri_params) %>%
  mutate(species = "ELRI") %>%
  merge(ELRI_climate, by.x = c("year_t1"), by.y = c("climate_year"))
ELVI_climate <- climate %>%
  mutate(census_month = 7, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>%
  filter(climate_year != 2006) %>%
  group_by(climate_year) %>%
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean))
ELVI <- yearly(elvi_params) %>%
  mutate(species = "ELVI") %>%
  merge(ELVI_climate, by.x = c("year_t1"), by.y = c("climate_year"))
FESU_climate <- climate %>% 
  mutate(census_month = 6, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>% 
  filter(climate_year != 2006) %>% 
  group_by(climate_year) %>% 
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean))
FESU <- yearly(fesu_params) %>% 
  mutate(species = "FESU") %>% 
  merge(FESU_climate, by.x = c("year_t1"), by.y = c("climate_year"))
LOAR_climate <- climate %>% 
  mutate(census_month = 7, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>% 
  filter(climate_year != 2006) %>% 
  group_by(climate_year) %>% 
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean))
LOAR <- yearly(loar_params) %>% 
  mutate(species = "LOAR") %>% 
  merge(LOAR_climate, by.x = c("year_t1"), by.y = c("climate_year"))
POAL_climate <- climate %>% 
  mutate(census_month = 5, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>% 
  filter(climate_year != 2006) %>% 
  group_by(climate_year) %>% 
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean))
POAL <- yearly(poal_params) %>% 
  mutate(species = "POAL") %>% 
  merge(POAL_climate, by.x = c("year_t1"), by.y = c("climate_year"))
POSY_climate <- climate %>% 
  mutate(census_month = 5, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>% 
  filter(climate_year != 2006) %>% 
  group_by(climate_year) %>% 
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean))
POSY <- yearly(posy_params) %>% 
  mutate(species = "POSY") %>% 
  merge(POSY_climate, by.x = c("year_t1"), by.y = c("climate_year"))

endo_climate <-  rbind(AGPE,ELRI,ELVI,FESU,LOAR,POAL,POSY) %>% 
  mutate(diff_lambda = yearly_eplus_rec - yearly_eminus_rec) %>% 
  melt(id.vars = c("species", "year_t1", 'Cumulative PPT (mm)', 'Mean Temp. (C˚)'))

# graph of difference in lambda for ppt
endo_climate %>% 
  filter(variable == "diff_lambda") %>% 
  ggplot() +
  geom_smooth(aes(x = `Cumulative PPT (mm)`, y = value),color = "grey", method = "glm", se = FALSE) +
  geom_point(aes(x = `Cumulative PPT (mm)`, y = value, color = variable, shape = species)) +
  facet_grid(cols = vars(species)) +
  theme_classic() + scale_color_manual(values = c("#ff7f00", "#6a3d9a")) +scale_shape_manual(values=c(0:6))+ labs(y = "Population Growth") + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))


endo_climate %>% 
  filter(variable == "diff_lambda") %>% 
  ggplot() +
  geom_smooth(aes(x = `Mean Temp. (C˚)`, y = value),color = "grey", method = "glm", se = FALSE) +
  geom_point(aes(x = `Mean Temp. (C˚)`, y = value, color = variable, shape = species)) +
  facet_grid(cols = vars(species)) +
  theme_classic() + scale_color_manual(values = c("#ff7f00", "#6a3d9a")) +scale_shape_manual(values=c(0:6))+ labs(y = "Population Growth") + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))


# Graphs of climate and ppt through time
endo_climate %>% 
  filter(species == "AGPE") %>% 
  ggplot()+
  geom_point(aes(x = year_t1, y = `Cumulative PPT (mm)`), color = "black") +
  theme_classic() + ylim(0,1400) + scale_color_manual(values = c("#ff7f00", "#6a3d9a"))+ scale_shape_manual(values=c(0:6))+ labs(y = "Cumulative Precipitation (mm)") + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))

endo_climate %>% 
  filter(species == "AGPE") %>% 
  ggplot()+
  geom_point(aes(x = year_t1, y = `Mean Temp. (C˚)`), color = "black") +
  theme_classic() + ylim(0,15) + scale_color_manual(values = c("#ff7f00", "#6a3d9a"))+ scale_shape_manual(values=c(0:6))+ labs(y = "Mean Temp. (C˚)") + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))




endo_climate %>% 
  filter(variable == "yearly_eminus_rec" | variable == "yearly_eplus_rec") %>% 
  ggplot() +
  geom_point(aes(x = year_t1, y = value, color = variable)) +
  facet_grid(rows = vars(species)) +
  theme_classic() + scale_color_manual(values = c("#ff7f00","#6a3d9a"))+ labs(y = "Population Growth") + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))

# Calculate effect size and sd for lambdas

endo <- endo_climate %>% 
  filter( variable == "yearly_eminus_rec" | variable == "yearly_eplus_rec" ) %>% 
  mutate(variable = recode(variable, yearly_eplus_rec = "E+", yearly_eminus_rec = "E-")) %>% 
  group_by(species,variable) %>% 
  summarize(mean_l = mean(value),
            sd_l = sd(value),
            lowCI = mean_l - 2*sd_l,
            highCI = mean_l + 2*sd_l)

# Plot of lambda distributions
endo_climate %>% 
  filter( variable == "yearly_eminus_rec" | variable == "yearly_eplus_rec" ) %>% 
  mutate(variable = recode(variable, yearly_eplus_rec = "E+", yearly_eminus_rec = "E-")) %>%
  ggplot(aes(x = variable, y = value, color = variable))+
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = .3)+
  facet_grid(rows = vars(species)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  coord_flip() + 
  xlab("Endophyte") + ylab("Mean Pop. Growth Rate") +
  theme_classic(base_size = 15) + scale_color_manual(values = c("#ff7f00","#6a3d9a"))+theme(axis.title = element_text(size = 30), axis.text = element_text(size = 20)) 



# plots of normal curves
xseq<-seq(0,2,.01)
AGPE_minus<-dnorm(xseq, mean(AGPE$yearly_eminus_rec),sd(AGPE$yearly_eminus_rec))
AGPE_plus<-dnorm(xseq, mean(AGPE$yearly_eplus_rec),sd(AGPE$yearly_eplus_rec))
ELRI_minus<-dnorm(xseq, mean(ELRI$yearly_eminus_rec),sd(ELRI$yearly_eminus_rec))
ELRI_plus<-dnorm(xseq, mean(ELRI$yearly_eplus_rec),sd(ELRI$yearly_eplus_rec))
ELVI_minus<-dnorm(xseq, mean(ELVI$yearly_eminus_rec),sd(ELVI$yearly_eminus_rec))
ELVI_plus<-dnorm(xseq, mean(ELVI$yearly_eplus_rec),sd(ELVI$yearly_eplus_rec))
FESU_minus<-dnorm(xseq, mean(FESU$yearly_eminus_rec),sd(FESU$yearly_eminus_rec))
FESU_plus<-dnorm(xseq, mean(FESU$yearly_eplus_rec),sd(FESU$yearly_eplus_rec))
LOAR_minus<-dnorm(xseq, mean(LOAR$yearly_eminus_rec),sd(LOAR$yearly_eminus_rec))
LOAR_plus<-dnorm(xseq, mean(LOAR$yearly_eplus_rec),sd(LOAR$yearly_eplus_rec))
POAL_minus<-dnorm(xseq, mean(POAL$yearly_eminus_rec),sd(POAL$yearly_eminus_rec))
POAL_plus<-dnorm(xseq, mean(POAL$yearly_eplus_rec),sd(POAL$yearly_eplus_rec))
POSY_minus<-dnorm(xseq, mean(POSY$yearly_eminus_rec),sd(POSY$yearly_eminus_rec))
POSY_plus<-dnorm(xseq, mean(POSY$yearly_eplus_rec),sd(POSY$yearly_eplus_rec))

par(mfrow=c(1,7), mar=c(3,2,4,2))

plot(xseq, AGPE_plus, col="#6a3d9a",xlab="", ylab = "Density", type="l",lwd=2, cex=2, main="AGPE", cex.axis=.8)
points(xseq, AGPE_minus, col="#ff7f00",xlab="", ylab="", type="l",lwd=2, cex=2, main="AGPE", cex.axis=.8)
plot(xseq, ELRI_plus, col="#6a3d9a",xlab="", ylab="", type="l",lwd=2, cex=2, main="ELRI", cex.axis=.8)
points(xseq, ELRI_minus, col="#ff7f00",xlab="", ylab="", type="l",lwd=2, cex=2, main="ELRI", cex.axis=.8)
plot(xseq, ELVI_minus, col="#ff7f00",xlab="", ylab="", type="l",lwd=2, cex=2, main="ELVI", cex.axis=.8)
points(xseq, ELVI_plus, col="#6a3d9a",xlab="", ylab="", type="l",lwd=2, cex=2, main="ELVI", cex.axis=.8)
plot(xseq, FESU_plus, col="#6a3d9a",xlab="", ylab="", type="l",lwd=2, cex=2, main="FESU", cex.axis=.8)
points(xseq, FESU_minus, col="#ff7f00",xlab="", ylab="", type="l",lwd=2, cex=2, main="FESU", cex.axis=.8)
plot(xseq, LOAR_plus, col="#6a3d9a",xlab="", ylab="", type="l",lwd=2, cex=2, main="LOAR", cex.axis=.8)
points(xseq, LOAR_minus, col="#ff7f00",xlab="", ylab="", type="l",lwd=2, cex=2, main="LOAR", cex.axis=.8)
plot(xseq, POAL_plus, col="#6a3d9a",xlab="", ylab="", type="l",lwd=2, cex=2, main="POAL", cex.axis=.8)
points(xseq, POAL_minus, col="#ff7f00",xlab="", ylab="", type="l",lwd=2, cex=2, main="POAL", cex.axis=.8)
plot(xseq, POSY_plus, col="#6a3d9a",xlab="", ylab="", type="l",lwd=2, cex=2, main="POSY", cex.axis=.8)
points(xseq, POSY_minus, col="#ff7f00",xlab="", ylab="", type="l",lwd=2, cex=2, main="POSY", cex.axis=.8)


# some mean and sd of the lambda values,
mean(yearly(agpe_mat)$yearly_eminus)
mean(yearly(agpe_mat)$yearly_eplus)
sd(yearly(agpe_mat)$yearly_eminus)
sd(yearly(agpe_mat)$yearly_eplus)
mean(yearly(agpe_mat)$yearly_eminus_rec)
mean(yearly(agpe_mat)$yearly_eplus_rec)
sd(yearly(agpe_mat)$yearly_eminus_rec)
sd(yearly(agpe_mat)$yearly_eplus_rec)











image(bigmatrix(agpe_params)$Fmat_eplus_y1)
image(bigmatrix(agpe_params)$Fmat_eminus)
image(bigmatrix(agpe_params)$Tmat_eplus)
image(bigmatrix(agpe_params)$Tmat_eminus)
image(bigmatrix(agpe_params)$MPMmat_eplus)
image(bigmatrix(agpe_params)$MPMmat_eminus)
lambda(bigmatrix(agpe_params)$MPMmat_eplus_rec)
lambda(bigmatrix(agpe_params)$MPMmat_eminus)

image(bigmatrix(fesu_params)$Fmat_eplus)
image(bigmatrix(fesu_params)$Fmat_eminus)
image(bigmatrix(fesu_params)$Tmat_eplus)
image(bigmatrix(fesu_params)$Tmat_eminus)
image(bigmatrix(fesu_params)$MPMmat_eplus)
image(bigmatrix(fesu_params)$MPMmat_eminus)
lambda(bigmatrix(fesu_params)$MPMmat_eplus)
lambda(bigmatrix(fesu_params)$MPMmat_eminus)

image(bigmatrix(loar_params)$Fmat_eplus)
image(bigmatrix(loar_params)$Fmat_eminus)
image(bigmatrix(loar_params)$Tmat_eplus)
image(bigmatrix(loar_params)$Tmat_eminus)
image(bigmatrix(loar_params)$MPMmat_eplus)
image(bigmatrix(loar_params)$MPMmat_eminus)
lambda(bigmatrix(loar_params)$MPMmat_eplus)
lambda(bigmatrix(loar_params)$MPMmat_eminus)

image(bigmatrix(poal_params)$Fmat_eplus)
image(bigmatrix(poal_params)$Fmat_eminus)
image(bigmatrix(poal_params)$Tmat_eplus)
image(bigmatrix(poal_params)$Tmat_eminus)
image(bigmatrix(poal_params)$MPMmat_eplus)
image(bigmatrix(poal_params)$MPMmat_eminus)
lambda(bigmatrix(poal_params)$MPMmat_eplus)
lambda(bigmatrix(poal_params)$MPMmat_eminus)

image(bigmatrix(posy_params)$Fmat_eplus)
image(bigmatrix(posy_params)$Fmat_eminus)
image(bigmatrix(posy_params)$Tmat_eplus)
image(bigmatrix(posy_params)$Tmat_eminus)
image(bigmatrix(posy_params)$MPMmat_eplus)
image(bigmatrix(posy_params)$MPMmat_eminus)
lambda(bigmatrix(posy_params)$MPMmat_eplus)
lambda(bigmatrix(posy_params)$MPMmat_eminus)




