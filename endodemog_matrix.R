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

getparams <- function(surv, grow, flw, fert, spike, seed, s_to_s, data){
params <- c()
params[1] <- lapply(rstan::extract(surv, pars = "beta[1]"), FUN = mean); names(params)[1] <- "surv_beta1"
params[2] <- lapply(rstan::extract(surv, pars = "beta[2]"), FUN = mean); names(params)[2] <- "surv_beta2"
params[3] <- lapply(rstan::extract(surv, pars = "beta[3]"), FUN = mean); names(params)[3] <- "surv_beta3"
params[4] <- lapply(rstan::extract(surv, pars = "beta[4]"), FUN = mean); names(params)[4] <- "surv_beta4"
# collect the growth parameters
params[5] <- lapply(rstan::extract(grow, pars = "beta[1]"), FUN = mean); names(params)[5] <- "grow_beta1"
params[6] <- lapply(rstan::extract(grow, pars = "beta[2]"), FUN = mean); names(params)[6] <- "grow_beta2"
params[7] <- lapply(rstan::extract(grow, pars = "beta[3]"), FUN = mean); names(params)[7] <- "grow_beta3"
params[8] <- lapply(rstan::extract(grow, pars = "beta[4]"), FUN = mean); names(params)[8] <- "grow_beta4"
# collect the flowering parameters
params[9] <- lapply(rstan::extract(flw, pars = "beta[1]"), FUN = mean); names(params)[9] <- "flw_beta1"
params[10] <- lapply(rstan::extract(flw, pars = "beta[2]"), FUN = mean); names(params)[10] <- "flw_beta2"
params[11] <- lapply(rstan::extract(flw, pars = "beta[3]"), FUN = mean); names(params)[11] <- "flw_beta3"
params[12] <- lapply(rstan::extract(flw, pars = "beta[4]"), FUN = mean); names(params)[12] <- "flw_beta4"
# collect the fertility parameters
params[13] <- lapply(rstan::extract(fert, pars = "beta[1]"), FUN = mean); names(params)[13] <- "fert_beta1"
params[14] <- lapply(rstan::extract(fert, pars = "beta[2]"), FUN = mean); names(params)[14] <- "fert_beta2"
params[15] <- lapply(rstan::extract(fert, pars = "beta[3]"), FUN = mean); names(params)[15] <- "fert_beta3"
params[16] <- lapply(rstan::extract(fert, pars = "beta[4]"), FUN = mean); names(params)[16] <- "fert_beta4"
# collect the spikelets/infl parameters
params[17] <- lapply(rstan::extract(spike, pars = "beta[1]"), FUN = mean); names(params)[17] <- "spike_beta1"
params[18] <- lapply(rstan::extract(spike, pars = "beta[2]"), FUN = mean); names(params)[18] <- "spike_beta2"
params[19] <- lapply(rstan::extract(spike, pars = "beta[3]"), FUN = mean); names(params)[19] <- "spike_beta3"
params[20] <- lapply(rstan::extract(spike, pars = "beta[4]"), FUN = mean); names(params)[20] <- "spike_beta4"
# collect the seed/spikelet parameters
params[21] <- lapply(rstan::extract(seed, pars = "mu_seed"), FUN = mean); names(params)[21] <- "mu_seed"
# collect the recruitment parameters
params[22] <- lapply(rstan::extract(s_to_s, pars = "beta[1]"), FUN = mean); names(params)[22] <- "s_to_s_beta1"
params[23] <- lapply(rstan::extract(s_to_s, pars = "beta[2]"), FUN = mean); names(params)[23] <- "s_to_s_beta2"
# min and max size
params[24] <- 1; names(params)[24]<-"min_size"
params[25] <- quantile(data$size_t1,0.95,na.rm=T); names(params)[25]<-"max_size"
# note that I define max size as the 95TH pctile of observed sized. The very max sizes observed often have very poor 
# replication, and I find that these few observations (and the corresponding vital rate predictions) can have a strong
# influence on the results. So this approach is conservative, but you can always experiment with this.
params[26] <- lapply(rstan::extract(grow, pars = "phi"), FUN = mean); names(params)[26] <- "grow_phi"
params[27] <- lapply(rstan::extract(fert, pars = "phi"), FUN = mean); names(params)[27] <- "fert_phi"
params[28] <- lapply(rstan::extract(spike, pars = "phi"), FUN = mean); names(params)[28] <- "spike_phi"
params <- unlist(params)
return(params)
}

agpe_params <- getparams(surv = survAGPE, grow = growAGPE, flw = flwAGPE, fert = fertAGPE, spike = spikeAGPE, seed = seedAGPE, s_to_s = s_to_sAGPE, data = AGPE_surv_data)
fesu_params <- getparams(surv = survFESU, grow = growFESU, flw = flwFESU, fert = fertFESU, spike = spikeFESU, seed = seedFESU, s_to_s = s_to_sFESU, data = FESU_surv_data)
loar_params <- getparams(surv = survLOAR, grow = growLOAR, flw = flwLOAR, fert = fertLOAR, spike = spikeLOAR, seed = seedLOAR, s_to_s = s_to_sLOAR, data = LOAR_surv_data)
poal_params <- getparams(surv = survPOAL, grow = growPOAL, flw = flwPOAL, fert = fertPOAL, spike = spikePOAL, seed = seedPOAL, s_to_s = s_to_sPOAL, data = POAL_surv_data)
posy_params <- getparams(surv = survPOSY, grow = growPOSY, flw = flwPOSY, fert = fertPOSY, spike = spikePOSY, seed = seedPOSY, s_to_s = s_to_sPOSY, data = POSY_surv_data)


# define functions that will be used to populate projection matrix
#SURVIVAL AT SIZE X.
# currently this is fitting as if the E- is the intercept

sx<-function(x,params){
  eminus_surv <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x))
  eplus_surv <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"])
  return(list(eminus_surv = eminus_surv, eplus_surv = eplus_surv))
}
# sx(x = 10, params = poal_params) # we can test out our functions for different sizes of x
gxy <- function(x,y,params){
  eminus_grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x)
  eminus_pr_grow <- dnbinom(x=y, mu = exp(eminus_grow.mean), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean), size = params["grow_phi"])))
  eplus_grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"]
  eplus_pr_grow <- dnbinom(x=y, mu = exp(eplus_grow.mean), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean), size = params["grow_phi"])))
  
return(list(eminus_pr_grow  = eminus_pr_grow, eplus_pr_grow = eplus_pr_grow))
}
# gxy(x = 10, y = 11, params = loar_params)$eplus_pr_grow
# I'll truncate this later...
# prob=dnbinom(1:n_post_draws, mu = exp(mu[i,j]), size = phi[i]/(1-dnbinom(0, mu = exp(mu[i,j]), size = phi[i]))))


#SURVIVAL*GROWTH
pxy<-function(x,y,params){
  eminuspxy <- sx(x,params)$eminus_surv * gxy(x,y,params)$eminus_pr_grow
  epluspxy <- sx(x,params)$eplus_surv * gxy(x,y,params)$eplus_pr_grow
  return(list(eminuspxy = eminuspxy, epluspxy = epluspxy))
}
# pxy(x = 2, y = 10, params = loar_params)

#FERTILITY--returns number of seedlings, which we will assume (for now) to be 1-tiller, produced by size X
fx<-function(x, params){
  p_flw <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x)) 
  fert.mean <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x))
  # p_fert <- exp(fert.mean)
  #p_fert <- dnbinom(x=y, prob = exp(fert.mean), size = params["fert_phi"])
  spike.mean <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x))
  # p_spike <- exp(spike.mean)
  #p_spike <-dnbinom(x=y, prob = exp(spike.mean), size = params["spike_phi"])
  seed.mean <- params["mu_seed"]
  p_rec <- invlogit(params["s_to_s_beta1"])
  seedlings <- p_flw * fert.mean * spike.mean * seed.mean * p_rec
  return(seedlings)
}

## note from Tom: we should think about adding a reproductive delay
# finally, here is the function that takes the parameter vector and assembles the matrix model from all of the pieces
bigmatrix<-function(params){   
  
  matdim<-params["max_size"]## matrix dimension
  y <- 1:params["max_size"]## size (tiller number) associated with each class
  # Fertility matrix
  Fmat<-matrix(0,matdim,matdim)
  # all seedlings get dumped into top row (1-tiller)
  Fmat[1,]<-fx(x = y, params=params) 
  
  # Growth/survival transition matrix
  Tmat<-matrix(0,matdim,matdim)
  Tmat<-t(outer(y,y,pxy,params=params))
  
  # Put it all together
  MPMmat<-Tmat+Fmat #sum the Tmat & Fmat to get the whole matrix
  
  return(list(MPMmat=MPMmat,Fmat=Fmat,Tmat=Tmat))
}

## population growth rate (eigenalaysis of the projection matrix)
# matrix <- bigmatrix(loar_params)
image(bigmatrix(agpe_params)$Fmat)
image(bigmatrix(agpe_params)$Tmat)
image(bigmatrix(agpe_params)$MPMmat)
lambda(bigmatrix(agpe_params)$MPMmat)

image(bigmatrix(fesu_params)$Fmat)
image(bigmatrix(fesu_params)$Tmat)
image(bigmatrix(fesu_params)$MPMmat)
lambda(bigmatrix(fesu_params)$MPMmat)

image(bigmatrix(loar_params)$Fmat)
image(bigmatrix(loar_params)$Tmat)
image(bigmatrix(loar_params)$MPMmat)
lambda(bigmatrix(loar_params)$MPMmat)

image(bigmatrix(poal_params)$Fmat)
image(bigmatrix(poal_params)$Tmat)
image(bigmatrix(poal_params)$MPMmat)
lambda(bigmatrix(poal_params)$MPMmat)

image(bigmatrix(posy_params)$Fmat)
image(bigmatrix(posy_params)$Tmat)
image(bigmatrix(posy_params)$MPMmat)
lambda(bigmatrix(posy_params)$MPMmat)
