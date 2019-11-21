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


#  Read in the model outputs ----------------------------------------------
source("endodemog_model_outputs.R")


# Assembling matrix model -------------------------------------------------

# LOAR
# collect all the parameters into one vector
# collect the survival parameters
loar_params <- c()
loar_params[1] <- rstan::extract(survLOAR, pars = "beta[1]"); names(loar_params)[1] <- "surv_beta1"
loar_params[2] <- rstan::extract(survLOAR, pars = "beta[2]"); names(loar_params)[2] <- "surv_beta2"
loar_params[3] <- rstan::extract(survLOAR, pars = "beta[3]"); names(loar_params)[3] <- "surv_beta3"
loar_params[4] <- rstan::extract(survLOAR, pars = "beta[4]"); names(loar_params)[4] <- "surv_beta4"
# collect the growth parameters
loar_params[5] <- rstan::extract(growLOAR, pars = "beta[1]"); names(loar_params)[5] <- "grow_beta1"
loar_params[6] <- rstan::extract(growLOAR, pars = "beta[2]"); names(loar_params)[6] <- "grow_beta2"
loar_params[7] <- rstan::extract(growLOAR, pars = "beta[3]"); names(loar_params)[7] <- "grow_beta3"
loar_params[8] <- rstan::extract(growLOAR, pars = "beta[4]"); names(loar_params)[8] <- "grow_beta4"
# collect the flowering parameters
loar_params[9] <- rstan::extract(flwLOAR, pars = "beta[1]"); names(loar_params)[9] <- "flw_beta1"
loar_params[10] <- rstan::extract(flwLOAR, pars = "beta[2]"); names(loar_params)[10] <- "flw_beta2"
loar_params[11] <- rstan::extract(flwLOAR, pars = "beta[3]"); names(loar_params)[11] <- "flw_beta3"
loar_params[12] <- rstan::extract(flwLOAR, pars = "beta[4]"); names(loar_params)[12] <- "flw_beta4"
# collect the fertility parameters
loar_params[13] <- rstan::extract(fertLOAR, pars = "beta[1]"); names(loar_params)[13] <- "fert_beta1"
loar_params[14] <- rstan::extract(fertLOAR, pars = "beta[2]"); names(loar_params)[14] <- "fert_beta2"
loar_params[15] <- rstan::extract(fertLOAR, pars = "beta[3]"); names(loar_params)[15] <- "fert_beta3"
loar_params[16] <- rstan::extract(fertLOAR, pars = "beta[4]"); names(loar_params)[16] <- "fert_beta4"
# collect the spikelets/infl parameters
loar_params[17] <- rstan::extract(spikeLOAR, pars = "beta[1]"); names(loar_params)[17] <- "spike_beta1"
loar_params[18] <- rstan::extract(spikeLOAR, pars = "beta[2]"); names(loar_params)[18] <- "spike_beta2"
loar_params[19] <- rstan::extract(spikeLOAR, pars = "beta[3]"); names(loar_params)[19] <- "spike_beta3"
loar_params[20] <- rstan::extract(spikeLOAR, pars = "beta[4]"); names(loar_params)[20] <- "spike_beta4"
# collect the seed/spikelet parameters
loar_params[21] <- rstan::extract(seedLOAR, pars = "mu_seed"); names(loar_params)[21] <- "mu_seed"
# collect the recruitment parameters
loar_params[22] <- rstan::extract(s_to_s_LOAR, pars = "beta[1]"); names(loar_params)[22] <- "s_to_s_beta1"
loar_params[23] <- rstan::extract(s_to_s_LOAR, pars = "beta[2]"); names(loar_params)[23] <- "s_to_s_beta2"
# min and max size
loar_params[24] <- 1; names(loar_params)[24]<-"min_size"
loar_params[25] <- quantile(LOAR_grow_data_list$size_t1,0.95,na.rm=T); names(loar_params)[25]<-"max_size"
# note that I define max size as the 95TH pctile of observed sized. The very max sizes observed often have very poor 
# replication, and I find that these few observations (and the corresponding vital rate predictions) can have a strong
# influence on the results. So this approach is conservative, but you can always experiment with this.
loar_params[26] <- rstan::extract(growLOAR, pars = "phi"); names(loar_params[26]) <- "grow_phi"
loar_params[27] <- rstan::extract(fertLOAR, pars = "phi"); names(loar_params[27]) <- "fert_phi"
loar_params[28] <- rstan::extract(spikeLOAR, pars = "phi"); names(loar_params[28]) <- "spike_phi"


# define functions that will be used to populate projection matrix
#SURVIVAL AT SIZE X.
# currently this is fitting as if the E- is the intercept

sx<-function(x,params){
  invlogit(params$surv_beta1 + params$surv_beta2*log(x))
}


gxy <- function(x,y,params){
  grow.mean <- params$grow_beta1 + params$grow_beta2*log(x)
  pr_grow <- sample(x=y, size=1, replace=T, prob=dnbinom(1:y, mu = exp(grow.mean), size = params$grow_phi/(1 - (dnbinom(0, mu = exp(grow.mean), size = params$grow_phi)) + sum(dnbinom(x = (params$max_size+1):500, mu = exp(grow.mean), size = params$grow_phi)))))
  return(pr_grow)
}

#SURVIVAL*GROWTH
pxy<-function(x,y,params){
  sx(x,params) * gxy(x,y,params)
}

#FERTILITY--returns number of seedlings, which we will assume (for now) to be 1-tiller, produced by size X
fx<-function(x,params){
  p_flw <- invlogit(params$flw_beta1 + params$flw_beta2*log(x))
  fert.mean <- params$fert_beta1 + params$fert_beta2*log(x)
  p_fert <- sample(x=y, size=1, replace=T, prob=dnbinom(1:y, mu = exp(fert.mean), size = params$fert_phi/(1 - (dnbinom(0, mu = exp(fert.mean), size = params$fert_phi)) + sum(dnbinom(x = (params$max_size+1):500, mu = exp(fert.mean), size = params$fert_phi)))))
  spike.mean <- params$spike_beta1 + params$spike_beta2*log(x)
  p_spike <- sample(x=y, size=1, replace=T, prob=dnbinom(1:y, mu = exp(spike.mean), size = params$spike_phi/(1 - (dnbinom(0, mu = exp(spike.mean), size = params$spike_phi)) + sum(dnbinom(x = (params$max_size+1):500, mu = exp(spike.mean), size = params$spike_phi)))))
  seed.mean <- params$mu_seed
  p_rec <- inv_logit(params$s_to_s_beta1)
  seedlings <- p_flw * p_fert * p_spike * seed.mean * p_rec
  return(seedlings)
}


# finally, here is the function that takes the parameter vector and assembles the matrix model from all of the pieces
bigmatrix<-function(params){   
  
  matdim<-params$max_size ## matrix dimension
  y<-1:params$max_size ## size (tiller number) associated with each class
  
  # Fertility matrix
  Fmat<-matrix(0,matdim,matdim)
  # all seedlings get dumped into top row (1-tiller)
  Fmat[1,]<-fx(x=y,params=params) 
  
  # Growth/survival transition matrix
  Tmat<-matrix(0,matdim,matdim)
  Tmat<-t(outer(y,y,pxy,params=params)) 
  
  # Put it all together
  MPMmat<-Tmat+Fmat #sum the Tmat & Fmat to get the whole matrix
  
  return(list(MPMmat=MPMmat,Fmat=Fmat,Tmat=Tmat))
}

## population growth rate (eigenalaysis of the projection matrix)
# matrix <- bigmatrix(loar_params)
lambda(bigmatrix(loar_params)$Tmat)


