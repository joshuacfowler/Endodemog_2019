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
tompath <- "C:/Users/tm9/Dropbox/"
# Growth data lists are generated in the endodemog_data_processing.R file
# within the section titled "Preparing datalists for Growth Kernel"
# source("endodemog_data_processing.R"

endo_demog <- read.csv("/Users/joshuacfowler/Dropbox/EndodemogData/fulldataplusmetadata/endo_demog_long.csv") 
loar <- endo_demog %>% 
  filter(species == "LOAR") %>% 
  mutate(log_size_t = log(size_t),
         log_size_t1 = log(size_t1),
         flow_t = as.integer(seed_t > 0))

# DATA ISSUES:
# I noticed that there are about 10 plants with a size of 0 tillers. They are:
filter(loar,size_t==0)
# I also noticed that there are a few instances where survival was > 1:
filter(loar,surv_t1 > 1)
# We'll need to cut these for the analysis.
loar <- loar %>% filter(size_t > 0,
                        surv_t1 <= 1)

# 
# #  Read in the model outputs ----------------------------------------------
# source("endodemog_model_outputs.R")

survLOAR <- read_rds(path = "C:/Users/tm9/Dropbox/EndodemogData/Model_Runs/endodemog_surv_LOAR.rds")

growLOAR <- read_rds(path = "C:/Users/tm9/Dropbox/EndodemogData/Model_Runs/endodemog_grow_LOAR.rds")
flwLOAR <- read_rds(path = "C:/Users/tm9/Dropbox/EndodemogData/Model_Runs/endodemog_flw_LOAR.rds")

fertLOAR <- read_rds(path = "C:/Users/tm9/Dropbox/EndodemogData/Model_Runs/endodemog_fert_LOAR.rds")

spikeLOAR <- read_rds(path = "C:/Users/tm9/Dropbox/EndodemogData/Model_Runs/endodemog_spike_LOAR.rds")
seedLOAR <- read_rds(path = "C:/Users/tm9/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_LOAR.rds")
s_to_sLOAR <- read_rds(path = "C:/Users/tm9/Dropbox/EndodemogData/Model_Runs/endodemog_s_to_s_LOAR.rds")

#############################################################################################
# Assembling matrix model -------------------------------------------------
#############################################################################################

# LOAR
# collect all the parameters into one vector
# collect the survival parameters
loar_params <- c()
loar_params[1] <- lapply(rstan::extract(survLOAR, pars = "beta[1]"), FUN = mean); names(loar_params)[1] <- "surv_beta1"
loar_params[2] <- lapply(rstan::extract(survLOAR, pars = "beta[2]"), FUN = mean); names(loar_params)[2] <- "surv_beta2"
loar_params[3] <- lapply(rstan::extract(survLOAR, pars = "beta[3]"), FUN = mean); names(loar_params)[3] <- "surv_beta3"
loar_params[4] <- lapply(rstan::extract(survLOAR, pars = "beta[4]"), FUN = mean); names(loar_params)[4] <- "surv_beta4"
# collect the growth parameters
loar_params[5] <- lapply(rstan::extract(growLOAR, pars = "beta[1]"), FUN = mean); names(loar_params)[5] <- "grow_beta1"
loar_params[6] <- lapply(rstan::extract(growLOAR, pars = "beta[2]"), FUN = mean); names(loar_params)[6] <- "grow_beta2"
loar_params[7] <- lapply(rstan::extract(growLOAR, pars = "beta[3]"), FUN = mean); names(loar_params)[7] <- "grow_beta3"
loar_params[8] <- lapply(rstan::extract(growLOAR, pars = "beta[4]"), FUN = mean); names(loar_params)[8] <- "grow_beta4"
# collect the flowering parameters
loar_params[9] <- lapply(rstan::extract(flwLOAR, pars = "beta[1]"), FUN = mean); names(loar_params)[9] <- "flw_beta1"
loar_params[10] <- lapply(rstan::extract(flwLOAR, pars = "beta[2]"), FUN = mean); names(loar_params)[10] <- "flw_beta2"
loar_params[11] <- lapply(rstan::extract(flwLOAR, pars = "beta[3]"), FUN = mean); names(loar_params)[11] <- "flw_beta3"
loar_params[12] <- lapply(rstan::extract(flwLOAR, pars = "beta[4]"), FUN = mean); names(loar_params)[12] <- "flw_beta4"
# collect the fertility parameters
loar_params[13] <- lapply(rstan::extract(fertLOAR, pars = "beta[1]"), FUN = mean); names(loar_params)[13] <- "fert_beta1"
loar_params[14] <- lapply(rstan::extract(fertLOAR, pars = "beta[2]"), FUN = mean); names(loar_params)[14] <- "fert_beta2"
loar_params[15] <- lapply(rstan::extract(fertLOAR, pars = "beta[3]"), FUN = mean); names(loar_params)[15] <- "fert_beta3"
loar_params[16] <- lapply(rstan::extract(fertLOAR, pars = "beta[4]"), FUN = mean); names(loar_params)[16] <- "fert_beta4"
# collect the spikelets/infl parameters
loar_params[17] <- lapply(rstan::extract(spikeLOAR, pars = "beta[1]"), FUN = mean); names(loar_params)[17] <- "spike_beta1"
loar_params[18] <- lapply(rstan::extract(spikeLOAR, pars = "beta[2]"), FUN = mean); names(loar_params)[18] <- "spike_beta2"
loar_params[19] <- lapply(rstan::extract(spikeLOAR, pars = "beta[3]"), FUN = mean); names(loar_params)[19] <- "spike_beta3"
loar_params[20] <- lapply(rstan::extract(spikeLOAR, pars = "beta[4]"), FUN = mean); names(loar_params)[20] <- "spike_beta4"
# collect the seed/spikelet parameters
loar_params[21] <- lapply(rstan::extract(seedLOAR, pars = "mu_seed"), FUN = mean); names(loar_params)[21] <- "mu_seed"
# collect the recruitment parameters
loar_params[22] <- lapply(rstan::extract(s_to_sLOAR, pars = "beta[1]"), FUN = mean); names(loar_params)[22] <- "s_to_s_beta1"
loar_params[23] <- lapply(rstan::extract(s_to_sLOAR, pars = "beta[2]"), FUN = mean); names(loar_params)[23] <- "s_to_s_beta2"
# min and max size
loar_params[24] <- 1; names(loar_params)[24]<-"min_size"
loar_params[25] <- quantile(loar$size_t1,0.95,na.rm=T); names(loar_params)[25]<-"max_size"
# note that I define max size as the 95TH pctile of observed sized. The very max sizes observed often have very poor 
# replication, and I find that these few observations (and the corresponding vital rate predictions) can have a strong
# influence on the results. So this approach is conservative, but you can always experiment with this.
loar_params[26] <- lapply(rstan::extract(growLOAR, pars = "phi"), FUN = mean); names(loar_params)[26] <- "grow_phi"
loar_params[27] <- lapply(rstan::extract(fertLOAR, pars = "phi"), FUN = mean); names(loar_params)[27] <- "fert_phi"
loar_params[28] <- lapply(rstan::extract(spikeLOAR, pars = "phi"), FUN = mean); names(loar_params)[28] <- "spike_phi"
loar_params <- unlist(loar_params)

# define functions that will be used to populate projection matrix
#SURVIVAL AT SIZE X.
# currently this is fitting as if the E- is the intercept

sx<-function(x,params){
  invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x))
}

gxy <- function(x,y,params){
  grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x)
  # pr_grow <- exp(grow.mean)
  pr_grow <- dnbinom(x=y, prob = exp(grow.mean), size = params["grow_phi"])
 return(pr_grow)
}

# I'll truncate this later...
# prob=dnbinom(1:n_post_draws, mu = exp(mu[i,j]), size = phi[i]/(1-dnbinom(0, mu = exp(mu[i,j]), size = phi[i]))))


#SURVIVAL*GROWTH
pxy<-function(x,y,params){
  sx(x,params) * gxy(x,y,params)
}

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
bigmatrix
image(bigmatrix(loar_params)$Fmat)
image(bigmatrix(loar_params)$Tmat)

lambda(bigmatrix(loar_params)$MPMmat)

y <- 1:33
grow.mean <- loar_params["grow_beta1"] + loar_params["grow_beta2"]*log(x)
spike.mean <- loar_params["spike_beta1"] + loar_params["spike_beta2"]*log(x)
fert.mean <- loar_params["fert_beta1"] + loar_params["fert_beta2"]*log(x)

g <- dnbinom(x = y, prob = grow.mean, size = loar_params["grow_phi"])
s <- dnbinom(x = y, prob = spike.mean, size = loar_params["spike_phi"])
f <- dnbinom(x = y, prob = fert.mean, size = loar_params["fert_phi"])

g
s
f


p_fert <- dnbinom(x=y, prob = exp(fert.mean), size = params["fert_phi"])

v <- gxy(x = 1:33, y = 1:33, params = loar_params)
v
z <- fx(x = 1:loar_params["max_size"], y = 1:loar_params["max_size"], params = loar_params)
z
