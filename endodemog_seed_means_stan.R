## Title: Grass endophyte population model with a bayesian framework
## Purpose: Creates seed production kernel written in STAN, 
## and does visualisation of posterior predictive checks
## Authors: Joshua and Tom
#############################################################

library(tidyverse)
library(rstan)
library(StanHeaders)
library(shinystan)
library(bayesplot)
library(devtools)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }


#############################################################################################
####### Data manipulation to prepare data as lists for Stan models------------------
#############################################################################################

# seed data lists are generated in the endodemog_data_processing.R file, 
# within the section titled "Preparing datalists for Seed Means Kernel"
source("endodemog_data_processing.R")

#########################################################################################################
# Stan model for mean of seed production per spikelet ------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <-10000
nb <- 5000
nc <- 3

# Stan model -------------
## here is the Stan model ##

sink("endodemog_seed_mean.stan")
cat("
    data { 
    int<lower=0> Nseed;                       // number of observations of seed/spikelet
    int<lower=0> K;                            // number of predictors
    int<lower=0,upper=1> endo_01[Nseed];          // Endophyte status
    real<lower=0> seed[Nseed];               // number of seeds per spikelet
    }
    
    parameters {
    vector[K] beta;
    real<lower=0> sigma_seed;
    }
    
    transformed parameters{
    vector[Nseed] mu_seed;
    for(n in 1:Nseed){
    mu_seed[n] = beta[1] + beta[2]*endo_01[n];
    }
    }
    
    model {

    // Priors
    beta ~ normal(0,100);
    // Likelihood
      seed ~ normal(mu_seed,sigma_seed);
    }

    ", fill = T)
sink()

stanmodel <- stanc("endodemog_seed_mean.stan")

## Run the model by calling stan()
## Save the outputs as rds files
smFESU <- stan(file = "endodemog_seed_mean.stan", data = FESU_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smFESU, file = "endodemog_seed_mean_FESU.rds")

smLOAR <- stan(file = "endodemog_seed_mean.stan", data = LOAR_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smLOAR, file = "endodemog_seed_mean_LOAR.rds")

smPOAL <- stan(file = "endodemog_seed_mean.stan", data = POAL_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOAL, file = "endodemog_seed_mean_POAL.rds")

smPOSY <- stan(file = "endodemog_seed_mean.stan", data = POSY_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOSY, file = "endodemog_seed_mean_POSY.rds")

# AGPE had data recorded as seed/spikelet already, but we are using this model to calculate seed/spikelet on average, so this is using the seed/spikelet info from AGPE
smAGPE<- stan(file = "endodemog_seed_mean.stan", data = AGPE_seed_data_list,
              iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smAGPE, file = "endodemog_seed_mean_AGPE.rds")


# ELRI and ELVI had data collected in a slightly different way. They recorded seeds/inflorescence, so this is our calculation for the mean seeds/infl
smELRI <- stan(file = "endodemog_seed_mean.stan", data = ELRI_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELRI, file = "endodemog_seed_mean_ELRI.rds")

smELVI <- stan(file = "endodemog_seed_mean.stan", data = ELVI_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELVI, file = "endodemog_seed_mean_ELVI.rds")


## to read in model output without rerunning models
smAGPE <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_AGPE.rds")
smELRI <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_ELRI.rds")
smELVI <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_ELVI.rds")
smFESU <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_FESU.rds")
smLOAR <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_LOAR.rds")
smPOAL <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_POAL.rds")
smPOSY <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_POSY.rds")


#########################################################################################################
###### Perform posterior predictive checks and assess model convergence-------------------------
#########################################################################################################
params <- c("beta", "sigma_seed")

##### POAL - survival
print(smPOAL)
# summary(smPOAL)

## plot traceplots of chains for select parameters
traceplot(smAGPE, pars = params)
traceplot(smELRI, pars = params)
traceplot(smELVI, pars = params)
traceplot(smFESU, pars = params)
traceplot(smLOAR, pars = params)
traceplot(smPOSY, pars = params)
traceplot(smPOSY, pars = params)


# Pull out the posteriors
post_survAGPE <- extract(smAGPE)
post_survELRI <- extract(smELRI)
post_survELVI <- extract(smELVI)
post_survFESU <- extract(smFESU)
post_survLOAR <- extract(smLOAR)
post_survPOAL <- extract(smPOAL)
post_survPOSY <- extract(smPOSY)


# This function extracts the posterior draws and generates replicate data for each given model
prediction <- function(data, fit, n_post_draws){
  post <- rstan::extract(fit)
  mu <- post$mu_seed
  sd <- post$sigma_seed
  yrep <- matrix(nrow = n_post_draws, ncol = data$N)
  for(n in 1:n_post_draws){
    yrep[n,] <- rnorm(n = data$N, mean = mu[n], sd = sd[n])
  }
  out <- list(yrep, mu)
  names(out) <- c("yrep", "mu")
  return(out)
}

# apply the function for each species
AGPE_seed_yrep <- prediction(data = AGPE_seed_data_list, fit = smAGPE, n_post_draws = 500)
ELRI_seed_yrep <- prediction(data = ELRI_seed_data_list, fit = smELRI, n_post_draws = 500)
ELVI_seed_yrep <- prediction(data = ELVI_seed_data_list, fit = smELVI, n_post_draws = 500)
FESU_seed_yrep <- prediction(data = FESU_seed_data_list, fit = smFESU, n_post_draws = 500)
LOAR_seed_yrep <- prediction(data = LOAR_seed_data_list, fit = smLOAR, n_post_draws = 500)
POAL_seed_yrep <- prediction(data = POAL_seed_data_list, fit = smPOAL, n_post_draws = 500)
POSY_seed_yrep <- prediction(data = POSY_seed_data_list, fit = smPOSY, n_post_draws = 500)


# overlay 100 replicates over the actual dataset
ppc_dens_overlay( y = AGPE_seed_data_list$seed, yrep = AGPE_seed_yrep$yrep[1:100,]) + ggtitle("AGPE")

ppc_dens_overlay( y = ELRI_seed_data_list$seed, yrep = ELRI_seed_yrep$yrep[1:100,]) + ggtitle("ELRI")

ppc_dens_overlay( y = ELVI_seed_data_list$seed, yrep = ELVI_seed_yrep$yrep[1:100,]) + ggtitle("ELVI")

ppc_dens_overlay( y = FESU_seed_data_list$seed, yrep = FESU_seed_yrep$yrep[1:100,]) + ggtitle("FESU")

ppc_dens_overlay( y = LOAR_seed_data_list$seed, yrep = LOAR_seed_yrep$yrep[1:100,]) + ggtitle("LOAR")

ppc_dens_overlay( y = POAL_seed_data_list$seed, yrep = POAL_seed_yrep$yrep[1:100,]) + ggtitle("POAL")

ppc_dens_overlay( y = POSY_seed_data_list$seed, yrep = POSY_seed_yrep$yrep[1:100,]) + ggtitle("POSY")




