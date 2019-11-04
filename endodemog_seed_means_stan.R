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
    real<lower=0> seed[Nseed];               // number of seeds per spikelet
    }
    
    
    parameters {
    real<lower=0> mu_seed;
    real<lower=0> sigma_seed;
    }
    
    model {

    // Priors
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
# saveRDS(smFESU, file = "endodemog_seed_mean_FESU.rds")

smLOAR <- stan(file = "endodemog_seed_mean.stan", data = LOAR_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smLOAR, file = "endodemog_seed_mean_LOAR.rds")

smPOAL <- stan(file = "endodemog_seed_mean.stan", data = POAL_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOAL, file = "endodemog_seed_mean_POAL.rds")

smPOSY <- stan(file = "endodemog_seed_mean.stan", data = POSY_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOSY, file = "endodemog_seed_mean_POSY.rds")

# AGPE had data recorded as seed/spikelet already, but we are using this model to calculate seed/spikelet on average, so this is using the seed/spikelet info from AGPE
smAGPE<- stan(file = "endodemog_seed_mean.stan", data = AGPE_seed_data_list,
              iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smAGPE, file = "endodemog_seed_mean_AGPE.rds")


# ELRI and ELVI had data collected in a slightly different way. They recorded seeds/inflorescence, so this is our calculation for the mean seeds/infl
smELRI <- stan(file = "endodemog_seed_mean.stan", data = ELRI_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELRI, file = "endodemog_seed_mean_ELRI.rds")

smELVI <- stan(file = "endodemog_seed_mean.stan", data = ELVI_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELVI, file = "endodemog_seed_mean_ELVI.rds")



########################################################################################################################################
##### Checking the mean from the output
########################################################################################################################################
post_seedAGPE <- extract(smAGPE)
# post_survELRI <- extract(smELRI)
# post_survELVI <- extract(smELVI)
post_seedFESU <- extract(smFESU)
post_seedLOAR <- extract(smLOAR)
post_seedPOAL <- extract(smPOAL)
post_seedPOSY <- extract(smPOSY)


seed_histogram <- function(fit,data){
  require(ggplot2)
  spp <- ggplot() + 
    geom_histogram(data = fit, aes(mu_seed), fill = "#ff7f77", alpha = .4, bins = 200)  +
    geom_vline(aes(xintercept = mean(data$seed))) +
    labs(x = "mean seed per inflorescence", y = "Posterior density",
         title = deparse(substitute(df))) + theme_classic()
  return(spp)
}
# survival
AGPE_seed <- seed_histogram(fit = as.data.frame(post_seedAGPE), data = AGPE_seed_data_list)
AGPE_seed
FESU_seed <- seed_histogram(fit = as.data.frame(post_seedFESU), data = FESU_seed_data_list)
FESU_seed
LOAR_seed <- seed_histogram(fit = as.data.frame(post_seedLOAR), data = LOAR_seed_data_list)
LOAR_seed
POAL_seed <- seed_histogram(fit = as.data.frame(post_seedPOAL), data = POAL_seed_data_list)
POAL_seed
POSY_seed <- seed_histogram(fit = as.data.frame(post_seedPOSY), data = POSY_seed_data_list)
POSY_seed


