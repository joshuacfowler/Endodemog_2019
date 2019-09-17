## Title: Grass endophyte population model with a bayesian framework
## Purpose: Creates spikelet production kernel written in STAN, 
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

# spikelet data lists are generated in the endodemog_data_processing.R file, 
# within the section titled "Preparing datalists for Spikelet Kernel"
source("endodemog_data_processing.R")
#########################################################################################################
# Stan model for spikelet production ------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <-1000
nb <- 500
nc <- 3

sink("endodemog_spike.stan")
cat("
    data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> K;                       // number of predictors
    int<lower=0> nYear;                       // number of years (used as index)
    int<lower=0> year_t[N];                      // year of observation
    int<lower=0> nEndo;                       // number of endo treatments
    int<lower=1, upper=2> endo_index[N];          // index for endophyte effect
    int<lower=0> nPlot;                         // number of plots
    int<lower=0> plot[N];                   // plot of observation
    real<lower=0> spike_t[N];      // plant size at time t+1 and target variable (response)
    vector<lower=0>[N] logsize_t;             // plant size at time t (predictor)
    int<lower=0,upper=1> endo_01[N];            // plant endophyte status (predictor)
    int<lower=0,upper=1> origin_01[N];          // plant origin status (predictor)
    }
    
    parameters {
    vector[K] beta;                     // predictor parameters

    vector[nYear] tau_year[nEndo];      // random year effect
    real<lower=0> sigma_e[nEndo];        //year variance by endophyte effect
    vector[nPlot] tau_plot;        // random plot effect
    real<lower=0> sigma_p;          // plot variance
    real<lower=0.001,upper=100> shape; // shape parameter for gamma distribution
    }
  
    model {
    vector[N] mu;
    
    // Linear Predictor
    for(n in 1:N){
       mu[n] = beta[1] + beta[2]*logsize_t[n] + beta[3]*endo_01[n] +beta[4]*origin_01[n]
       + beta[5]*logsize_t[n]*endo_01[n] 
       + tau_year[endo_index[n],year_t[n]]
       + tau_plot[plot[n]];
    }
    
    // Priors
    beta ~ normal(0,100);
    tau_plot ~ normal(0,sigma_p);   // prior for plot random effects
    to_vector(tau_year[1]) ~ normal(0,sigma_e[1]);   // prior for E- year random effects
    to_vector(tau_year[2]) ~ normal(0,sigma_e[2]);   // prior for E+ year random effects
      
    // Likelihood 
    for(i in 1:N)
      spike_t[i] ~ gamma(shape, (shape / exp(mu[i])));
    }
    
    generated quantities{
    }
  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_spike.stan")



## Run the model by calling stan()
## Save the outputs as rds files

smAGPE<- stan(file = "endodemog_spike.stan", data = AGPE_spike_data_list,
              iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smAGPE, file = "endodemog_spike_AGPE_withplot.rds")

smELRI <- stan(file = "endodemog_spike.stan", data = ELRI_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELRI, file = "endodemog_spike_ELRI_withplot.rds")

smELVI <- stan(file = "endodemog_spike.stan", data = ELVI_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELVI, file = "endodemog_spike_ELVI_withplot.rds")

smFESU <- stan(file = "endodemog_spike.stan", data = FESU_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smFESU, file = "endodemog_spike_FESU_withplot.rds")

smLOAR <- stan(file = "endodemog_spike.stan", data = LOAR_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smLOAR, file = "endodemog_spike_LOAR_withplot.rds")

smPOAL <- stan(file = "endodemog_spike.stan", data = POAL_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOAL, file = "endodemog_spike_POAL_withplot.rds")

smPOSY <- stan(file = "endodemog_spike.stan", data = POSY_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOSY, file = "endodemog_spike_POSY_withplot.rds")

