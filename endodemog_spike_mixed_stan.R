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
library(countreg)

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
ni <-5000
nb <- 2000
nc <- 3

sink("endodemog_spike.stan")
cat("
    data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> K;                       // number of predictors
    int<lower=0> lowerlimit;                         //lower limit for truncated negative binomial
    int<lower=0> nYear;                       // number of years (used as index)
    int<lower=0> year_t[N];                      // year of observation
    int<lower=0> nEndo;                       // number of endo treatments
    int<lower=1, upper=2> endo_index[N];          // index for endophyte effect
    int<lower=0> nPlot;                         // number of plots
    int<lower=0> plot[N];                   // plot of observation
    int<lower=0> spike_t[N];      // plant size at time t+1 and target variable (response)
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
    real<lower=0> phi;            // negative binomial dispersion parameter
    }
    
   transformed parameters{
    real mu[N];                         // Linear Predictor
 
    for(n in 1:N){
       mu[n] = beta[1] + beta[2]*logsize_t[n] + beta[3]*endo_01[n] + beta[4]*origin_01[n]
       + tau_year[endo_index[n],year_t[n]]
       + tau_plot[plot[n]];
    }
    }
    
    model {
    
    // Priors
    beta ~ normal(0,100);
    tau_plot ~ normal(0,sigma_p);   // prior for plot random effects
    to_vector(tau_year[1]) ~ normal(0,sigma_e[1]);   // prior for E- year random effects
    to_vector(tau_year[2]) ~ normal(0,sigma_e[2]);   // prior for E+ year random effects
    phi ~ cauchy(0., 5.);


    // Likelihood 
        for(n in 1:N){
      spike_t[n] ~ neg_binomial_2_log(mu[n],phi);
    }
    }
    
    generated quantities{
    }
  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_spike.stan")

# I fit this as a neg. binom. because there are some zeros in the data. The POAL for that individ has a notes that says that the flower was browsed so there was a flowering tiller but no spikelets.
# There are also zeros for some ELVI and several FESU. I think some of this could be from seeds that were collected in the field and counted later, so there could be some mix ups.

## Run the model by calling stan()
## Save the outputs as rds files

smAGPE<- stan(file = "endodemog_spike.stan", data = AGPE_spike_data_list,
              iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smAGPE, file = "endodemog_spike_AGPE.rds")

smELRI <- stan(file = "endodemog_spike.stan", data = ELRI_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELRI, file = "endodemog_spike_ELRI.rds")

smELVI <- stan(file = "endodemog_spike.stan", data = ELVI_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELVI, file = "endodemog_spike_ELVI.rds")

smFESU <- stan(file = "endodemog_spike.stan", data = FESU_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smFESU, file = "endodemog_spike_FESU.rds")

smLOAR <- stan(file = "endodemog_spike.stan", data = LOAR_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smLOAR, file = "endodemog_spike_LOAR.rds")

smPOAL <- stan(file = "endodemog_spike.stan", data = POAL_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOAL, file = "endodemog_spike_POAL.rds")

smPOSY <- stan(file = "endodemog_spike.stan", data = POSY_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOSY, file = "endodemog_spike_POSY.rds")

## to read in model output without rerunning models

smAGPE <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_AGPE.rds")
smELRI <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_ELRI.rds")
smELVI <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_ELVI.rds")
smFESU <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_FESU.rds")
smLOAR <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_LOAR.rds")
smPOAL <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_POAL.rds")
smPOSY <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_POSY.rds")



#########################################################################################################
###### Perform posterior predictive checks and assess model convergence-------------------------
#########################################################################################################
params <- c("beta[1]", "beta[2]", "tau_year[1,4]", "sigma_e[1]", "sigma_e[2]")


##### POAL - growth
print(smPOAL)

## plot traceplots of chains for select parameters
traceplot(smAGPE, pars = params)
traceplot(smELRI, pars = params)
traceplot(smELVI, pars = params)
traceplot(smFESU, pars = params)
traceplot(smLOAR, pars = params)
traceplot(smPOSY, pars = params)
traceplot(smPOSY, pars = params)


# Pull out the posteriors
post_spikeAGPE <- extract(smAGPE)
post_spikeELRI <- extract(smELRI)
post_spikeELVI <- extract(smELVI)
post_spikeFESU <- extract(smFESU)
post_spikeLOAR <- extract(smLOAR)
post_spikePOAL <- extract(smPOAL)
post_spikePOSY <- extract(smPOSY)



# This function generates replicate data for each given data point using the model matrix within the datalist for the random effects
prediction <- function(data, fit, n_post_draws){
  post <- extract(fit)
  mu <- post$mu
  phi <- post$phi
  yrep <- matrix(nrow = n_post_draws, ncol = data$N)
  for(n in 1:n_post_draws){
    yrep[n,] <- rnbinom(n = data$N, size = phi[n], mu = exp(mu[n,]))
  }
  out <- list(yrep, mu)
  names(out) <- c("yrep", "mu")
  return(out)
}


# apply the function for each species
AGPE_spike_yrep <- prediction(data = AGPE_spike_data_list, fit = smAGPE, n_post_draws = 500)
ELRI_spike_yrep <- prediction(data = ELRI_spike_data_list, fit = smELRI, n_post_draws = 500)
ELVI_spike_yrep <- prediction(data = ELVI_spike_data_list, fit = smELVI, n_post_draws = 500)
FESU_spike_yrep <- prediction(data = FESU_spike_data_list, fit = smFESU, n_post_draws = 500)
LOAR_spike_yrep <- prediction(data = LOAR_spike_data_list, fit = smLOAR, n_post_draws = 500)
POAL_spike_yrep <- prediction(data = POAL_spike_data_list, fit = smPOAL, n_post_draws = 500)
POSY_spike_yrep <- prediction(data = POSY_spike_data_list, fit = smPOSY, n_post_draws = 500)



# overlay 100 replicates over the actual dataset
ppc_dens_overlay( y = AGPE_spike_data_list$spike_t, yrep = AGPE_spike_yrep$yrep[1:100,])+ xlab("prob. of y") + ggtitle("AGPE")

ppc_dens_overlay( y = ELRI_spike_data_list$spike_t, yrep = ELRI_spike_yrep$yrep[1:100,])+ xlab("prob. of y") + ggtitle("ELRI")

ppc_dens_overlay( y = ELVI_spike_data_list$spike_t, yrep = ELVI_spike_yrep$yrep[1:100,]) + xlab("prob. of y") + ggtitle("ELVI")

ppc_dens_overlay( y = FESU_spike_data_list$spike_t, yrep = FESU_spike_yrep$yrep[1:100,]) + xlab("prob. of y") + ggtitle("FESU")

ppc_dens_overlay( y = LOAR_spike_data_list$spike_t, yrep = LOAR_spike_yrep$yrep[1:100,]) + xlab("prob. of y") + ggtitle("LOAR")

ppc_dens_overlay( y = POAL_spike_data_list$spike_t, yrep = POAL_spike_yrep$yrep[1:100,]) + xlab("prob. of y") + ggtitle("POAL")

ppc_dens_overlay( y = POSY_spike_data_list$spike_t, yrep = POSY_spike_yrep$yrep[1:100,]) + xlab("prob. of y") + ggtitle("POSY")








