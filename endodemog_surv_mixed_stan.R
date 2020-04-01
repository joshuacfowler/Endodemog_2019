## Title: Grass endophyte population model with a bayesian framework
## Purpose: Creates survival kernel written in STAN with mixed effects, 
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

# survival data lists are generated in the endodemog_data_processing.R file, 
# within the section titled "Preparing datalists for Survival Kernel"
source("endodemog_data_processing.R")

#########################################################################################################
# GLMM for Surv ~ size_t + Endo + Origin + size_t*Endo with year and plot random effects------------------------------
#########################################################################################################
## run this code recommended to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)
## MCMC settings
ni <- 1000
nb <- 500
nc <- 3

# Stan model -------------
## here is the Stan model

sink("endodemog_surv.stan")
cat("
    data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> K;                       // number of predictors
    int<lower=0> nYear;                       // number of years
    int<lower=0, upper=11> year_t[N];         // year of observation
    int<lower=0> nEndo;                       // number of endo treatments
    int<lower=1, upper=2> endo_index[N];          // index for endophyte effect
    int<lower=0> nPlot;                         // number of plots
    int<lower=0> plot[N];                   // plot of observation
    int<lower=0, upper=1> surv_t1[N];      // plant survival at time t+1 and target variable (response)
    vector<lower=0>[N] logsize_t;             // plant size at time t (predictor)
    int<lower=0,upper=1> endo_01[N];            // plant endophyte status (predictor)
    int<lower=0,upper=1> origin_01[N];          // plant origin status (predictor)
    }

    parameters {
    vector[K] beta;                     // predictor parameters
    vector[nYear] tau_year[nEndo];      // random year effect
    real<lower=0> sigma_e[nEndo];        //year variance by endophyte effect
    vector[nPlot] tau_plot;        // random plot effect
    real<lower=0> sigma_p;          // plot variance effect
    }
    
    transformed parameters {
    real mu[N];                           // Linear Predictor

       for(n in 1:N){
    mu[n] = beta[1] + beta[2]*logsize_t[n] + beta[3]*endo_01[n] +beta[4]*origin_01[n]
    + tau_year[endo_index[n],year_t[n]] 
    + tau_plot[plot[n]];
    }
    }
    
    model {
    // Priors
    beta ~ normal(0,100);      // prior for predictor intercepts
    tau_plot ~ normal(0,sigma_p);   // prior for plot random effects
    to_vector(tau_year[1]) ~ normal(0,sigma_e[1]);   // prior for E- year random effects
    to_vector(tau_year[2]) ~ normal(0,sigma_e[2]);   // prior for E+ year random effects
    
    // Likelihood
      surv_t1 ~ bernoulli_logit(mu);
    }
    
    generated quantities{
    }
  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_surv.stan")



## Run the model by calling stan()
## and save the output to .rds files so that they can be called laters

smAGPE <- stan(file = "endodemog_surv.stan", data = AGPE_surv_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smAGPE, file = "endodemog_surv_AGPE.rds")

smELRI <- stan(file = "endodemog_surv.stan", data = ELRI_surv_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELRI, file = "endodemog_surv_ELRI.rds")

smELVI <- stan(file = "endodemog_surv.stan", data = ELVI_surv_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELVI, file = "endodemog_surv_ELVI.rds")

smFESU <- stan(file = "endodemog_surv.stan", data = FESU_surv_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smFESU, file = "endodemog_surv_FESU.rds")

smLOAR <- stan(file = "endodemog_surv.stan", data = LOAR_surv_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smLOAR, file = "endodemog_surv_LOAR.rds")

smPOAL <- stan(file = "endodemog_surv.stan", data = POAL_surv_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOAL, file = "endodemog_surv_POAL.rds")

smPOSY <- stan(file = "endodemog_surv.stan", data = POSY_surv_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOSY, file = "endodemog_surv_POSY.rds")

print(sm)
summary(sm)
print(sm, pars = "sigma_e")





## to read in model output without rerunning models
smAGPE <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_AGPE.rds")
smELRI <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_ELRI.rds")
smELVI <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_ELVI.rds")
smFESU <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_FESU.rds")
smLOAR <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_LOAR.rds")
smPOAL <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_POAL.rds")
smPOSY <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_surv_POSY.rds")


#########################################################################################################
###### Perform posterior predictive checks and assess model convergence-------------------------
#########################################################################################################
params <- c("beta", "tau_year[1,1]", "tau_plot[1]", "sigma_e", "sigma_p")

##### POAL - survival
print(smPOAL)
# summary(smPOAL)

## plot traceplots of chains for select parameters
traceplot(smAGPE, pars = params)
traceplot(smELRI, pars = params)
traceplot(smELVI, pars = params)
traceplot(smFESU, pars = params)
traceplot(smLOAR, pars = params)
traceplot(smPOAL, pars = params)
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
  post <- extract(fit)
  mu <- post$mu
  yrep <- matrix(nrow = n_post_draws, ncol = data$N)
  for(n in 1:n_post_draws){
    yrep[n,] <- rbinom(n = data$N, size = 1, prob = invlogit(mu[n,]))
  }
  out <- list(yrep, mu)
  names(out) <- c("yrep", "mu")
  return(out)
}

# apply the function for each species
AGPE_surv_yrep <- prediction(data = AGPE_surv_data_list, fit = smAGPE, n_post_draws = 500)
ELRI_surv_yrep <- prediction(data = ELRI_surv_data_list, fit = smELRI, n_post_draws = 500)
ELVI_surv_yrep <- prediction(data = ELVI_surv_data_list, fit = smELVI, n_post_draws = 500)
FESU_surv_yrep <- prediction(data = FESU_surv_data_list, fit = smFESU, n_post_draws = 500)
LOAR_surv_yrep <- prediction(data = LOAR_surv_data_list, fit = smLOAR, n_post_draws = 500)
POAL_surv_yrep <- prediction(data = POAL_surv_data_list, fit = smPOAL, n_post_draws = 500)
POSY_surv_yrep <- prediction(data = POSY_surv_data_list, fit = smPOSY, n_post_draws = 500)


# overlay 100 replicates over the actual dataset
ppc_dens_overlay( y = AGPE_surv_data_list$surv_t1, yrep = AGPE_surv_yrep$yrep[1:100,]) + ggtitle("AGPE")

ppc_dens_overlay( y = ELRI_surv_data_list$surv_t1, yrep = ELRI_surv_yrep$yrep[1:500,]) + ggtitle("ELRI")

ppc_dens_overlay( y = ELVI_surv_data_list$surv_t1, yrep = ELVI_surv_yrep$yrep[1:100,]) + ggtitle("ELVI")

ppc_dens_overlay( y = FESU_surv_data_list$surv_t1, yrep = FESU_surv_yrep$yrep[1:100,]) + ggtitle("FESU")

ppc_dens_overlay( y = LOAR_surv_data_list$surv_t1, yrep = LOAR_surv_yrep$yrep[1:100,]) + ggtitle("LOAR")

ppc_dens_overlay( y = POAL_surv_data_list$surv_t1, yrep = POAL_surv_yrep$yrep[1:100,]) + ggtitle("POAL")

ppc_dens_overlay( y = POSY_surv_data_list$surv_t1, yrep = POSY_surv_yrep$yrep[1:100,]) + ggtitle("POSY")




# Pairs plots to diagnose sampling

pairs(smPOAL, pars = "beta")


