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
ni <-10000
nb <- 5000
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
    real<lower=0> reciprocal_phi;            // inverse dispersion parameter
    }
    
   transformed parameters{
    real<lower=0> phi;                    // negative binomial dispersion parameter
    real mu[N];                         // Linear Predictor
    
    phi = 1. / reciprocal_phi;

    
        for(n in 1:N){
       mu[n] = beta[1] + beta[2]*logsize_t[n] + beta[3]*endo_01[n] +beta[4]*origin_01[n]
       + beta[5]*logsize_t[n]*endo_01[n] 
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
    reciprocal_phi ~ cauchy(0., 5.);


    // Likelihood 
      spike_t ~ neg_binomial_2_log(mu, phi);
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
# saveRDS(smAGPE, file = "endodemog_spike_AGPE_withplot.rds")

smELRI <- stan(file = "endodemog_spike.stan", data = ELRI_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELRI, file = "endodemog_spike_ELRI_withplot.rds")

smELVI <- stan(file = "endodemog_spike.stan", data = ELVI_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELVI, file = "endodemog_spike_ELVI_withplot.rds")

smFESU <- stan(file = "endodemog_spike.stan", data = FESU_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smFESU, file = "endodemog_spike_FESU_withplot.rds")

smLOAR <- stan(file = "endodemog_spike.stan", data = LOAR_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smLOAR, file = "endodemog_spike_LOAR_withplot.rds")

smPOAL <- stan(file = "endodemog_spike.stan", data = POAL_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOAL, file = "endodemog_spike_POAL_withplot.rds")

smPOSY <- stan(file = "endodemog_spike.stan", data = POSY_spike_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOSY, file = "endodemog_spike_POSY_withplot.rds")



#########################################################################################################
###### Perform posterior predictive checks and assess model convergence-------------------------
#########################################################################################################
params <- c("beta[1]", "beta[2]", "tau_year[1,1]", "sigma_e[1]", "sigma_e[2]")


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
prediction<- function(x, fit, reps) {
  post <- extract(fit)
  beta_post <- post$beta
  tau_plot_post <- post$tau_plot
  tau_year_post <- post$tau_year
  dim(tau_year_post) <- c(15000,22)
  lin_comb <- matrix(nrow = x$N, ncol = reps)
  yrep <- matrix(nrow = x$N, ncol = reps)
  for(n in 1:reps){
    lin_comb[,n] <- sample(beta_post[,1], size = x$N)+ x$logsize_t*sample(beta_post[,2], size = x$N) 
    + x$endo_01*sample(beta_post[,3], size = x$N) + x$origin_01*sample(beta_post[,4], size = x$N) 
    + x$logsize_t*x$endo_01*sample(beta_post[,5], size = x$N) 
    + x$plot_Xs[]*sample(tau_plot_post[], size = x$N) 
    + x$yearendo_Xs[]*sample(tau_year_post[], size = x$N)
    
    prob <- exp(lin_comb)
    yrep[,n] <- rnbinom(prob = prob[,n], n = x$N, size = mean(post$phi))
    print(paste("rep", n, "of", reps))
    
  }
  out <- list(yrep, prob, lin_comb) 
  names(out) <- c("yrep", "prob", "lin_comb")
  
  return(out)
}



# apply the function for each species
AGPE_spike_yrep <- prediction(AGPE_spike_data_list, smAGPE, 500)
ELRI_spike_yrep <- prediction(ELRI_spike_data_list, smELRI, 500)
ELVI_spike_yrep <- prediction(ELVI_spike_data_list, smELVI, 500)
FESU_spike_yrep <- prediction(FESU_spike_data_list, smFESU, 500)
LOAR_spike_yrep <- prediction(LOAR_spike_data_list, smLOAR, 500)
POAL_spike_yrep <- prediction(POAL_spike_data_list, smPOAL, 500)
POSY_spike_yrep <- prediction(POSY_spike_data_list, smPOSY, 500)



# overlay 100 replicates over the actual dataset
mspike_yrep_AGPE <- t(AGPE_spike_yrep$yrep)
ppc_dens_overlay( y = AGPE_spike_data_list$spike_t, yrep = mspike_yrep_AGPE[1:100,])+ xlab("prob. of y") + ggtitle("AGPE")

mspike_yrep_ELRI <- t(ELRI_spike_yrep$yrep)
ppc_dens_overlay( y = ELRI_spike_data_list$spike_t, yrep = mspike_yrep_ELRI[1:100,])+ xlab("prob. of y") + ggtitle("ELRI")

mspike_yrep_ELVI <- t(ELVI_spike_yrep$yrep)
ppc_dens_overlay( y = ELVI_spike_data_list$spike_t, yrep = mspike_yrep_ELVI[1:100,]) + xlab("prob. of y") + ggtitle("ELVI")

mspike_yrep_FESU <- t(FESU_spike_yrep$yrep)
ppc_dens_overlay( y = FESU_spike_data_list$spike_t, yrep = mspike_yrep_FESU[1:100,]) + xlab("prob. of y") + ggtitle("FESU")

mspike_yrep_LOAR <- t(LOAR_spike_yrep$yrep)
ppc_dens_overlay( y = LOAR_spike_data_list$spike_t, yrep = mspike_yrep_LOAR[1:100,]) + xlab("prob. of y") + ggtitle("LOAR")

mspike_yrep_POAL <- t(POAL_spike_yrep$yrep)
ppc_dens_overlay( y = POAL_spike_data_list$spike_t, yrep = mspike_yrep_POAL[1:100,]) + xlab("prob. of y") + ggtitle("POAL")

mspike_yrep_POSY <- t(POSY_spike_yrep$yrep)
ppc_dens_overlay( y = POSY_spike_data_list$spike_t, yrep = mspike_yrep_POSY[1:100,]) + xlab("prob. of y") + ggtitle("POSY")








