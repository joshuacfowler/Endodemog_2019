## Title: Grass endophyte population model with a bayesian framework
## Purpose: Fits all vital rate models, written in STAN
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
# GLMM for all Vital Rates ~ size_t + Endo + Origin with year and plot random effects------------------------------
#########################################################################################################
## run this code recommended to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)
## MCMC settings
ni <- 10000
nb <- 5000
nc <- 3


# Stan model -------------
## here is the Stan model

sink("endodemog_all_vr.stan")
cat("
    data { 
    // indices
    int<lower=0> nYear;                       // number of years
    int<lower=0> nPlot;                         // number of plots
    int<lower=0> nEndo;                       // number of endo treatments

    // surv data
    int<lower=0> N_s;                       // number of observations for surv model
    int<lower=0> K_s;                       // number of predictors for surv model
    int<lower=0, upper=11> year_t_s[N_s];         // year of observation for surv model
    int<lower=1, upper=2> endo_index_s[N_s];          // index for endophyte effect for surv model
    int<lower=0> plot_s[N_s];                   // plot of observation for surv model
    int<lower=0, upper=1> surv_t1[N_s];      // plant survival at time t+1
    vector<lower=0>[N_s] logsize_t_s;             // plant size at time t for surv model
    int<lower=0,upper=1> endo_01_s[N_s];            // plant endophyte status for surv model
    int<lower=0,upper=1> origin_01_s[N_s];          // plant origin status for surv model
    
    // growth data
    int<lower=0> N_g;                       // number of observations for growth model
    int<lower=0> K_g;                       // number of predictors for growth model
    int<lower=0> lowerlimit_g;                         //lower limit for truncated negative binomial
    int<lower=0, upper=11> year_t_g[N_g];         // year of observation for growth model
    int<lower=1, upper=2> endo_index_g[N_g];          // index for endophyte effect for growth model
    int<lower=0> plot_g[N_g];                   // plot of observation for growth model
    int<lower=lowerlimit_g> size_t1[N_g];      // plant size at time t+1
    vector<lower=0>[N_g] logsize_t_g;             // plant size at time t for growth model
    int<lower=0,upper=1> endo_01_g[N_g];            // plant endophyte status for growth model
    int<lower=0,upper=1> origin_01_g[N_g];          // plant origin status for growth model
    
    // flowering data
    int<lower=0> N_fl;                       // number of observations
    int<lower=0> K_fl;                       // number of predictors
    int<lower=0> year_t_fl[N_fl];                      // year of observation
    int<lower=0> plot_fl[N_fl];                      // plot of observation
    int<lower=0, upper=1> flw_t[N_fl];      // flowering status at time t
    vector<lower=-1>[N_fl] logsize_t_fl;                  // log of plant size at time t 
    int<lower=0, upper=1> endo_01_fl[N_fl];            // endophyte status 
    int<lower=1, upper=2> endo_index_fl[N_fl];        // index for endophyte effect
    int<lower=0, upper=1> origin_01_fl[N_fl];            // origin status
    
    // # of flw tiller data
    int<lower=0> N_ft;                       // number of observations
    int<lower=0> K_ft;                       // number of predictors
    int<lower=0> lowerlimit_ft;                         //lower limit for truncated negative binomial
    int<lower=0> year_t_ft[N_ft];                      // year of observation
    int<lower=1, upper=2> endo_index_ft[N_ft];          // index for endophyte effect
    int<lower=0> plot_ft[N_ft];                   // plot of observation
    int<lower=lowerlimit_ft> flw_count_t[N_ft];      // plant size at time t and target variable 
    vector<lower=0>[N_ft] logsize_t_ft;             // plant size at time t
    int<lower=0,upper=1> endo_01_ft[N_ft];            // plant endophyte status 
    int<lower=0,upper=1> origin_01_ft[N_ft];          // plant origin status 
  
    // spikelet/infl data
    int<lower=0> N_sp;                       // number of observations
    int<lower=0> K_sp;                       // number of predictors
    int<lower=0> year_t_sp[N_sp];                      // year of observation
    int<lower=1, upper=2> endo_index_sp[N_sp];          // index for endophyte effect
    int<lower=0> plot_sp[N_sp];                   // plot of observation
    int<lower=0> spike_t[N_sp];              // no. of spike per infl at time t
    vector<lower=0>[N_sp] logsize_t_sp;             // plant size at time t 
    int<lower=0,upper=1> endo_01_sp[N_sp];            // plant endophyte status
    int<lower=0,upper=1> origin_01_sp[N_sp];          // plant origin status 
    
    // seed/spiklet data
    int<lower=0> N_se;                       // number of observations of seed/spikelet
    int<lower=0> K_se;                    // number of predictors
    real<lower=0> seed[N_se];               // number of seeds per spikelet
    int<lower=0,upper=1> endo_01_se[N_se];            // plant endophyte status

    }

    parameters {
    // surv params
    vector[K_s] beta_s;                     // predictor parameters
    vector[nYear] tau_year_s[nEndo];      // random year effect
    real<lower=0> sigma_e_s[nEndo];        //year variance by endophyte effect
    vector[nPlot] tau_plot_s;        // random plot effect
    real<lower=0> sigma_p_s;          // plot variance effect
    
    // growth params
    vector[K_g] beta_g;                     // predictor parameters
    vector[nYear] tau_year_g[nEndo];      // random year effect
    real<lower=0> sigma_e_g[nEndo];        //year variance by endophyte effect
    vector[nPlot] tau_plot_g;        // random plot effect
    real<lower=0> sigma_p_g;          // plot variance
    real<lower=0> phi_g;            // dispersion parameter
    
    // flower params
    vector[K_fl] beta_fl;                     // predictor parameters
    vector[nYear] tau_year_fl[nEndo];      // random year effect
    real<lower=0> sigma_e_fl[nEndo];        //year variance by endophyte effect
    vector[nPlot] tau_plot_fl;        // random plot effect
    real<lower=0> sigma_p_fl;         // plot variance
    
    // fertility params
    vector[K_ft] beta_ft;                     // predictor parameters
    vector[nYear] tau_year_ft[nEndo];      // random year effect
    real<lower=0> sigma_e_ft[nEndo];        //year variance by endophyte effect
    vector[nPlot] tau_plot_ft;        // random plot effect
    real<lower=0> sigma_p_ft;          // plot variance
    real<lower=0> phi_ft;            // negative binomial dispersion parameter
    
    // spike/inf params
    vector[K_sp] beta_sp;                     // predictor parameters
    vector[nYear] tau_year_sp[nEndo];      // random year effect
    real<lower=0> sigma_e_sp[nEndo];        //year variance by endophyte effect
    vector[nPlot] tau_plot_sp;        // random plot effect
    real<lower=0> sigma_p_sp;          // plot variance
    real<lower=0> phi_sp;            // negative binomial dispersion parameter
    
    // seed/spike params
    vector[K_se] beta_se;                     // predictor parameters
    real<lower=0> sigma_se;         // seed per spikelet variance
    }
    
    transformed parameters {
    real mu_s[N_s];                           // surv Linear Predictor
    real mu_g[N_g];                           // growth Linear Predictor
    real mu_fl[N_fl];                           // flowering Linear Predictor
    real mu_ft[N_ft];                           // # of flowering tillers Linear Predictor
    real mu_sp[N_sp];                           // # of spikelets/infl Linear Predictor
    real mu_se[N_se];                           // mean seed per spikelet


       for(n in 1:N_s){
    mu_s[n] = beta_s[1] + beta_s[2]*logsize_t_s[n] + beta_s[3]*endo_01_s[n] +beta_s[4]*origin_01_s[n]
    + tau_year_s[endo_index_s[n],year_t_s[n]] 
    + tau_plot_s[plot_s[n]];
       }
       for(n in 1:N_g){
    mu_g[n] = beta_g[1] + beta_g[2]*logsize_t_g[n] + beta_g[3]*endo_01_g[n] +beta_g[4]*origin_01_g[n]
    + tau_year_g[endo_index_g[n],year_t_g[n]] 
    + tau_plot_g[plot_g[n]];
       }
       for(n in 1:N_fl){
    mu_fl[n] = beta_fl[1] + beta_fl[2]*logsize_t_fl[n] + beta_fl[3]*endo_01_fl[n] +beta_fl[4]*origin_01_fl[n]
    + tau_year_fl[endo_index_fl[n],year_t_fl[n]] 
    + tau_plot_fl[plot_fl[n]];
       }
       for(n in 1:N_ft){
    mu_ft[n] = beta_ft[1] + beta_ft[2]*logsize_t_ft[n] + beta_ft[3]*endo_01_ft[n] +beta_ft[4]*origin_01_ft[n]
    + tau_year_ft[endo_index_ft[n],year_t_ft[n]] 
    + tau_plot_ft[plot_ft[n]];
       }
       for(n in 1:N_sp){
    mu_sp[n] = beta_sp[1] + beta_sp[2]*logsize_t_sp[n] + beta_sp[3]*endo_01_sp[n] +beta_sp[4]*origin_01_sp[n]
    + tau_year_sp[endo_index_sp[n],year_t_sp[n]] 
    + tau_plot_sp[plot_sp[n]];
       }
       for(n in 1:N_se){
    mu_se[n] = beta_se[1] + beta_se[2]*endo_01_se[n];
       }
    }
           
    model {
    // Priors
    // surv priors
    beta_s ~ normal(0,100);      // prior for predictor intercepts
    tau_plot_s ~ normal(0,sigma_p_s);   // prior for plot random effects
    to_vector(tau_year_s[1]) ~ normal(0,sigma_e_s[1]);   // prior for E- year random effects
    to_vector(tau_year_s[2]) ~ normal(0,sigma_e_s[2]);   // prior for E+ year random effects
    // growth priors
    beta_g ~ normal(0,100);      // prior for predictor intercepts
    tau_plot_g ~ normal(0,sigma_p_g);   // prior for plot random effects
    to_vector(tau_year_g[1]) ~ normal(0,sigma_e_g[1]);   // prior for E- year random effects
    to_vector(tau_year_g[2]) ~ normal(0,sigma_e_g[2]);   // prior for E+ year random effects
    phi_g ~ cauchy(0., 5.);
    // flw priors
    beta_fl ~ normal(0,100);      // prior for predictor intercepts
    tau_plot_fl ~ normal(0,sigma_p_fl);   // prior for plot random effects
    to_vector(tau_year_fl[1]) ~ normal(0,sigma_e_fl[1]);   // prior for E- year random effects
    to_vector(tau_year_fl[2]) ~ normal(0,sigma_e_fl[2]);   // prior for E+ year random effects
    // # of flw priors
    beta_ft ~ normal(0,100);      // prior for predictor intercepts
    tau_plot_ft ~ normal(0,sigma_p_ft);   // prior for plot random effects
    to_vector(tau_year_ft[1]) ~ normal(0,sigma_e_ft[1]);   // prior for E- year random effects
    to_vector(tau_year_ft[2]) ~ normal(0,sigma_e_ft[2]);   // prior for E+ year random effects
    phi_ft ~ cauchy(0., 5.);
    // # of spike priors
    beta_sp ~ normal(0,100);
    tau_plot_sp ~ normal(0,sigma_p_sp);   // prior for plot random effects
    to_vector(tau_year_sp[1]) ~ normal(0,sigma_e_sp[1]);   // prior for E- year random effects
    to_vector(tau_year_sp[2]) ~ normal(0,sigma_e_sp[2]);   // prior for E+ year random effects
    phi_sp ~ cauchy(0., 5.);

    // seed mean priors
    beta_se ~ normal(0,100);

    // Likelihoods
      // surv
      surv_t1 ~ bernoulli_logit(mu_s);
      // growth
        for(n in 1:N_g){
      size_t1[n] ~ neg_binomial_2_log(mu_g[n],phi_g);
      target += -log1m(neg_binomial_2_log_lpmf(lowerlimit_g | mu_g[n], phi_g)); // manually adjusting computation of likelihood because T[,] truncation syntax doesn't compile for neg binomial
        }
      //flowering  
      flw_t ~ bernoulli_logit(mu_fl);
      // # of flow
        for(n in 1:N_ft){
      flw_count_t[n] ~ neg_binomial_2_log(mu_ft[n],phi_ft);
      target += -log1m(neg_binomial_2_log_lpmf(lowerlimit_ft | mu_ft[n], phi_ft)); // manually adjusting computation of likelihood because T[,] truncation syntax doesn't compile for neg binomial
        }
      // # of spikelet/infl  
      spike_t ~ neg_binomial_2_log(mu_sp,phi_sp);
      // mean seed/spike
      seed ~ normal(mu_se,sigma_se);
    }
    
    generated quantities{
    }
  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_all_vr.stan")


## Run the model by calling stan()
## and save the output to .rds files so that they can be called laters

smAGPE <- stan(file = "endodemog_all_vr.stan", data = AGPE_all_vr_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smAGPE, file = "endodemog_all_vr_AGPE.rds")

smELRI <- stan(file = "endodemog_all_vr.stan", data = ELRI_all_vr_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELRI, file = "endodemog_all_vr_ELRI.rds")

smELVI <- stan(file = "endodemog_all_vr.stan", data = ELVI_all_vr_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELVI, file = "endodemog_all_vr_ELVI.rds")

smFESU <- stan(file = "endodemog_all_vr.stan", data = FESU_all_vr_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smFESU, file = "endodemog_all_vr_FESU.rds")

smLOAR <- stan(file = "endodemog_all_vr.stan", data = LOAR_all_vr_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smLOAR, file = "endodemog_all_vr_LOAR.rds")

smPOAL <- stan(file = "endodemog_all_vr.stan", data = POAL_all_vr_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOAL, file = "endodemog_all_vr_POAL.rds")

smPOSY <- stan(file = "endodemog_all_vr.stan", data = POSY_all_vr_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOSY, file = "endodemog_all_vr_POSY.rds")

print(sm)
summary(sm)
print(sm, pars = "sigma_e")





## to read in model output without rerunning models
smAGPE <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_all_vr_AGPE.rds")
smELRI <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_all_vr_ELRI.rds")
smELVI <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_all_vr_ELVI.rds")
smFESU <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_all_vr_FESU.rds")
smLOAR <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_all_vr_LOAR.rds")
smPOAL <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_all_vr_POAL.rds")
smPOSY <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_all_vr_POSY.rds")



