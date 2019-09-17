## Title: Grass endophyte population model with a bayesian framework
## Purpose: Creates growth kernel written in STAN with mixed effects, 
## and does visualisation of posterior predictive checks
## Authors: Joshua and Tom
#############################################################

library(tidyverse)
library(rstan)
library(StanHeaders)
library(shinystan)
library(bayesplot)
library(devtools)


#############################################################################################
####### Data manipulation to prepare data as lists for Stan models------------------
#############################################################################################

# Growth data lists are generated in the endodemog_data_processing.R file
# within the section titled "Preparing datalists for Growth Kernel"
source("endodemog_data_processing.R")

#########################################################################################################
# GLMM for size_t1 ~ size_t + Endo + Origin + size_t*Endo with year and plot random effects------------------------------
#########################################################################################################
## run this code recommended to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <-10000
nb <- 5000
nc <- 3

# Stan model -------------
## here is the Stan model with a model matrix and species effects ##



sink("endodemog_grow.stan")
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
    int<lower=lowerlimit> size_t1[N];      // plant size at time t+1 and target variable (response)
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
    phi = 1. / reciprocal_phi;
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
    beta ~ normal(0,100);      // prior for predictor intercepts
    tau_plot ~ normal(0,sigma_p);   // prior for plot random effects
    to_vector(tau_year[1]) ~ normal(0,sigma_e[1]);   // prior for E- year random effects
    to_vector(tau_year[2]) ~ normal(0,sigma_e[2]);   // prior for E+ year random effects
    reciprocal_phi ~ cauchy(0., 5.);


    // Likelihood
    for(n in 1:N){
      size_t1[n] ~ neg_binomial_2_log(mu[n],phi);
      target += -log1m(neg_binomial_2_log_lpmf(lowerlimit | mu[n], phi)); // manually adjusting computation of likelihood because T[,] truncation syntax doesn't compile for neg binomial
    }
    }
    
   generated quantities{
      vector[N] mu;
    
    // for posterior predictive check
    for(n in 1:N){
      mu[n] = beta[1] + beta[2]*logsize_t[n] + beta[3]*endo_01[n] +beta[4]*origin_01[n]
      + beta[5]*logsize_t[n]*endo_01[n] 
      + tau_year[endo_index[n],year_t[n]]
      + tau_plot[plot[n]];
    }
   }
  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_grow.stan")



## Run the model by calling stan()
## and save the output to .rds files so that they can be called laters

smAGPE <- stan(file = "endodemog_grow.stan", data = AGPE_grow_data_list,
           iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smAGPE, file = "endodemog_grow_AGPE.rds")

smELRI <- stan(file = "endodemog_grow.stan", data = ELRI_grow_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELRI, file = "endodemog_grow_ELRI.rds")

smELVI <- stan(file = "endodemog_grow.stan", data = ELVI_grow_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELVI, file = "endodemog_grow_ELVI.rds")

smFESU <- stan(file = "endodemog_grow.stan", data = FESU_grow_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smFESU, file = "endodemog_grow_FESU.rds")

smLOAR <- stan(file = "endodemog_grow.stan", data = LOAR_grow_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smLOAR, file = "endodemog_grow_LOAR.rds")

smPOAL <- stan(file = "endodemog_grow.stan", data = POAL_grow_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOAL, file = "endodemog_grow_POAL.rds")

smPOSY <- stan(file = "endodemog_grow.stan", data = POSY_grow_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOSY, file = "endodemog_grow_POSY.rds")





print(sm)
summary(sm)
print(sm, pars = "tau_year")


## to read in model output without rerunning models
smPOAL <- readRDS(file = "model_run_MAR7/endodemog_grow_POAL_withplot.rds")
smPOSY <- readRDS(file = "model_run_MAR7/endodemog_grow_POSY_withplot.rds")
smLOAR <- readRDS(file = "model_run_MAR7/endodemog_grow_LOAR_withplot.rds")
smELVI <- readRDS(file = "model_run_MAR7/endodemog_grow_ELVI_withplot.rds")
smELRI <- readRDS(file = "model_run_MAR7/endodemog_grow_ELRI_withplot.rds")
smFESU <- readRDS(file = "model_run_MAR7/endodemog_grow_FESU_withplot.rds")
smAGPE <- readRDS(file = "model_run_MAR7/endodemog_grow_AGPE_withplot.rds")




#########################################################################################################
###### Perform posterior predictive checks and assess model convergence-------------------------
#########################################################################################################
params <- c("beta[1]", "beta[2]", "tau_year[1,1]", "sigma_e[1]", "sigma_e[2]")


##### POAL - growth
print(smPOAL)

#extract posteriors into dataframe
posterior <- as.data.frame(smPOAL)

yrep <- gen_est(POAL_surv_data_list)



y_rep <- apply(posterior, MARGIN = 1:2, FUN = function(draw) {
  rnorm(n, mean = draw[grepl("mu", names(draw))], sd = "sigma")
})


# Extract Entire posterior (for all parameters) - 3 chains
posterior <- extract(smPOAL, inc_warmup = FALSE, permute = FALSE)

# Generate new data

y_rep <- apply(posterior, MARGIN = 1:2, FUN = function(draw) {
  rnorm(n, mean = draw[grepl("^mu", names(draw))], sd = draw["sigma"])
})


dim(y_rep) # 16 replicates by 5000 iterations by 3 chains
error_rep <- y - y_rep # y = 16 replications
mean(error_rep^2) # far greater than sigma2 (92 vs true = 25) - also differs from estimated sigma2 of ~40

# alternative way to extract mu MCMC estimates
mu_rep <- apply(posterior, MARGIN = 1:2, FUN = function(draw) {
  mu1 = draw[grepl("^mu", names(draw))]
})

# Fit statistics that are done within WinBUGS
residual <- y - mu_rep
sq <- residual^2
sq.new <- (y_rep - mu_rep)^2

for(i in 1:5000){
  fit[i] <- sum(sq[ , i, 1]) # Not sure how to get this to loop right for 3 chains - might be easier with permuted extract function
}

for(i in 1:5000) fit.new[i] <- sum(sq.new[, i, 1])

# Posterior predictive distributions and bayesian p-values
test <- ifelse(test = (fit.new - fit) > 0, yes=1, no = 0) # Test whether new data set more extreme
(bpvalue <- mean(test)) # Bayesian p-value = 0.0358 indicates poor fit (should be ~0.5)

plot(fit, fit.new) 
abline(0, 1) # Also shows poor fit - bias of idealized (new) fit data being lower than fit data. The SSQ are also MUCH larger than those found in WinBUGS. Seems like I did something incorrectly since the summary output was the same for WinBUGS and Stan.
