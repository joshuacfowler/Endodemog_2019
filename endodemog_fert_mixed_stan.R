## Title: Grass endophyte population model with a bayesian framework
## Purpose: Creates fertility (flowering tiller #) kernel written in STAN with mixed effects, 
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
# Fertility data lists are generated in the endodemog_data_processing.R file
# within the section titled "Preparing datalists for Fertility Kernel"
source("endodemog_data_processing.R")

#########################################################################################################
# GLMM for flw_t1~ size +Endo + Origin  with year random effects------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <-1000
nb <- 500
nc <- 3

# Stan model -------------
## here is the Stan model with a model matrix and species effects ##



sink("endodemog_fert.stan")
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
    int<lower=lowerlimit> flw_t[N];      // plant size at time t+1 and target variable (response)
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
    // Linear Predictor
    real mu[N];
    
    for(n in 1:N){
      mu[n] = beta[1] + beta[2]*logsize_t[n] + beta[3]*endo_01[n] +beta[4]*origin_01[n]
      + beta[5]*logsize_t[n]*endo_01[n] 
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
    phi ~ cauchy(0., 5.);

    // Likelihood
    for(n in 1:N){
      flw_t[n] ~ neg_binomial_2_log(mu[n],phi);
      target += -log1m(neg_binomial_2_log_lpmf(0 | mu[n], phi)); // manually adjusting computation of likelihood because T[,] truncation syntax doesn't compile for neg binomial
    }
    }
    
   generated quantities{
   }
  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_fert.stan")

## Run the model by calling stan()
## Save the outputs as rds files

smAGPE<- stan(file = "endodemog_fert.stan", data = AGPE_fert_data_list,
              iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smAGPE, file = "endodemog_fert_AGPE_withplot.rds")

smELRI <- stan(file = "endodemog_fert.stan", data = ELRI_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELRI, file = "endodemog_fert_ELRI_withplot.rds")

smELVI <- stan(file = "endodemog_fert.stan", data = ELVI_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELVI, file = "endodemog_fert_ELVI_withplot.rds")

smFESU <- stan(file = "endodemog_fert.stan", data = FESU_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smFESU, file = "endodemog_fert_FESU_withplot.rds")

smLOAR <- stan(file = "endodemog_fert.stan", data = LOAR_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smLOAR, file = "endodemog_fert_LOAR_withplot.rds")

smPOAL <- stan(file = "endodemog_fert.stan", data = POAL_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOAL, file = "endodemog_fert_POAL_withplot.rds")

smPOSY <- stan(file = "endodemog_fert.stan", data = POSY_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOSY, file = "endodemog_fert_POSY_withplot.rds")












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
prediction <- function(data, fit, n_post_draws){
  post <- extract(fit)
  mu <- post$mu
  phi <- post$phi
  yrep <- matrix(nrow = n_post_draws, ncol = data$N)
  for(n in 1:n_post_draws){
    yrep[n,] <- rztnbinom(n = data$N, size = phi[n], mu = exp(mu[n,]))
  }
  out <- list(yrep, mu)
  names(out) <- c("yrep", "mu")
  return(out)
}


# apply the function for each species
AGPE_fert_yrep <- prediction(data = AGPE_fert_data_list, fit = smAGPE, n_post_draws = 500)
ELRI_fert_yrep <- prediction(data = ELRI_fert_data_list, fit = smELRI, n_post_draws = 500)
ELVI_fert_yrep <- prediction(data = ELVI_fert_data_list, fit = smELVI, n_post_draws = 500)
FESU_fert_yrep <- prediction(data = FESU_fert_data_list, fit = smFESU, n_post_draws = 500)
LOAR_fert_yrep <- prediction(data = LOAR_fert_data_list, fit = smLOAR, n_post_draws = 500)
POAL_fert_yrep <- prediction(data = POAL_fert_data_list, fit = smPOAL, n_post_draws = 500)
POSY_fert_yrep <- prediction(data = POSY_fert_data_list, fit = smPOSY, n_post_draws = 500)



# overlay 100 replicates over the actual dataset
ppc_dens_overlay( y = AGPE_fert_data_list$flw_t, yrep = AGPE_fert_yrep$yrep[1:100,]) + xlim(0,40) + xlab("prob. of y") + ggtitle("AGPE")

ppc_dens_overlay( y = ELRI_fert_data_list$flw_t, yrep = ELRI_fert_yrep$yrep[1:100,])+ xlab("prob. of y") + ggtitle("ELRI")

ppc_dens_overlay( y = ELVI_fert_data_list$flw_t, yrep = ELVI_fert_yrep$yrep[1:100,]) + xlab("prob. of y") + ggtitle("ELVI")

ppc_dens_overlay( y = FESU_fert_data_list$flw_t, yrep = FESU_fert_yrep$yrep[1:100,]) + xlab("prob. of y") + ggtitle("FESU")

ppc_dens_overlay( y = LOAR_fert_data_list$flw_t, yrep = LOAR_fert_yrep$yrep[1:100,])+ xlim(0,20) + xlab("prob. of y") + ggtitle("LOAR")

ppc_dens_overlay( y = POAL_fert_data_list$flw_t, yrep = POAL_fert_yrep$yrep[1:100,]) + xlim(0,40) + xlab("prob. of y") + ggtitle("POAL")

ppc_dens_overlay( y = POSY_fert_data_list$flw_t, yrep = POSY_fert_yrep$yrep[1:100,]) + xlab("prob. of y") + ggtitle("POSY")





















# Working on a function to get fertility estimates to use with the seed means estimates to get estimates of seed production. 



get_est <- function(post, data, size){
  df <- as_data_frame(cbind(data$year_t, data$endo_index, data$plot, data$flw_t, data$logsize_t, data$origin_01))
  df <- setNames(object = df, nm = c("year_t", "endo_index", "plot", "flw_t", "logsize_t", "origin_01"))
  fixedeff <- sample(post$`beta[1]`,size = size)+ sample(post$`beta[2]`,size = size)*data$logsize_t + sample(post$`beta[3]`,size = size)*data$endo_01 + sample(post$`beta[4]`,size = size)*data$origin_01 + sample(post$`beta[5]`,size = size)*data$logsize_t*data$endo_01
  year_rand_eff <- ifelse(df$year_t == 1 & df$endo_index == 1, sample(post$`tau_year[1,1]`,size = size),ifelse(df$year_t == 1 & df$endo_index == 2, sample(post$`tau_year[2,1]`,size = size),ifelse(df$year_t == 2 & df$endo_index == 1, sample(post$`tau_year[1,2]`,size = size),ifelse(df$year_t == 2 & df$endo_index == 2, sample(post$`tau_year[2,2]`,size = size),ifelse(df$year_t == 3 & df$endo_index == 1, sample(post$`tau_year[1,3]`,size = size),ifelse(df$year_t == 3 & df$endo_index == 2, sample(post$`tau_year[2,3]`,size = size),ifelse(df$year_t == 4 & df$endo_index == 1, sample(post$`tau_year[1,4]`,size = size),ifelse(df$year_t == 4 & df$endo_index == 2, sample(post$`tau_year[2,4]`,size = size),ifelse(df$year_t == 5 & df$endo_index == 1, sample(post$`tau_year[1,5]`,size = size),ifelse(df$year_t == 5 & df$endo_index == 2, sample(post$`tau_year[2,5]`,size = size),ifelse(df$year_t == 6 & df$endo_index == 1, sample(post$`tau_year[1,6]`,size = size),ifelse(df$year_t == 6 & df$endo_index == 2, sample(post$`tau_year[2,6]`,size = size),ifelse(df$year_t == 7 & df$endo_index == 1, sample(post$`tau_year[1,7]`,size = size),ifelse(df$year_t == 7 & df$endo_index == 2, sample(post$`tau_year[2,7]`,size = size),ifelse(df$year_t ==8 & df$endo_index == 1, sample(post$`tau_year[1,8]`,size = size),ifelse(df$year_t == 8 & df$endo_index == 2, sample(post$`tau_year[2,8]`,size = size),ifelse(df$year_t == 9 & df$endo_index == 1, sample(post$`tau_year[1,9]`,size = size),ifelse(df$year_t == 9 & df$endo_index ==2, sample(post$`tau_year[2,9]`,size = size),ifelse(df$year_t == 10 & df$endo_index == 1, sample(post$`tau_year[1,10]`,size = size),ifelse(df$year_t ==10 & df$endo_index == 2, sample(post$`tau_year[2,10]`,size = size),ifelse(df$year_t == 11 & df$endo_index == 1, sample(post$`tau_year[1,11]`,size = size),ifelse(df$year_t == 11 & df$endo_index == 2, sample(post$`tau_year[2,11]`,size = size),ifelse(df$year_t == 12 & df$endo_index == 1, sample(post$`tau_year[1,12]`,size = size),ifelse(df$year_t == 12 & df$endo_index == 2, sample(post$`tau_year[2,12]`,size = size),NA))))))))))))))))))))))))
  plot_rand_eff <- ifelse(df$plot == 1, sample(post$`tau_plot[1]`,size = size), ifelse(df$plot == 2, sample(post$`tau_plot[2]`,size = size), ifelse(df$plot == 3, sample(post$`tau_plot[3]`,size = size), ifelse(df$plot == 4, sample(post$`tau_plot[4]`,size = size), ifelse(df$plot == 5, sample(post$`tau_plot[5]`,size = size), ifelse(df$plot == 6, sample(post$`tau_plot[6]`,size = size), ifelse(df$plot == 7, sample(post$`tau_plot[7]`,size = size), ifelse(df$plot == 8, sample(post$`tau_plot[8]`,size = size), ifelse(df$plot == 9, sample(post$`tau_plot[9]`,size = size), ifelse(df$plot == 10, sample(post$`tau_plot[10]`,size = size), ifelse(df$plot == 11, sample(post$`tau_plot[11]`,size = size), ifelse(df$plot == 12, sample(post$`tau_plot[12]`,size = size), ifelse(df$plot == 13, sample(post$`tau_plot[13]`,size = size), ifelse(df$plot == 14, sample(post$`tau_plot[14]`,size = size), ifelse(df$plot == 15, sample(post$`tau_plot[15]`,size = size), ifelse(df$plot == 16, sample(post$`tau_plot[16]`,size = size), ifelse(df$plot == 17, sample(post$`tau_plot[17]`,size = size), ifelse(df$plot == 18, sample(post$`tau_plot[18]`,size = size), ifelse(df$plot == 19, sample(post$`tau_plot[19]`,size = size), ifelse(df$plot == 20, sample(post$`tau_plot20]`,size = size) ,NA))))))))))))))))))))
  df$est<- as.integer(exp(fixedeff+year_rand_eff+plot_rand_eff))
  
  return(df)
}

AGPE_post <- as.data.frame(smAGPE)
AGPE_data_for_est <- LTREB_data %>% 
  mutate(flw_stat_t1 = as.integer(flw_stat_t1)) %>%
  mutate(flw_stat_t = as.integer(flw_stat_t)) %>% 
  mutate(flw_t = as.integer(flw_t)) %>% 
  mutate(flw_t1 = as.integer(flw_t1)) %>%
  filter(!is.na(logsize_t)) %>% 
  filter(logsize_t >= 0) %>% 
  filter(!is.na(endo_01)) %>% 
  filter(species == "AGPE") %>% 
  mutate(year_t = as.integer(as.factor(as.character(year_t_index)))) %>% 
  mutate(plot = as.integer(as.factor(as.integer(as.character(plot_fixed)))))

AGPE_est <-get_est(post = AGPE_post, data = AGPE_data_for_est, size = length(AGPE_data_for_est$flw_t))

View(AGPE_est)

plot(AGPE_est$logsize_t, AGPE_est$est)
points(AGPE_est$logsize_t, AGPE_est$flw_t, col = "blue")

ELRI_post <- as.data.frame(smELRI)

ELRI_est <- get_est(post = ELRI_post, data = ELRI_fert_data_list)
View(ELRI_est)

plot(ELRI_est$logsize_t, ELRI_est$est)
points(ELRI_est$logsize_t, ELRI_est$flw_t, col = "blue")

# Here is a simpler model for those species with little data


sink("endodemog_fert_woyearplot.stan")
cat("
    data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> K;                       // number of predictors
    int<lower=0> lowerlimit;                         //lower limit for truncated negative binomial
    int<lower=0> nEndo;                       // number of endo treatments
    int<lower=1, upper=2> endo_index[N];          // index for endophyte effect
    int<lower=lowerlimit> flw_t[N];      // plant size at time t+1 and target variable (response)
    vector<lower=0>[N] logsize_t;             // plant size at time t (predictor)
    int<lower=0,upper=1> endo_01[N];            // plant endophyte status (predictor)
    int<lower=0,upper=1> origin_01[N];          // plant origin status (predictor)
    }
    
    parameters {
    vector[3] beta;                     // predictor parameters

    real<lower=0> phi;            // neg. binomial disperion parameter
    }
    
    model {
    vector[N] mu;
    
    // Linear Predictor
    for(n in 1:N){
      mu[n] = beta[1] + beta[2]*logsize_t[n];
    }
    // Priors
    beta ~ normal(0,10);      // prior for predictor intercepts
    phi ~ gamma(2,0.1);
    // Likelihood
    for(n in 1:N){
      flw_t[n] ~ neg_binomial_2_log(mu[n],phi);
      target += -log1m(neg_binomial_2_log_lpmf(lowerlimit | mu[n], phi)); // manually adjusting computation of likelihood because T[,] truncation syntax doesn't compile for neg binomial
    }
    }
    
   generated quantities{
   }
  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_fert_woyearplot.stan")

smELRI2 <- stan(file = "endodemog_fert_woyearplot.stan", data = ELRI_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELRI, file = "endodemog_fert_ELRI.rds")

smELVI2 <- stan(file = "endodemog_fert_woyearplot.stan", data = ELVI_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELVI, file = "endodemog_fert_ELVI.rds")


smFESU2 <- stan(file = "endodemog_fert_woyearplot.stan", data = FESU_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smLOAR, file = "endodemog_fert_LOAR.rds")


## to read in model output without rerunning models
smPOAL <- readRDS(file = "endodemog_fert_POAL_withplot.rds")
smPOSY <- readRDS(file = "endodemog_fert_POSY_withplot.rds")
smLOAR <- readRDS(file = "endodemog_fert_LOAR_withplot.rds")
smELVI <- readRDS(file = "endodemog_fert_ELVI_withplot.rds")
smELRI <- readRDS(file = "endodemog_fert_ELRI_withplot.rds")
smFESU <- readRDS(file = "endodemog_fert_FESU_withplot.rds")
smAGPE <- readRDS(file = "endodemog_fert_AGPE_withplot.rds")


# Temporary section to look at the problems with the LOAR, ELVI, ELRI data

#########################################################################################################
###### Perform posterior predictive checks and assess model convergence-------------------------
#########################################################################################################
params = c("beta[1]", "beta[2]", "tau_year[1,1]", "sigma_e[1]", "sigma_e[2]")
##### POAL - flowering
print(smPOAL)
# summary(smPOAL)

## plot traceplots of chains for select parameters
traceplot(smPOAL, pars = c("beta[2]", "sigma_e[1]"))

## Plotting residuals
flw_t1 <- as.vector(POAL_flw_dat$flw_t1)
yrep <- as.matrix(smPOAL, pars = "yrep")
mu <- as.matrix(smPOAL, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(flw_t1 - p)
yrep_resid <- abs(yrep - p[1,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "POAL flw Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## plot of neff ratios
ratios_smPOAL <- neff_ratio(smPOAL)
mcmc_neff(ratios_smPOAL, size =3)

## Overlay plot of yrep vs flw_t1
ppc_dens_overlay(flw_t1, yrep[1:500, ])


## Density plot of postieror distribution for select parameters
stan_dens(smPOAL, pars = params)

## Option to analyse in shiny
# will probably want to choose certain parameters
shiny <- as.shinystan(smPOAL)
launch_shinystan(shiny)



##### POSY - flowering
print(smPOSY)
# summary(smPOSY)

## plot traceplots of chains for select parameters
traceplot(smPOSY, pars = c("beta[2]", "sigma_e[1]"))

## Plotting residuals
flw_t1 <- as.vector(POSY_data$flw_stat_t1)
yrep <- as.matrix(smPOSY, pars = "yrep")
mu <- as.matrix(smPOSY, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(flw_t1 - p[])
yrep_resid <- abs(yrep[] - p[1,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "POSY flw Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## plot of neff ratios
ratios_smPOSY <- neff_ratio(smPOSY)
mcmc_neff(ratios_smPOSY, size =3)

## Overlay plot of yrep vs flw_t1 
ppc_dens_overlay(flw_t1[], yrep[1:500, ])

## Density plot of postieror distribution for select parameters
stan_dens(smPOSY, pars = c("beta[1]", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Option to analyse in shiny
# will probably want to choose certain parameters
shiny <- as.shinystan(smPOSY)
launch_shinystan(shiny)


##### LOAR - flowering
print(smLOAR)
# summary(smLOAR)

## plot traceplots of chains for select parameters
traceplot(smLOAR, pars = params)

## Plotting residuals
flw_t1 <- as.vector(LOAR_data$flw_stat_t1)
yrep <- as.matrix(smLOAR, pars = "yrep")
mu <- as.matrix(smLOAR, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(flw_t1 - p)
yrep_resid <- abs(yrep - p[1,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "LOAR flw Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## plot of neff ratios
ratios_smLOAR <- neff_ratio(smLOAR)
mcmc_neff(ratios_smLOAR, size =3)

## Overlay plot of yrep vs flw_t1 
ppc_dens_overlay(flw_t1, yrep[1:500, ])

## Density plot of postieror distribution for select parameters
stan_dens(smLOAR, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Option to analyse in shiny
# will probably want to choose certain parameters
shiny <- as.shinystan(smLOAR)
launch_shinystan(shiny)



##### ELVI - flowering
print(smELVI)
# summary(smELVI)

## plot traceplots of chains for select parameters
traceplot(smELVI, pars = c("beta[2]"))

## Plotting residuals
flw_t1 <- as.vector(ELVI_flw_dat$flw_t1)
yrep <- as.matrix(smELVI, pars = "yrep")
mu <- as.matrix(smELVI, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(flw_t1 - p)
yrep_resid <- abs(yrep - p[1,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "ELVI flw Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## plot of neff ratios
ratios_smELVI <- neff_ratio(smELVI)
mcmc_neff(ratios_smELVI, size =3)

## Overlay plot of yrep vs flw_t1 
ppc_dens_overlay(flw_t1, yrep[1:500, ])

## Density plot of postieror distribution for select parameters
stan_dens(smELVI, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Option to analyse in shiny
# will probably want to choose certain parameters
shiny <- as.shinystan(smELVI)
launch_shinystan(shiny)




##### ELRI - flowering
print(smELRI)
# summary(smELRI)

## plot traceplots of chains for select parameters
traceplot(smELRI, pars = params)

## Plotting residuals
flw_t1 <- as.vector(ELRI_flw_dat$flw_t1)
yrep <- as.matrix(smELRI, pars = "yrep")
mu <- as.matrix(smELRI, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(flw_t1 - p)
yrep_resid <- abs(yrep - p[i,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "ELRI flw Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## plot of neff ratios
ratios_smELRI <- neff_ratio(smELRI)
mcmc_neff(ratios_smELRI, size =3)

## Overlay plot of yrep vs flw_t1 
ppc_dens_overlay(flw_t1, yrep[1:500, ])

## Density plot of postieror distribution for select parameters
stan_dens(smELRI, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Option to analyse in shiny
# will probably want to choose certain parameters
shiny <- as.shinystan(smELRI)
launch_shinystan(shiny)


##### FESU - flowering
print(smFESU)
# summary(smFESU)

## plot traceplots of chains for select parameters
traceplot(smFESU, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Plotting residuals
flw_t1 <- as.vector(FESU_flw_dat$flw_t1)
yrep <- as.matrix(smFESU, pars = "yrep")
mu <- as.matrix(smFESU, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(flw_t1 - p)
yrep_resid <- abs(yrep - p[i,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "FESU flw Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## plot of neff ratios
ratios_smFESU <- neff_ratio(smFESU)
mcmc_neff(ratios_smFESU, size =3)

## Overlay plot of yrep vs flw_t1 
ppc_dens_overlay(flw_t1, yrep[1:500, ])

## Density plot of postieror distribution for select parameters
stan_dens(smFESU, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Option to analyse in shiny
# will probably want to choose certain parameters
shiny <- as.shinystan(smFESU)
launch_shinystan(shiny)


##### AGPE - flowering
print(smAGPE)
# summary(smAGPE)

## plot traceplots of chains for select parameters
traceplot(smAGPE, pars = params)

## Plotting residuals
flw_t1 <- as.vector(AGPE_flw_dat$flw_t1)
yrep <- as.matrix(smAGPE, pars = "yrep")
mu <- as.matrix(smAGPE, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(flw_t1 - p)
yrep_resid <- abs(yrep - p[1,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "AGPE flw Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## plot of neff ratios
ratios_smAGPE <- neff_ratio(smAGPE)
mcmc_neff(ratios_smAGPE, size =3)

## Overlay plot of yrep vs flw_t1 
ppc_dens_overlay(flw_t1, yrep[1:500, ])

## Density plot of postieror distribution for select parameters
stan_dens(smAGPE, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Option to analyse in shiny
# will probably want to choose certain parameters
shiny <- as.shinystan(smAGPE)
launch_shinystan(shiny)

