## Title: Grass endophyte population model with a bayesian framework
## Purpose: Creates flowering (flowering tiller status) kernel written in STAN with mixed effects, 
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

# FLower data lists are generated in the endodemog_data_processing.R file
# within the section titled "Preparing datalists for Flower Kernel"
source("endodemog_data_processing.R")

#########################################################################################################
# GLMM for Surv~ size +Endo + Origin  with year random effects-------------------------
########################################################################################3
## here is the Stan model ##
## run this to optimize computer system settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <- 10000
nb <- 5000
nc <- 3


sink("endodemog_flw.stan")
cat("
    data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> K;                       // number of predictors
    
    int<lower=0> nYear;                       // number of years (used as index)
    int<lower=0> year_t[N];                      // year of observation
    int<lower=0> nPlot;                       // number of plots (used as index)
    int<lower=0> plot[N];                      // plot of observation

    int<lower=0, upper=1> flw_t[N];      // flowering status at time t+1 and target variable (response)
    vector<lower=-1>[N] logsize_t;                  // log of plant size at time t (predictor)
    int<lower=0> nEndo;                           // number of endo treatments
    int<lower=0, upper=1> endo_01[N];            // endophyte status (predictor)
    int<lower=1, upper=2> endo_index[N];        // index for endophyte effect
    int<lower=0, upper=1> origin_01[N];            // origin status(predictor)
    }
    
    parameters {
    vector[K] beta;                     // predictor parameters
    vector[nYear] tau_year[nEndo];      // random year effect
    real<lower=0> sigma_e[nEndo];        //year variance by endophyte effect
    vector[nPlot] tau_plot;        // random plot effect
    real<lower=0> sigma_p;         // plot variance
    }

    transformed parameters{
    real mu[N];                             //Linear Predictor
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
    
    // Likelihood
    flw_t ~ bernoulli_logit(mu);
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

stanmodel <- stanc("endodemog_flw.stan")


## Run the model by calling stan()
## Save the outputs as rds files

smAGPE<- stan(file = "endodemog_flw.stan", data = AGPE_flw_data_list,
              iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smAGPE, file = "endodemog_flw_AGPE.rds")

smELRI <- stan(file = "endodemog_flw.stan", data = ELRI_flw_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELRI, file = "endodemog_flw_ELRI.rds")

smELVI <- stan(file = "endodemog_flw.stan", data = ELVI_flw_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELVI, file = "endodemog_flw_ELVI.rds")

smFESU <- stan(file = "endodemog_flw.stan", data = FESU_flw_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smFESU, file = "endodemog_flw_FESU.rds")

smLOAR <- stan(file = "endodemog_flw.stan", data = LOAR_flw_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smLOAR, file = "endodemog_flw_LOAR.rds")

smPOAL <- stan(file = "endodemog_flw.stan", data = POAL_flw_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOAL, file = "endodemog_flw_POAL.rds")

smPOSY <- stan(file = "endodemog_flw.stan", data = POSY_flw_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOSY, file = "endodemog_flw_POSY.rds")






## to read in model output without rerunning models

smAGPE <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_AGPE.rds")
smELRI <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_ELRI.rds")
smELVI <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_ELVI.rds")
smFESU <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_FESU.rds")
smLOAR <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_LOAR.rds")
smPOAL <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_POAL.rds")
smPOSY <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_flw_POSY.rds")



#########################################################################################################
###### Perform posterior predictive checks and assess model convergence-------------------------
#########################################################################################################
params = c("beta[1]", "beta[2]", "tau_year[1,1]", "sigma_e[1]", "sigma_e[2]")

##### POAL - flowering
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
post_flwAGPE <- extract(smAGPE)
post_flwELRI <- extract(smELRI)
post_flwELVI <- extract(smELVI)
post_flwFESU <- extract(smFESU)
post_flwLOAR <- extract(smLOAR)
post_flwPOAL <- extract(smPOAL)
post_flwPOSY <- extract(smPOSY)




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
    
    prob <- invlogit(lin_comb)
    yrep[,n] <- rbinom(x$N, 1, prob[,n])
    print(paste("rep", n, "of", reps))
    
  }
  out <- list(yrep, prob, lin_comb) 
  names(out) <- c("yrep", "prob", "lin_comb")
  
  return(out)
}

# apply the function for each species
AGPE_flw_yrep <- prediction(AGPE_flw_data_list, smAGPE, 500)
ELRI_flw_yrep <- prediction(ELRI_flw_data_list, smELRI, 500)
ELVI_flw_yrep <- prediction(ELVI_flw_data_list, smELVI, 500)
FESU_flw_yrep <- prediction(FESU_flw_data_list, smFESU, 500)
LOAR_flw_yrep <- prediction(LOAR_flw_data_list, smLOAR, 500)
POAL_flw_yrep <- prediction(POAL_flw_data_list, smPOAL, 500)
POSY_flw_yrep <- prediction(POSY_flw_data_list, smPOSY, 500)


# overlay 100 replicates over the actual dataset
mflw_yrep_AGPE <- t(AGPE_flw_yrep$yrep)
ppc_dens_overlay( y = AGPE_flw_data_list$flw_t, yrep = mflw_yrep_AGPE[1:100,])+ xlab("prob. of y") + ggtitle("AGPE")

mflw_yrep_ELRI <- t(ELRI_flw_yrep$yrep)
ppc_dens_overlay( y = ELRI_flw_data_list$flw_t, yrep = mflw_yrep_ELRI[1:100,])+ xlab("prob. of y") + ggtitle("ELRI")

mflw_yrep_ELVI <- t(ELVI_flw_yrep$yrep)
ppc_dens_overlay( y = ELVI_flw_data_list$flw_t, yrep = mflw_yrep_ELVI[1:100,]) + xlab("prob. of y") + ggtitle("ELVI")

mflw_yrep_FESU <- t(FESU_flw_yrep$yrep)
ppc_dens_overlay( y = FESU_flw_data_list$flw_t, yrep = mflw_yrep_FESU[1:100,]) + xlab("prob. of y") + ggtitle("FESU")

mflw_yrep_LOAR <- t(LOAR_flw_yrep$yrep)
ppc_dens_overlay( y = LOAR_flw_data_list$flw_t, yrep = mflw_yrep_LOAR[1:100,]) + xlab("prob. of y") + ggtitle("LOAR")

mflw_yrep_POAL <- t(POAL_flw_yrep$yrep)
ppc_dens_overlay( y = POAL_flw_data_list$flw_t, yrep = mflw_yrep_POAL[1:100,]) + xlab("prob. of y") + ggtitle("POAL")

mflw_yrep_POSY <- t(POSY_flw_yrep$yrep)
ppc_dens_overlay( y = POSY_flw_data_list$flw_t, yrep = mflw_yrep_POSY[1:100,]) + xlab("prob. of y") + ggtitle("POSY")





































## plot traceplots of chains for select parameters
traceplot(smPOAL, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma_0[1]", "sigma_0[2]"))

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
traceplot(smPOSY, pars = c("beta[2]"))

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
traceplot(smELRI, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

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





