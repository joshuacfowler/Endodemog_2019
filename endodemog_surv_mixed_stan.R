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
ni <- 2000
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
    
    model {
        // Linear Predictor
    vector[N] mu;
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
saveRDS(smAGPE, file = "endodemog_surv_AGPE.rds")

smELRI <- stan(file = "endodemog_surv.stan", data = ELRI_surv_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELRI, file = "endodemog_surv_ELRI.rds")

smELVI <- stan(file = "endodemog_surv.stan", data = ELVI_surv_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELVI, file = "endodemog_surv_ELVI.rds")

smFESU <- stan(file = "endodemog_surv.stan", data = FESU_surv_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smFESU, file = "endodemog_surv_FESU.rds")

smLOAR <- stan(file = "endodemog_surv.stan", data = LOAR_surv_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smLOAR, file = "endodemog_surv_LOAR.rds")

smPOAL <- stan(file = "endodemog_surv.stan", data = POAL_surv_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOAL, file = "endodemog_surv_POAL.rds")

smPOSY <- stan(file = "endodemog_surv.stan", data = POSY_surv_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOSY, file = "endodemog_surv_POSY.rds")

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
AGPE_surv_yrep <- prediction(AGPE_surv_data_list, smAGPE, 500)
ELRI_surv_yrep <- prediction(ELRI_surv_data_list, smELRI, 500)
ELVI_surv_yrep <- prediction(ELVI_surv_data_list, smELVI, 500)
FESU_surv_yrep <- prediction(FESU_surv_data_list, smFESU, 500)
LOAR_surv_yrep <- prediction(LOAR_surv_data_list, smLOAR, 500)
POAL_surv_yrep <- prediction(POAL_surv_data_list, smPOAL, 500)
POSY_surv_yrep <- prediction(POSY_surv_data_list, smPOSY, 500)


# overlay 100 replicates over the actual dataset
msurv_yrep_AGPE <- t(AGPE_surv_yrep$yrep)
ppc_dens_overlay( y = AGPE_surv_data_list$surv_t1, yrep = msurv_yrep_AGPE[1:100,])+ xlab("prob. of y") + ggtitle("AGPE")

msurv_yrep_ELRI <- t(ELRI_surv_yrep$yrep)
ppc_dens_overlay( y = ELRI_surv_data_list$surv_t1, yrep = msurv_yrep_ELRI[1:100,])+ xlab("prob. of y") + ggtitle("ELRI")

msurv_yrep_ELVI <- t(ELVI_surv_yrep$yrep)
ppc_dens_overlay( y = ELVI_surv_data_list$surv_t1, yrep = msurv_yrep_ELVI[1:100,]) + xlab("prob. of y") + ggtitle("ELVI")

msurv_yrep_FESU <- t(FESU_surv_yrep$yrep)
ppc_dens_overlay( y = FESU_surv_data_list$surv_t1, yrep = msurv_yrep_FESU[1:100,]) + xlab("prob. of y") + ggtitle("FESU")

msurv_yrep_LOAR <- t(LOAR_surv_yrep$yrep)
ppc_dens_overlay( y = LOAR_surv_data_list$surv_t1, yrep = msurv_yrep_LOAR[1:100,]) + xlab("prob. of y") + ggtitle("LOAR")

msurv_yrep_POAL <- t(POAL_surv_yrep$yrep)
ppc_dens_overlay( y = POAL_surv_data_list$surv_t1, yrep = msurv_yrep_POAL[1:100,]) + xlab("prob. of y") + ggtitle("POAL")

msurv_yrep_POSY <- t(POSY_surv_yrep$yrep)
ppc_dens_overlay( y = POSY_surv_data_list$surv_t1, yrep = msurv_yrep_POSY[1:100,]) + xlab("prob. of y") + ggtitle("POSY")























resid_calc <- function(y,yrep){
y_resid <- y$surv_t1 - yrep$prob
yrep_resid <- yrep$yrep - yrep$prob  

y_resid_abs <-  abs(y_resid)
yrep_resid_abs <- abs(yrep_resid)

fit_col <- colSums(y_resid_abs)
fit_yrep_col <- colSums(yrep_resid_abs)

fit_row <- rowSums(y_resid_abs)
fit_yrep_row <- rowSums(yrep_resid_abs)

out <- list(y_resid, yrep_resid, y_resid_abs, yrep_resid_abs, fit_col, fit_yrep_col, fit_row, fit_yrep_row)
names(out) <- c("y_resid", "yrep_resid", "y_resid_abs", "yrep_resid_abs", "fit_col", "fit_yrep_col", "fit_row", "fit_yrep_row")

return(out)
}

resid <- resid_calc(POAL_surv_data_list, POAL_surv_yrep)
col_resid <- cbind(resid$yrep_resid_abs, resid$y_resid_abs)


col <- cbind(resid$fit_yrep_col, resid$fit_col)
plot(col)

row <- cbind(resid$fit_yrep_row, resid$fit_row)

plot(row)
abline(a = 0, b = 1, col = "blue", lwd = 2)

plot(resid$fit,resid$fit_yrep)
abline(a = 0, b = 1, col = "blue", lwd = 2)

yrep <- t(POAL_surv_yrep$yrep)
ppc_dens_overlay( y = POAL_surv_data_list$surv_t1, yrep = yrep[1:500,])



yresid <- resid[[1]]
yrep_resid <- resid[[2]]


fit <- colSums(yresid)
fit_yrep <- colSums(yrep_resid)

plot(x = fit_yrep, y = fit)
abline(a = 0, b = 1, col = "blue", lwd = 2)





ivalues <- length(beta1_post)
nvalues <- length(POAL_surv_data$logsize_t)
xdummy <- seq(min(POAL_surv_data$logsize_t), max(POAL_surv_data$logsize_t), length.out = nvalues)

surv_fits <- function(df_samples){
  require(tidyverse)
  
  nvalues <- length(df_raw$logsize_t)
  xdummy <- seq(min(df_raw$logsize_t), max(df_raw$logsize_t), length.out = nvalues)
  
  #endophyte negative fits  
  ydummy_eminus <- as.vector(invlogit(mean(df_samples["beta[1]"[i]]) + mean(df_samples["beta[2]"])*xdummy + mean(df_samples["beta[3]"])*0 + mean(df_samples["beta[4]"])*mean(df_raw$origin_01) + mean(df_samples["beta[5]"])*xdummy*0))

}
post_survPOAL <- extract(smPOAL)
# generate credible intervals
POALpost_list <- as.list(post_survPOAL)

    for(i in 1:7500){
    POAL_ydummy_eminus <- invlogit((post_survPOAL$beta[i,1]) + (post_survPOAL$beta[i,2])*2 + (post_survPOAL$beta[i,3])*0 + (post_survPOAL$beta[i,4])*mean(POAL_surv_data$origin_01) + (post_survPOAL$beta[i,5])*xdummy*0)
  }

POAL_ydummy_eminus <- invlogit((post_survPOAL$beta[1:7500,1]) + (post_survPOAL$beta[1:7500,2])*2 + (post_survPOAL$beta[1:7500,3])*0 + (post_survPOAL$beta[1:7500,4])*mean(POAL_surv_data$origin_01) + (post_survPOAL$beta[1:7500,5])*2*0)


POAL_ydummy_eminus <- 


"The \"R\" project for statistical computing"

POALpost$mu[1]
mu <-  beta1_post + beta2_post*xdummy + beta3_post*0 + beta4_post*origin_01 + beta5_post*xdummy*0


warnings()
mu_plus <- 
tauyear_post <- POALpost$tau_year
tauplot_post <- POALpost$tau_plot
s <- rnorm(mean = mean(tauplot_post), sd = sd(tauplot_post), n = 1000000)
mean(s)
mean(tauplot_post)
str(POAL_surv_data_list)

N <- as.integer(POAL_surv_data_list$N)
logsize_t <- POAL_surv_data_list$logsize_t
origin_01 <- POAL_surv_data_list$origin_01
endo_01 <- POAL_surv_data_list$endo_01
surv_t1 <- POAL_surv_data_list$surv_t1
# Function for replicating y based on x
lin_pred <-    replicate(1000, sample(beta1_post, size = N)+ sample(beta2_post,size = N)*logsize_t + sample(beta3_post,size = N)*endo_01 + sample(beta4_post,size = N)*origin_01 + sample(beta5_post,size = N)*logsize_t*endo_01 +sample(rnorm(mean = mean(tauyear_post), sd = sd(tauyear_post), n = 1000000), size = N)+ sample(rnorm(mean = mean(tauplot_post), sd = sd(tauplot_post), n = 1000000), size = N))

prob <- as.list(as.data.frame(invlogit(lin_pred)))


yrep <- list()
yrep <- lapply(prob,(rbinom(3241, 1, prob = prob[1:1000])))
  yrep <- 



# generate residuals

y_resid <-  abs(surv_t1 - prob)

yrep_resid <- abs(yrep - prob[])


fit <- colSums(y_resid)
fit_yrep <- colSums(yrep_resid)

plot(x = fit_yrep, y = fit)
abline(a = 0, b = 1, col = "blue", lwd = 2)



plot(y = yrep_resid[,1], x = logsize_t)
ppc_dens_overlay(yrep = yrep[,1:10],y = surv_t1)


length(surv_t1)
# Run the function on x_test
set.seed(56)
y_pred_r <- gen_quant_r(x_test)
# Accuracy
mean(y_pred_r == y_test)
# generate yrep from samples of parameters
estimate_mu_rep <- function(x, post, length){
  beta1 <- sample(post$beta[,1], size = length)
  beta2 <- sample(post$beta[,2], size = length)
  beta3 <- sample(post$beta[,3], size = length)
  beta4 <- sample(post$beta[,4], size = length)
  beta5 <- sample(post$beta[,5], size = length)
  lin_pred <- c(1:length)
  for(i in 1:length){
    lin_pred[i] <-    beta1[i] + beta2[i]*x$logsize_t[i] + beta3[i]*x$endo_01[i] + beta4[i]*x$origin_01[i] + beta5[i]*x$logsize_t[i]*x$endo_01[i]
  }
  mu <- invlogit(lin_pred)
  
  return(mu)
}
estimate_y_rep <- function(x, post, length){
  beta1 <- sample(post$beta[,1], size = length)
  beta2 <- sample(post$beta[,2], size = length)
  beta3 <- sample(post$beta[,3], size = length)
  beta4 <- sample(post$beta[,4], size = length)
  beta5 <- sample(post$beta[,5], size = length)
  lin_pred <- c(1:length)
for(i in 1:length){
  lin_pred[i] <-    beta1[i] + beta2[i]*x$logsize_t[i] + beta3[i]*x$endo_01[i] + beta4[i]*x$origin_01[i] + beta5[i]*x$logsize_t[i]*x$endo_01[i]
}
  prob <- invlogit(lin_pred)
  estimate <- rbinom(length, 1, prob)

  return(estimate)
}

yrep <- estimate_y_rep(x = POAL_surv_data_list, post = POALpost, length = length(POAL_surv_data_list$surv_t1))
yrep <- replicate(100, estimate_y_rep(x = POAL_surv_data_list, post = POALpost, length = length(POAL_surv_data_list$surv_t1)))

murep <- estimate_mu_rep(x = POAL_surv_data_list, post = POALpost, length = length(POAL_surv_data_list$surv_t1))
murep <- replicate(100, estimate_mu_rep(x = POAL_surv_data_list, post = POALpost, length = length(POAL_surv_data_list$surv_t1)))

## Plotting residuals
surv_t1 <- as.vector(POAL_surv_data_list$surv_t1)
hist(surv_t1)
hist(yrep[,1])
y_resid <-  abs(surv_t1 - murep)
yrep_resid <- abs(yrep - murep[,1])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = yrep_resid[,1], y = fit, main = "POAL surv Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

plot(x = POAL_surv_data_list$logsize_t, y = surv_t1, main = "POAL surv Residual plot")
plot(x = POAL_surv_data_list$logsize_t, y = yrep[,1], main = "POAL surv Residual plot")
plot(x = POAL_surv_data_list$logsize_t, y = yrep[,2], main = "POAL surv Residual plot")
mean(yrep)
ppc_dens_overlay(yrep[1:3241], surv_t1)

## plot of neff ratios
ratios_smPOAL <- neff_ratio(smPOAL)
mcmc_neff(ratios_smPOAL, size =3)

## Overlay plot of yrep vs surv_t1 
ppc_dens_overlay(surv_t1, yrep[,1])

## Density plot of postieror distribution for select parameters
stan_dens(smPOAL, pars = params)

## Option to analyse in shiny
# will probably want to choose certain parameters
shiny <- as.shinystan(smPOAL)
launch_shinystan(shiny)


##### POSY - survival
print(smPOSY)
# summary(smPOSY)

## plot traceplots of chains for select parameters
traceplot(smPOSY, pars = params)

## Plotting residuals
surv_t1 <- as.vector(POSY_surv_dat$surv_t1)
yrep <- as.matrix(smPOSY, pars = "yrep")
mu <- as.matrix(smPOSY, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(surv_t1 - p)
yrep_resid <- abs(yrep - p[1,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "POSY surv Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## plot of neff ratios
ratios_smPOSY <- neff_ratio(smPOSY)
mcmc_neff(ratios_smPOSY, size =3)

## Overlay plot of yrep vs surv_t1 
ppc_dens_overlay(POAL_surv_data_list$surv_t1, resid$yrep[,1:10])

## Density plot of postieror distribution for select parameters
stan_dens(smPOSY, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Option to analyse in shiny
# will probably want to choose certain parameters
shiny <- as.shinystan(smPOSY)
launch_shinystan(shiny)


##### LOAR - survival
print(smLOAR)
# summary(smLOAR)

## plot traceplots of chains for select parameters
traceplot(smLOAR, pars = params)

## Plotting residuals
surv_t1 <- as.vector(LOAR_surv_dat$surv_t1)
yrep <- as.matrix(smLOAR, pars = "yrep")
mu <- as.matrix(smLOAR, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(surv_t1 - p)
yrep_resid <- abs(yrep - p[1,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "LOAR surv Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## plot of neff ratios
ratios_smLOAR <- neff_ratio(smLOAR)
mcmc_neff(ratios_smLOAR, size =3)

## Overlay plot of yrep vs surv_t1 
ppc_dens_overlay(surv_t1, yrep[1:500, ])

## Density plot of postieror distribution for select parameters
stan_dens(smLOAR, pars = params)

## Option to analyse in shiny
# will probably want to choose certain parameters
shiny <- as.shinystan(smLOAR)
launch_shinystan(shiny)



##### ELVI - survival
print(smELVI)
# summary(smELVI)

## plot traceplots of chains for select parameters
traceplot(smELVI, pars = params)

## Plotting residuals
surv_t1 <- as.vector(ELVI_surv_dat$surv_t1)
yrep <- as.matrix(smELVI, pars = "yrep")
mu <- as.matrix(smELVI, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(surv_t1 - p)
yrep_resid <- abs(yrep - p[1,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "ELVI surv Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## plot of neff ratios
ratios_smELVI <- neff_ratio(smELVI)
mcmc_neff(ratios_smELVI, size =3)

## Overlay plot of yrep vs surv_t1 
ppc_dens_overlay(surv_t1, yrep[1:500, ])

## Density plot of postieror distribution for select parameters
stan_dens(smELVI, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Option to analyse in shiny
# will probably want to choose certain parameters
shiny <- as.shinystan(smELVI)
launch_shinystan(shiny)




##### ELRI - survival
print(smELRI)
# summary(smELRI)

## plot traceplots of chains for select parameters
traceplot(smELRI, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Plotting residuals
surv_t1 <- as.vector(ELRI_surv_dat$surv_t1)
yrep <- as.matrix(smELRI, pars = "yrep")
mu <- as.matrix(smELRI, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(surv_t1 - p)
yrep_resid <- abs(yrep - p[1,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "ELRI surv Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## plot of neff ratios
ratios_smELRI <- neff_ratio(smELRI)
mcmc_neff(ratios_smELRI, size =3)

## Overlay plot of yrep vs surv_t1 
ppc_dens_overlay(surv_t1, yrep[1:500, ])

## Density plot of postieror distribution for select parameters
stan_dens(smELRI, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Option to analyse in shiny
# will probably want to choose certain parameters
shiny <- as.shinystan(smELRI)
launch_shinystan(shiny)




##### FESU - survival
print(smFESU)
# summary(smFESU)

## plot traceplots of chains for select parameters
traceplot(smFESU, pars = c("beta[1]", "beta[2]", "tau_plot[1]", "sigma_e[1]", "sigma_e[2]"))

## Plotting residuals
surv_t1 <- as.vector(FESU_data$surv_t1)
yrep <- as.matrix(smFESU, pars = "yrep")
mu <- as.matrix(smFESU, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(surv_t1 - p)
yrep_resid <- abs(yrep - p[1,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "FESU surv with plot 2 Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## plot of neff ratios
ratios_smFESU <- neff_ratio(smFESU)
mcmc_neff(ratios_smFESU, size =3)

## Overlay plot of yrep vs surv_t1 
ppc_dens_overlay(surv_t1, yrep[1:500, ])

## Density plot of postieror distribution for select parameters
stan_dens(smFESU, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Option to analyse in shiny
# will probably want to choose certain parameters
shiny <- as.shinystan(smFESU)
launch_shinystan(shiny)


##### AGPE - survival
print(smAGPE)
# summary(smAGPE)

## plot traceplots of chains for select parameters
traceplot(smAGPE, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Plotting residuals
surv_t1 <- as.vector(AGPE_surv_dat$surv_t1)
yrep <- as.matrix(smAGPE, pars = "yrep")
mu <- as.matrix(smAGPE, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(surv_t1 - p)
yrep_resid <- abs(yrep - p[1,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "AGPE surv Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## plot of neff ratios
ratios_smAGPE <- neff_ratio(smAGPE)
mcmc_neff(ratios_smAGPE, size =3)

## Overlay plot of yrep vs surv_t1 
ppc_dens_overlay(surv_t1, yrep[1:500, ])

## Density plot of postieror distribution for select parameters
stan_dens(smAGPE, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Option to analyse in shiny
# will probably want to choose certain parameters
shiny <- as.shinystan(smAGPE)
launch_shinystan(shiny)










