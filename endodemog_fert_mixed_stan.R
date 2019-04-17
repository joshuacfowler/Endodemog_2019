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
## Load in full data frame
# LTREB_endodemog <- 
# read.csv(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_long.csv")


#"C:/Users/MillerLab/Desktop/Endodemog-master/endo_demog_long.csv"
#"/Users/joshuacfowler/Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_long.csv")
str(LTREB_endodemog)
dim(LTREB_endodemog)


# ####### Run Endodemog data processing script to get this data frame, which we will clean up!
# This is pulled from the raw data files for flw tiller status, and then merged with LTREB_endodemog
LTREB_flw


## Clean up the main data frame for NA's, other small data entry errors, and change standardize the coding for variables.
LTREB_data <- LTREB_flw %>% 
  rename(endo = Endo) %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(size_t1, logsize_t1 = log(size_t1)) %>%  
  mutate(surv_t1 = as.integer(recode(surv_t1, "0" = 0, "1" =1, "2" = 1, "4" = 1))) %>% 
  mutate(endo_01 = as.integer(case_when(endo == "0" | endo == "minus" ~ 0,
                                        endo == "1"| endo =="plus" ~ 1))) %>% 
  mutate(endo_index = as.integer(as.factor(endo_01+1)))  %>% 
  mutate(species_index = as.integer(recode_factor(species,                   
                                                  "AGPE" = 1, "ELRI" = 2, "ELVI" = 3, 
                                                  "FESU" = 4, "LOAR" = 5, "POAL" = 6, 
                                                  "POSY" = 7))) %>% 
  mutate(year_t_index = as.integer(recode(year_t, 
                                          '2007' = 1, '2008' = 2, '2009' = 3, 
                                          '2010' = 4, '2011' = 5, '2012' = 6, 
                                          '2013' = 7, '2014' = 8, '2015' = 9, 
                                          '2016' = 10, '2017' = 11))) %>%             
  mutate(year_t1_index = as.integer(recode(year_t1, 
                                           '2008' = 2, '2009' = 3, '2010' = 4, 
                                           '2011' = 5, '2012' = 6, '2013' = 7, 
                                           '2014' = 8, '2015' = 9, '2016' = 10, 
                                           '2017' = 11, '2018' = 12))) %>%               
  mutate(origin_01 = as.integer(case_when(origin == "O" ~ 0, 
                                          origin == "R" ~ 1, 
                                          origin != "R" | origin != "O" ~ 1))) %>%   
  mutate(plot_fixed = (case_when(plot != "R" ~ plot, 
                                 plot == "R" ~ origin)))                       

dim(LTREB_data)

# NA's in survival come from mostly 2017 recruits.
LTREB_data1 <- LTREB_data %>%
  mutate(flw_stat_t1 = as.integer(flw_stat_t1)) %>%
  mutate(flw_stat_t = as.integer(flw_stat_t)) %>% 
  mutate(flw_t = as.integer(flw_t)) %>% 
  mutate(flw_t1 = as.integer(flw_t1)) %>% 
  filter(!is.na(flw_t1)) %>% 
  filter(flw_t1>0) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(logsize_t >= 0) %>% 
  filter(!is.na(endo_01))


dim(LTREB_data1)

#########################################################################################################
# Creating individual species data lists to be passed to the model------------------------------
#########################################################################################################
# Split up the main dataframe by species and recode plot to be used as an index for each species
AGPE_data <- LTREB_data1 %>% 
  filter(species == "AGPE") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
ELRI_data <- LTREB_data1 %>% 
  filter(species == "ELRI") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
ELVI_data <- LTREB_data1 %>% 
  filter(species == "ELVI") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
FESU_data <- LTREB_data1 %>% 
  filter(species == "FESU") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed))))) %>% 
  mutate(year_t_index = as.integer(recode(year_t, 
                                          '2007' = 0, '2008' = 1, '2009' = 2, 
                                          '2010' = 3, '2011' = 4, '2012' = 5, 
                                          '2013' = 6, '2014' = 7, '2015' = 8, 
                                          '2016' = 9, '2017' = 10)))
LOAR_data <- LTREB_data1 %>% 
  filter(species == "LOAR") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed))))) %>% 
  mutate(year_t_index = as.integer(recode(year_t, 
                                          '2007' = 0, '2008' = 0, '2009' = 0, 
                                          '2010' = 1, '2011' = 0, '2012' = 2, 
                                          '2013' = 3, '2014' = 4, '2015' = 5, 
                                          '2016' = 6, '2017' = 7)))
POAL_data <- LTREB_data1 %>% 
  filter(species == "POAL") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed))))) %>%
  mutate(year_t_index = as.integer(recode(year_t, 
                                          '2007' = 1, '2008' = 2, '2009' = 3, 
                                          '2010' = 4, '2011' = 5, '2012' = 6, 
                                          '2013' = 7, '2014' = 0, '2015' = 8, 
                                          '2016' = 9, '2017' = 10)))
POSY_data <- LTREB_data1 %>% 
  filter(species == "POSY") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))


# Create model matrices for each species
AGPE_for_matrix <- model.frame(flw_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = AGPE_data)
AGPE_Xs <- model.matrix(flw_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =AGPE_for_matrix)

ELRI_for_matrix <- model.frame(flw_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = ELRI_data)
ELRI_Xs <- model.matrix(flw_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =ELRI_for_matrix)

ELVI_for_matrix <- model.frame(flw_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = ELVI_data)
ELVI_Xs <- model.matrix(flw_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =ELVI_for_matrix)

FESU_for_matrix <- model.frame(flw_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = FESU_data)
FESU_Xs <- model.matrix(flw_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =FESU_for_matrix)

LOAR_for_matrix <- model.frame(flw_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = LOAR_data)
LOAR_Xs <- model.matrix(flw_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =LOAR_for_matrix)

POAL_for_matrix <- model.frame(flw_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = POAL_data)
POAL_Xs <- model.matrix(flw_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =POAL_for_matrix)

POSY_for_matrix <- model.frame(flw_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = POSY_data)
POSY_Xs <- model.matrix(flw_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =POSY_for_matrix)



# Create data lists to be used for the Stan model
AGPE_fert_data_list <- list(flw_t = AGPE_data$flw_t,
                           Xs = AGPE_Xs,
                           logsize_t = AGPE_data$logsize_t,
                           origin_01 = AGPE_data$origin_01,
                           endo_01 = AGPE_data$endo_01,
                           endo_index = AGPE_data$endo_index,
                           year_t = AGPE_data$year_t_index,
                           plot = AGPE_data$plot_index,
                           N = nrow(AGPE_data),
                           K = ncol(AGPE_Xs),
                           lowerlimit = as.integer(min(AGPE_data$flw_t)),
                           nYear = length(unique(AGPE_data$year_t_index)),
                           nPlot = length(unique(AGPE_data$plot_index)),
                           nEndo =   length(unique(AGPE_data$endo_01)))
str(AGPE_fert_data_list)

ELRI_fert_data_list <- list(flw_t = ELRI_data$flw_t,
                           Xs = ELRI_Xs,
                           logsize_t = ELRI_data$logsize_t,
                           origin_01 = ELRI_data$origin_01,
                           endo_01 = ELRI_data$endo_01,
                           endo_index = ELRI_data$endo_index,
                           year_t = ELRI_data$year_t_index,
                           plot = ELRI_data$plot_index,
                           N = nrow(ELRI_data),
                           K = ncol(ELRI_Xs),
                           lowerlimit = as.integer(min(ELRI_data$flw_t)),
                           nYear = length(unique(ELRI_data$year_t_index)),
                           nPlot = length(unique(ELRI_data$plot_index)),
                           nEndo =   length(unique(ELRI_data$endo_01)))
str(ELRI_fert_data_list)

ELVI_fert_data_list <- list(flw_t = ELVI_data$flw_t,
                           Xs = ELVI_Xs,
                           logsize_t = ELVI_data$logsize_t,
                           origin_01 = ELVI_data$origin_01,
                           endo_01 = ELVI_data$endo_01,
                           endo_index = ELVI_data$endo_index,
                           year_t = ELVI_data$year_t_index,
                           plot = ELVI_data$plot_index,
                           N = nrow(ELVI_data),
                           K = ncol(ELVI_Xs),
                           lowerlimit = as.integer(min(ELVI_data$flw_t)),
                           nYear = length(unique(ELVI_data$year_t_index)),
                           nPlot = length(unique(ELVI_data$plot_index)),
                           nEndo =   length(unique(ELVI_data$endo_01)))
str(ELVI_fert_data_list)

FESU_fert_data_list <- list(flw_t = FESU_data$flw_t,
                           Xs = FESU_Xs,
                           logsize_t = FESU_data$logsize_t,
                           origin_01 = FESU_data$origin_01,
                           endo_01 = FESU_data$endo_01,
                           endo_index = FESU_data$endo_index,
                           year_t = FESU_data$year_t_index,
                           plot = FESU_data$plot_index,
                           N = nrow(FESU_data),
                           K = ncol(FESU_Xs),
                           lowerlimit = as.integer(min(FESU_data$flw_t)),
                           nYear = length(unique(FESU_data$year_t_index)),
                           nPlot = length(unique(FESU_data$plot_index)),
                           nEndo =   length(unique(FESU_data$endo_01)))
str(FESU_fert_data_list)

LOAR_fert_data_list <- list(flw_t = LOAR_data$flw_t,
                           Xs = LOAR_Xs,
                           logsize_t = LOAR_data$logsize_t,
                           origin_01 = LOAR_data$origin_01,
                           endo_01 = LOAR_data$endo_01,
                           endo_index = LOAR_data$endo_index,
                           year_t = LOAR_data$year_t_index,
                           plot = LOAR_data$plot_index,
                           N = nrow(LOAR_data),
                           K = ncol(LOAR_Xs),
                           lowerlimit = as.integer(min(LOAR_data$flw_t)),
                           nYear = length(unique(LOAR_data$year_t_index)),
                           nPlot = length(unique(LOAR_data$plot_index)),
                           nEndo =   length(unique(LOAR_data$endo_01)))
str(LOAR_fert_data_list)
POAL_fert_data_list <- list(flw_t = POAL_data$flw_t,
                           Xs = POAL_Xs,
                           logsize_t = POAL_data$logsize_t,
                           origin_01 = POAL_data$origin_01,
                           endo_01 = POAL_data$endo_01,
                           endo_index = POAL_data$endo_index,
                           year_t = POAL_data$year_t_index,
                           plot = POAL_data$plot_index,
                           N = nrow(POAL_data),
                           K = ncol(POAL_Xs),
                           lowerlimit = as.integer(min(POAL_data$flw_t)),
                           nYear = length(unique(POAL_data$year_t_index)),
                           nPlot = length(unique(POAL_data$plot_index)),
                           nEndo =   length(unique(POAL_data$endo_01)))
str(POAL_fert_data_list)

POSY_fert_data_list <- list(flw_t = POSY_data$flw_t,
                           Xs = POSY_Xs,
                           logsize_t = POSY_data$logsize_t,
                           origin_01 = POSY_data$origin_01,
                           endo_01 = POSY_data$endo_01,
                           endo_index = POSY_data$endo_index,
                           year_t = POSY_data$year_t_index,
                           plot = POSY_data$plot_index,
                           N = nrow(POSY_data),
                           K = ncol(POSY_Xs),
                           lowerlimit = as.integer(min(POSY_data$flw_t)),
                           nYear = length(unique(POSY_data$year_t_index)),
                           nPlot = length(unique(POSY_data$plot_index)),
                           nEndo =   length(unique(POSY_data$endo_01)))
str(POSY_fert_data_list)


#########################################################################################################
# GLMM for flw_t1~ size +Endo + Origin  with year random effects------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <-500
nb <- 250
nc <- 1

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
    //vector[nPlot] tau_plot;        // random plot effect
    //real<lower=0> sigma_p;          // plot variance
    real<lower=0> phi;
    }
    
    model {
    vector[N] mu;
    
    // Linear Predictor
    for(n in 1:N){
      mu[n] = beta[1] + beta[2]*logsize_t[n] + beta[3]*endo_01[n] +beta[4]*origin_01[n]
      + beta[5]*logsize_t[n]*endo_01[n] 
      + tau_year[endo_index[n],year_t[n]];
      //+ tau_plot[plot[n]];
    }
    // Priors
    beta ~ normal(0,100);      // prior for predictor intercepts
    //tau_plot ~ normal(0,sigma_p);   // prior for plot random effects
    to_vector(tau_year[1]) ~ normal(0,sigma_e[1]);   // prior for E- year random effects
    to_vector(tau_year[2]) ~ normal(0,sigma_e[2]);   // prior for E+ year random effects
    

    // Likelihood
    for(n in 1:N){
      flw_t[n] ~ neg_binomial_2_log(mu[n],phi);
      target += -log1m(neg_binomial_2_log_lpmf(lowerlimit | mu[n], phi)); // manually adjusting computation of likelihood because T[,] truncation syntax doesn't compile for neg binomial
    }
    }
    
   generated quantities{
      vector[N] mu;
    
    // for posterior predictive check
    for(n in 1:N){
      mu[n] = beta[1] + beta[2]*logsize_t[n] + beta[3]*endo_01[n] +beta[4]*origin_01[n]
      + beta[5]*logsize_t[n]*endo_01[n] 
      + tau_year[endo_index[n],year_t[n]];
     // + tau_plot[plot[n]];
    }
   }
  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_fert.stan")

## Run the model by calling stan()
## Save the outputs as rds files

smAGPE<- stan(file = "endodemog_fert.stan", data = AGPE_fert_data_list,
              iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smAGPE, file = "endodemog_fert_AGPE.rds")

smELRI <- stan(file = "endodemog_fert.stan", data = ELRI_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELRI, file = "endodemog_fert_ELRI.rds")

smELVI <- stan(file = "endodemog_fert.stan", data = ELVI_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELVI, file = "endodemog_fert_ELVI.rds")

smFESU <- stan(file = "endodemog_fert.stan", data = FESU_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smFESU, file = "endodemog_fert_FESU.rds")

smLOAR <- stan(file = "endodemog_fert.stan", data = LOAR_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smLOAR, file = "endodemog_fert_LOAR.rds")

smPOAL <- stan(file = "endodemog_fert.stan", data = POAL_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOAL, file = "endodemog_fert_POAL_withplot.rds")

smPOSY <- stan(file = "endodemog_fert.stan", data = POSY_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOSY, file = "endodemog_fert_POSY_withplot.rds")



# Here is a simpler model for those species with little data


sink("endodemog_fert_woyearplot.stan")
cat("
    data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> K;                       // number of predictors
    int<lower=0> lowerlimit;                         //lower limit for truncated negative binomial
    int<lower=0> nEndo;                       // number of endo treatments
    int<lower=1, upper=2> endo_index[N];          // index for endophyte effect
    int<lower=lowerlimit> flw_t1[N];      // plant size at time t+1 and target variable (response)
    vector<lower=0>[N] logsize_t;             // plant size at time t (predictor)
    int<lower=0,upper=1> endo_01[N];            // plant endophyte status (predictor)
    int<lower=0,upper=1> origin_01[N];          // plant origin status (predictor)
    }
    
    parameters {
    vector[3] beta;                     // predictor parameters

    real<lower=0> phi;
    }
    
    model {
    vector[N] mu;
    
    // Linear Predictor
    for(n in 1:N){
      mu[n] = beta[1] + beta[2]*logsize_t[n];
    }
    // Priors
    beta ~ normal(0,100);      // prior for predictor intercepts

    // Likelihood
    for(n in 1:N){
      flw_t1[n] ~ neg_binomial_2_log(mu[n],phi);
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

smELVI <- stan(file = "endodemog_fert_woyearplot.stan", data = ELVI_fert_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELVI, file = "endodemog_fert_ELVI.rds")


smLOAR <- stan(file = "endodemog_fert_woyearplot.stan", data = LOAR_fert_data_list,
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

