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

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

#############################################################################################
####### Data manipulation to prepare data as lists for Stan models------------------
#############################################################################################
## Load in full data frame
LTREB_endodemog <- 
  read.csv(file = "endo_demog_long.csv")


str(LTREB_endodemog)
dim(LTREB_endodemog)


## Clean up the main data frame for NA's, other small data entry errors, and change standardize the coding for variables.
LTREB_data <- LTREB_endodemog %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(size_t1, logsize_t1 = log(size_t1)) %>%  
  mutate(surv_t1 = as.integer(recode(surv_t1, "0" = 0, "1" =1, "2" = 1, "4" = 1))) %>% 
  mutate(species_0 = (recode_factor(species,                   
                                    "AGPE" = 1, "ELRI" = 2, "ELVI" = 3, 
                                    "FESU" = 4, "LOAR" = 5, "POAL" = 6, 
                                    "POSY" = 7))) %>%        
  mutate(species_index = as.integer(species_0)) %>% 
  mutate(year_t_index = as.integer(recode_factor(year_t, 
                                                 '2007' = 1, '2008' = 2, '2009' = 3, 
                                                 '2010' = 4, '2011' = 5, '2012' = 6, 
                                                 '2013' = 7, '2014' = 8, '2015' = 9, 
                                                 '2016' = 10, '2017' = 11))) %>%             
  mutate(year_t1_index = as.integer(recode_factor(year_t1, 
                                                  '2008' = 2, '2009' = 3, '2010' = 4, 
                                                  '2011' = 5, '2012' = 6, '2013' = 7, 
                                                  '2014' = 8, '2015' = 9, '2016' = 10, 
                                                  '2017' = 11, '2018' = 12))) %>%               
  mutate(origin_01 = as.integer(case_when(origin == "O" ~ 0, 
                                          origin == "R" ~ 1, 
                                          origin != "R" | origin != "O" ~ 1))) %>%   
  mutate(plot_fixed = (case_when(plot != "R" ~ plot, 
                                 plot == "R" ~ origin))) %>%                         
  mutate(plot_index = as.integer(as.factor(plot_fixed))) %>%                         
  mutate(endo_01 = case_when(endo == "0" | endo == "minus" ~ 0,
                             endo == "1"| endo =="plus" ~ 1)) %>% 
  mutate(endo_index = as.integer(as.factor(endo_01+1)))                              

dim(LTREB_data)
unique(LTREB_data$logsize_t)

# NA's in survival come from mostly 2017 recruits.
LTREB_data1 <- LTREB_data %>%
  filter(!is.na(logsize_t1)) %>%
  filter(!is.na(logsize_t)) %>% 
  filter(logsize_t >= 0) %>% 
  filter(!is.na(endo_01))

dim(LTREB_data1)
LTREB_for_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01 + species_0)^2 + origin_01 
                                , data = LTREB_data1)
Xs <- model.matrix(size_t1~ (logsize_t + endo_01 + species_0)^2 + origin_01 
                   , data =LTREB_for_matrix)


# Create data list for Stan model
LTREB_grow_data_list <- list(size_t1 = LTREB_data1$size_t1,
                             Xs = Xs,    
                             endo_index = LTREB_data1$endo_index,
                             species_index = LTREB_data1$species_index,
                             year_t = LTREB_data1$year_t_index, 
                             N = nrow(LTREB_data1), 
                             K = ncol(Xs), 
                             nyear = length(unique(LTREB_data1$year_t_index)), 
                             nEndo =   length(unique(LTREB_data1$endo_01)),
                             nSpp = length(unique(LTREB_data1$species_index)))


str(LTREB_grow_data_list)

LTREB_sample <- sample_n(LTREB_data1, 1000)
sample_for_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01 + species_0) + origin_01 
                                 , data = LTREB_sample)
Xs <- model.matrix(size_t1 ~ (logsize_t + endo_01 + species_0) + origin_01 
                   , data =sample_for_matrix)


# Create sample data list for Stan model
sample_grow_data_list <- list(size_t1 = LTREB_sample$size_t1,
                              Xs = Xs,    
                              endo_index = LTREB_sample$endo_index,
                              species_index = LTREB_sample$species_index,
                              year_t = LTREB_sample$year_t_index, 
                              N = nrow(LTREB_sample), 
                              K = ncol(Xs), 
                              nyear = length(unique(LTREB_sample$year_t_index)), 
                              nEndo =   length(unique(LTREB_sample$endo_01)),
                              nSpp = length(unique(LTREB_sample$species_index)))


str(sample_grow_data_list)

#########################################################################################################
# GLMM for Surv~ size +Endo + Origin  with year random effects------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <- 100
nb <- 50
nc <- 1

# Stan model -------------
## here is the Stan model with a model matrix and species effects ##



sink("endodemog_grow_matrix.stan")
cat("
    data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> K;                       // number of predictors
    
    int<lower=0> nyear;                       // number of years (used as index)
    int<lower=0> year_t[N];                      // year of observation
    int<lower=0> nEndo;                       // number of endo treatments
    int<lower=0> nSpp;                         // number of species
    int<lower=1, upper=2> endo_index[N];       // index for endophyte effect
    int<lower=1, upper=7> species_index[N];  // index for species effect
    int<lower=0> size_t1[N];      // plant survival at time t+1 and target variable (response)
    matrix[N,K] Xs;                  //  predictor matrix - surv_t1~logsize_t+endo+origin+logsize_t*endo
    }
    
    parameters {
    vector[K] beta;                     // predictor parameters
    //matrix[nEndo,nyear] tau_year[nSpp];      // mean random year effect
    //matrix<lower=0>[nSpp,nEndo] sigma_0;        //year variance intercept
    }
 
 
    model {
    vector[N] mu;
    
    // Linear Predictor
    for(n in 1:N){
    mu = Xs*beta;
    //+ tau_year[species_index[n],endo_index[n],year_t[n]];
    }
    // Priors
    beta ~ normal(0,1e6);      // prior for predictor intercepts
    //to_matrix(tau_year) ~ normal(0,100);
    
    //sigma_0 ~ gamma(0,1/10);
   // for(s in 1:nSpp){
    //for(e in 1:nSpp){
   // tau_year[s,e] ~ normal(0,sigma_0[s,e]);  // prior for year random effects
   // }}
    // Likelihood
      size_t1 ~ binomial_logit(N,mu);
    }
    
   //generated quantities{
    //int yrep[N];
    //vector[N] mu;
    
    // for posterior predictive check
    //for(n in 1:N){
     // mu[n] = Xs[n]*beta
    //   + tau_year[endo_index[n], year_t[n]];
      
    //  yrep[n] = bernoulli_logit_rng(mu[n]);
    //}
    
   // }
  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_grow_matrix.stan")



## Run the model by calling stan()
## and save the output to .rds files so that they can be called laters

sm <- stan(file = "endodemog_grow_matrix.stan", data = sample_grow_data_list,
           iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)

print(sm)
summary(sm)
print(sm, pars = "tau_year")
