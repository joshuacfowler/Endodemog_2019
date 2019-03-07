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

# LTREB_for_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01 + species)^3 + origin_01 
#                                 , data = LTREB_data1)
# Xs <- model.matrix(size_t1~ (logsize_t + endo_01 + species)^3 + origin_01 
#                    , data =LTREB_for_matrix)
# 
# 
# # Create data list for Stan model
# LTREB_grow_data_list <- list(size_t1 = LTREB_data1$size_t1,
#                              Xs = Xs,    
#                              spp_endo_index = LTREB_data1$spp_endo_index,
#                              year_t = LTREB_data1$year_t_index, 
#                              N = nrow(LTREB_data1), 
#                              K = ncol(Xs), 
#                              lowerlimit = min(LTREB_data1$size_t1),
#                              nyear = length(unique(LTREB_data1$year_t_index)), 
#                              nEndo =   length(unique(LTREB_data1$endo_01)),
#                              nSpp = length(unique(LTREB_data1$species_index)))
# 
# 
# str(LTREB_grow_data_list)
# 
# LTREB_sample <- sample_n(LTREB_data1, 1000)
# sample_for_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01 + species)^3 + origin_01 
#                                  , data = LTREB_sample)
# Xs <- model.matrix(size_t1 ~ (logsize_t + endo_01 + species)^3 + origin_01 
#                    , data =sample_for_matrix)
# 
# 
# # Create sample data list for Stan model
# sample_grow_data_list <- list(size_t1 = LTREB_sample$size_t1,
#                               Xs = Xs,    
#                               spp_endo_index = LTREB_sample$spp_endo_index,
#                               year_t = LTREB_sample$year_t_index, 
#                               N = nrow(LTREB_sample), 
#                               K = ncol(Xs), 
#                               lowerlimit = min(LTREB_sample$size_t1),
#                               nyear = length(unique(LTREB_sample$year_t_index)), 
#                               nEndo =   length(unique(LTREB_sample$endo_01)),
#                               nSpp = length(unique(LTREB_sample$species_index)))
# 
# 
# str(sample_grow_data_list)


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
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
LOAR_data <- LTREB_data1 %>% 
  filter(species == "LOAR") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
POAL_data <- LTREB_data1 %>% 
  filter(species == "POAL") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))
POSY_data <- LTREB_data1 %>% 
  filter(species == "POSY") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot_fixed)))))


# Create model matrices for each species
AGPE_for_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01 + origin_01)^3
                               , data = AGPE_data)
AGPE_Xs <- model.matrix(size_t1 ~ (logsize_t + endo_01 + origin_01)^3
                        , data =AGPE_for_matrix)

ELRI_for_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = ELRI_data)
ELRI_Xs <- model.matrix(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =ELRI_for_matrix)

ELVI_for_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = ELVI_data)
ELVI_Xs <- model.matrix(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =ELVI_for_matrix)

FESU_for_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = FESU_data)
FESU_Xs <- model.matrix(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =FESU_for_matrix)

LOAR_for_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = LOAR_data)
LOAR_Xs <- model.matrix(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =LOAR_for_matrix)

POAL_for_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = POAL_data)
POAL_Xs <- model.matrix(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =POAL_for_matrix)

POSY_for_matrix <- model.frame(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                               , data = POSY_data)
POSY_Xs <- model.matrix(size_t1 ~ (logsize_t + endo_01)^2 + origin_01
                        , data =POSY_for_matrix)
=======

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

>>>>>>> 3b195afb74218a5ee8d92add3f8b76cff775ad42

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
    int<lower=1, upper=2> endo_index[N];          // index for endophyte effect
    //int<lower=0> nPlot;                         // number of plots
    //int<lower=0> plot[N];                   // plot of observation
    int<lower=lowerlimit> size_t1[N];      // plant survival at time t+1 and target variable (response)
    matrix[N,K] Xs;                  //  predictor matrix - surv_t1~logsize_t+endo+origin+logsize_t*endo
    }
    
    parameters {
    vector[K] beta;                     // predictor parameters

    real tau_year[nEndo,nYear];      // random year effect
    real<lower=0> sigma_e[nEndo];        //year variance by endophyte effect
    //vector[nPlot] tau_plot;        // random plot effect
    real<lower=0> phi;
    }
    
    model {
    vector[N] mu;
    
    // Linear Predictor
    for(n in 1:N){
    mu = Xs*beta + tau_year[endo_index[n],year_t[n]]; // + tau_plot[plot[n]];
    }
    // Priors
    beta ~ normal(0,100);      // prior for predictor intercepts
    //to_vector(tau_plot) ~ normal(0,100);   // prior for plot random effects
    for(e in 1:nEndo){
    for(y in 1:nYear){
    tau_year[e,y]~ normal(0,sigma_e[e]); // prior for year random effects
    }}

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
# saveRDS(smAGPE, file = "endodemog_grow_AGPE.rds")

smELRI <- stan(file = "endodemog_grow.stan", data = ELRI_grow_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELRI, file = "endodemog_grow_ELRI.rds")

smELVI <- stan(file = "endodemog_grow.stan", data = ELVI_grow_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smSLVI, file = "endodemog_grow_ELVI.rds")
smFESU <- stan(file = "endodemog_grow.stan", data = FESU_grow_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smFESU, file = "endodemog_grow_FESU.rds")

smLOAR <- stan(file = "endodemog_grow.stan", data = LOAR_grow_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smLOAR, file = "endodemog_grow_LOAR.rds")

smPOAL <- stan(file = "endodemog_grow.stan", data = POAL_grow_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOAL, file = "endodemog_grow_POAL.rds")

smPOSY <- stan(file = "endodemog_grow.stan", data = POSY_grow_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOSY, file = "endodemog_grow_POSY.rds")





print(sm)
summary(sm)
print(sm, pars = "tau_year")
