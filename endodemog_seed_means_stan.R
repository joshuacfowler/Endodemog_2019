## Title: Grass endophyte population model with a bayesian framework
## Purpose: Creates seed production kernel written in STAN, 
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

LTREB_repro

## Clean up the main data frame for NA's, other small data entry errors, and change standardize the coding for variables.
LTREB_data <- LTREB_repro %>% 
  rename(endo = Endo) %>% 
  select(-contains("Birth Year")) %>% 
  mutate(endo_01 = as.integer(case_when(endo == "0" | endo == "minus" ~ 0,
                                        endo == "1"| endo =="plus" ~ 1))) %>% 
  mutate(endo_index = as.integer(as.factor(endo_01+1))) %>% 
  mutate(species_index = as.integer(recode_factor(species,                   
                                                  "AGPE" = 1, "ELRI" = 2, "ELVI" = 3, 
                                                  "FESU" = 4, "LOAR" = 5, "POAL" = 6, 
                                                  "POSY" = 7))) %>% 
  mutate(year_index = as.integer(recode(year, 
                                           '2008' = 2, '2009' = 3, '2010' = 4, 
                                           '2011' = 5, '2012' = 6, '2013' = 7, 
                                           '2014' = 8, '2015' = 9, '2016' = 10, 
                                           '2017' = 11, '2018' = 12)))  %>% 
  filter(!is.na(plot))

dim(LTREB_data)

LTREB_data1 <- LTREB_data %>%
  filter(!is.na(flw)) %>% 
  filter(flw>0) %>% 
  mutate(seedperspike = case_when(species == "ELVI" | species == "ELRI" ~ NA_real_,
                                  species == "POAL" | species == "POSY" | species == "FESU" | species == "LOAR" ~ seed/spikelets, # for these species, they are recorded as seeds/infl (collected for a few years) and spikelets/infl (collected for all year), so we are calculating avg seed/spike to be able to multiply byspiklets/infl in years without seed counts 
                                  species == "AGPE" ~ seed)) %>% 
  mutate(seedperinf = case_when(species == "ELVI" | species == "ELRI" ~ seed,
                                species == "POAL" | species == "POSY" | species == "FESU" | species == "LOAR" ~ NA_real_,
                                species == "AGPE"~ NA_real_)) %>% 
  group_by(plot, tag, species, endo_01,endo_index, year, year_index) %>% 
  summarize(seedperspike = mean(seedperspike, na.rm = TRUE),
            seedperinf = mean(seedperinf, na.rm = TRUE),
            spikeperinf = mean(spikelets),
            flw = mean(flw),
            samplesize = n())
  
  

dim(LTREB_data1)

View(LTREB_data1)
#########################################################################################################
# Creating individual species data lists to be passed to the model------------------------------
#########################################################################################################
# Split up the main dataframe by species and recode plot to be used as an index for each species
AGPE_seed_data <- LTREB_data1 %>% 
  filter(species == "AGPE") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot)))))
ELRI_seed_data <- LTREB_data1 %>% 
  filter(species == "ELRI") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot)))))
ELVI_seed_data <- LTREB_data1 %>% 
  filter(species == "ELVI") %>%  
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot)))))
FESU_seed_data <- LTREB_data1 %>% 
  filter(species == "FESU") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot))))) 
LOAR_seed_data <- LTREB_data1 %>% 
  filter(species == "LOAR") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot))))) 
POAL_seed_data <- LTREB_data1 %>% 
  filter(species == "POAL") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot))))) 
POSY_seed_data <- LTREB_data1 %>% 
  filter(species == "POSY") %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot)))))



# Create data lists to be used for the Stan model
AGPE_seed_data_list <- list(seed = na.omit(AGPE_seed_data$seedperspike),
                            spike = na.omit(AGPE_seed_data$spikeperinf),
                            endo_01 = AGPE_seed_data$endo_01,
                            endo_index = AGPE_seed_data$endo_index,
                            year = AGPE_seed_data$year_index,
                            plot = AGPE_seed_data$plot_index,
                            Nseed = length(na.omit(AGPE_seed_data$seedperspike)),
                            Nspike = length(na.omit(AGPE_seed_data$spikeperinf)),
                            K = 2,
                            nYear = length(unique(AGPE_seed_data$year_index)),
                            nPlot = length(unique(AGPE_seed_data$plot_index)),
                            nEndo =   length(unique(AGPE_seed_data$endo_01)))
str(AGPE_seed_data_list)

ELRI_seed_data_list <- list(seed = na.omit(ELRI_seed_data$seedperinf),
                            endo_01 = na.omit(ELRI_seed_data$endo_01),
                            endo_index = ELRI_seed_data$endo_index,
                            year = ELRI_seed_data$year_index,
                            plot = ELRI_seed_data$plot_index,
                            N = nrow(ELRI_seed_data),
                            K = 2,
                            nYear = length(unique(ELRI_seed_data$year_index)),
                            nPlot = length(unique(ELRI_seed_data$plot_index)),
                            nEndo =   length(unique(ELRI_seed_data$endo_01)))
str(ELRI_seed_data_list)

ELVI_seed_data_list <- list(seed = na.omit(ELVI_seed_data$seedperinf),
                            endo_01 = na.omit(ELVI_seed_data$endo_01),
                            endo_index = ELVI_seed_data$endo_index,
                            year = ELVI_seed_data$year_index,
                            plot = ELVI_seed_data$plot_index,
                            N = nrow(ELVI_seed_data),
                            K = 2,
                            nYear = length(unique(ELVI_seed_data$year_index)),
                            nPlot = length(unique(ELVI_seed_data$plot_index)),
                            nEndo =   length(unique(ELVI_seed_data$endo_01)))
str(ELVI_seed_data_list)

FESU_seed_data_list <- list(seed = na.omit(FESU_seed_data$seedperspike),
                            spike = na.omit(FESU_seed_data$spikeperinf),
                            endo_01 = na.omit(FESU_seed_data$endo_01),
                            endo_index = FESU_seed_data$endo_index,
                            year = FESU_seed_data$year_index,
                            plot = FESU_seed_data$plot_index,
                            Nseed = length(na.omit(FESU_seed_data$seedperspike)),
                            Nspike = length(na.omit(FESU_seed_data$spikeperinf)),
                            K = 2,
                            nYear = length(unique(FESU_seed_data$year_index)),
                            nPlot = length(unique(FESU_seed_data$plot_index)),
                            nEndo =   length(unique(FESU_seed_data$endo_01)))
str(FESU_seed_data_list)

LOAR_seed_data_list <- list(seed = na.omit(LOAR_seed_data$seedperspike),
                            spike = na.omit(LOAR_seed_data$spikeperinf),
                            endo_01 = na.omit(LOAR_seed_data$endo_01),
                            endo_index = LOAR_seed_data$endo_index,
                            year = LOAR_seed_data$year_index,
                            plot = LOAR_seed_data$plot_index,
                            Nseed = length(na.omit(LOAR_seed_data$seedperspike)),
                            Nspike = length(na.omit(LOAR_seed_data$spikeperinf)),
                            K = 2,
                            nYear = length(unique(LOAR_seed_data$year_index)),
                            nPlot = length(unique(LOAR_seed_data$plot_index)),
                            nEndo =   length(unique(LOAR_seed_data$endo_01)))
str(LOAR_seed_data_list)

POAL_seed_data_list <- list(seed = na.omit(POAL_seed_data$seedperspike),
                            spike = na.omit(POAL_seed_data$spikeperinf),
                            endo_01 = na.omit(POAL_seed_data$endo_01),
                            endo_index = POAL_seed_data$endo_index,
                            year = POAL_seed_data$year_index,
                            plot = POAL_seed_data$plot_index,
                            Nseed = length(na.omit(POAL_seed_data$seedperspike)),
                            Nspike = length(na.omit(POAL_seed_data$spikeperinf)),
                            K = 2,
                            nYear = length(unique(POAL_seed_data$year_index)),
                            nPlot = length(unique(POAL_seed_data$plot_index)),
                            nEndo =   length(unique(POAL_seed_data$endo_01)))
str(POAL_seed_data_list)

POSY_seed_data_list <- list(seed = na.omit(POSY_seed_data$seedperspike),
                            spike = na.omit(POSY_seed_data$spikeperinf),
                            endo_01 = na.omit(POSY_seed_data$endo_01),
                            endo_index = POSY_seed_data$endo_index,
                            year = POSY_seed_data$year_index,
                            plot = POSY_seed_data$plot_index,
                            Nseed = length(na.omit(POSY_seed_data$seedperspike)),
                            Nspike = length(na.omit(POSY_seed_data$spikeperinf)),
                            K = 2,
                            nYear = length(unique(POSY_seed_data$year_index)),
                            nPlot = length(unique(POSY_seed_data$plot_index)),
                            nEndo =   length(unique(POSY_seed_data$endo_01)))
str(POSY_seed_data_list)


#########################################################################################################
# Stan model for mean of seed production ------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <-2000
nb <- 1000
nc <- 3

# Stan model -------------
## here is the Stan model with a model matrix and species effects ##

sink("endodemog_seed_mean.stan")
cat("
    data { 
    int<lower=0> Nseed;                       // number of observations of seed/spikelet
    int<lower=0> Nspike;                       // number of observations of spikelets/inf
    int<lower=0> K;                       // number of predictors
    real<lower=0> seed[Nseed];               // number of seeds per spikelet or seeds per inflorescence
    real<lower=0> spike[Nspike];
    }
    
    
    parameters {
    real<lower=0> mu_seed;
    real<lower=0> sigma_seed;
    real<lower=0> mu_spike;
    real<lower=0> sigma_spike;
    }
    
    model {

    // Priors
    // Likelihood
      seed ~ normal(mu_seed,sigma_seed);
      spike ~ normal(mu_spike, sigma_spike);
    }

  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_seed_mean.stan")

## Run the model by calling stan()
## Save the outputs as rds files
smFESU <- stan(file = "endodemog_seed_mean.stan", data = FESU_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smFESU, file = "endodemog_seed_mean_FESU.rds")

smLOAR <- stan(file = "endodemog_seed_mean.stan", data = LOAR_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smLOAR, file = "endodemog_seed_mean_LOAR.rds")

smPOAL <- stan(file = "endodemog_seed_mean.stan", data = POAL_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOAL, file = "endodemog_seed_mean_POAL_withplot.rds")

smPOSY <- stan(file = "endodemog_seed_mean.stan", data = POSY_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOSY, file = "endodemog_seed_mean_POSY_withplot.rds")

# AGPE had data recorded as seed/spikelet already, but we are using this model to calculate seed/spikelet on average, so this is using the seed/spikelet info from AGPE
smAGPE<- stan(file = "endodemog_seed_mean.stan", data = AGPE_seed_data_list,
              iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smAGPE, file = "endodemog_seed_mean_AGPE.rds")


# ELRI and ELVI had data collected in a slightly different way. They recorded seeds/inflorescence


smELRI <- stan(file = "endodemog_seed_mean.stan", data = ELRI_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELRI, file = "endodemog_seed_mean_ELRI.rds")

smELVI <- stan(file = "endodemog_seed_mean.stan", data = ELVI_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELVI, file = "endodemog_seed_mean_ELVI.rds")


est_seed <- function(df,x){
v <- sample(df$mu, size = nrow(x), replace = TRUE)*x$spikeperinf*x$flw
return(v)
}

# calculate POAL seed estimate
POAL <- as.data.frame(smPOAL)
POAL_est_data <- LTREB_data1 %>% 
  filter(species == "POAL")


POAL_est_data$seed_est <- est_seed(POAL, POAL_est_data)
View(POAL_est_data)
