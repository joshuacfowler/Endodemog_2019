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
  mutate(seedperspike = case_when(species == "AGPE" ~ seed, 
                                  species == "ELVI" | species == "ELRI" ~ seed/1,
                                  species == "POAL" | species == "POSY" | species == "FESU" ~ seed/spikelets)) %>% 
  group_by(plot, tag, species, endo_01, year) %>% 
  summarize(seedperspike = mean(seedperspike, na.rm = TRUE))
  
  

dim(LTREB_data1)

# View(LTREB_data1)
#########################################################################################################
# Creating individual species data lists to be passed to the model------------------------------
#########################################################################################################
# Split up the main dataframe by species and recode plot to be used as an index for each species
AGPE_seed_data <- LTREB_data1 %>% 
  filter(species == "AGPE") %>% 
  filter(!is.na(seed)) %>% 
  filter(spikelets > 0) %>% 
  group_by(tag) %>% 
  summarize(seedperspike = seed/spikelets)
  mutate(seedperspike = seed)
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot)))))
ELRI_seed_data <- LTREB_data1 %>% 
  filter(species == "ELRI") %>% 
  filter(!is.na(seed)) %>% 
  group_by(tag) %>% 
  summarize(seedperspike = mean(seed/1)) %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot)))))
ELVI_seed_data <- LTREB_data1 %>% 
  filter(species == "ELVI") %>%  
  filter(!is.na(seed)) %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot)))))
FESU_seed_data <- LTREB_data1 %>% 
  filter(species == "FESU") %>% 
  filter(!is.na(seed)) %>% 
  filter(spikelets > 0) %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot))))) 
LOAR_seed_data <- LTREB_data1 %>% 
  filter(species == "LOAR") %>% 
  filter(!is.na(seed)) %>% 
  filter(spikelets > 0) %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot))))) 
POAL_seed_data <- LTREB_data1 %>% 
  filter(species == "POAL") %>% 
  filter(!is.na(seed)) %>% 
  filter(spikelets > 0) %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot))))) 
POSY_seed_data <- LTREB_data1 %>% 
  filter(species == "POSY") %>% 
  filter(!is.na(seed)) %>% 
  filter(spikelets > 0) %>% 
  mutate(plot_index = as.integer(as.factor(as.integer(as.character(plot)))))



# Create data lists to be used for the Stan model
AGPE_seed_data_list <- list(seed = AGPE_seed_data$seed,
                            spike = AGPE_seed_data$spikelets,
                            flw = AGPE_seed_data$flw,
                            endo_01 = AGPE_seed_data$endo_01,
                            endo_index = AGPE_seed_data$endo_index,
                            year = AGPE_seed_data$year_index,
                            plot = AGPE_seed_data$plot_index,
                            N = nrow(AGPE_seed_data),
                            K = 2,
                            nYear = length(unique(AGPE_seed_data$year_index)),
                            nPlot = length(unique(AGPE_seed_data$plot_index)),
                            nEndo =   length(unique(AGPE_seed_data$endo_01)))
str(AGPE_seed_data_list)

ELRI_seed_data_list <- list(seed = ELRI_seed_data$seed,
                            spike = ELRI_seed_data$spikelets,
                            flw = ELRI_seed_data$flw,
                            endo_01 = ELRI_seed_data$endo_01,
                            endo_index = ELRI_seed_data$endo_index,
                            year = ELRI_seed_data$year_index,
                            plot = ELRI_seed_data$plot_index,
                            N = nrow(ELRI_seed_data),
                            K = 2,
                            nYear = length(unique(ELRI_seed_data$year_index)),
                            nPlot = length(unique(ELRI_seed_data$plot_index)),
                            nEndo =   length(unique(ELRI_seed_data$endo_01)))
str(ELRI_seed_data_list)

ELVI_seed_data_list <- list(seed = ELVI_seed_data$seed,
                            spike = ELVI_seed_data$spikelets,
                            flw = ELVI_seed_data$flw,
                            endo_01 = ELVI_seed_data$endo_01,
                            endo_index = ELVI_seed_data$endo_index,
                            year = ELVI_seed_data$year_index,
                            plot = ELVI_seed_data$plot_index,
                            N = nrow(ELVI_seed_data),
                            K = 2,
                            nYear = length(unique(ELVI_seed_data$year_index)),
                            nPlot = length(unique(ELVI_seed_data$plot_index)),
                            nEndo =   length(unique(ELVI_seed_data$endo_01)))
str(ELVI_seed_data_list)

FESU_seed_data_list <- list(seed = FESU_seed_data$seed,
                            spike = FESU_seed_data$spikelets,
                            flw = FESU_seed_data$flw,
                            endo_01 = FESU_seed_data$endo_01,
                            endo_index = FESU_seed_data$endo_index,
                            year = FESU_seed_data$year_index,
                            plot = FESU_seed_data$plot_index,
                            N = nrow(FESU_seed_data),
                            K = 2,
                            nYear = length(unique(FESU_seed_data$year_index)),
                            nPlot = length(unique(FESU_seed_data$plot_index)),
                            nEndo =   length(unique(FESU_seed_data$endo_01)))
str(FESU_seed_data_list)

LOAR_seed_data_list <- list(seed = LOAR_seed_data$seed,
                            spike = LOAR_seed_data$spikelets,
                            flw = LOAR_seed_data$flw,
                            endo_01 = LOAR_seed_data$endo_01,
                            endo_index = LOAR_seed_data$endo_index,
                            year = LOAR_seed_data$year_index,
                            plot = LOAR_seed_data$plot_index,
                            N = nrow(LOAR_seed_data),
                            K = 2,
                            nYear = length(unique(LOAR_seed_data$year_index)),
                            nPlot = length(unique(LOAR_seed_data$plot_index)),
                            nEndo =   length(unique(LOAR_seed_data$endo_01)))
str(LOAR_seed_data_list)

POAL_seed_data_list <- list(seed = POAL_seed_data$seed,
                            spike = POAL_seed_data$spikelets,
                            flw = POAL_seed_data$flw,
                            endo_01 = POAL_seed_data$endo_01,
                            endo_index = POAL_seed_data$endo_index,
                            year = POAL_seed_data$year_index,
                            plot = POAL_seed_data$plot_index,
                            N = nrow(POAL_seed_data),
                            K = 2,
                            nYear = length(unique(POAL_seed_data$year_index)),
                            nPlot = length(unique(POAL_seed_data$plot_index)),
                            nEndo =   length(unique(POAL_seed_data$endo_01)))
str(POAL_seed_data_list)

POSY_seed_data_list <- list(seed = POSY_seed_data$seed,
                            spike = POSY_seed_data$spikelets,
                            flw = POSY_seed_data$flw,
                            endo_01 = POSY_seed_data$endo_01,
                            endo_index = POSY_seed_data$endo_index,
                            year = POSY_seed_data$year_index,
                            plot = POSY_seed_data$plot_index,
                            N = nrow(POSY_seed_data),
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
    int<lower=0> N;                       // number of observations
    int<lower=0> K;                       // number of predictors
    int<lower=0> nEndo;                       // number of endo treatments
    int<lower=1, upper=2> endo_index[N];          // index for endophyte effect
    real<lower=0> seed[N];               // number of seeds per inflorescence
    real<lower=0> spike[N];                  // number of spikelets per inflorescence
    }
    
    transformed data{
    vector<lower=0>[N] seedperspike; // number of seeds per spikelet 
    for(i in 1:N){
    seedperspike[i] = seed[i]/spike[i];
    }
    }
    
    parameters {
    real<lower=0> mu;
    real<lower=0> sigma_0;
    }
    
    model {

    // Priors
    mu ~ normal(100,100);      // prior for mu

    // Likelihood
      seedperspike ~ normal(mu,sigma_0);
    }

  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_seed_mean.stan")

## Run the model by calling stan()
## Save the outputs as rds files
Pdf <- POAL_seed_data %>% 
  mutate(seed_spike = seed/spikelets) %>% 
  summarize(mean_seed_spike = mean(seed_spike),
            sample_size = n()) 
  summarize(mean = mean(mean_seed_spike),
            sample_size = n())
View(Pdf)


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


# ELRI and ELVI had data collected in a slightly different way




smELRI <- stan(file = "endodemog_seed_mean_elymus.stan", data = ELRI_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELRI, file = "endodemog_seed_mean_ELRI.rds")

smELVI <- stan(file = "endodemog_seed_mean_elymus.stan", data = ELVI_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELVI, file = "endodemog_seed_mean_ELVI.rds")


# AGPE has seed data in a different format
smAGPE<- stan(file = "endodemog_seed_mean.stan", data = AGPE_seed_data_list,
              iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smAGPE, file = "endodemog_seed_mean_AGPE.rds")



