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

# View(LTREB_data1)
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
                            spike = na.omit(ELRI_seed_data$spikeperinf),
                            endo_01 = na.omit(ELRI_seed_data$endo_01),
                            endo_index = ELRI_seed_data$endo_index,
                            year = ELRI_seed_data$year_index,
                            plot = ELRI_seed_data$plot_index,
                            Nseed = length(na.omit(ELRI_seed_data$seedperinf)),
                            Nspike = length(na.omit(ELRI_seed_data$spikeperinf)),
                            N = nrow(ELRI_seed_data),
                            K = 2,
                            nYear = length(unique(ELRI_seed_data$year_index)),
                            nPlot = length(unique(ELRI_seed_data$plot_index)),
                            nEndo =   length(unique(ELRI_seed_data$endo_01)))
str(ELRI_seed_data_list)

ELVI_seed_data_list <- list(seed = na.omit(ELVI_seed_data$seedperinf),
                            spike = na.omit(ELVI_seed_data$spikeperinf),
                            endo_01 = na.omit(ELVI_seed_data$endo_01),
                            endo_index = ELVI_seed_data$endo_index,
                            year = ELVI_seed_data$year_index,
                            plot = ELVI_seed_data$plot_index,
                            Nseed = length(na.omit(ELVI_seed_data$seedperinf)),
                            Nspike = length(na.omit(ELVI_seed_data$spikeperinf)),
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
    real<lower=0> seed[Nseed];               // number of seeds per spikelet
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
sink("endodemog_seed_mean_elymus.stan")
cat("
    data { 
    int<lower=0> Nseed;                       // number of observations of seed/infl
    int<lower=0> K;                       // number of predictors
    real<lower=0> seed[Nseed];               // number of  seeds per inflorescence
    }
    
    
    parameters {
    real<lower=0> mu_seed;
    real<lower=0> sigma_seed;
    }
    
    model {

    // Priors
    // Likelihood
      seed ~ normal(mu_seed,sigma_seed);
    }

  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_seed_mean_elymus.stan")


smELRI <- stan(file = "endodemog_seed_mean_elymus.stan", data = ELRI_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELRI, file = "endodemog_seed_mean_ELRI.rds")

smELVI <- stan(file = "endodemog_seed_mean_elymus.stan", data = ELVI_seed_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELVI, file = "endodemog_seed_mean_ELVI.rds")

########################################################################################################################################
##### Using the seed means to get seed production estimates 
########################################################################################################################################
# For the non-Elymus species, this function uses the posteriors generated by the model and the measured reproductive data to estimate seed production (currently I'm generating this in three different ways because I'm not sure what the best way to incorporate the measurements is)
est_seed <- function(post,data){
set.seed(5-1-2019)
est <- data
est$seedperspike_est <-sample(post$mu_seed, size = nrow(data))
est$spikeperinf_est <- sample(post$mu_spike, size = nrow(data))
est$seed_a<- as.integer(data$seedperspike*sample(post$mu_spike, size = nrow(data))*data$flw)
est$seed_b <- as.integer(sample(post$mu_seed, size = nrow(data))*data$spikeperinf*data$flw)
est$seed <- as.integer(sample(post$mu_seed, size = nrow(data))*sample(post$mu_spike, size = nrow(data))*data$flw) # this is the seed estimate I'm using for now.
est$seed_norma <- as.integer(data$seedperspike*rnorm(mean = mean(post$mu_spike), sd = mean(post$sigma_spike), n = nrow(data))*data$flw)
est$seed_normb <- as.integer(rnorm(mean = mean(post$mu_seed),sd = mean(post$sigma_seed), n = nrow(data))*data$spikeperinf*data$flw)
est$seed_normc <- as.integer(rnorm(mean = mean(post$mu_seed), sd = mean(post$sigma_seed), n = nrow(data))*rnorm(mean = mean(post$mu_spike), sd = mean(post$sigma_spike), n= nrow(data))*data$flw)

return(est)
}


# calculate POAL seed estimate
POAL_post <- as.data.frame(smPOAL)
POAL_seed_est <- est_seed(post = POAL_post, data = POAL_seed_data)
  select(-seed_a, -seed_b, -seedperspike_est, -spikeperinf_est, -seed_norma, -seed_normb, -seed_normc)

# calculate POSY seed estimate
POSY_post <- as.data.frame(smPOSY)
POSY_seed_est <- est_seed(post = POSY_post, data = POSY_seed_data)%>% 
  select(-seed_a, -seed_b, -seedperspike_est, -spikeperinf_est, -seed_norma, -seed_normb, -seed_normc)


# calculate FESU seed estimate
FESU_post <- as.data.frame(smFESU)
FESU_seed_est <- est_seed(post = FESU_post, data = FESU_seed_data)%>% 
  select(-seed_a, -seed_b, -seedperspike_est, -spikeperinf_est, -seed_norma, -seed_normb, -seed_normc)


# calculate LOAR seed estimate
LOAR_post <- as.data.frame(smLOAR)
LOAR_seed_est <- est_seed(post = LOAR_post, data = LOAR_seed_data)%>% 
  select(-seed_a, -seed_b, -seedperspike_est, -spikeperinf_est, -seed_norma, -seed_normb, -seed_normc)


# calculate AGPE seed estimate
AGPE_post <- as.data.frame(smAGPE)
AGPE_seed_est <- est_seed(post = AGPE_post, data = AGPE_seed_data)%>% 
  select(-seed_a, -seed_b, -seedperspike_est, -spikeperinf_est, -seed_norma, -seed_normb, -seed_normc)



# For the Elymus species, this function uses the posteriors generated by the model and the measured reproductive data to estimate seed production (currently I'm generating this in three different ways because I'm not sure what the best way to incorporate the measurements is)
est_seed_Elymus <- function(post,data){
  set.seed(5-1-2019)
  est <- data
  est$seedperinf_est <-sample(post$mu_seed, size = nrow(data))
  est$seedperinf_normest <- rnorm(mean = mean(post$mu_seed), sd = mean(post$sigma_seed), n = nrow(data))
  est$seed_a <- as.integer(data$seedperinf*data$flw)
  est$seed <- as.integer(sample(post$mu_seed, size = nrow(data))*data$flw) #This is the seed estimate I'm using for now
  est$seed_normb <- as.integer(rnorm(mean = mean(post$mu_seed), sd = mean(post$sigma_seed), n= nrow(data))*data$flw)
  return(est)
}

# calculate ELRI seed estimate
ELRI_post <- as.data.frame(smELRI)
ELRI_seed_est <- est_seed_Elymus(post = ELRI_post, data = ELRI_seed_data)%>% 
  select(-seed_a, -seedperinf_est, -seedperinf_normest,  -seed_normb)


# calculate ELVI seed estimate
ELVI_post <- as.data.frame(smELVI)
ELVI_seed_est <- est_seed_Elymus(post = ELVI_post, data = ELVI_seed_data)%>% 
  select(-seed_a, -seedperinf_est, -seedperinf_normest,  -seed_normb)





# Now I am going to merge the seed estimates for the different species into one dataframe. This will be merged with endo_demog_long for the seed to seedling transition rate.
seed_est_data <- ELVI_seed_est %>% 
  merge(ELRI_seed_est, by = c("plot","tag", "species", "endo_01", "endo_index", "year", "year_index", "seedperspike", "seedperinf", "spikeperinf", "flw", "samplesize", "plot_index", "seed" ),all = TRUE) %>% 
  merge(AGPE_seed_est, by = c("plot","tag", "species", "endo_01", "endo_index", "year", "year_index", "seedperspike", "seedperinf", "spikeperinf", "flw", "samplesize", "plot_index", "seed" ),all = TRUE) %>% 
  merge(FESU_seed_est, by = c("plot","tag", "species", "endo_01", "endo_index", "year", "year_index", "seedperspike", "seedperinf", "spikeperinf", "flw", "samplesize", "plot_index", "seed" ),all = TRUE) %>% 
  merge(LOAR_seed_est, by = c("plot","tag", "species", "endo_01", "endo_index", "year", "year_index", "seedperspike", "seedperinf", "spikeperinf", "flw", "samplesize", "plot_index", "seed" ),all = TRUE) %>% 
  merge(POAL_seed_est, by = c("plot","tag", "species", "endo_01", "endo_index", "year", "year_index", "seedperspike", "seedperinf", "spikeperinf", "flw", "samplesize", "plot_index", "seed" ),all = TRUE) %>% 
  merge(POSY_seed_est, by = c("plot","tag", "species", "endo_01", "endo_index", "year", "year_index", "seedperspike", "seedperinf", "spikeperinf", "flw", "samplesize", "plot_index", "seed" ),all = TRUE)
  

# View(seed_est_data)
# getting a dataframe with seed_t and seed_t1
seed_est_t1 <-seed_est_data %>%
  filter(year!= min(year)) %>% 
  rename(year_t1 = year, seed_t1 = seed, flw_t1 = flw) %>%  
  mutate(year_t = year_t1 - 1) %>% 
  select(-samplesize, -seedperinf, -seedperspike, -spikeperinf, -year_index)

seed_est_t <- seed_est_data %>%
  filter(year != max(year)) %>% 
  rename(year_t = year, seed_t = seed, flw_t = flw) %>% 
  mutate(year_t1 = year_t + 1) %>% 
  select(-samplesize, -seedperinf, -seedperspike, -spikeperinf, -year_index)


seed_estmerge <- seed_est_t1 %>% 
  full_join(seed_est_t, by = c("plot", "tag", "endo_01", "endo_index", "year_t", "year_t1", "species", "plot_index"),
            all.x = all, all.y = all)




# I'm going to merge this with the main data frame
LTREB_endodemog <-  read.csv(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_long.csv")

LTREB_data_seed_est <- LTREB_endodemog %>% 
  rename(tag = id) %>% 
  select(-seed_t, -seed_t1, -flw_t1) %>% 
  merge(seed_estmerge, by = c("plot", "tag", "species", "year_t", "year_t1"), all = TRUE)
# View(LTREB_data_seed_est)


hist(AGPE_seed_est$seed)
hist(ELVI_seed_est$seed)
hist(ELRI_seed_est$seed)
hist(LOAR_seed_est$seed)
hist(POAL_seed_est$seed)
hist(POSY_seed_est$seed)
hist(FESU_seed_est$seed)
