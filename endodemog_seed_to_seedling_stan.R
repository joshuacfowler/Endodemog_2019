## Title: Grass endophyte population model with a bayesian framework
## Purpose: Creates seed to seedling kernel written in STAN with mixed effects, 
## and does visualisation of posterior predictive checks
## Authors: Joshua and Tom
#############################################################

library(tidyverse)
library(rstan)
library(StanHeaders)
library(shinystan)
library(bayesplot)
library(devtools)
library(countreg)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }
#############################################################################################
####### Data manipulation to prepare data as lists for Stan models------------------
#############################################################################################

# Using seed_means and spike data to generate seed estimates for all plots and years
# read in the model outputs for seed means and for spike/inflorescence
# ELRI and ELVI had data recorded as seed/infl
sm_seedAGPE <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_AGPE.rds")
sm_seedELRI <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_ELRI.rds")
sm_seedELVI <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_ELVI.rds")
sm_seedFESU <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_FESU.rds")
sm_seedLOAR <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_LOAR.rds")
sm_seedPOAL <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_POAL.rds")
sm_seedPOSY <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_seed_mean_POSY.rds")

sm_spikeAGPE <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_AGPE.rds")
sm_spikeELRI <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_ELRI.rds")
sm_spikeELVI <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_ELVI.rds")
sm_spikeFESU <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_FESU.rds")
sm_spikeLOAR <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_LOAR.rds")
sm_spikePOAL <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_POAL.rds")
sm_spikePOSY <- readRDS(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Model_Runs/endodemog_spike_POSY.rds")


# Surv data lists are generated in the endodemog_data_processing.R file
# within the section titled "Preparing datalists for Growth Kernel"
source("endodemog_data_processing.R")


# I'm gonna merge the seed estimates with the surv data which includes more than just the plants that are flowering but also includes recruits.
est_seed <- function(sm_seed,sm_spike,data){
  post_seed <- rstan::extract(sm_seed)
  post_spike <- rstan::extract(sm_spike)
  data$fixedeffect <- sample(post_spike$beta[,1], size = length(data$FLW_STAT_T)) + 
    sample(post_spike$beta[,2], size = length(data$FLW_STAT_T))*data$logsize_t +
    sample(post_spike$beta[,3], size = length(data$FLW_STAT_T))*data$endo_01 +
    sample(post_spike$beta[,4], size = length(data$FLW_STAT_T))*data$origin_01 
  data$plot_effect <- ifelse(data$plot_index == 1, sample(post_spike$tau_plot[,1], size = length(data$FLW_STAT_T)), ifelse(data$plot_index == 2, sample(post_spike$tau_plot[,2], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 3, sample(post_spike$tau_plot[,3], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 4, sample(post_spike$tau_plot[,4], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 5, sample(post_spike$tau_plot[,5], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 6, sample(post_spike$tau_plot[,6], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 7, sample(post_spike$tau_plot[,7], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 8, sample(post_spike$tau_plot[,8], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 9, sample(post_spike$tau_plot[,9], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 10, sample(post_spike$tau_plot[,10], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 11, sample(post_spike$tau_plot[,11], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 12, sample(post_spike$tau_plot[,12], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 13, sample(post_spike$tau_plot[,13], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 14, sample(post_spike$tau_plot[,14], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 15, sample(post_spike$tau_plot[,15], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 16, sample(post_spike$tau_plot[,16], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 17, sample(post_spike$tau_plot[,17], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 18, sample(post_spike$tau_plot[,18], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 19, sample(post_spike$tau_plot[,19], size = length(data$FLW_STAT_T)),ifelse(data$plot_index == 20, sample(post_spike$tau_plot[,20], size = length(data$FLW_STAT_T)), NA))))))))))))))))))))
  data$yearendo_effect <- ifelse(data$year_t_index == 1 & data$endo_index == 1, sample(post_spike$tau_year[,1,1], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 1 & data$endo_index == 2, sample(post_spike$tau_year[,2,1], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 2 & data$endo_index == 1, sample(post_spike$tau_year[,1,2], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 2 & data$endo_index == 2, sample(post_spike$tau_year[,2,2], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 3 & data$endo_index == 1, sample(post_spike$tau_year[,1,3], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 3 & data$endo_index == 2, sample(post_spike$tau_year[,2,3], size = length(data$FLW_STAT_T)) , ifelse(data$year_t_index == 4 & data$endo_index == 1, sample(post_spike$tau_year[,1,4], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 4 & data$endo_index == 2, sample(post_spike$tau_year[,2,4], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 5 & data$endo_index == 1, sample(post_spike$tau_year[,1,5], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 5 & data$endo_index == 2, sample(post_spike$tau_year[,2,5], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 6 & data$endo_index == 1, sample(post_spike$tau_year[,1,6], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 6 & data$endo_index == 2, sample(post_spike$tau_year[,2,6], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 7 & data$endo_index == 1, sample(post_spike$tau_year[,1,7], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 7 & data$endo_index == 2, sample(post_spike$tau_year[,2,7], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 8 & data$endo_index == 1, sample(post_spike$tau_year[,1,8], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 8 & data$endo_index == 2, sample(post_spike$tau_year[,2,8], size = length(data$FLW_STAT_T)) ,ifelse(data$year_t_index == 9 & data$endo_index == 1, sample(post_spike$tau_year[,1,9], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 9 & data$endo_index == 2, sample(post_spike$tau_year[,2,9], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 10 & data$endo_index == 1, sample(post_spike$tau_year[,1,10], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 10 & data$endo_index == 2, sample(post_spike$tau_year[,2,10], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 11 & data$endo_index == 1, sample(post_spike$tau_year[,1,11], size = length(data$FLW_STAT_T)),ifelse(data$year_t_index == 11 & data$endo_index == 2, sample(post_spike$tau_year[,2,11], size = length(data$FLW_STAT_T)),NA))))))))))))))))))))))
  data$seedperspike_mean <- sample(post_seed$beta[,1], size = length(data$FLW_STAT_T)) + sample(post_seed$beta[,2], size = length(data$FLW_STAT_T))*data$endo_01
  
  data$linpred <-  data$fixedeffect+data$plot_effect +data$yearendo_effect
  data$spikeperinf_prob = exp(data$linpred)
  
  data$spikeperinf_pred <- rztnbinom(n = nrow(data), size = mean(post_spike$phi), mu = data$spikeperinf_prob)
  
  
  s_to_s_data <- data %>% 
    mutate(seed_est = as.integer(FLW_STAT_T*FLW_COUNT_T*spikeperinf_pred*seedperspike_mean)) %>% 
    group_by(species, plot_fixed, plot_index, year_t, year_t_index, year_t1, year_t1_index, endo_01, endo_index) %>%
    summarize(tot_seed_t = as.integer(round(sum(seed_est, na.rm = TRUE))),
              tot_recruit_t1 = length((origin_01 == "1"& year_t == birth)),
              samplesize = n()) %>% 
    filter(tot_seed_t>tot_recruit_t1)

  return(s_to_s_data)
}


AGPE_s_to_s_data <- est_seed(sm_seed = sm_seedAGPE, sm_spike = sm_spikeAGPE, data = AGPE_surv_data)
write_csv(AGPE_s_to_s_data, "~/Dropbox/EndodemogData/Fulldataplusmetadata/AGPE_s_to_s_data.csv")
ELRI_s_to_s_data <- est_seed(sm_seed = sm_seedELRI, sm_spike = sm_spikeELRI, data = ELRI_surv_data)
write_csv(ELRI_s_to_s_data, "~/Dropbox/Endodemogdata/fulldataplusmetadata/ELRI_s_to_s_data.csv")
ELVI_s_to_s_data <- est_seed(sm_seed = sm_seedELVI, sm_spike = sm_spikeELVI, data = ELVI_surv_data)
write_csv(ELVI_s_to_s_data, "~/Dropbox/Endodemogdata/fulldataplusmetadata/ELVI_s_to_s_data.csv")
FESU_s_to_s_data <- est_seed(sm_seed = sm_seedFESU, sm_spike = sm_spikeFESU, data = FESU_surv_data)
write_csv(FESU_s_to_s_data, "~/Dropbox/Endodemogdata/fulldataplusmetadata/FESU_s_to_s_data.csv")
LOAR_s_to_s_data <- est_seed(sm_seed = sm_seedLOAR, sm_spike = sm_spikeLOAR, data = LOAR_surv_data)
write_csv(LOAR_s_to_s_data, "~/Dropbox/Endodemogdata/fulldataplusmetadata/LOAR_s_to_s_data.csv")
POAL_s_to_s_data <- est_seed(sm_seed = sm_seedPOAL, sm_spike = sm_spikePOAL, data = POAL_surv_data)
write_csv(POAL_s_to_s_data, "~/Dropbox/Endodemogdata/fulldataplusmetadata/POAL_s_to_s_data.csv")
POSY_s_to_s_data <- est_seed(sm_seed = sm_seedPOSY, sm_spike = sm_spikePOSY, data = POSY_surv_data)
write_csv(POSY_s_to_s_data, "~/Dropbox/Endodemogdata/fulldataplusmetadata/POSY_s_to_s_data.csv")

# Create data lists to be used for the Stan model
AGPE_s_to_s_data_list <- list(tot_recruit_t1 = AGPE_s_to_s_data$tot_recruit_t1,
                              tot_seed_t = AGPE_s_to_s_data$tot_seed_t,
                              endo_01 = AGPE_s_to_s_data$endo_01,
                              endo_index = AGPE_s_to_s_data$endo_index,
                              year_t = AGPE_s_to_s_data$year_t_index,
                              plot = AGPE_s_to_s_data$plot_index,
                              N = nrow(AGPE_s_to_s_data),
                              K = 3L,
                              lowerlimit = as.integer(min(AGPE_s_to_s_data$tot_recruit_t1)),
                              nYear = 11L,
                              nPlot = 10L,
                              nEndo =   2L)
str(AGPE_s_to_s_data_list)

ELRI_s_to_s_data_list <- list(tot_recruit_t1 = ELRI_s_to_s_data$tot_recruit_t1,
                              tot_seed_t = ELRI_s_to_s_data$tot_seed_t,
                              endo_01 = ELRI_s_to_s_data$endo_01,
                              endo_index = ELRI_s_to_s_data$endo_index,
                              year_t = ELRI_s_to_s_data$year_t_index,
                              plot = ELRI_s_to_s_data$plot_index,
                              N = nrow(ELRI_s_to_s_data),
                              K = 3L,
                              lowerlimit = as.integer(min(ELRI_s_to_s_data$tot_recruit_t1)),
                              nYear = 11L,
                              nPlot = 10L,
                              nEndo =   2L)
str(ELRI_s_to_s_data_list)

ELVI_s_to_s_data_list <- list(tot_recruit_t1 = ELVI_s_to_s_data$tot_recruit_t1,
                              tot_seed_t = ELVI_s_to_s_data$tot_seed_t,
                              endo_01 = ELVI_s_to_s_data$endo_01,
                              endo_index = ELVI_s_to_s_data$endo_index,
                              year_t = ELVI_s_to_s_data$year_t_index,
                              plot = ELVI_s_to_s_data$plot_index,
                              N = nrow(ELVI_s_to_s_data),
                              K = 3L,
                              lowerlimit = as.integer(min(ELVI_s_to_s_data$tot_recruit_t1)),
                              nYear = 11L,
                              nPlot = 10L,
                              nEndo =   2L)
str(ELVI_s_to_s_data_list)

FESU_s_to_s_data_list <- list(tot_recruit_t1 = FESU_s_to_s_data$tot_recruit_t1,
                              tot_seed_t = FESU_s_to_s_data$tot_seed_t,
                              endo_01 = FESU_s_to_s_data$endo_01,
                              endo_index = FESU_s_to_s_data$endo_index,
                              year_t = FESU_s_to_s_data$year_t_index,
                              plot = FESU_s_to_s_data$plot_index,
                              N = nrow(FESU_s_to_s_data),
                              K = 3L,
                              lowerlimit = as.integer(min(FESU_s_to_s_data$tot_recruit_t1)),
                              nYear = 11L,
                              nPlot = 10L,
                              nEndo =   2L)
str(FESU_s_to_s_data_list)

LOAR_s_to_s_data_list <- list(tot_recruit_t1 = LOAR_s_to_s_data$tot_recruit_t1,
                              tot_seed_t = LOAR_s_to_s_data$tot_seed_t,
                              endo_01 = LOAR_s_to_s_data$endo_01,
                              endo_index = LOAR_s_to_s_data$endo_index,
                              year_t = LOAR_s_to_s_data$year_t_index,
                              plot = LOAR_s_to_s_data$plot_index,
                              N = nrow(LOAR_s_to_s_data),
                              K = 3L,
                              lowerlimit = as.integer(min(LOAR_s_to_s_data$tot_recruit_t1)),
                              nYear = 11L,
                              nPlot = 10L,
                              nEndo =   2L)
str(LOAR_s_to_s_data_list)

POAL_s_to_s_data_list <- list(tot_recruit_t1 = POAL_s_to_s_data$tot_recruit_t1,
                              tot_seed_t = POAL_s_to_s_data$tot_seed_t,
                              endo_01 = POAL_s_to_s_data$endo_01,
                              endo_index = POAL_s_to_s_data$endo_index,
                              year_t = POAL_s_to_s_data$year_t_index,
                              plot = POAL_s_to_s_data$plot_index,
                              N = nrow(POAL_s_to_s_data),
                              K = 3L,
                              lowerlimit = as.integer(min(POAL_s_to_s_data$tot_recruit_t1)),
                              nYear = 11L,
                              nPlot = 18L,
                              nEndo =   2L)
str(POAL_s_to_s_data_list)
POSY_s_to_s_data_list <- list(tot_recruit_t1 = POSY_s_to_s_data$tot_recruit_t1,
                              tot_seed_t = POSY_s_to_s_data$tot_seed_t,
                              endo_01 = POSY_s_to_s_data$endo_01,
                              endo_index = POSY_s_to_s_data$endo_index,
                              year_t = POSY_s_to_s_data$year_t_index,
                              plot = POSY_s_to_s_data$plot_index,
                              N = nrow(POSY_s_to_s_data),
                              K = 3L,
                              lowerlimit = as.integer(min(POSY_s_to_s_data$tot_recruit_t1)),
                              nYear = 11L,
                              nPlot = 20L,
                              nEndo =   2L)
str(POSY_s_to_s_data_list)



#########################################################################################################
# GLMM for tot_recruit_t1~ tot_seed_t1 +Endo  with year and plot random effects------------------------------
#########################################################################################################
## run this code recommended to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)
## MCMC settings
ni <- 10000
nb <- 5000
nc <- 3

# Stan model -------------
## here is the Stan model

sink("endodemog_s_to_s.stan")
cat("
    data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> K;                       // number of predictors
    int<lower=0> lowerlimit;                         //lower limit for truncated negative binomial
    int<lower=0> nYear;                       // number of years
    int<lower=0, upper=11> year_t[N];         // year of observation
    int<lower=0> nEndo;                       // number of endo treatments
    int<lower=1, upper=2> endo_index[N];          // index for endophyte effect
    int<lower=0> nPlot;                         // number of plots
    int<lower=0> plot[N];                   // plot of observation
    int<lower=0> tot_recruit_t1[N];      // total recruits into the plot (response)
    int<lower=0> tot_seed_t[N];             // total seeds in plot (predictor)
    int<lower=0,upper=1> endo_01[N];            // plant endophyte status (predictor)
    }

    parameters {
    vector[2] beta;                     // predictor parameters
    vector[nYear] tau_year[nEndo];      // random year effect
    real<lower=0> sigma_e[nEndo];        //year variance by endophyte effect
    vector[nPlot] tau_plot;        // random plot effect
    real<lower=0> sigma_p;          // plot variance effect
    }
    
    transformed parameters{
    // Linear Predictor
    vector[N] mu;
    
    for(n in 1:N){
      mu[n] = beta[1] + beta[2]*endo_01[n]
      + tau_year[endo_index[n],year_t[n]]
      + tau_plot[plot[n]];
    }
    
    }

    model {
  
    // Priors
    beta ~ normal(0,1);      // prior for predictor intercepts
    tau_plot ~ normal(0,sigma_p);   // prior for plot random effects
    to_vector(tau_year[1]) ~ normal(0,sigma_e[1]);   // prior for E- year random effects
    to_vector(tau_year[2]) ~ normal(0,sigma_e[2]);   // prior for E+ year random effects
    
    // Likelihood
    tot_recruit_t1 ~ binomial_logit(tot_seed_t, mu);
    }
    
    generated quantities{
    }
  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_s_to_s.stan")



## Run the model by calling stan()
## and save the output to .rds files so that they can be called laters

smAGPE <- stan(file = "endodemog_s_to_s.stan", data = AGPE_s_to_s_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE, control = list(adapt_delta = 0.99, max_treedepth = 20))
# saveRDS(smAGPE, file = "endodemog_s_to_s_AGPE.rds")

smELRI <- stan(file = "endodemog_s_to_s.stan", data = ELRI_s_to_s_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE, control = list(adapt_delta = 0.99, max_treedepth = 20))
# saveRDS(smELRI, file = "endodemog_s_to_s_ELRI.rds")

smELVI <- stan(file = "endodemog_s_to_s.stan", data = ELVI_s_to_s_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE, control = list(adapt_delta = 0.99, max_treedepth = 20))
# saveRDS(smELVI, file = "endodemog_s_to_s_ELVI.rds")

smFESU <- stan(file = "endodemog_s_to_s.stan", data = FESU_s_to_s_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE, control = list(adapt_delta = 0.99, max_treedepth = 20))
# saveRDS(smFESU, file = "endodemog_s_to_s_FESU.rds")

smLOAR <- stan(file = "endodemog_s_to_s.stan", data = LOAR_s_to_s_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE, control = list(adapt_delta = 0.99, max_treedepth = 20))
# saveRDS(smLOAR, file = "endodemog_s_to_s_LOAR.rds")

smPOAL <- stan(file = "endodemog_s_to_s.stan", data = POAL_s_to_s_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE, control = list(adapt_delta = 0.99, max_treedepth = 20))
# saveRDS(smPOAL, file = "endodemog_s_to_s_POAL.rds")

smPOSY <- stan(file = "endodemog_s_to_s.stan", data = POSY_s_to_s_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE, control = list(adapt_delta = 0.99, max_treedepth = 20))
# saveRDS(smPOSY, file = "endodemog_s_to_s_POSY.rds")

print(sm)
summary(sm)
print(sm, pars = "sigma_e")





## to read in model output without rerunning models
smPOAL <- readRDS(file = "model_run_MAR7/endodemog_s_to_s_POAL.rds")
smPOSY <- readRDS(file = "model_run_MAR7/endodemog_s_to_s_POSY.rds")
smLOAR <- readRDS(file = "model_run_MAR7/endodemog_s_to_s_LOAR.rds")
smELVI <- readRDS(file = "model_run_MAR7/endodemog_s_to_s_ELVI.rds")
smELRI <- readRDS(file = "model_run_MAR7/endodemog_s_to_s_ELRI.rds")
smFESU <- readRDS(file = "model_run_MAR7/endodemog_s_to_s_FESU.rds")
smAGPE <- readRDS(file = "model_run_MAR7/endodemog_s_to_s_AGPE.rds")


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


# This function extracts the posterior draws and generates replicate data for each given model
prediction <- function(data, fit, n_post_draws){
  post <- rstan::extract(fit)
  mu <- post$mu
  yrep <- matrix(nrow = n_post_draws, ncol = data$N)
  for(n in 1:n_post_draws){
    yrep[n,] <- rbinom(n = data$N, size = data$tot_seed_t, prob = invlogit(mu[n,]))
  }
  out <- list(yrep, mu)
  names(out) <- c("yrep", "mu")
  return(out)
}

# apply the function for each species
AGPE_s_to_s_yrep <- prediction(data = AGPE_s_to_s_data_list, fit = smAGPE, n_post_draws = 500)
ELRI_s_to_s_yrep <- prediction(data = ELRI_s_to_s_data_list, fit = smELRI, n_post_draws = 500)
ELVI_s_to_s_yrep <- prediction(data = ELVI_s_to_s_data_list, fit = smELVI, n_post_draws = 500)
FESU_s_to_s_yrep <- prediction(data = FESU_s_to_s_data_list, fit = smFESU, n_post_draws = 500)
LOAR_s_to_s_yrep <- prediction(data = LOAR_s_to_s_data_list, fit = smLOAR, n_post_draws = 500)
POAL_s_to_s_yrep <- prediction(data = POAL_s_to_s_data_list, fit = smPOAL, n_post_draws = 500)
POSY_s_to_s_yrep <- prediction(data = POSY_s_to_s_data_list, fit = smPOSY, n_post_draws = 500)


# overlay 100 replicates over the actual dataset
ppc_dens_overlay( y = AGPE_s_to_s_data_list$tot_recruit_t1, yrep = AGPE_s_to_s_yrep$yrep[1:100,]) + ggtitle("AGPE")

ppc_dens_overlay( y = ELRI_s_to_s_data_list$tot_recruit_t1, yrep = ELRI_s_to_s_yrep$yrep[1:100,]) + ggtitle("ELRI")

ppc_dens_overlay( y = ELVI_s_to_s_data_list$tot_recruit_t1, yrep = ELVI_s_to_s_yrep$yrep[1:100,]) + ggtitle("ELVI")

ppc_dens_overlay( y = FESU_s_to_s_data_list$tot_recruit_t1, yrep = FESU_s_to_s_yrep$yrep[1:100,]) + ggtitle("FESU")

ppc_dens_overlay( y = LOAR_s_to_s_data_list$tot_recruit_t1, yrep = LOAR_s_to_s_yrep$yrep[1:100,]) + ggtitle("LOAR")

ppc_dens_overlay( y = POAL_s_to_s_data_list$tot_recruit_t1, yrep = POAL_s_to_s_yrep$yrep[1:100,]) + ggtitle("POAL")

ppc_dens_overlay( y = POSY_s_to_s_data_list$tot_recruit_t1, yrep = POSY_s_to_s_yrep$yrep[1:100,]) + ggtitle("POSY")



