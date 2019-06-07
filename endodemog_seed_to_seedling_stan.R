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

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }
#############################################################################################
####### Data manipulation to prepare data as lists for Stan models------------------
#############################################################################################
## Load in full data frame
LTREB_endodemog <-read.csv(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_long.csv")
LTREB_endodemog <- LTREB_data_seed_est

#"C:/Users/MillerLab/Desktop/Endodemog-master/endo_demog_long.csv"
#"/Users/joshuacfowler/Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_long.csv")
str(LTREB_endodemog)
dim(LTREB_endodemog)

## Clean up the main data frame for NA's, other small data entry errors, and change standardize the coding for variables.
LTREB_data <- LTREB_endodemog %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(size_t1, logsize_t1 = log(size_t1)) %>%  
  mutate(surv_t1 = as.integer(recode(surv_t1, "0" = 0, "1" =1, "2" = 1, "4" = 1))) %>% 
  mutate(endo_01 = as.integer(case_when(endo == "0" | endo == "minus" ~ 0,
                                        endo == "1"| endo =="plus" ~ 1))) %>% 
  filter(!is.na(endo_01)) %>% 
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

LTREB_data1 <- LTREB_data %>%
  group_by(species, plot_fixed, year_t, year_t_index, year_t1, year_t1_index, endo_01, endo_index) %>% 
  summarize(tot_seed_t = as.integer(round(sum(seed_t, na.rm = TRUE))),
            tot_recruit_t1 = length((origin_01 == "1"& year_t == birth)),
            samplesize = n()) %>% 
  filter(tot_seed_t>tot_recruit_t1)

dim(LTREB_data1)
# View(LTREB_data1)

#########################################################################################################
# Split up the main dataframe by species and recode plot to be used as an index for each species
AGPE_data <- LTREB_data1 %>% 
  filter(species == "AGPE") %>% 
  mutate(year_t_index_new = case_when(year_t_index ==2~1, year_t_index ==3~2, year_t_index ==4~3, year_t_index ==5~4, year_t_index ==6~5, year_t_index ==7~6, year_t_index ==8~7, year_t_index ==9~8)) %>% 
  mutate(plot_index = case_when(plot_fixed == 111 ~ 1, plot_fixed ==112 ~2, plot_fixed ==113 ~ 3, plot_fixed ==114 ~ 4, plot_fixed ==115 ~ 5,plot_fixed == 116 ~ 6, plot_fixed ==117~ 7, plot_fixed ==118 ~ 8, plot_fixed ==119 ~ 9,plot_fixed == 120 ~ 10))
ELRI_data <- LTREB_data1 %>% 
  filter(species == "ELRI") %>% 
  mutate(year_t_index_new = case_when(year_t_index ==2~1, year_t_index ==3~2, year_t_index ==4~3, year_t_index ==5~4, year_t_index ==6~5, year_t_index ==7~6, year_t_index ==8~7, year_t_index ==9~8)) %>% 
  mutate(plot_index = case_when(plot_fixed == 101 ~ 1, plot_fixed ==102 ~2, plot_fixed ==103 ~ 3, plot_fixed ==104 ~ 4, plot_fixed ==105 ~ 5,plot_fixed == 106 ~ 6, plot_fixed ==107~ 7, plot_fixed ==108 ~ 8, plot_fixed ==109 ~ 9,plot_fixed == 110 ~ 10)) 
ELVI_data <- LTREB_data1 %>% 
  filter(species == "ELVI") %>% 
  mutate(year_t_index_new = case_when(year_t_index ==3~1, year_t_index ==4~2, year_t_index ==5~3, year_t_index ==6~4, year_t_index ==7~5, year_t_index ==8~6, year_t_index ==9~7)) %>% 
  mutate(plot_index = case_when(plot_fixed == 91 ~ 1, plot_fixed ==92 ~2, plot_fixed ==93 ~ 3, plot_fixed ==94 ~ 4, plot_fixed ==95 ~ 5,plot_fixed == 96 ~ 6, plot_fixed ==97~ 7, plot_fixed ==98 ~ 8, plot_fixed ==99 ~ 9,plot_fixed == 100 ~ 10, plot_fixed == 101 ~ 11))
FESU_data <- LTREB_data1 %>% 
  filter(species == "FESU") %>% 
  mutate(year_t_index_new = case_when(year_t_index ==3~1, year_t_index ==4~2, year_t_index ==5~3, year_t_index ==6~4, year_t_index ==7~5, year_t_index ==8~6, year_t_index ==9~7)) %>% 
  mutate(plot_index = case_when(plot_fixed == 121 ~ 1, plot_fixed ==122 ~2, plot_fixed ==123 ~ 3, plot_fixed ==124 ~ 4, plot_fixed ==125 ~ 5,plot_fixed == 126 ~ 6, plot_fixed ==127~ 7, plot_fixed ==128 ~ 8, plot_fixed ==129 ~ 9,plot_fixed == 130 ~ 10)) 
LOAR_data <- LTREB_data1 %>% 
  filter(species == "LOAR") %>% 
  mutate(year_t_index_new = case_when(year_t_index ==2~1, year_t_index ==3~2, year_t_index ==4~3, year_t_index ==5~4, year_t_index ==6~5, year_t_index ==7~6, year_t_index ==8~7, year_t_index ==9~8)) %>% 
  mutate(plot_index = case_when(plot_fixed == 31 ~ 1, plot_fixed ==32 ~2, plot_fixed ==33 ~ 3, plot_fixed ==34 ~ 4, plot_fixed ==35 ~ 5,plot_fixed == 36 ~ 6, plot_fixed ==37~ 7, plot_fixed ==38 ~ 8, plot_fixed ==39 ~ 9,plot_fixed == 40 ~ 10, plot_fixed == 41 ~11,  plot_fixed == 42 ~12, plot_fixed == 43 ~13, plot_fixed == 44 ~14,)) 
POAL_data <- LTREB_data1 %>% 
  filter(species == "POAL") %>% 
  mutate(year_t_index_new = case_when(year_t_index ==2~1, year_t_index ==3~2, year_t_index ==4~3, year_t_index ==5~4, year_t_index ==6~5, year_t_index ==7~6, year_t_index ==8~7, year_t_index ==9~8)) %>% 
  mutate(plot_index = case_when(plot_fixed == 3 ~ 1, plot_fixed ==4 ~2, plot_fixed ==8 ~ 3, plot_fixed ==9 ~ 4, plot_fixed ==10 ~ 5,plot_fixed == 11 ~ 6, plot_fixed ==15~ 7, plot_fixed ==16~ 8, plot_fixed ==17~ 9, plot_fixed ==151 ~ 10, plot_fixed ==152 ~ 11,plot_fixed == 153 ~ 12, plot_fixed == 154 ~13,  plot_fixed == 155 ~14, plot_fixed == 156 ~15, plot_fixed == 157 ~16,plot_fixed == 158 ~17, plot_fixed == 19 ~18, )) 
POSY_data <- LTREB_data1 %>% 
  filter(species == "POSY") %>% 
  mutate(year_t_index_new = case_when(year_t_index ==2~1, year_t_index ==3~2, year_t_index ==4~3, year_t_index ==5~4, year_t_index ==6~5, year_t_index ==7~6, year_t_index ==8~7, year_t_index ==9~8)) %>% 
  mutate(plot_index = case_when(plot_fixed == 1 ~ 1, plot_fixed ==2 ~2, plot_fixed ==5 ~ 3, plot_fixed ==6 ~ 4, plot_fixed ==7 ~ 5,plot_fixed == 12 ~ 6, plot_fixed ==13~ 7, plot_fixed ==14 ~ 8, plot_fixed ==18 ~ 9,plot_fixed == 20 ~ 10, plot_fixed == 141 ~ 11, plot_fixed == 142 ~12,  plot_fixed == 143 ~13, plot_fixed == 144 ~14, plot_fixed == 145 ~15,plot_fixed == 146 ~16,plot_fixed == 147 ~17,plot_fixed == 148 ~18,plot_fixed == 149 ~19,plot_fixed == 150 ~20,))



# Create data lists to be used for the Stan model
AGPE_s_to_s_data_list <- list(tot_recruit_t1 = AGPE_data$tot_recruit_t1,
                            tot_seed_t = AGPE_data$tot_seed_t,
                            endo_01 = AGPE_data$endo_01,
                            endo_index = AGPE_data$endo_index,
                            year_t = AGPE_data$year_t_index_new,
                            plot = AGPE_data$plot_index,
                            N = nrow(AGPE_data),
                            K = 3L,
                            lowerlimit = as.integer(min(AGPE_data$tot_recruit_t1)),
                            nYear = length(unique(AGPE_data$year_t_index)),
                            nPlot = length(unique(AGPE_data$plot_index)),
                            nEndo =   length(unique(AGPE_data$endo_01)))
str(AGPE_s_to_s_data_list)

ELRI_s_to_s_data_list <- list(tot_recruit_t1 = ELRI_data$tot_recruit_t1,
                            tot_seed_t = ELRI_data$tot_seed_t,
                            endo_01 = ELRI_data$endo_01,
                            endo_index = ELRI_data$endo_index,
                            year_t = ELRI_data$year_t_index_new,
                            plot = ELRI_data$plot_index,
                            N = nrow(ELRI_data),
                            K = 3L,
                            lowerlimit = 1L,
                            nYear = length(unique(ELRI_data$year_t_index)),
                            nPlot = length(unique(ELRI_data$plot_index)),
                            nEndo =   length(unique(ELRI_data$endo_01)))
str(ELRI_s_to_s_data_list)

ELVI_s_to_s_data_list <- list(tot_recruit_t1 = ELVI_data$tot_recruit_t1,
                            tot_seed_t = ELVI_data$tot_seed_t,
                            endo_01 = ELVI_data$endo_01,
                            endo_index = ELVI_data$endo_index,
                            year_t = ELVI_data$year_t_index_new,
                            plot = ELVI_data$plot_index,
                            N = nrow(ELVI_data),
                            K = 3L,
                            lowerlimit = as.integer(min(ELVI_data$tot_recruit_t1)),
                            nYear = length(unique(ELVI_data$year_t_index)),
                            nPlot = length(unique(ELVI_data$plot_index)),
                            nEndo =   length(unique(ELVI_data$endo_01)))
str(ELVI_s_to_s_data_list)

FESU_s_to_s_data_list <- list(tot_recruit_t1 = FESU_data$tot_recruit_t1,
                              tot_seed_t = FESU_data$tot_seed_t,
                            endo_01 = FESU_data$endo_01,
                            endo_index = FESU_data$endo_index,
                            year_t = FESU_data$year_t_index_new,
                            plot = FESU_data$plot_index,
                            N = nrow(FESU_data),
                            K = 3L,
                            lowerlimit = as.integer(min(FESU_data$tot_recruit_t1)),
                            nYear = length(unique(FESU_data$year_t_index)),
                            nPlot = length(unique(FESU_data$plot_index)),
                            nEndo =   length(unique(FESU_data$endo_01)))
str(FESU_s_to_s_data_list)

LOAR_s_to_s_data_list <- list(tot_recruit_t1 = LOAR_data$tot_recruit_t1,
                              tot_seed_t = LOAR_data$tot_seed_t,
                            endo_01 = LOAR_data$endo_01,
                            endo_index = LOAR_data$endo_index,
                            year_t = LOAR_data$year_t_index_new,
                            plot = LOAR_data$plot_index,
                            N = nrow(LOAR_data),
                            K = 3L,
                            lowerlimit = as.integer(min(LOAR_data$tot_recruit_t1)),
                            nYear = length(unique(LOAR_data$year_t_index)),
                            nPlot = length(unique(LOAR_data$plot_index)),
                            nEndo =   length(unique(LOAR_data$endo_01)))
str(LOAR_s_to_s_data_list)

POAL_s_to_s_data_list <- list(tot_recruit_t1 = POAL_data$tot_recruit_t1,
                              tot_seed_t = POAL_data$tot_seed_t,
                            endo_01 = POAL_data$endo_01,
                            endo_index = POAL_data$endo_index,
                            year_t = POAL_data$year_t_index_new,
                            plot = POAL_data$plot_index,
                            N = nrow(POAL_data),
                            K = 3L,
                            lowerlimit = as.integer(min(POAL_data$tot_recruit_t1)),
                            nYear = length(unique(POAL_data$year_t_index)),
                            nPlot = length(unique(POAL_data$plot_index)),
                            nEndo =   length(unique(POAL_data$endo_01)))
str(POAL_s_to_s_data_list)

POSY_s_to_s_data_list <- list(tot_recruit_t1 = POSY_data$tot_recruit_t1,
                              tot_seed_t = POSY_data$tot_seed_t,
                            endo_01 = POSY_data$endo_01,
                            endo_index = POSY_data$endo_index,
                            year_t = POSY_data$year_t_index_new,
                            plot = POSY_data$plot_index,
                            N = nrow(POSY_data),
                            K = 3L,
                            lowerlimit = as.integer(min(POSY_data$tot_recruit_t1)),
                            nYear = length(unique(POSY_data$year_t_index)),
                            nPlot = length(unique(POSY_data$plot_index)),
                            nEndo =   length(unique(POSY_data$endo_01)))
str(POSY_s_to_s_data_list)


#########################################################################################################
# GLMM for tot_recruit_t1~ tot_seed_t1 +Endo  with year and plot random effects------------------------------
#########################################################################################################
## run this code recommended to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)
## MCMC settings
ni <- 2000
nb <- 1000
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
    real<lower=0> phi;
    }

    model {
    
    vector[N] mu;
    
    // Linear Predictor
    for(n in 1:N){
    mu[n] = beta[1] + beta[2]*endo_01[n]
    + tau_year[endo_index[n],year_t[n]]
    + tau_plot[plot[n]];
    }
    
    // Priors
    beta ~ normal(0,100);      // prior for predictor intercepts
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
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smAGPE, file = "endodemog_s_to_s_AGPE.rds")

smELRI <- stan(file = "endodemog_s_to_s.stan", data = ELRI_s_to_s_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELRI, file = "endodemog_s_to_s_ELRI.rds")

smELVI <- stan(file = "endodemog_s_to_s.stan", data = ELVI_s_to_s_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smELVI, file = "endodemog_s_to_s_ELVI.rds")

smFESU <- stan(file = "endodemog_s_to_s.stan", data = FESU_s_to_s_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smFESU, file = "endodemog_s_to_s_FESU.rds")

smLOAR <- stan(file = "endodemog_s_to_s.stan", data = LOAR_s_to_s_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smLOAR, file = "endodemog_s_to_s_LOAR.rds")

smPOAL <- stan(file = "endodemog_s_to_s.stan", data = POAL_s_to_s_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOAL, file = "endodemog_s_to_s_POAL.rds")

smPOSY <- stan(file = "endodemog_s_to_s.stan", data = POSY_s_to_s_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
# saveRDS(smPOSY, file = "endodemog_s_to_s_POSY.rds")

print(sm)
summary(sm)
print(sm, pars = "sigma_e")





## to read in model output without rerunning models
smPOAL <- readRDS(file = "model_run_MAR7/endodemog_surv_POAL.rds")
smPOSY <- readRDS(file = "model_run_MAR7/endodemog_surv_POSY.rds")
smLOAR <- readRDS(file = "model_run_MAR7/endodemog_surv_LOAR.rds")
smELVI <- readRDS(file = "model_run_MAR7/endodemog_surv_ELVI.rds")
smELRI <- readRDS(file = "model_run_MAR7/endodemog_surv_ELRI.rds")
smFESU <- readRDS(file = "model_run_MAR7/endodemog_surv_FESU.rds")
smAGPE <- readRDS(file = "model_run_MAR7/endodemog_surv_AGPE.rds")



