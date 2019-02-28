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
## Load in full data frame
LTREB_endodemog <- 
  read.csv(file = "/Users/joshuacfowler/Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_long.csv")


#"C:/Users/MillerLab/Desktop/Endodemog-master/endo_demog_long.csv"
#"/Users/joshuacfowler/Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_long.csv")
str(LTREB_endodemog)
dim(LTREB_endodemog)


## Clean up the main data frame for NA's, other small data entry errors, and change standardize the coding for variables.
LTREB_data <- LTREB_endodemog %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(size_t1, logsize_t1 = log(size_t1)) %>%  
  mutate(surv_t1 = as.integer(recode(surv_t1, "0" = 0, "1" =1, "2" = 1, "4" = 1))) %>% 
  mutate(endo_01 = case_when(endo == "0" | endo == "minus" ~ 0,
                             endo == "1"| endo =="plus" ~ 1)) %>% 
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
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(logsize_t >= 0) %>% 
  filter(!is.na(endo_01))
  
dim(LTREB_data1)

# LTREB_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01 + species)^2 + origin_01 
                                 # , data = LTREB_data1)
# Xs <- model.matrix(surv_t1 ~ (logsize_t + endo_01 + species)^2 + origin_01 
                                 # , data =LTREB_for_matrix)

# 
# # Create data list for Stan model
# # # LTREB_surv_data_list <- list(surv_t1 = LTREB_data1$surv_t1,
# #                              Xs = Xs,    
# #                              endo_index = LTREB_data1$endo_index,
# #                              species_index = LTREB_data1$species_index,
# #                              spp_endo_index = LTREB_data1$spp_endo_index,
# #                              spp_year_index = LTREB_data1$spp_year_index,
# #                              year_t = LTREB_data1$year_t_index, 
# #                              N = nrow(LTREB_data1), 
# #                              K = ncol(Xs), 
# #                              nyear = length(unique(LTREB_data1$year_t_index)), 
# #                              nEndo =   length(unique(LTREB_data1$endo_01)),
# #                              nSpp = length(unique(LTREB_data1$species_index)))
# # /*** /*
# 
# str(LTREB_surv_data_list)
# 
# LTREB_sample <- sample_n(LTREB_data1, 1000)
# sample_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01 + species)^3 + origin_01 
#                                 , data = LTREB_sample)
# Xs <- model.matrix(surv_t1 ~ (logsize_t + endo_01 + species)^3 + origin_01 
#                    , data =sample_for_matrix)
# 
# 
# # Create sample data list for Stan model
# sample_surv_data_list <- list(surv_t1 = LTREB_sample$surv_t1,
#                              Xs = Xs,    
#                              endo_index = LTREB_sample$endo_index,
#                              species_index = LTREB_sample$species_index,
#                              spp_endo_index = LTREB_sample$spp_endo_index,
#                              spp_year_index = LTREB_sample$spp_year_index,
#                              year_t = LTREB_sample$year_t_index, 
#                              N = nrow(LTREB_sample), 
#                              K = ncol(Xs), 
#                              nyear = length(unique(LTREB_sample$year_t_index)), 
#                              nEndo =   length(unique(LTREB_sample$endo_01)),
#                              nSpp = length(unique(LTREB_sample$species_index)))
# 
# 
# str(sample_surv_data_list)

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
AGPE_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01 + origin_01)^3
, data = AGPE_data)
AGPE_Xs <- model.matrix(surv_t1 ~ (logsize_t + endo_01 + origin_01)^3
, data =AGPE_for_matrix)

ELRI_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01 + origin_01)^3
                               , data = ELRI_data)
ELRI_Xs <- model.matrix(surv_t1 ~ (logsize_t + endo_01 + origin_01)^3
                        , data =ELRI_for_matrix)

ELVI_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01 + origin_01)^3
                               , data = ELVI_data)
ELVI_Xs <- model.matrix(surv_t1 ~ (logsize_t + endo_01 + origin_01)^3
                        , data =ELVI_for_matrix)

FESU_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01 + origin_01)^3
                               , data = FESU_data)
FESU_Xs <- model.matrix(surv_t1 ~ (logsize_t + endo_01 + origin_01)^3
                        , data =FESU_for_matrix)

LOAR_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01 + origin_01)^3
                               , data = LOAR_data)
LOAR_Xs <- model.matrix(surv_t1 ~ (logsize_t + endo_01 + origin_01)^3
                        , data =LOAR_for_matrix)

POAL_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01 + origin_01)^3
                               , data = POAL_data)
POAL_Xs <- model.matrix(surv_t1 ~ (logsize_t + endo_01 + origin_01)^3
                        , data =POAL_for_matrix)

POSY_for_matrix <- model.frame(surv_t1 ~ (logsize_t + endo_01 + origin_01)^3
                               , data = POSY_data)
POSY_Xs <- model.matrix(surv_t1 ~ (logsize_t + endo_01 + origin_01)^3
                        , data =POSY_for_matrix)



# Create data lists to be used for the Stan model
AGPE_surv_data_list <- list(surv_t1 = AGPE_data$surv_t1,
                             Xs = AGPE_Xs,
                             endo_index = AGPE_data$endo_index,
                             year_t = AGPE_data$year_t_index,
                             plot = AGPE_data$plot_index,
                             N = nrow(AGPE_data),
                             K = ncol(AGPE_Xs),
                             nYear = length(unique(AGPE_data$year_t_index)),
                             nPlot = length(unique(AGPE_data$plot_index)),
                             nEndo =   length(unique(AGPE_data$endo_01)))
str(AGPE_surv_data_list)

ELRI_surv_data_list <- list(surv_t1 = ELRI_data$surv_t1,
                            Xs = ELRI_Xs,
                            endo_index = ELRI_data$endo_index,
                            year_t = ELRI_data$year_t_index,
                            plot = ELRI_data$plot_index,
                            N = nrow(ELRI_data),
                            K = ncol(ELRI_Xs),
                            nYear = length(unique(ELRI_data$year_t_index)),
                            nPlot = length(unique(ELRI_data$plot_index)),
                            nEndo =   length(unique(ELRI_data$endo_01)))
str(ELRI_surv_data_list)

ELVI_surv_data_list <- list(surv_t1 = ELVI_data$surv_t1,
                            Xs = ELVI_Xs,
                            endo_index = ELVI_data$endo_index,
                            year_t = ELVI_data$year_t_index,
                            plot = ELVI_data$plot_index,
                            N = nrow(ELVI_data),
                            K = ncol(ELVI_Xs),
                            nYear = length(unique(ELVI_data$year_t_index)),
                            nPlot = length(unique(ELVI_data$plot_index)),
                            nEndo =   length(unique(ELVI_data$endo_01)))
str(ELVI_surv_data_list)

FESU_surv_data_list <- list(surv_t1 = FESU_data$surv_t1,
                            Xs = FESU_Xs,
                            endo_index = FESU_data$endo_index,
                            year_t = FESU_data$year_t_index,
                            plot = FESU_data$plot_index,
                            N = nrow(FESU_data),
                            K = ncol(FESU_Xs),
                            nYear = length(unique(FESU_data$year_t_index)),
                            nPlot = length(unique(FESU_data$plot_index)),
                            nEndo =   length(unique(FESU_data$endo_01)))
str(FESU_surv_data_list)

LOAR_surv_data_list <- list(surv_t1 = LOAR_data$surv_t1,
                            Xs = LOAR_Xs,
                            endo_index = LOAR_data$endo_index,
                            year_t = LOAR_data$year_t_index,
                            plot = LOAR_data$plot_index,
                            N = nrow(LOAR_data),
                            K = ncol(LOAR_Xs),
                            nYear = length(unique(LOAR_data$year_t_index)),
                            nPlot = length(unique(LOAR_data$plot_index)),
                            nEndo =   length(unique(LOAR_data$endo_01)))
str(LOAR_surv_data_list)

POAL_surv_data_list <- list(surv_t1 = POAL_data$surv_t1,
                            Xs = POAL_Xs,
                            endo_index = POAL_data$endo_index,
                            year_t = POAL_data$year_t_index,
                            plot = POAL_data$plot_index,
                            N = nrow(POAL_data),
                            K = ncol(POAL_Xs),
                            nYear = length(unique(POAL_data$year_t_index)),
                            nPlot = length(unique(POAL_data$plot_index)),
                            nEndo =   length(unique(POAL_data$endo_01)))
str(POAL_surv_data_list)

POSY_surv_data_list <- list(surv_t1 = POSY_data$surv_t1,
                            Xs = POSY_Xs,
                            endo_index = POSY_data$endo_index,
                            year_t = POSY_data$year_t_index,
                            plot = POSY_data$plot_index,
                            N = nrow(POSY_data),
                            K = ncol(POSY_Xs),
                            nYear = length(unique(POSY_data$year_t_index)),
                            nPlot = length(unique(POSY_data$plot_index)),
                            nEndo =   length(unique(POSY_data$endo_01)))
str(POSY_surv_data_list)

#########################################################################################################
# GLMM for Surv~ size +Endo + Origin  with year and plot random effects------------------------------
#########################################################################################################
## run this code recommended to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <- 100
nb <- 50
nc <- 1

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
    int<lower=1, upper=14> endo_index[N];          // index for endophyte effect
    int<lower=0> nPlot;                         // number of plots
    int<lower=0> plot[N];                   // plot of observation
    int<lower=0, upper=1> surv_t1[N];      // plant survival at time t+1 and target variable (response)
    matrix[N,K] Xs;                     //  predictor matrix 
    }
 
    parameters {
    vector[K] beta;                     // predictor parameters

    vector[nYear] tau_year[nEndo];      // random year effect
    real<lower=0> sigma_e[nEndo];        //year variance by endophyte effect
    vector[nPlot] tau_plot;        // random plot effect
    }

    model {
    
    vector[N] mu;
    
    // Linear Predictor
    for(n in 1:N){
    mu = Xs*beta
    + tau_year[endo_index[n],year_t[n]] + tau_plot[plot[n]];
    }
    
    // Priors
    beta ~ normal(0,100);      // prior for predictor intercepts
    sigma_e ~ gamma(2,.1);
    tau_plot ~ normal(0,100);
    for(e in 1:nEndo){
    to_vector(tau_year[e][]) ~ normal(0,sigma_e[e]); // prior for year random effects
    }
   
    // Likelihood
      surv_t1 ~ bernoulli_logit(mu);
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

stanmodel <- stanc("endodemog_surv.stan")



## Run the model by calling stan()
## and save the output to .rds files so that they can be called laters

sm <- stan(file = "endodemog_surv.stan", data = AGPE_surv_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)

print(sm)
summary(sm)
print(sm, pars = "sigma_0")







## Create dataframes with individual species datasets
POAL_data <- LTREB_endodemog %>% 
  filter(species == "POAL")
dim(POAL_data)

POSY_data <- LTREB_endodemog %>% 
  filter(species == "POSY")
dim(POSY_data)

LOAR_data <- LTREB_endodemog %>% 
  filter(species == "LOAR")
dim(LOAR_data)

ELVI_data <- LTREB_endodemog %>% 
  filter(species == "ELVI")
dim(ELVI_data)

ELRI_data <- LTREB_endodemog %>% 
  filter(species == "ELRI")
dim(ELRI_data)

FESU_data <- LTREB_endodemog %>% 
  filter(species == "FESU")
dim(FESU_data)

AGPE_data <- LTREB_endodemog %>% 
  filter(species == "AGPE")
dim(AGPE_data)

## Data manipulation for POAL to prepare as list for Stan
## create new column with log(size) and recode years and plot (including fixing values entered as R where they should be plot 16)
POAL_data1 <- POAL_data %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(year_t_index = as.factor(recode(year_t, '2007'=1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11))) %>% 
  mutate(year_t1_index = as.factor(recode(year_t1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11, '2018'=12))) %>% 
  mutate(year_t1_index = as.factor(recode(year_t1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11, '2018'=12))) %>% 
  mutate(plot_index = as.factor(recode(plot, 'R'=8, '3'=1, '4'=2, '8'=3, '9'=4, '10'=5, '11'=6, '15'=7, '16'=8, '17'=9, '19'=10, '151'=11, '152'=12, '153'=13, '154'=14, '155'=15, '156'=16, '157'=17, '158'=18)))
dim(POAL_data1)

POAL_data1 <- POAL_data1 %>% 
  filter(!is.na(surv_t1)) 
dim(POAL_data1)

POAL_surv_dat <- POAL_data1 %>% 
  select(surv_t1) %>% 
  filter(!is.na(surv_t1))
POAL_size_dat <- POAL_data1 %>% 
  select(logsize_t) %>% 
  filter(!is.na(logsize_t))
POAL_year_dat <- POAL_data1 %>% 
  select(year_t_index) %>% 
  filter(!is.na(year_t_index)) %>% 
  mutate(year_t_index = as.integer(year_t_index))
POAL_endo_dat <- POAL_data1 %>% 
  select(endo) %>% 
  filter(!is.na(endo)) %>% 
  mutate(endo1 = as.integer(as.integer(endo)-4)) %>% 
  mutate(endo_index = as.integer(endo1+1))   ## recoding endo to 0 or 1
POAL_origin_dat <- POAL_data1 %>% 
  select(origin) %>% 
  filter(!is.na(origin)) %>% 
  mutate(origin1 = as.integer(case_when(origin == "O" ~ 0,
                                        origin == "R" ~ 1,
                                        origin == 16 ~ 1)))
# Plot data not included as a random effect due to poor estimation
POAL_plot_dat <- POAL_data1 %>%
  select(plot, plot_index) %>% 
  filter(!is.na(plot_index)) %>% 
  mutate(plot_index = as.integer(plot_index))

dim(POAL_surv_dat)
dim(POAL_size_dat)
dim(POAL_year_dat)
dim(POAL_endo_dat)
dim(POAL_origin_dat)
dim(POAL_plot_dat)

# Data list for POAL
POAL_data_list <- list(surv_t1 = POAL_surv_dat$surv_t1, 
                       logsize_t = POAL_size_dat$logsize_t, 
                       endo = POAL_endo_dat$endo1, 
                       endo_index = POAL_endo_dat$endo_index,
                       origin = POAL_origin_dat$origin1,
                       year_t = POAL_year_dat$year_t_index, 
                       N = 3241L, K = 4L, nyear = 11L, nEndo = 2L)
str(POAL_data_list) 


## Data manipulation for POSY to prepare as list for Stan
## create new column with log(size) and recode years and origin (including fixing values entered as plot numbers where they should be recruit origin status)
POSY_data1 <- POSY_data %>% 
  filter(size_t !=0) %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(year_t_index = as.factor(recode(year_t, '2007'=1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11))) %>% 
  mutate(year_t1_index = as.factor(recode(year_t1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11, '2018'=12))) %>% 
  mutate(origin = as.factor(recode(origin, 'O' = 0, 'R' = 1, '1' = 1, '2'=1, '2'=1, '5'=1, '7'=1, '12'=1, '13'=1, '14'=1, '18'=1, '20'=1)))
dim(POSY_data1)


POSY_data1 <- POSY_data1 %>%
  filter(!is.na(size_t)) %>% 
  filter(!is.na(surv_t1)) 

dim(POSY_data1)
str(POSY_data1)
# View(POSY_data1)

POSY_surv_dat <- POSY_data1 %>% 
  select(surv_t1) %>% 
  filter(!is.na(surv_t1))
POSY_size_dat <- POSY_data1 %>% 
  select(logsize_t) %>% 
  filter(!is.na(logsize_t))
POSY_year_dat <- POSY_data1 %>% 
  select(year_t_index) %>% 
  filter(!is.na(year_t_index)) %>% 
  mutate(year_t_index = as.integer(year_t_index))
POSY_endo_dat <- POSY_data1 %>% 
  select(endo) %>% 
  filter(!is.na(endo)) %>% 
  mutate(endo1 = as.integer(as.integer(endo)-4)) %>% 
  mutate(endo_index = as.integer(endo1+1))   ## recoding endo to 0 or 1
POSY_origin_dat <- POSY_data1 %>% 
  select(origin) %>% 
  filter(!is.na(origin)) %>% 
  mutate(origin1 = as.integer(as.integer(origin)-1))

dim(POSY_surv_dat)
dim(POSY_size_dat)
dim(POSY_year_dat)
dim(POSY_endo_dat)
dim(POSY_origin_dat)

# Data list for POSY
POSY_data_list <- list(surv_t1 = POSY_surv_dat$surv_t1, 
                       logsize_t = POSY_size_dat$logsize_t, 
                       endo = POSY_endo_dat$endo1, 
                       endo_index = POSY_endo_dat$endo_index,
                       origin = POSY_origin_dat$origin1,
                       year_t = POSY_year_dat$year_t_index, 
                       N = 12913L, K = 4L, nyear = 11L, nEndo = 2L)
str(POSY_data_list) 

# Data manipulation for LOAR to prepare as list for Stan
## create new column with log(size) and recode years, plot, and origin 
LOAR_data1 <- LOAR_data %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(year_t_index = as.factor(recode(year_t, '2007'=1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11))) %>% 
  mutate(year_t1_index = as.factor(recode(year_t1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11, '2018'=12))) %>% 
  mutate(plot_index = as.factor(recode(plot, '31'=1, '32'=2, '33'=3, '34'=4, '35'=5, '36'=6, '37'=7, '38'=8, '39'=9, '40'=10, '41'=11, '42'=12, '43'=13, '44'=14))) %>% 
  mutate(origin = as.factor(recode(origin, 'O'=0, 'R'=1)))
dim(LOAR_data1)


LOAR_data1 <- LOAR_data1 %>% 
  filter(surv_t1 != 4) %>% 
  filter(surv_t1 != 2) %>% 
  filter(size_t !=0) %>% 
  filter(!is.na(surv_t1)) %>%
  filter(!is.na(size_t)) %>% 
  filter(endo != "")

dim(LOAR_data1)

LOAR_surv_dat <- LOAR_data1 %>% 
  select(surv_t1) %>% 
  filter(!is.na(surv_t1))
LOAR_size_dat <- LOAR_data1 %>% 
  select(logsize_t) %>% 
  filter(!is.na(logsize_t))
LOAR_year_dat <- LOAR_data1 %>% 
  select(year_t_index) %>% 
  filter(!is.na(year_t_index)) %>% 
  mutate(year_t_index = as.integer(year_t_index))
LOAR_endo_dat <- LOAR_data1 %>% 
  select(endo) %>% 
  filter(!is.na(endo)) %>% 
  filter(endo != "") %>% 
  mutate(endo1 = as.integer(as.integer(endo)-2)) %>%      ## recoding endo to 1 or 2 for the endo index
  mutate(endo_index = as.integer(endo1+1))   
LOAR_origin_dat <- LOAR_data1 %>% 
  select(origin) %>% 
  filter(!is.na(origin)) %>% 
  mutate(origin = as.integer(as.integer(origin)-1))

# Plot data not included as a random effect due to poor estimation
LOAR_plot_dat <- LOAR_data1 %>%
  select(plot, plot_index) %>% 
  filter(!is.na(plot_index)) %>% 
  mutate(plot_index = as.integer(plot_index))

dim(LOAR_surv_dat)
dim(LOAR_size_dat)
dim(LOAR_year_dat)
dim(LOAR_endo_dat)
dim(LOAR_origin_dat)
dim(LOAR_plot_dat)

# Data list for LOAR
LOAR_data_list <- list(surv_t1 = LOAR_surv_dat$surv_t1, 
                       logsize_t = LOAR_size_dat$logsize_t, 
                       endo = LOAR_endo_dat$endo1, 
                       endo_index = LOAR_endo_dat$endo_index,
                       origin = LOAR_origin_dat$origin,
                       year_t = LOAR_year_dat$year_t_index, 
                       N = 1759L, K = 4L, nyear = 11L, nEndo = 2L)
str(LOAR_data_list) 



# Data manipulation for ELVI to prepare as list for Stan
## create new column with log(size) and recode years, plot and origin 
ELVI_data1 <- ELVI_data %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(year_t_index = as.factor(recode(year_t, '2007'=1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11))) %>% 
  mutate(year_t1_index = as.factor(recode(year_t1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11, '2018'=12))) %>% 
  mutate(plot_index = as.factor(recode(plot, '91'=1, '92'=2, '93'=3, '94'=4, '95'=5, '96'=6, '97'=7, '98'=8, '99'=9, '100'=10, '101'=11)))
dim(ELVI_data1)


ELVI_data1 <- ELVI_data1 %>% 
  filter(!is.na(surv_t1)) 

dim(ELVI_data1)

# View(ELVI_data1)
ELVI_surv_dat <- ELVI_data1 %>% 
  select(surv_t1) %>% 
  filter(!is.na(surv_t1))
ELVI_size_dat <- ELVI_data1 %>% 
  select(logsize_t) %>% 
  filter(!is.na(logsize_t))
ELVI_year_dat <- ELVI_data1 %>% 
  select(year_t_index) %>% 
  filter(!is.na(year_t_index)) %>% 
  mutate(year_t_index = as.integer(year_t_index))
ELVI_endo_dat <- ELVI_data1 %>% 
  select(endo) %>% 
  filter(!is.na(endo)) %>% 
  mutate(endo1 = as.integer(as.integer(endo)-4)) %>% 
  mutate(endo_index = as.integer(endo1+1))   ## recoding endo to 0 or 1
ELVI_origin_dat <- ELVI_data1 %>% 
  select(origin) %>% 
  filter(!is.na(origin)) %>% 
  mutate(origin1 = as.integer(case_when(origin == "O" ~ 0,
                                        origin == "R" ~ 1)))
# Plot data not included as a random effect due to poor estimation
ELVI_plot_dat <- ELVI_data1 %>%
  select(plot, plot_index) %>% 
  filter(!is.na(plot_index)) %>% 
  mutate(plot_index = as.integer(plot_index))

dim(ELVI_surv_dat)
dim(ELVI_size_dat)
dim(ELVI_year_dat)
dim(ELVI_endo_dat)
dim(ELVI_origin_dat)

# Data list for ELVI
ELVI_data_list <- list(surv_t1 = ELVI_surv_dat$surv_t1, 
                       logsize_t = ELVI_size_dat$logsize_t, 
                       endo = ELVI_endo_dat$endo1, 
                       endo_index = ELVI_endo_dat$endo_index,
                       origin = ELVI_origin_dat$origin1,
                       year_t = ELVI_year_dat$year_t_index, 
                       N = 1490L, K = 4L, nyear = 11L, nEndo = 2L)
str(ELVI_data_list) 

# Data manipulation for ELRI to prepare as list for Stan
## create new column with log(size) and recode years and plot 
ELRI_data1 <- ELRI_data %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(year_t_index = as.factor(recode(year_t, '2007'=1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11))) %>% 
  mutate(year_t1_index = as.factor(recode(year_t1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11, '2018'=12))) %>% 
  mutate(plot_index = as.factor(recode(plot, '101'=1, '102'=2, '103'=3, '104'=4, '105'=5, '106'=6, '107'=7, '108'=8, '109'=9, '110'=10)))
dim(ELRI_data1)


ELRI_data1 <- ELRI_data1 %>% 
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(size_t))

dim(ELRI_data1)
str(ELRI_data1)
# View(ELRI_data1)
ELRI_surv_dat <- ELRI_data1 %>% 
  select(surv_t1) %>% 
  filter(!is.na(surv_t1))
ELRI_size_dat <- ELRI_data1 %>% 
  select(logsize_t) %>% 
  filter(!is.na(logsize_t))
ELRI_year_dat <- ELRI_data1 %>% 
  select(year_t_index) %>% 
  filter(!is.na(year_t_index)) %>% 
  mutate(year_t_index = as.integer(year_t_index))
ELRI_endo_dat <- ELRI_data1 %>% 
  select(endo) %>% 
  filter(!is.na(endo)) %>% 
  mutate(endo1 = as.integer(as.integer(endo)-4)) %>% 
  mutate(endo_index = as.integer(endo1+1))   ## recoding endo to 0 or 1
ELRI_origin_dat <- ELRI_data1 %>% 
  select(origin) %>% 
  filter(!is.na(origin)) %>% 
  mutate(origin1 = as.integer(case_when(origin == "O" ~ 0,
                                        origin == "R" ~ 1)))
# Plot data not included as a random effect due to poor estimation
ELRI_plot_dat <- ELRI_data1 %>%
  select(plot, plot_index) %>% 
  filter(!is.na(plot_index)) %>% 
  mutate(plot_index = as.integer(plot_index))

dim(ELRI_surv_dat)
dim(ELRI_size_dat)
dim(ELRI_year_dat)
dim(ELRI_endo_dat)
dim(ELRI_origin_dat)
dim(ELRI_plot_dat)

# Data list for ELRI
ELRI_data_list <- list(surv_t1 = ELRI_surv_dat$surv_t1, 
                       logsize_t = ELRI_size_dat$logsize_t, 
                       endo = ELRI_endo_dat$endo1, 
                       endo_index = ELRI_endo_dat$endo_index,
                       origin = ELRI_origin_dat$origin1,
                       year_t = ELRI_year_dat$year_t_index, 
                       N = 1635L, K = 4L, nyear = 11L, nEndo = 2L)
str(ELRI_data_list) 

# Data manipulation for FESU to prepare as list for Stan
## create new column with log(size) and recode years 
FESU_data1 <- FESU_data %>%
  filter(size_t !=0) %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(year_t_index = as.factor(recode(year_t, '2007'=1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11))) %>% 
  mutate(year_t1_index = as.factor(recode(year_t1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11, '2018'=12)))  
dim(FESU_data1)


FESU_data1 <- FESU_data1 %>% 
  filter(surv_t1 != 2) %>% 
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(size_t))
dim(FESU_data1)

FESU_surv_dat <- FESU_data1 %>% 
  select(surv_t1) %>% 
  filter(!is.na(surv_t1))
FESU_size_dat <- FESU_data1 %>% 
  select(logsize_t) %>% 
  filter(!is.na(logsize_t))
FESU_year_dat <- FESU_data1 %>% 
  select(year_t_index) %>% 
  filter(!is.na(year_t_index)) %>% 
  mutate(year_t_index = as.integer(year_t_index))
FESU_endo_dat <- FESU_data1 %>% 
  select(endo) %>% 
  filter(!is.na(endo)) %>% 
  mutate(endo1 = as.integer(as.integer(endo)-4)) %>% 
  mutate(endo_index = as.integer(endo1+1))   ## recoding endo to 0 or 1
FESU_origin_dat <- FESU_data1 %>% 
  select(origin, plot) %>% 
  filter(!is.na(origin)) %>% 
  mutate(origin1 = as.integer(case_when(origin == "O" ~ 0,
                                        origin == "R" ~ 1,
                                        origin == 122 ~ 1,
                                        origin == 124 ~ 1,
                                        origin == 126 ~ 1,
                                        origin == 127 ~ 1,
                                        origin == 128 ~ 1,
                                        origin == 129 ~ 1,
                                        origin == 130 ~ 1)))

dim(FESU_surv_dat)
dim(FESU_size_dat)
dim(FESU_year_dat)
dim(FESU_endo_dat)
dim(FESU_origin_dat)

# Data list for FESU
FESU_data_list <- list(surv_t1 = FESU_surv_dat$surv_t1, 
                       logsize_t = FESU_size_dat$logsize_t, 
                       endo = FESU_endo_dat$endo1, 
                       endo_index = FESU_endo_dat$endo_index,
                       origin = FESU_origin_dat$origin1,
                       year_t = FESU_year_dat$year_t_index, 
                       N = 4626L, K = 4L, nyear = 11L, nEndo = 2L)
str(FESU_data_list) 

# Data manipulation for AGPE to prepare as list for Stan
## create new column with log(size) and recode years and plot
AGPE_data1 <- AGPE_data %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(year_t_index = as.factor(recode(year_t, '2007'=1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11))) %>% 
  mutate(year_t1_index = as.factor(recode(year_t1, '2008'=2, '2009'=3, '2010'=4, '2011'=5, '2012'=6, '2013'=7, '2014'=8, '2015'=9, '2016'=10, '2017'=11, '2018'=12))) %>% 
  mutate(plot_index = as.factor(recode(plot, '111'=1, '112'=2, '113'=3, '114'=4, '115'=5, '116'=6, '117'=7, '118'=8, '119'=9, '120'=10)))
dim(AGPE_data1)


AGPE_data1 <- AGPE_data1 %>% 
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(size_t))

dim(AGPE_data1)
str(AGPE_data1)
# View(AGPE_data1)
AGPE_surv_dat <- AGPE_data1 %>% 
  select(surv_t1) %>% 
  filter(!is.na(surv_t1))
AGPE_size_dat <- AGPE_data1 %>% 
  select(logsize_t) %>% 
  filter(!is.na(logsize_t))
AGPE_year_dat <- AGPE_data1 %>% 
  select(year_t_index) %>% 
  filter(!is.na(year_t_index)) %>% 
  mutate(year_t_index = as.integer(year_t_index))
AGPE_endo_dat <- AGPE_data1 %>% 
  select(endo) %>% 
  filter(!is.na(endo)) %>% 
  mutate(endo1 = as.integer(case_when(endo == "minus" ~ 0,
                                     endo == "plus" ~ 1,
                                     endo == "0" ~ 0,
                                     endo == "1" ~ 1))) %>% 
  mutate(endo_index = as.integer(endo1+1))   ## recoding endo to 0 or 1 and 1 or 2 for indexing
AGPE_origin_dat <- AGPE_data1 %>% 
  select(origin) %>% 
  filter(!is.na(origin)) %>% 
  mutate(origin1 = as.integer(case_when(origin == "O" ~ 0,
                                        origin == "R" ~ 1)))
# Plot data not included as a random effect due to poor estimation
AGPE_plot_dat <- AGPE_data1 %>%
  select(plot, plot_index) %>% 
  filter(!is.na(plot_index)) %>% 
  mutate(plot_index = as.integer(plot_index))

dim(AGPE_surv_dat)
dim(AGPE_size_dat)
dim(AGPE_year_dat)
dim(AGPE_endo_dat)
dim(AGPE_origin_dat)

# Data list for AGPE
AGPE_data_list <- list(surv_t1 = AGPE_surv_dat$surv_t1, 
                       logsize_t = AGPE_size_dat$logsize_t, 
                       endo = AGPE_endo_dat$endo1, 
                       endo_index = AGPE_endo_dat$endo_index,
                       origin = AGPE_origin_dat$origin1,
                       year_t = AGPE_year_dat$year_t_index, 
                       N = 3108L, K = 4L, nyear = 11L, nEndo = 2L)
str(AGPE_data_list) 

#########################################################################################################
# GLMM for Surv~ size +Endo + Origin  with year random effects------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <- 5000
nb <- 2500
nc <- 4

## here is the Stan model ##
sink("endodemog_surv_full.stan")
cat("
    data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> K;                       // number of predictors
    
    int<lower=0> nyear;                       // number of years (used as index)
    int<lower=0> year_t[N];                      // year of observation
    int<lower=0> nEndo;                       // number of endo treatments
    
    int<lower=0, upper=1> surv_t1[N];      // plant survival at time t+1 and target variable (response)
    vector<lower=-1>[N] logsize_t;                  // log of plant size at time t (predictor)
    int<lower=0, upper=1> endo[N];            // endophyte status 
    int<lower=1, upper=2> endo_index[N];       // index for endophyte effect
    int<lower=0, upper=1> origin[N];            // origin status
    }
    
    parameters {
    real alpha;                  // intercept
    vector[K] beta;              // predictor parameters
    matrix[nEndo, nyear] tau_year;      // random year effect
      
    real<lower=0> sigma_0[nEndo];        //year variance intercept
    }
    

    model {
    vector[N] mu;
   
       // Linear Predictor
    for(n in 1:N){
       mu[n] = alpha + beta[1]*logsize_t[n] + beta[2]*endo[n] + beta[3]*origin[n] 
       + beta[4]*logsize_t[n]*endo[n] +
       + tau_year[endo_index[n], year_t[n]];
    }
    // Priors
    alpha ~ normal(0,1e6);      // prior for fixed intercept
    beta ~ normal(0,1e6);      // prior for predictor intercepts
    for(n in 1:nyear){
    for(e in 1:nEndo){
    tau_year[e,n] ~ normal(0,sigma_0[e]);   // prior for year random effects
    }}
    // Likelihood
      surv_t1 ~ bernoulli_logit(mu);
    }
    
   generated quantities{
    int yrep[N];
    vector[N] mu;
    
    // for posterior predictive check
    for(n in 1:N){
      mu[n] = alpha + beta[1]*logsize_t[n] + beta[2]*endo[n] + beta[3]*origin[n] 
      +beta[4]*logsize_t[n]*endo[n]
      + tau_year[endo_index[n], year_t[n]];
      
      yrep[n] = bernoulli_logit_rng(mu[n]);
    }
    
    }
  
    ", fill = T)
sink()

stanmodel <- stanc("endodemog_surv_full.stan")

## Run the model by calling stan()
## and save the output to .rds files so that they can be called laters

smPOAL <- stan(file = "endodemog_surv_full.stan", data = POAL_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOAL, file = "endodemog_surv_full_POAL.rds")

smPOSY <- stan(file = "endodemog_surv_full.stan", data = POSY_data_list,
           iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smPOSY, file = "endodemog_surv_full_POSY.rds")

smLOAR <- stan(file = "endodemog_surv_full.stan", data = LOAR_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smLOAR, file = "endodemog_surv_full_LOAR.rds")

smELVI <- stan(file = "endodemog_surv_full.stan", data = ELVI_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELVI, file = "endodemog_surv_full_ELVI.rds")

smELRI <- stan(file = "endodemog_surv_full.stan", data = ELRI_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smELRI, file = "endodemog_surv_full_ELRI.rds")

smFESU <- stan(file = "endodemog_surv_full.stan", data = FESU_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smFESU, file = "endodemog_surv_full_FESU.rds")

smAGPE<- stan(file = "endodemog_surv_full.stan", data = AGPE_data_list,
               iter = ni, warmup = nb, chains = nc, save_warmup = FALSE)
saveRDS(smAGPE, file = "endodemog_surv_full_AGPE.rds")

## to read in model output without rerunning models
smPOAL <- readRDS(file = "endodemog_surv_full_POAL.rds")
smPOSY <- readRDS(file = "endodemog_surv_full_POSY.rds")
smLOAR <- readRDS(file = "endodemog_surv_full_LOAR.rds")
smELVI <- readRDS(file = "endodemog_surv_full_ELVI.rds")
smELRI <- readRDS(file = "endodemog_surv_full_ELRI.rds")
smFESU <- readRDS(file = "endodemog_surv_full_FESU.rds")
smAGPE <- readRDS(file = "endodemog_surv_full_AGPE.rds")


#########################################################################################################
###### Perform posterior predictive checks and assess model convergence-------------------------
#########################################################################################################
params <- c("alpha", "beta[1]", "beta[2]", "tau_year[1,1]", "sigma_0[1]", "sigma_0[2]")

##### POAL - survival
print(smPOAL)
# summary(smPOAL)

## plot traceplots of chains for select parameters
traceplot(smPOAL, pars = params)

## Plotting residuals
surv_t1 <- as.vector(POAL_surv_dat$surv_t1)
yrep <- as.matrix(smPOAL, pars = "yrep")
mu <- as.matrix(smPOAL, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(surv_t1 - p)
yrep_resid <- abs(yrep - p[1,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "POAL surv Residual plot")
abline(a = 0, b = 1, col = "blue", lwd = 2)

## plot of neff ratios
ratios_smPOAL <- neff_ratio(smPOAL)
mcmc_neff(ratios_smPOAL, size =3)

## Overlay plot of yrep vs surv_t1 
ppc_dens_overlay(surv_t1, yrep[1:500, ])

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
ppc_dens_overlay(surv_t1, yrep[1:500, ])

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
traceplot(smFESU, pars = c("alpha", "beta[2]", "tau_year[1,1]", "sigma[1]", "sigma[2]"))

## Plotting residuals
surv_t1 <- as.vector(FESU_surv_dat$surv_t1)
yrep <- as.matrix(smFESU, pars = "yrep")
mu <- as.matrix(smFESU, pars = "mu")
p <- invlogit(mu)

y_resid <-  abs(surv_t1 - p)
yrep_resid <- abs(yrep - p[1,])

fit <- as.matrix(rowSums(y_resid))
fit_yrep <- as.matrix(rowSums(yrep_resid))

plot(x = fit, y = fit_yrep, main = "FESU surv Residual plot")
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










