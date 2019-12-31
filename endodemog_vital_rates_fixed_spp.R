library(tidyverse)
library(rstan)
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )
library(bayesplot)
library(scales)
invlogit<-function(x){exp(x)/(1+exp(x))}


joshpath <- "/Users/joshuacfowler/Dropbox/EndodemogData/"
tompath <- "C:/Users/tm9/Dropbox/EndodemogData/"
#tompath <- "C:/Users/tm634/Dropbox/EndodemogData/"

LTREB_full <- read_csv(paste0(tompath,"Fulldataplusmetadata/LTREB_full.csv"))
LTREB_data_forsurv <- LTREB_full %>%
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(logsize_t >= 0) %>% 
  filter(!is.na(endo_01))
## note that I am fitting flowering model to size_t1 / flowering_t1 (Josh did size_t/flowering_t). We gain 2018 data my way.
LTREB_data_forflw <- LTREB_full %>% 
  filter(!is.na(FLW_STAT_T1)) %>% 
  filter(!is.na(logsize_t1)) %>% 
  filter(!is.na(endo_01))
LTREB_data_forgrow <- LTREB_full %>%
  filter(!is.na(logsize_t)) %>% 
  filter(!is.na(size_t1)) %>% 
  filter(!is.na(endo_01)) %>% 
  mutate(size_t1 = as.integer(size_t1))
## again, here I am using size_t1 and flw_count_t1
LTREB_data_forfert <- LTREB_full %>% 
  filter(!is.na(FLW_COUNT_T1)) %>% 
  filter(FLW_COUNT_T1 > 0) %>% 
  filter(!is.na(logsize_t1))
LTREB_data_forspike <- LTREB_full %>%
  dplyr::select(-FLW_COUNT_T1, -FLW_STAT_T1, -SPIKE_A_T1, -SPIKE_B_T1, -SPIKE_C_T1, -SPIKE_D_T1, -endo_status_from_check, -plot_endo_for_check, -endo_mismatch, -dist_a, -dist_b) %>% 
  filter(!is.na(FLW_STAT_T)) %>% 
  filter(FLW_STAT_T>0) %>% 
  melt(id.var = c("plot_fixed" ,            "pos"         ,           "id",
                  "species"       ,         "species_index"  ,        "endo_01",
                  "endo_index"  ,           "origin_01"       ,       "birth" ,
                  "year_t1"         ,       "year_t1_index"       ,   "surv_t1" ,
                  "size_t1"         ,       "logsize_t1"       ,
                  "year_t",
                  "year_t_index"     ,      "size_t"           ,      "logsize_t"  ,
                  "FLW_COUNT_T"      ,      "FLW_STAT_T"),
       value.name = "spike_count_t") %>% 
  rename(spikelet_id = variable) %>% 
  filter(!is.na(spike_count_t), spike_count_t > 0) %>% 
  mutate(spike_count_t = as.integer(spike_count_t))

ggplot(LTREB_data_forspike)+
  geom_histogram(aes(x=spike_count_t))+
  facet_grid(year_t~species)
## I don't think there are enough data to fit year variances
## so I am just going to fit fixed effects of size and endo
rm(LTREB_full)

## this is good, there are no plot numbers duplicated across species
table(LTREB_data_forsurv$plot_fixed,LTREB_data_forsurv$species)

## note that this is how species are indexing:
cbind(unique(LTREB_data_forsurv$species),
as.integer(as.numeric(as.factor(unique(LTREB_data_forsurv$species)))))

table(LTREB_data_forfert$FLW_COUNT_T1,LTREB_data_forfert$species,LTREB_data_forfert$year_t1)

surv_dat <- list(nYear = length(unique(LTREB_data_forsurv$year_t - (min(LTREB_data_forsurv$year_t)-1))),
                 nPlot = max(LTREB_data_forsurv$plot_fixed),
                 nSpp = length(unique(LTREB_data_forsurv$species)),
                 nEndo=length(unique(LTREB_data_forsurv$endo_01)),
                 N = nrow(LTREB_data_forsurv),
                 year_t = LTREB_data_forsurv$year_t - (min(LTREB_data_forsurv$year_t)-1),
                 plot = LTREB_data_forsurv$plot_fixed,
                 spp = as.integer(as.numeric(as.factor(LTREB_data_forsurv$species))),
                 y = LTREB_data_forsurv$surv_t1,
                 logsize_t = LTREB_data_forsurv$logsize_t,
                 endo_01 = LTREB_data_forsurv$endo_01,
                 origin_01 = LTREB_data_forsurv$origin_01);rm(LTREB_data_forsurv)
flow_dat <- list(nYear = length(unique(LTREB_data_forflw$year_t - (min(LTREB_data_forflw$year_t)-1))),
                 nPlot = max(LTREB_data_forflw$plot_fixed),
                 nSpp = length(unique(LTREB_data_forflw$species)),
                 nEndo=length(unique(LTREB_data_forflw$endo_01)),
                 N = nrow(LTREB_data_forflw),
                 year_t = LTREB_data_forflw$year_t - (min(LTREB_data_forflw$year_t)-1),
                 plot = LTREB_data_forflw$plot_fixed,
                 spp = as.integer(as.numeric(as.factor(LTREB_data_forflw$species))),
                 y = LTREB_data_forflw$FLW_STAT_T1,
                 logsize_t = LTREB_data_forflw$logsize_t1,
                 endo_01 = LTREB_data_forflw$endo_01,
                 origin_01 = LTREB_data_forflw$origin_01);rm(LTREB_data_forflw)
grow_dat <- list(nYear = length(unique(LTREB_data_forgrow$year_t - (min(LTREB_data_forgrow$year_t)-1))),
                 nPlot = max(LTREB_data_forgrow$plot_fixed),
                 nSpp = length(unique(LTREB_data_forgrow$species)),
                 nEndo=length(unique(LTREB_data_forgrow$endo_01)),
                 N = nrow(LTREB_data_forgrow),
                 year_t = LTREB_data_forgrow$year_t - (min(LTREB_data_forgrow$year_t)-1),
                 plot = LTREB_data_forgrow$plot_fixed,
                 spp = as.integer(as.numeric(as.factor(LTREB_data_forgrow$species))),
                 y = LTREB_data_forgrow$size_t1,
                 logsize_t = LTREB_data_forgrow$logsize_t,
                 endo_01 = LTREB_data_forgrow$endo_01,
                 origin_01 = LTREB_data_forgrow$origin_01);rm(LTREB_data_forgrow)
fert_dat <- list(nYear = length(unique(LTREB_data_forfert$year_t - (min(LTREB_data_forfert$year_t)-1))),
                 nPlot = max(LTREB_data_forfert$plot_fixed),
                 nSpp = length(unique(LTREB_data_forfert$species)),
                 #nEndo=length(unique(LTREB_data_forfert$endo_01)),
                 N = nrow(LTREB_data_forfert),
                 year_t = LTREB_data_forfert$year_t - (min(LTREB_data_forfert$year_t)-1),
                 plot = LTREB_data_forfert$plot_fixed,
                 spp = as.integer(as.numeric(as.factor(LTREB_data_forfert$species))),
                 y = LTREB_data_forfert$FLW_COUNT_T1,
                 logsize_t = LTREB_data_forfert$logsize_t1,
                 #endo_01 = LTREB_data_forfert$endo_01,
                 origin_01 = LTREB_data_forfert$origin_01);rm(LTREB_data_forfert)

spike_dat <- list(nYear = length(unique(LTREB_data_forspike$year_t - (min(LTREB_data_forspike$year_t)-1))),
                  nPlot = max(LTREB_data_forspike$plot_fixed),
                  nSpp = length(unique(LTREB_data_forspike$species)),
                  nEndo=length(unique(LTREB_data_forspike$endo_01)),
                  N = nrow(LTREB_data_forspike),
                  year_t = LTREB_data_forspike$year_t - (min(LTREB_data_forspike$year_t)-1),
                  plot = LTREB_data_forspike$plot_fixed,
                  spp = as.integer(as.numeric(as.factor(LTREB_data_forspike$species))),
                  y = LTREB_data_forspike$spike_count_t,
                  logsize_t = LTREB_data_forspike$logsize_t,
                  endo_01 = LTREB_data_forspike$endo_01,
                  origin_01 = LTREB_data_forspike$origin_01);rm(LTREB_data_forspike)

sim_pars <- list(
  warmup = 1000, 
  iter = 5000, 
  thin = 3, 
  chains = 1
)

flow_fit <- stan(
  file = 'survival_flowering_fixed_spp.stan',
  data = flow_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(flow_fit,paste0(tompath,"Fulldataplusmetadata/SppRFX/flow_fit_fixed.rds"))

surv_fit <- stan(
  file = 'survival_flowering_fixed_spp.stan',
  data = surv_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(surv_fit,paste0(tompath,"Fulldataplusmetadata/SppRFX/surv_fit_fixed.rds"))

grow_fit <- stan(
  file = 'growth_fertility_fixed_spp.stan',
  data = grow_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(grow_fit,paste0(tompath,"Fulldataplusmetadata/SppRFX/grow_fit_fixed.rds"))

fert_fit <- stan(
  file = 'growth_fertility_fixed_spp.stan',
  data = fert_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(fert_fit,paste0(tompath,"Fulldataplusmetadata/SppRFX/fert_fit_fixed.rds"))

## the stan file here is fitting endo effects on mean only
fert_fit_endo_mean <- stan(
  file = 'fertility_fixed_spp_endo_mean.stan', 
  data = fert_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(fert_fit_endo_mean,paste0(tompath,"Fulldataplusmetadata/SppRFX/fert_fit_fixed_endo_mean.rds"))

## the stan file here is fitting no endo effects
fert_fit_no_endo <- stan(
  file = 'fertility_fixed_spp_no_endo.stan', 
  data = fert_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(fert_fit_no_endo,paste0(tompath,"Fulldataplusmetadata/SppRFX/fert_fit_no_endo.rds"))

## wow, fuck this shit
library(lme4)
fert_glmm <- glmer(fert_dat$y ~ as.factor(fert_dat$spp)-1 +fert_dat$logsize_t + fert_dat$origin_01 + 
                     (1|fert_dat$plot) + (1|fert_dat$year_t),
                   family="poisson")
summary(fert_glmm)
write_rds(fert_glmm,paste0(tompath,"Fulldataplusmetadata/SppRFX/fert_glmm.rds"))

spike_fit_endo_mean <- stan(
  file = 'spikelets_fixed_spp_endo_mean.stan', 
  data = spike_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(spike_fit_endo_mean,paste0(tompath,"Fulldataplusmetadata/SppRFX/spike_fit_fixed_endo_mean.rds"))

spike_fit_no_endo <- stan(
  file = 'spikelets_fixed_spp_no_endo.stan', 
  data = spike_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(spike_fit_no_endo,paste0(tompath,"Fulldataplusmetadata/SppRFX/spike_fit_no_endo.rds"))


# Diagnostics and results -------------------------------------------------
#survival
surv_fit <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/surv_fit_fixed.rds"))
predS <- rstan::extract(surv_fit, pars = c("p"))$p
n_post_draws <- 100
post_draws <- sample.int(dim(predS)[1], n_post_draws)
y_s_sim <- matrix(NA,n_post_draws,length(surv_dat$y))
for(i in 1:n_post_draws){
  y_s_sim[i,] <- rbinom(n=length(surv_dat$y), size=1, prob = invlogit(predS[post_draws[i],]))
}
ppc_dens_overlay(surv_dat$y, y_s_sim)

#growth
grow_fit <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/grow_fit_fixed.rds"))
predG <- rstan::extract(grow_fit, pars = c("lambda"))$lambda
n_post_draws <- 100
post_draws <- sample.int(dim(predG)[1], n_post_draws)
y_g_sim <- matrix(NA,n_post_draws,length(grow_dat$y))

#flowering
flow_fit <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/flow_fit_fixed.rds"))
predF <- rstan::extract(flow_fit, pars = c("p"))$p

## fertility
fert_fit <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/fert_fit_fixed.rds"))

mcmc_trace(fert_fit,pars=c("betaendo[1]","betaendo[2]","betaendo[3]"
                           ,"betaendo[4]","betaendo[5]","betaendo[6]"
                           ,"betaendo[7]"))

mcmc_trace(fert_fit,pars=c("betasize[1]","betasize[2]","betasize[3]"
                           ,"betasize[4]","betasize[5]","betasize[6]"
                           ,"betasize[7]"))

mcmc_trace(fert_fit,pars=c("sigmaendo[1]","sigmaendo[2]","sigmaendo[3]"
                           ,"sigmaendo[4]","sigmaendo[5]","sigmaendo[6]"
                           ,"sigmaendo[7]"))
mcmc_trace(fert_fit,pars=c("sigma0[1]","sigma0[2]","sigma0[3]"
                           ,"sigma0[4]","sigma0[5]","sigma0[6]"
                           ,"sigma0[7]"))

fert_fit_endo_mean <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/fert_fit_fixed_endo_mean.rds"))
mcmc_trace(fert_fit_endo_mean,pars=c("betaendo[1]","betaendo[2]","betaendo[3]"
                           ,"betaendo[4]","betaendo[5]","betaendo[6]"
                           ,"betaendo[7]"))

## spikelets
spike_fit_no_endo <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/spike_fit_no_endo.rds"))
mcmc_trace(spike_fit_no_endo,pars=c("betasize[1]","betasize[2]","betasize[3]"
                           ,"betasize[4]","betasize[5]","betasize[6]"
                           ,"betasize[7]"))
predSpike <- rstan::extract(spike_fit_no_endo, pars = c("lambda"))$lambda
n_post_draws <- 100
post_draws <- sample.int(dim(predSpike)[1], n_post_draws)
y_spike_sim <- matrix(NA,n_post_draws,length(spike_dat$y))
for(i in 1:n_post_draws){
  y_spike_sim[i,] <- rpois(n=length(spike_dat$y), lambda = exp(predSpike[post_draws[i],]))
}
ppc_dens_overlay(spike_dat$y, y_spike_sim)


