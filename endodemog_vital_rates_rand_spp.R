library(tidyverse)
library(rstan)
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )
library(bayesplot)
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

## this is good, there are no plot numbers duplicated across species
table(LTREB_data_forsurv$plot_fixed,LTREB_data_forsurv$species)

## note that this is how species are indexing:
cbind(unique(LTREB_data_forsurv$species),
as.integer(as.numeric(as.factor(unique(LTREB_data_forsurv$species)))))

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
                 origin_01 = LTREB_data_forsurv$origin_01)
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
                 origin_01 = LTREB_data_forflw$origin_01)
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
                 origin_01 = LTREB_data_forgrow$origin_01)
fert_dat <- list(nYear = length(unique(LTREB_data_forfert$year_t - (min(LTREB_data_forfert$year_t)-1))),
                 nPlot = max(LTREB_data_forfert$plot_fixed),
                 nSpp = length(unique(LTREB_data_forfert$species)),
                 nEndo=length(unique(LTREB_data_forfert$endo_01)),
                 N = nrow(LTREB_data_forfert),
                 year_t = LTREB_data_forfert$year_t - (min(LTREB_data_forfert$year_t)-1),
                 plot = LTREB_data_forfert$plot_fixed,
                 spp = as.integer(as.numeric(as.factor(LTREB_data_forfert$species))),
                 y = LTREB_data_forfert$FLW_COUNT_T1,
                 logsize_t = LTREB_data_forfert$logsize_t1,
                 endo_01 = LTREB_data_forfert$endo_01,
                 origin_01 = LTREB_data_forfert$origin_01)

sim_pars <- list(
  warmup = 1000, 
  iter = 5000, 
  thin = 3, 
  chains = 1
)

surv_fit <- stan(
  file = 'survival_flowering_random_spp.stan',
  data = surv_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
#write_rds(surv_fit,paste0(tompath,"Fulldataplusmetadata/SppRFX/surv_fit.rds"))

flow_fit <- stan(
  file = 'survival_flowering_random_spp.stan',
  data = flow_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(flow_fit,paste0(tompath,"Fulldataplusmetadata/SppRFX/flow_fit.rds"))

grow_fit <- stan(
  file = 'growth_fertility_random_spp.stan',
  data = grow_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(grow_fit,paste0(tompath,"Fulldataplusmetadata/SppRFX/grow_fit.rds"))

fert_fit <- stan(
  file = 'growth_fertility_random_spp.stan',
  data = fert_dat,
  warmup = sim_pars$warmup,
  iter = sim_pars$iter,
  thin = sim_pars$thin,
  chains = sim_pars$chains )
write_rds(fert_fit,paste0(tompath,"Fulldataplusmetadata/SppRFX/fert_fit.rds"))


# Diagnostics and results -------------------------------------------------
surv_fit <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/surv_fit.rds"))
mcmc_dens_overlay(surv_fit,pars = c("mu_betaendo","mu_sigmaendo"))
mcmc_dens_overlay(surv_fit,pars = c("mu_betaendo","mu_sigmaendo"))
mcmc_rhat(rhat(surv_fit,pars = c("mu_betaendo","mu_sigmaendo")))
mcmc_dens_overlay(surv_fit,pars = c("betaendo"))

predS <- rstan::extract(surv_fit, pars = c("p"))$p
n_post_draws <- 100
post_draws <- sample.int(dim(predS)[1], n_post_draws)
y_s_sim <- matrix(NA,n_post_draws,length(surv_dat$y))
for(i in 1:n_post_draws){
  y_s_sim[i,] <- rbinom(n=length(surv_dat$y), size=1, prob = invlogit(predS[post_draws[i],]))
}
ppc_dens_overlay(surv_dat$y, y_s_sim)

mcmc_trace(surv_fit,par=c("beta0"))
mcmc_trace(surv_fit,par=c("betasize"))
mcmc_trace(surv_fit,par=c("betaendo"))

betaendo_surv<-rstan::extract(surv_fit, pars = c("mu_betaendo","betaendo[1]","betaendo[2]",
                                                 "betaendo[3]","betaendo[4]","betaendo[5]",
                                                 "betaendo[6]","betaendo[7]"))
betaendo_surv_mean <- lapply(betaendo_surv,"mean")
betaendo_surv_quant <- as.matrix(data.frame(lapply(betaendo_surv,"quantile",probs=c(0.05,0.25,0.75,0.95))))

sigmaendo_surv<-rstan::extract(surv_fit, pars = c("mu_sigmaendo","sigmaendo[1]","sigmaendo[2]",
                                                 "sigmaendo[3]","sigmaendo[4]","sigmaendo[5]",
                                                 "sigmaendo[6]","sigmaendo[7]"))

sigmaendo_surv<-rstan::extract(surv_fit, pars = c("sigma0[1]","sigma0[2]",
                                                  "sigma0[3]","sigma0[4]","sigma0[5]",
                                                  "sigma0[6]","sigma0[7]"))
sigmaendo_surv_mean <- lapply(sigmaendo_surv,"mean")
sigmaendo_surv_quant <- as.matrix(data.frame(lapply(sigmaendo_surv,"quantile",probs=c(0.05,0.25,0.75,0.95))))

par(mfrow=c(2,1))
plot(rep(0,8),1:8,type="n")
abline(v=0,lty=2,col="gray")
points(betaendo_surv_mean,1:8,cex=c(4,rep(2,7)),pch=16)
arrows(betaendo_surv_quant[1,],1:8,betaendo_surv_quant[4,],1:8,length=0,lwd=2)
arrows(betaendo_surv_quant[2,],1:8,betaendo_surv_quant[3,],1:8,length=0,lwd=4)

plot(rep(0,8),1:8,type="n")
abline(v=0,lty=2,col="gray")
points(sigmaendo_surv_mean,1:8,cex=c(4,rep(2,7)),pch=16)
arrows(sigmaendo_surv_quant[1,],1:8,sigmaendo_surv_quant[4,],1:8,length=0,lwd=2)
arrows(sigmaendo_surv_quant[2,],1:8,sigmaendo_surv_quant[3,],1:8,length=0,lwd=4)
