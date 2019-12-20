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


# Visualize endo effects --------------------------------------------------
# survival
betaendo_surv<-rstan::extract(surv_fit, pars = c("betaendo[1]","betaendo[2]",
                                                 "betaendo[3]","betaendo[4]","betaendo[5]",
                                                 "betaendo[6]","betaendo[7]"))
betaendo_surv$betaendo_mean <- (betaendo_surv[[1]] + betaendo_surv[[2]] + betaendo_surv[[3]] + betaendo_surv[[4]]
                                + betaendo_surv[[5]] + betaendo_surv[[6]] + betaendo_surv[[7]])/7

betaendo_surv_mean <- lapply(betaendo_surv,"mean")
betaendo_surv_quant <- as.matrix(data.frame(lapply(betaendo_surv,"quantile",probs=c(0.05,0.25,0.75,0.95))))

sigmaendo_surv<-rstan::extract(surv_fit, pars = c("sigmaendo[1]","sigmaendo[2]",
                                                 "sigmaendo[3]","sigmaendo[4]","sigmaendo[5]",
                                                 "sigmaendo[6]","sigmaendo[7]"))
sigmaendo_surv$sigmaendo_mean <- (sigmaendo_surv[[1]] + sigmaendo_surv[[2]] + sigmaendo_surv[[3]] + sigmaendo_surv[[4]]
                                + sigmaendo_surv[[5]] + sigmaendo_surv[[6]] + sigmaendo_surv[[7]])/7
sigmaendo_surv_mean <- lapply(sigmaendo_surv,"mean")
sigmaendo_surv_quant <- as.matrix(data.frame(lapply(sigmaendo_surv,"quantile",probs=c(0.05,0.25,0.75,0.95))))

## growth
betaendo_grow<-rstan::extract(grow_fit, pars = c("betaendo[1]","betaendo[2]",
                                                 "betaendo[3]","betaendo[4]","betaendo[5]",
                                                 "betaendo[6]","betaendo[7]"))
betaendo_grow$betaendo_mean <- (betaendo_grow[[1]] + betaendo_grow[[2]] + betaendo_grow[[3]] + betaendo_grow[[4]]
                                + betaendo_grow[[5]] + betaendo_grow[[6]] + betaendo_grow[[7]])/7

betaendo_grow_mean <- lapply(betaendo_grow,"mean")
betaendo_grow_quant <- as.matrix(data.frame(lapply(betaendo_grow,"quantile",probs=c(0.05,0.25,0.75,0.95))))

sigmaendo_grow<-rstan::extract(grow_fit, pars = c("sigmaendo[1]","sigmaendo[2]",
                                                  "sigmaendo[3]","sigmaendo[4]","sigmaendo[5]",
                                                  "sigmaendo[6]","sigmaendo[7]"))
sigmaendo_grow$sigmaendo_mean <- (sigmaendo_grow[[1]] + sigmaendo_grow[[2]] + sigmaendo_grow[[3]] + sigmaendo_grow[[4]]
                                  + sigmaendo_grow[[5]] + sigmaendo_grow[[6]] + sigmaendo_grow[[7]])/7
sigmaendo_grow_mean <- lapply(sigmaendo_grow,"mean")
sigmaendo_grow_quant <- as.matrix(data.frame(lapply(sigmaendo_grow,"quantile",probs=c(0.05,0.25,0.75,0.95))))

## flowering
betaendo_flow<-rstan::extract(flow_fit, pars = c("betaendo[1]","betaendo[2]",
                                                 "betaendo[3]","betaendo[4]","betaendo[5]",
                                                 "betaendo[6]","betaendo[7]"))
betaendo_flow$betaendo_mean <- (betaendo_flow[[1]] + betaendo_flow[[2]] + betaendo_flow[[3]] + betaendo_flow[[4]]
                                + betaendo_flow[[5]] + betaendo_flow[[6]] + betaendo_flow[[7]])/7

betaendo_flow_mean <- lapply(betaendo_flow,"mean")
betaendo_flow_quant <- as.matrix(data.frame(lapply(betaendo_flow,"quantile",probs=c(0.05,0.25,0.75,0.95))))

sigmaendo_flow<-rstan::extract(flow_fit, pars = c("sigmaendo[1]","sigmaendo[2]",
                                                  "sigmaendo[3]","sigmaendo[4]","sigmaendo[5]",
                                                  "sigmaendo[6]","sigmaendo[7]"))
sigmaendo_flow$sigmaendo_mean <- (sigmaendo_flow[[1]] + sigmaendo_flow[[2]] + sigmaendo_flow[[3]] + sigmaendo_flow[[4]]
                                  + sigmaendo_flow[[5]] + sigmaendo_flow[[6]] + sigmaendo_flow[[7]])/7
sigmaendo_flow_mean <- lapply(sigmaendo_flow,"mean")
sigmaendo_flow_quant <- as.matrix(data.frame(lapply(sigmaendo_flow,"quantile",probs=c(0.05,0.25,0.75,0.95))))

## make a nice figure
spp_names <- c(data.frame(cbind(unique(LTREB_data_forsurv$species),
                 as.integer(as.numeric(as.factor(unique(LTREB_data_forsurv$species)))))
) %>% 
  arrange(X2) %>% select(X1)); spp_names<- c(as.character(spp_names$X1),"Mean")

spp_cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d","black")
spp_alpha <- 0.75


## start with LOAR and FESU survival as example
win.graph()
par(mfrow=c(2,1),mar=c(5,5,1,1))
plot(rep(0,7),1:7,type="n",xlim=c(-2,2),axes=F,ylab=" ",xlab="Endophyte effect on mean survival",cex.lab=1.4)
axis(side=1)
arrows(0,1,0,8,lty=2,col="gray",length=0)
arrows(betaendo_surv_quant[1,5],5,betaendo_surv_quant[4,5],5,length=0,lwd=2,col=alpha(spp_cols[5],spp_alpha))
arrows(betaendo_surv_quant[2,5],5,betaendo_surv_quant[3,5],5,length=0,lwd=8,col=alpha(spp_cols[5],spp_alpha))
points(betaendo_surv_mean[5],5,cex=3,pch=16,col=alpha(spp_cols[5],spp_alpha))
#axis(side=2,at=5,labels=spp_names[5],las=1,cex.axis=1.5,tick=F)

plot(rep(0,7),1:7,type="n",xlim=c(-2,2),axes=F,ylab=" ",xlab="Endophyte effect on survival variance",cex.lab=1.4)
axis(side=1)
abline(v=0,lty=2,col="gray")#;box()
arrows(sigmaendo_surv_quant[1,5],5,sigmaendo_surv_quant[4,5],5,length=0,lwd=2,col=alpha(spp_cols[5],spp_alpha))
arrows(sigmaendo_surv_quant[2,5],5,sigmaendo_surv_quant[3,5],5,length=0,lwd=8,col=alpha(spp_cols[5],spp_alpha))
points(sigmaendo_surv_mean[5],5,cex=3,pch=16,col=alpha(spp_cols[5],spp_alpha))
#axis(side=2,at=5,labels=spp_names[5],las=1,cex.axis=1.5,tick=F)

## now all of them
win.graph()
par(mfrow=c(2,3),mar=c(5,5,1,1))
## survival
plot(rep(0,8),1:8,type="n",xlim=c(-2,2),axes=F,ylab=" ",xlab="Endophyte effect on mean survival",cex.lab=1.6)
axis(side=2,at=1:8,labels=spp_names[1:8],las=1,cex.axis=1.6)
axis(side=1)
#abline(v=0,lty=2,col="gray")#;box()
arrows(0,1,0,8,lty=2,col="gray",length=0)
arrows(betaendo_surv_quant[1,],1:8,betaendo_surv_quant[4,],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(betaendo_surv_quant[2,],1:8,betaendo_surv_quant[3,],1:8,length=0,lwd=8,col=alpha(spp_cols,spp_alpha))
points(betaendo_surv_mean,1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))

# growth
plot(rep(0,8),1:8,type="n",xlim=c(-1,1),axes=F,ylab=" ",xlab="Endophyte effect on mean growth",cex.lab=1.5)
axis(side=2,at=1:8,labels=spp_names[1:8],las=1,cex.axis=1.6)
axis(side=1)
arrows(0,1,0,8,lty=2,col="gray",length=0)
arrows(betaendo_grow_quant[1,],1:8,betaendo_grow_quant[4,],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(betaendo_grow_quant[2,],1:8,betaendo_grow_quant[3,],1:8,length=0,lwd=8,col=alpha(spp_cols,spp_alpha))
points(betaendo_grow_mean,1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))

# flowering
plot(rep(0,8),1:8,type="n",xlim=c(-3,3),axes=F,ylab=" ",xlab="Endophyte effect on mean flowering",cex.lab=1.5)
axis(side=2,at=1:8,labels=spp_names[1:8],las=1,cex.axis=1.6)
axis(side=1)
arrows(0,1,0,8,lty=2,col="gray",length=0)
arrows(betaendo_flow_quant[1,],1:8,betaendo_flow_quant[4,],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(betaendo_flow_quant[2,],1:8,betaendo_flow_quant[3,],1:8,length=0,lwd=8,col=alpha(spp_cols,spp_alpha))
points(betaendo_flow_mean,1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))

## same for sigmas
plot(rep(0,8),1:8,type="n",xlim=c(-1.5,1.5),axes=F,ylab=" ",xlab="Endophyte effect on survival variance",cex.lab=1.5)
axis(side=2,at=1:8,labels=spp_names[1:8],las=1,cex.axis=1.6)
axis(side=1)
abline(v=0,lty=2,col="gray")#;box()
arrows(sigmaendo_surv_quant[1,],1:8,sigmaendo_surv_quant[4,],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(sigmaendo_surv_quant[2,],1:8,sigmaendo_surv_quant[3,],1:8,length=0,lwd=8,col=alpha(spp_cols,spp_alpha))
points(sigmaendo_surv_mean,1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))

plot(rep(0,8),1:8,type="n",xlim=c(-2,2),axes=F,ylab=" ",xlab="Endophyte effect on growth variance",cex.lab=1.5)
axis(side=2,at=1:8,labels=spp_names[1:8],las=1,cex.axis=1.6)
axis(side=1)
abline(v=0,lty=2,col="gray")#;box()
arrows(sigmaendo_grow_quant[1,],1:8,sigmaendo_grow_quant[4,],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(sigmaendo_grow_quant[2,],1:8,sigmaendo_grow_quant[3,],1:8,length=0,lwd=8,col=alpha(spp_cols,spp_alpha))
points(sigmaendo_grow_mean,1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))

plot(rep(0,8),1:8,type="n",xlim=c(-1,1),axes=F,ylab=" ",xlab="Endophyte effect on flowering variance",cex.lab=1.5)
axis(side=2,at=1:8,labels=spp_names[1:8],las=1,cex.axis=1.6)
axis(side=1)
abline(v=0,lty=2,col="gray")#;box()
arrows(sigmaendo_flow_quant[1,],1:8,sigmaendo_flow_quant[4,],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(sigmaendo_flow_quant[2,],1:8,sigmaendo_flow_quant[3,],1:8,length=0,lwd=8,col=alpha(spp_cols,spp_alpha))
points(sigmaendo_flow_mean,1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))


# fertility -- may not show w the others bc endo effects on sigma not converging
plot(rep(0,7),1:7,type="n",xlim=c(-3,3),axes=F,ylab=" ",xlab="Endophyte effect",cex.lab=1.4)
axis(side=2,at=1:7,labels=spp_names[1:7],las=1)
axis(side=1)
arrows(0,1,0,8,lty=2,col="gray",length=0)
arrows(betaendo_fert_quant[1,],1:8,betaendo_fert_quant[4,],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(betaendo_fert_quant[2,],1:8,betaendo_fert_quant[3,],1:8,length=0,lwd=8,col=alpha(spp_cols,spp_alpha))
points(betaendo_fert_mean,1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))

plot(rep(0,8),1:8,type="n",xlim=c(-10,10),axes=F,ylab=" ",xlab="Endophyte effect")
axis(side=2,at=1:8,labels=spp_names,las=1)
axis(side=1)
abline(v=0,lty=2,col="gray")#;box()
arrows(sigmaendo_fert_quant[1,],1:8,sigmaendo_fert_quant[4,],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(sigmaendo_fert_quant[2,],1:8,sigmaendo_fert_quant[3,],1:8,length=0,lwd=8,col=alpha(spp_cols,spp_alpha))
points(sigmaendo_fert_mean,1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))


# Examples of mean and variance effects -----------------------------------


## cherry-pick nice examples of endo mean and variance effects
## LOAR and FESU survival
mean_surv <- LTREB_data_forsurv %>% 
  mutate(size_bin = as.integer(cut_interval(logsize_t,12))) %>% 
  filter(species=="LOAR" | species=="FESU") %>% 
  group_by(species,size_bin,endo_01) %>% 
  summarise(mean_size = mean(logsize_t),
            mean_surv = mean(surv_t1),
            n_surv = n()) %>% 
  mutate(endo_pch=ifelse(endo_01==0,1,16))

x_seq <- seq(min(LTREB_data_forsurv$logsize_t),max(LTREB_data_forsurv$logsize_t),length.out = 50)
surv_mean<-lapply(rstan::extract(surv_fit, pars = c("beta0[4]","beta0[5]",
                                             "betasize[4]","betasize[5]",
                                             "betaendo[4]","betaendo[5]",
                                             "tau_year[4,1,1]","tau_year[4,1,2]","tau_year[4,1,3]","tau_year[4,1,4]",
                                             "tau_year[4,1,5]","tau_year[4,1,6]","tau_year[4,1,7]","tau_year[4,1,8]",
                                             "tau_year[4,1,9]","tau_year[4,1,10]","tau_year[4,1,11]",
                                             "tau_year[4,2,1]","tau_year[4,2,2]","tau_year[4,2,3]","tau_year[4,2,4]",
                                             "tau_year[4,2,5]","tau_year[4,2,6]","tau_year[4,2,7]","tau_year[4,2,8]",
                                             "tau_year[4,2,9]","tau_year[4,2,10]","tau_year[4,2,11]",
                                             "tau_year[5,1,1]","tau_year[5,1,2]","tau_year[5,1,3]","tau_year[5,1,4]",
                                             "tau_year[5,1,5]","tau_year[5,1,6]","tau_year[5,1,7]","tau_year[5,1,8]",
                                             "tau_year[5,1,9]","tau_year[5,1,10]","tau_year[5,1,11]",
                                             "tau_year[5,2,1]","tau_year[5,2,2]","tau_year[5,2,3]","tau_year[5,2,4]",
                                             "tau_year[5,2,5]","tau_year[5,2,6]","tau_year[5,2,7]","tau_year[5,2,8]",
                                             "tau_year[5,2,9]","tau_year[5,2,10]","tau_year[5,2,11]"
                                             
                                             )),mean)

LOAR_rfx <- matrix(unlist(surv_mean[7:28]),nrow=11,ncol=2,byrow = T)
FESU_rfx <- matrix(unlist(surv_mean[29:50]),nrow=11,ncol=2,byrow = T)

win.graph()
par(mfrow=c(1,1),mar=c(5,5,1,1))
plot(mean_surv$mean_size[mean_surv$species=="LOAR"],mean_surv$mean_surv[mean_surv$species=="LOAR"],
     pch=mean_surv$endo_pch[mean_surv$species=="LOAR"],ylim=c(0,1),col=alpha(spp_cols[5],.5),xlab="log Size",ylab="Pr.(Survival)",cex.lab=1.4,
     lwd=2,cex=2 + 3*(mean_surv$n_surv[mean_surv$species=="LOAR"]/max(mean_surv$n_surv[mean_surv$species=="LOAR"])))
lines(x_seq,invlogit(surv_mean$`beta0[5]`+surv_mean$`betasize[5]`*x_seq),lwd=4,lty=2,col=spp_cols[5])
lines(x_seq,invlogit(surv_mean$`beta0[5]`+surv_mean$`betasize[5]`*x_seq + surv_mean$`betaendo[5]`),lwd=4,lty=1,col=spp_cols[5])
legend("bottomright",legend=c("E+","E-"),lty=1:2,lwd=3,pch=c(16,1),
       col=spp_cols[5],cex=1.8,bty="n")

plot(mean_surv$mean_size[mean_surv$species=="FESU"],mean_surv$mean_surv[mean_surv$species=="FESU"],
     pch=mean_surv$endo_pch[mean_surv$species=="FESU"],ylim=c(0,1),col=alpha(spp_cols[4],.5),xlab="log Size",ylab="Pr.(Survival)",cex.lab=1.4,
     lwd=2,cex=2 + 3*(mean_surv$n_surv[mean_surv$species=="FESU"]/max(mean_surv$n_surv[mean_surv$species=="FESU"])))
lines(x_seq,invlogit(surv_mean$`beta0[4]`+surv_mean$`betasize[4]`*x_seq),lwd=4,col=spp_cols[4],lty=2)
lines(x_seq,invlogit(surv_mean$`beta0[4]`+surv_mean$`betasize[4]`*x_seq + surv_mean$`betaendo[4]`),col=spp_cols[4],lwd=4,lty=1)
legend("bottomright",legend=c("E+","E-"),lty=1:2,lwd=3,pch=c(16,1),
       col=spp_cols[4],cex=1.8,bty="n")

##now year-specific
year_surv <- LTREB_data_forsurv %>% 
  mutate(size_bin = as.integer(cut_interval(logsize_t,12))) %>% 
  filter(species=="LOAR" | species=="FESU") %>% 
  group_by(year_t,species,size_bin,endo_01) %>% 
  summarise(mean_size = mean(logsize_t),
            mean_surv = mean(surv_t1),
            n_surv = n()) %>% 
  mutate(endo_pch=ifelse(endo_01==0,1,16))
year_surv_years <- unique(year_surv$year_t)
  
plot(mean_surv$mean_size[mean_surv$species=="FESU"],mean_surv$mean_surv[mean_surv$species=="FESU"],
     pch=mean_surv$endo_pch[mean_surv$species=="FESU"],ylim=c(0,1),col=alpha(spp_cols[4],.5),xlab="log Size",ylab="Pr.(Survival)",cex.lab=1.4,
     lwd=2,cex=2 + 3*(mean_surv$n_surv[mean_surv$species=="FESU"]/max(mean_surv$n_surv[mean_surv$species=="FESU"])))
for(t in 1:length(year_surv_years)){
  lines(x_seq,invlogit(surv_mean$`beta0[4]`+surv_mean$`betasize[4]`*x_seq + 
                         FESU_rfx[t,1]),lwd=4,col=spp_cols[4],lty=2)
  lines(x_seq,invlogit(surv_mean$`beta0[4]`+surv_mean$`betasize[4]`*x_seq + 
                         FESU_rfx[t,2]),lwd=4,col=spp_cols[4],lty=1)
}

win.graph()
par(mfrow=c(2,1),mar=c(5,5,1,1))
plot(mean_surv$mean_size[mean_surv$species=="LOAR" & mean_surv$endo_01==0],mean_surv$mean_surv[mean_surv$species=="LOAR" & mean_surv$endo_01==0],
     pch=mean_surv$endo_pch[mean_surv$species=="LOAR" & mean_surv$endo_01==0],ylim=c(0,1),col=alpha(spp_cols[5],.5),xlab="log Size",ylab="Pr.(Survival)",cex.lab=1.4,
     lwd=2,cex=2 + 3*(mean_surv$n_surv[mean_surv$species=="LOAR" & mean_surv$endo_01==0]/max(mean_surv$n_surv[mean_surv$species=="LOAR" & mean_surv$endo_01==0])))
for(t in 1:length(year_surv_years)){
  lines(x_seq,invlogit(surv_mean$`beta0[5]`+surv_mean$`betasize[4]`*x_seq + 
                         LOAR_rfx[t,1]),lwd=4,col=spp_cols[5],lty=2)
}

plot(mean_surv$mean_size[mean_surv$species=="LOAR" & mean_surv$endo_01==1],mean_surv$mean_surv[mean_surv$species=="LOAR" & mean_surv$endo_01==1],
     pch=mean_surv$endo_pch[mean_surv$species=="LOAR" & mean_surv$endo_01==1],ylim=c(0,1),col=alpha(spp_cols[5],.5),xlab="log Size",ylab="Pr.(Survival)",cex.lab=1.4,
     lwd=2,cex=2 + 3*(mean_surv$n_surv[mean_surv$species=="LOAR" & mean_surv$endo_01==1]/max(mean_surv$n_surv[mean_surv$species=="LOAR" & mean_surv$endo_01==1])))
for(t in 1:length(year_surv_years)){
  lines(x_seq,invlogit(surv_mean$`beta0[5]`+surv_mean$`betasize[4]`*x_seq + 
                         LOAR_rfx[t,2]),lwd=4,col=spp_cols[5],lty=1)
}
# climate data ------------------------------------------------------------
climate <- read_csv(file = "C:/Users/tm9/Dropbox/EndodemogData/PRISMClimateData_BrownCo.csv") %>% 
  mutate(year = year(Date), month = month(Date), day = day(Date)) %>% 
  rename(ppt = `ppt (mm)`, tmean = `tmean (degrees C)`) %>% 
  mutate(site_lat = 39.235900000000, site_long = -86.218100000000)

LOAR_climate <- climate %>% 
  mutate(census_month = 7, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>% 
  filter(climate_year != 2006) %>% 
  group_by(climate_year) %>% 
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean)) %>% 
  filter(climate_year <= max(LTREB_full$year_t))

FESU_lm_Em <- lm(FESU_rfx[,1] ~ LOAR_climate$`Cumulative PPT (mm)`)
FESU_lm_Ep <- lm(FESU_rfx[,2] ~ LOAR_climate$`Cumulative PPT (mm)`)

win.graph()
plot(LOAR_climate$`Cumulative PPT (mm)`,FESU_rfx[,2],col=alpha(spp_cols[4],.5),cex.lab=1.4,
     cex=3,lwd=2,xlab="Annual precipitation (mm)",ylab="Survival (standardized)")
points(LOAR_climate$`Cumulative PPT (mm)`,FESU_rfx[,1],pch=16,col=alpha(spp_cols[4],.5),cex=3,lwd=2)


abline(coef(FESU_lm_Em),lty=2)
abline(coef(FESU_lm_Ep))

LOAR_lm_Em <- lm(LOAR_rfx[,1] ~ LOAR_climate$`Cumulative PPT (mm)`)
LOAR_lm_Ep <- lm(LOAR_rfx[,2] ~ LOAR_climate$`Cumulative PPT (mm)`)

plot(LOAR_climate$`Cumulative PPT (mm)`,LOAR_rfx[,1])
points(LOAR_climate$`Cumulative PPT (mm)`,LOAR_rfx[,2],pch=16)
abline(coef(LOAR_lm_Em),lty=2)
abline(coef(LOAR_lm_Ep))

plot(LOAR_climate$`Mean Temp. (C°)`,LOAR_rfx[,1])
points(LOAR_climate$`Mean Temp. (C°)`,LOAR_rfx[,2],pch=16)

## check out AGPE growth
AGPE_grow_rfx<-lapply(rstan::extract(grow_fit, pars = c("tau_year[1,1,1]","tau_year[1,1,2]","tau_year[1,1,3]","tau_year[1,1,4]",
                                                    "tau_year[1,1,5]","tau_year[1,1,6]","tau_year[1,1,7]","tau_year[1,1,8]",
                                                    "tau_year[1,1,9]","tau_year[1,1,10]","tau_year[1,1,11]",
                                                    "tau_year[1,2,1]","tau_year[1,2,2]","tau_year[1,2,3]","tau_year[1,2,4]",
                                                    "tau_year[1,2,5]","tau_year[1,2,6]","tau_year[1,2,7]","tau_year[1,2,8]",
                                                    "tau_year[1,2,9]","tau_year[1,2,10]","tau_year[1,2,11]"
                                                    
)),mean)
AGPE_grow_rfx <- matrix(unlist(AGPE_grow_rfx[1:22]),nrow=11,ncol=2,byrow = T)
plot(LOAR_climate$`Cumulative PPT (mm)`,AGPE_grow_rfx[,1])
points(LOAR_climate$`Cumulative PPT (mm)`,AGPE_grow_rfx[,2],pch=16)