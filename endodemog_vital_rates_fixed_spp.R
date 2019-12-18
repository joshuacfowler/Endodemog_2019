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
                 nEndo=length(unique(LTREB_data_forfert$endo_01)),
                 N = nrow(LTREB_data_forfert),
                 year_t = LTREB_data_forfert$year_t - (min(LTREB_data_forfert$year_t)-1),
                 plot = LTREB_data_forfert$plot_fixed,
                 spp = as.integer(as.numeric(as.factor(LTREB_data_forfert$species))),
                 y = LTREB_data_forfert$FLW_COUNT_T1,
                 logsize_t = LTREB_data_forfert$logsize_t1,
                 endo_01 = LTREB_data_forfert$endo_01,
                 origin_01 = LTREB_data_forfert$origin_01);rm(LTREB_data_forfert)

sim_pars <- list(
  warmup = 2000, 
  iter = 15000, 
  thin = 3, 
  chains = 3
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

## fertility
betaendo_fert<-rstan::extract(fert_fit, pars = c("betaendo[1]","betaendo[2]",
                                                 "betaendo[3]","betaendo[4]","betaendo[5]",
                                                 "betaendo[6]","betaendo[7]"))
betaendo_fert$betaendo_mean <- (betaendo_fert[[1]] + betaendo_fert[[2]] + betaendo_fert[[3]] + betaendo_fert[[4]]
                                + betaendo_fert[[5]] + betaendo_fert[[6]] + betaendo_fert[[7]])/7

betaendo_fert_mean <- lapply(betaendo_fert,"mean")
betaendo_fert_quant <- as.matrix(data.frame(lapply(betaendo_fert,"quantile",probs=c(0.05,0.25,0.75,0.95))))

sigmaendo_fert<-rstan::extract(fert_fit, pars = c("sigmaendo[1]","sigmaendo[2]",
                                                  "sigmaendo[3]","sigmaendo[4]","sigmaendo[5]",
                                                  "sigmaendo[6]","sigmaendo[7]"))
sigmaendo_fert$sigmaendo_mean <- (sigmaendo_fert[[1]] + sigmaendo_fert[[2]] + sigmaendo_fert[[3]] + sigmaendo_fert[[4]]
                                  + sigmaendo_fert[[5]] + sigmaendo_fert[[6]] + sigmaendo_fert[[7]])/7
sigmaendo_fert_mean <- lapply(sigmaendo_fert,"mean")
sigmaendo_fert_quant <- as.matrix(data.frame(lapply(sigmaendo_fert,"quantile",probs=c(0.05,0.25,0.75,0.95))))

## make a nice figure
spp_names <- c(data.frame(cbind(unique(LTREB_data_forsurv$species),
                 as.integer(as.numeric(as.factor(unique(LTREB_data_forsurv$species)))))
) %>% 
  arrange(X2) %>% select(X1)); spp_names<- c(as.character(spp_names$X1),"Mean")

spp_cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d","black")
spp_alpha <- 0.75


par(mfrow=c(2,3),xpd=T)
## survival
plot(rep(0,7),1:7,type="n",xlim=c(-2,2),axes=F,ylab=" ",xlab="Endophyte effect",cex.lab=1.4)
axis(side=2,at=1:7,labels=spp_names[1:7],las=1)
axis(side=1)
#abline(v=0,lty=2,col="gray")#;box()
arrows(0,1,0,8,lty=2,col="gray",length=0)
arrows(betaendo_surv_quant[1,],1:8,betaendo_surv_quant[4,],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(betaendo_surv_quant[2,],1:8,betaendo_surv_quant[3,],1:8,length=0,lwd=8,col=alpha(spp_cols,spp_alpha))
points(betaendo_surv_mean,1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))

# growth
plot(rep(0,7),1:7,type="n",xlim=c(-1,1),axes=F,ylab=" ",xlab="Endophyte effect",cex.lab=1.4)
axis(side=2,at=1:7,labels=spp_names[1:7],las=1)
axis(side=1)
arrows(0,1,0,8,lty=2,col="gray",length=0)
arrows(betaendo_grow_quant[1,],1:8,betaendo_grow_quant[4,],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(betaendo_grow_quant[2,],1:8,betaendo_grow_quant[3,],1:8,length=0,lwd=8,col=alpha(spp_cols,spp_alpha))
points(betaendo_grow_mean,1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))

# flowering
plot(rep(0,7),1:7,type="n",xlim=c(-3,3),axes=F,ylab=" ",xlab="Endophyte effect",cex.lab=1.4)
axis(side=2,at=1:7,labels=spp_names[1:7],las=1)
axis(side=1)
arrows(0,1,0,8,lty=2,col="gray",length=0)
arrows(betaendo_flow_quant[1,],1:8,betaendo_flow_quant[4,],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(betaendo_flow_quant[2,],1:8,betaendo_flow_quant[3,],1:8,length=0,lwd=8,col=alpha(spp_cols,spp_alpha))
points(betaendo_flow_mean,1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))

## same for sigmas
plot(rep(0,8),1:8,type="n",xlim=c(-1.5,1.5),axes=F,ylab=" ",xlab="Endophyte effect")
axis(side=2,at=1:8,labels=spp_names,las=1)
axis(side=1)
abline(v=0,lty=2,col="gray")#;box()
arrows(sigmaendo_surv_quant[1,],1:8,sigmaendo_surv_quant[4,],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(sigmaendo_surv_quant[2,],1:8,sigmaendo_surv_quant[3,],1:8,length=0,lwd=8,col=alpha(spp_cols,spp_alpha))
points(sigmaendo_surv_mean,1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))

plot(rep(0,8),1:8,type="n",xlim=c(-2,2),axes=F,ylab=" ",xlab="Endophyte effect")
axis(side=2,at=1:8,labels=spp_names,las=1)
axis(side=1)
abline(v=0,lty=2,col="gray")#;box()
arrows(sigmaendo_grow_quant[1,],1:8,sigmaendo_grow_quant[4,],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(sigmaendo_grow_quant[2,],1:8,sigmaendo_grow_quant[3,],1:8,length=0,lwd=8,col=alpha(spp_cols,spp_alpha))
points(sigmaendo_grow_mean,1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))

plot(rep(0,8),1:8,type="n",xlim=c(-1,1),axes=F,ylab=" ",xlab="Endophyte effect")
axis(side=2,at=1:8,labels=spp_names,las=1)
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


## cherry-pick nice examples of endo mean and variance effects
## LOAR and FESU survival
mean_surv <- LTREB_data_forsurv %>% 
  mutate(size_bin = as.integer(cut_interval(logsize_t,9))) %>% 
  filter(species=="LOAR" | species=="FESU") %>% 
  group_by(species,size_bin,endo_01) %>% 
  summarise(mean_size = mean(logsize_t),
            mean_surv = mean(surv_t1),
            n_surv = n()) %>% 
  mutate(endo_pch=ifelse(endo_01==0,1,16))

plot(mean_surv$mean_size[mean_surv$species=="LOAR"],mean_surv$mean_surv[mean_surv$species=="LOAR"],
     pch=mean_surv$endo_pch[mean_surv$species=="LOAR"],ylim=c(0,1),
     cex=2 + 3*(mean_surv$n_surv[mean_surv$species=="LOAR"]/max(mean_surv$n_surv[mean_surv$species=="LOAR"])))

plot(mean_surv$mean_size[mean_surv$species=="FESU"],mean_surv$mean_surv[mean_surv$species=="FESU"],
     pch=mean_surv$endo_pch[mean_surv$species=="FESU"],ylim=c(0,1),
     cex=2 + 3*(mean_surv$n_surv[mean_surv$species=="FESU"]/max(mean_surv$n_surv[mean_surv$species=="FESU"])))
