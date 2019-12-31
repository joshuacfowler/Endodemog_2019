## purpose use the fitted vital rate models to compute the stochastic growth rate and
## decompose endophyte effects on mean versus variance
source("endodemog_stochastic_MPM_source.R")

# Read in data and stan objects -------------------------------------------
joshpath <- "/Users/joshuacfowler/Dropbox/EndodemogData/"
tompath <- "C:/Users/tm9/Dropbox/EndodemogData/"
#tompath <- "C:/Users/tm634/Dropbox/EndodemogData/"

# Read in data and stan objects -------------------------------------------
## big demography data frame
LTREB_full <- read_csv(paste0(tompath,"Fulldataplusmetadata/LTREB_full.csv"))
## seeds per spikelet data frame
seed_par <- read_csv(paste0(tompath,"Fulldataplusmetadata/LTREB_repro1.csv")) %>% 
  mutate(seed_per_spike = seed/spikelets) %>% 
  mutate(SEED_PER_SPIKE= case_when(species != "AGPE" ~ seed_per_spike,
                                   species == "AGPE" & tillerid_fixed == "multitillermean" ~ seed, # AGPE has some of its seeds data already recorded as seed/spike
                                   species == "AGPE" & tillerid_fixed != "multitillermean" ~ seed_per_spike)) %>% 
  filter(!is.na(SEED_PER_SPIKE)) %>% 
  select(species,SEED_PER_SPIKE) %>% group_by(species) %>% 
  summarise(mean_seeds = mean(SEED_PER_SPIKE))
## I am a little concerned about AGPE mean = 0.427. SHould be 1!
seed_par$mean_seeds[seed_par$species=="AGPE"] <- 1

## recruitment data frame
recruit_par <- bind_rows(read_csv(paste0(tompath,"Fulldataplusmetadata/AGPE_s_to_s_data.csv")),
                              read_csv(paste0(tompath,"Fulldataplusmetadata/POAL_s_to_s_data.csv")),
                              read_csv(paste0(tompath,"Fulldataplusmetadata/POSY_s_to_s_data.csv")),
                              read_csv(paste0(tompath,"Fulldataplusmetadata/FESU_s_to_s_data.csv")),
                              read_csv(paste0(tompath,"Fulldataplusmetadata/LOAR_s_to_s_data.csv")),
                              read_csv(paste0(tompath,"Fulldataplusmetadata/ELVI_s_to_s_data.csv")),
                              read_csv(paste0(tompath,"Fulldataplusmetadata/ELRI_s_to_s_data.csv"))
) %>% 
  mutate(recruitment = tot_recruit_t1 /tot_seed_t) %>% 
  group_by(species) %>% 
  summarise(mean_rec = mean(recruitment))

# find max sizes
max_size <- LTREB_full %>% 
  select(species,size_t) %>% 
  filter(!is.na(size_t)) %>% 
  group_by(species) %>% 
  summarise(max_size = quantile(size_t,probs=0.975))

# stan fits for surv, grow, flow, fert, spikelets
surv_fit <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/surv_fit_fixed.rds"))
grow_fit <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/grow_fit_fixed.rds"))
flow_fit <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/flow_fit_fixed.rds"))
fert_fit <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/fert_glmm.rds"))#read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/fert_fit_fixed_endo_mean.rds"))
spike_fit <- read_rds(paste0(tompath,"Fulldataplusmetadata/SppRFX/spike_fit_no_endo.rds"))

surv_par <- extract(surv_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,tau_year))
grow_par <- extract(grow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,phi,tau_year))
flow_par <- extract(flow_fit, pars = quote_bare(beta0,betasize,betaendo,betaorigin,tau_year))
fert_par <- fixef(fert_fit)
spike_par <- extract(spike_fit, pars = quote_bare(beta0,betasize,betaorigin,tau_year))

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
spp_names <- c(data.frame(cbind(unique(LTREB_full$species),
                                as.integer(as.numeric(as.factor(unique(LTREB_full$species)))))
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

focal_spp <- "FESU"
focal_index <- which(spp_names==focal_spp)
focal_month <- c(9,7,7,5,7,5,5,NA)[focal_index]
cut_n <- 10

mean_surv <- LTREB_full %>% 
  mutate(size_bin = as.integer(cut_interval(logsize_t,cut_n))) %>% 
  filter(species==focal_spp) %>% 
  group_by(species,size_bin,endo_01) %>% 
  summarise(mean_size = mean(logsize_t),
            mean_surv = mean(surv_t1),
            n_surv = n()) %>% 
  mutate(endo_pch=ifelse(endo_01==0,1,16))

x_seq <- seq(min(LTREB_full$logsize_t,na.rm=T),max(LTREB_full$logsize_t,na.rm=T),length.out = 50)
surv_mean<-unlist(lapply(rstan::extract(surv_fit, pars = c(paste0("beta0[",focal_index,"]"),
                                                  paste0("betasize[",focal_index,"]"),
                                                  paste0("betaendo[",focal_index,"]"),
                                                  paste0("tau_year[",focal_index,",1,",1:11,"]"),
                                                  paste0("tau_year[",focal_index,",2,",1:11,"]")
                                                    
)),mean))

win.graph()
par(mfrow=c(3,1),mar=c(5,5,1,1))
plot(mean_surv$mean_size,mean_surv$mean_surv,
     pch=mean_surv$endo_pch,ylim=c(0,1),col=alpha(spp_cols[focal_index],.5),xlab="log Size",ylab="Pr.(Survival)",cex.lab=1.4,
     lwd=2,cex=2 + 3*(mean_surv$n_surv/max(mean_surv$n_surv)))
lines(x_seq,invlogit(surv_mean[1]+surv_mean[2]*x_seq),lwd=4,lty=2,col=spp_cols[focal_index])
lines(x_seq,invlogit(surv_mean[1]+surv_mean[2]*x_seq + surv_mean[3]),lwd=4,lty=1,col=spp_cols[focal_index])
legend("bottomright",legend=c("E+","E-"),lty=1:2,lwd=3,pch=c(16,1),
       col=spp_cols[focal_index],cex=1.8,bty="n")


##now year-specific
year_surv <- LTREB_full %>% 
  mutate(size_bin = as.integer(cut_interval(logsize_t,cut_n))) %>% 
  filter(species==focal_spp) %>% 
  group_by(year_t,species,size_bin,endo_01) %>% 
  summarise(mean_size = mean(logsize_t),
            mean_surv = mean(surv_t1),
            n_surv = n()) %>% 
  mutate(endo_pch=ifelse(endo_01==0,1,16))
year_surv_years <- unique(year_surv$year_t)

focal_rfx <- matrix(surv_mean[4:25],nrow=11,ncol=2,byrow = F)


plot(mean_surv$mean_size[mean_surv$endo_01==0],mean_surv$mean_surv[mean_surv$endo_01==0],
     pch=mean_surv$endo_pch[mean_surv$endo_01==0],ylim=c(0,1),col=alpha(spp_cols[focal_index],.5),xlab="log Size",ylab="Pr.(Survival)",cex.lab=1.4,
     lwd=2,cex=2 + 3*(mean_surv$n_surv[mean_surv$endo_01==0]/max(mean_surv$n_surv[mean_surv$endo_01==0],na.rm=T)))
for(t in 1:length(year_surv_years)){
  lines(x_seq,invlogit(surv_mean[1]+surv_mean[2]*x_seq + 
                         focal_rfx[t,1]),lwd=4,col=spp_cols[focal_index],lty=2)
}

plot(mean_surv$mean_size[mean_surv$endo_01==1],mean_surv$mean_surv[mean_surv$endo_01==1],
     pch=mean_surv$endo_pch[mean_surv$endo_01==1],ylim=c(0,1),col=alpha(spp_cols[focal_index],.5),xlab="log Size",ylab="Pr.(Survival)",cex.lab=1.4,
     lwd=2,cex=2 + 3*(mean_surv$n_surv[mean_surv$endo_01==1]/max(mean_surv$n_surv[mean_surv$endo_01==1],na.rm=T)))
for(t in 1:length(year_surv_years)){
  lines(x_seq,invlogit(surv_mean[1]+surv_mean[2]*x_seq + surv_mean[3]+
                         focal_rfx[t,2]),lwd=4,col=spp_cols[focal_index],lty=1)
}
# climate data ------------------------------------------------------------
climate <- read_csv(file = "C:/Users/tm9/Dropbox/EndodemogData/PRISMClimateData_BrownCo.csv") %>% 
  #mutate(year = year(Date), month = month(Date), day = day(Date)) %>% 
  mutate(year = as.numeric(str_sub(Date,1,4)), 
         month = as.numeric(str_sub(Date,6,7)), 
         day = as.numeric(str_sub(Date,9,10)))%>% 
  rename(ppt = `ppt (mm)`, tmean = `tmean (degrees C)`) %>% 
  mutate(site_lat = 39.235900000000, site_long = -86.218100000000) %>% 
  mutate(census_month = focal_month, 
         climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>% 
  filter(climate_year != 2006) %>% 
  group_by(climate_year) %>% 
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean)) %>% 
  filter(climate_year <= max(LTREB_full$year_t))


par(mfrow=c(2,1),mar=c(5,5,1,1))

plot(climate$`Cumulative PPT (mm)`,focal_rfx[,1],col=alpha(spp_cols[focal_index],.5),cex.lab=1.4,
     cex=3,lwd=2,xlab="Annual precipitation (mm)",ylab="Survival (standardized)",pch=1,
     ylim=c(min(focal_rfx),max(focal_rfx)))
points(climate$`Cumulative PPT (mm)`,focal_rfx[,2],pch=16,col=alpha(spp_cols[focal_index],.5),cex=3)
abline(coef(lm(focal_rfx[,1] ~ climate$`Cumulative PPT (mm)`)),lty=2,lwd=2,col=spp_cols[focal_index])
abline(coef(lm(focal_rfx[,2] ~ climate$`Cumulative PPT (mm)`)),lwd=2,col=spp_cols[focal_index])

plot(climate$`Mean Temp. (C°)`,focal_rfx[,1],col=alpha(spp_cols[focal_index],.5),cex.lab=1.4,
     cex=3,lwd=2,xlab="Mean temp (C)",ylab="Survival (standardized)",pch=1,
     ylim=c(min(focal_rfx),max(focal_rfx)))
points(climate$`Mean Temp. (C°)`,focal_rfx[,2],pch=16,col=alpha(spp_cols[focal_index],.5),cex=3)
abline(coef(lm(focal_rfx[,1] ~ climate$`Mean Temp. (C°)`)),lty=2)
abline(coef(lm(focal_rfx[,2] ~ climate$`Mean Temp. (C°)`)))


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

## try out the make_params function
## note that this is how species are indexing:
cbind(unique(LTREB_full$species),
      as.integer(as.numeric(as.factor(unique(LTREB_full$species)))))

## draw a posterior sample
i <- 765
## get "mean" parameters
test_params <- make_params(species=4,
            endo_mean=1,
            endo_var=1,
            draw=i,
            max_size=max_size,
            rfx=T,
            year=1,
            surv_par=surv_par,
            grow_par=grow_par,
            flow_par=flow_par,
            fert_par=fert_par,
            spike_par=spike_par,
            seed_par=seed_par,
            recruit_par=recruit_par)
lambda(bigmatrix(test_params)$MPMmat)

# MPM stuff ---------------------------------------------------------------

## estimate endo effects on mean lambda and variance in lambda  
n_draws <- 500
post_draws <- sample.int(5000,size=n_draws)
lambda_mean <- array(dim = c(8,2,n_draws))
for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
        lambda_mean[s,e,i] <- lambda(bigmatrix(make_params(species=s,
                                                         endo_mean=(e-1),
                                                         endo_var=(e-1),
                                                         draw=post_draws[i],
                                                         max_size=max_size,
                                                         rfx=F,
                                                         surv_par=surv_par,
                                                         grow_par=grow_par,
                                                         flow_par=flow_par,
                                                         fert_par=fert_par,
                                                         spike_par=spike_par,
                                                         seed_par=seed_par,
                                                         recruit_par=recruit_par))$MPMmat)
    }
    lambda_mean[8,e,i] <- mean(lambda_mean[1:7,e,i])
  }
}

lambda_mean_diff <- matrix(NA,8,5)
for(s in 1:8){
  lambda_mean_diff[s,1] = mean(lambda_mean[s,2,] - lambda_mean[s,1,])
  lambda_mean_diff[s,2:5] = quantile(lambda_mean[s,2,] - lambda_mean[s,1,],probs=c(0.05,0.25,0.75,0.95))
}

## now do variance in lambda
lambda_hold <- array(dim = c(10,7,2,n_draws))
lambda_var <- array(dim = c(8,2,n_draws))
for(i in 1:length(post_draws)){
  for(e in 1:2){
    for(s in 1:7){
      for(y in 1:10){
        
      lambda_hold[y,s,e,i] <- lambda(bigmatrix(make_params(species=s,
                                                         endo_mean=(e-1),
                                                         endo_var=(e-1),
                                                         draw=post_draws[i],
                                                         max_size=max_size,
                                                         rfx=T,
                                                         year=y,
                                                         surv_par=surv_par,
                                                         grow_par=grow_par,
                                                         flow_par=flow_par,
                                                         fert_par=fert_par,
                                                         spike_par=spike_par,
                                                         seed_par=seed_par,
                                                         recruit_par=recruit_par))$MPMmat)
      }
      lambda_var[s,e,i] <- sd(lambda_hold[,s,e,i])
    }
    lambda_var[8,e,i] <- mean(lambda_var[1:7,e,i])
  }
}

lambda_var_diff <- matrix(NA,8,5)
for(s in 1:8){
  lambda_var_diff[s,1] = mean(lambda_var[s,2,] - lambda_var[s,1,])
  lambda_var_diff[s,2:5] = quantile(lambda_var[s,2,] - lambda_var[s,1,],probs=c(0.05,0.25,0.75,0.95))
}

spp_names <- c(data.frame(cbind(unique(LTREB_full$species),
                                as.integer(as.numeric(as.factor(unique(LTREB_full$species)))))
) %>% 
  arrange(X2) %>% select(X1)); spp_names<- c(as.character(spp_names$X1),"Mean")
spp_cols <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02","#a6761d","black")
spp_alpha <- 0.75

win.graph()
par(mfrow=c(1,2),mar=c(5,5,1,1))
plot(rep(0,8),1:8,type="n",xlim=c(-2,2),axes=F,ylab=" ",xlab=expression(paste("Endophyte effect on ",bar(lambda))),cex.lab=1.4)
axis(side=1)
axis(side=2,at=1:8,labels=spp_names[1:8],las=1,cex.axis=1.6)
arrows(0,1,0,8,lty=2,col="gray",length=0)
arrows(lambda_mean_diff[,2],1:8,lambda_mean_diff[,5],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(lambda_mean_diff[,3],1:8,lambda_mean_diff[,4],1:8,length=0,lwd=6,col=alpha(spp_cols,spp_alpha))
points(lambda_mean_diff[,1],1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))

plot(rep(0,8),1:8,type="n",xlim=c(-0.5,0.5),axes=F,ylab=" ",xlab=expression(paste("Endophyte effect on Var(",lambda,")")),cex.lab=1.4)
axis(side=1)
axis(side=2,at=1:8,labels=spp_names[1:8],las=1,cex.axis=1.6)
arrows(0,1,0,8,lty=2,col="gray",length=0)
arrows(lambda_var_diff[,2],1:8,lambda_var_diff[,5],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(lambda_var_diff[,3],1:8,lambda_var_diff[,4],1:8,length=0,lwd=6,col=alpha(spp_cols,spp_alpha))
points(lambda_var_diff[,1],1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))

## make a list of year-specific transition matrices

lambdaS_out <- array(dim = c(8,4,n_draws))
for(i in 1:length(post_draws)){
for(s in 1:7){
  eminus_list <- eplus_list <- eplus__mean_only_list <- eplus__var_only_list <- list()
    for(y in 1:10){
    eminus_list[[y]] <- bigmatrix(make_params(species=s,endo_mean=0,endo_var=0,draw=post_draws[i],max_size=max_size,rfx=T,
                                              year=y,surv_par=surv_par,grow_par=grow_par,flow_par=flow_par,fert_par=fert_par,
                                              spike_par=spike_par,seed_par=seed_par,recruit_par=recruit_par))$MPMmat
    eplus_list[[y]] <- bigmatrix(make_params(species=s,endo_mean=1,endo_var=1,draw=post_draws[i],max_size=max_size,rfx=T,
                                              year=y,surv_par=surv_par,grow_par=grow_par,flow_par=flow_par,fert_par=fert_par,
                                              spike_par=spike_par,seed_par=seed_par,recruit_par=recruit_par))$MPMmat
    eplus__mean_only_list[[y]] <- bigmatrix(make_params(species=s,endo_mean=1,endo_var=0,draw=post_draws[i],max_size=max_size,rfx=T,
                                             year=y,surv_par=surv_par,grow_par=grow_par,flow_par=flow_par,fert_par=fert_par,
                                             spike_par=spike_par,seed_par=seed_par,recruit_par=recruit_par))$MPMmat
    eplus__var_only_list[[y]] <- bigmatrix(make_params(species=s,endo_mean=0,endo_var=1,draw=post_draws[i],max_size=max_size,rfx=T,
                                                        year=y,surv_par=surv_par,grow_par=grow_par,flow_par=flow_par,fert_par=fert_par,
                                                        spike_par=spike_par,seed_par=seed_par,recruit_par=recruit_par))$MPMmat
  }
  lambdaS_out[s,1,i] <- lambdaSim(eminus_list)$lambdaS
  lambdaS_out[s,2,i] <- lambdaSim(eplus__mean_only_list)$lambdaS
  lambdaS_out[s,3,i] <- lambdaSim(eplus__var_only_list)$lambdaS
  lambdaS_out[s,4,i] <- lambdaSim(eplus_list)$lambdaS
}
  lambdaS_out[8,1,i] <- mean(lambdaS_out[1:7,1,i])
  lambdaS_out[8,2,i] <- mean(lambdaS_out[1:7,2,i])
  lambdaS_out[8,3,i] <- mean(lambdaS_out[1:7,3,i])
  lambdaS_out[8,4,i] <- mean(lambdaS_out[1:7,4,i])
}

lambdaS_diff <- lambdaS_diff_mean_only <- lambdaS_diff_var_only <- matrix(NA,8,7)
for(s in 1:8){
  lambdaS_diff[s,1] = mean(lambdaS_out[s,4,] - lambdaS_out[s,1,])
  lambdaS_diff[s,2:7] = quantile(lambdaS_out[s,4,] - lambdaS_out[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
  
  lambdaS_diff_mean_only[s,1] = mean(lambdaS_out[s,2,] - lambdaS_out[s,1,])
  lambdaS_diff_mean_only[s,2:7] = quantile(lambdaS_out[s,2,] - lambdaS_out[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
  
  lambdaS_diff_var_only[s,1] = mean(lambdaS_out[s,3,] - lambdaS_out[s,1,])
  lambdaS_diff_var_only[s,2:7] = quantile(lambdaS_out[s,3,] - lambdaS_out[s,1,],probs=c(0.05,0.125,0.25,0.75,0.875,0.95))
}


win.graph()
par(mfrow=c(4,2),mar=c(5,5,4,1),xpd=T)
for(s in 1:8){
  plot(rep(0,3),1:3,type="n",xlim=c(-1.2,1.2),axes=F,ylab=" ",xlab=expression(paste("Endophyte effect on ",lambda,"S")),cex.lab=1.4)
  axis(side=1)
  abline(v=0,lty=2,col="gray")
  arrows(lambdaS_diff[s,2],1,lambdaS_diff[s,7],1,length=0,lwd=2,col=alpha(spp_cols[s],spp_alpha))
  arrows(lambdaS_diff[s,3],1,lambdaS_diff[s,6],1,length=0,lwd=4,col=alpha(spp_cols[s],spp_alpha))
  arrows(lambdaS_diff[s,4],1,lambdaS_diff[s,5],1,length=0,lwd=6,col=alpha(spp_cols[s],spp_alpha))
  points(lambdaS_diff[s,1],1,cex=3,pch=16,col=alpha(spp_cols[s],spp_alpha))
  
  arrows(lambdaS_diff_mean_only[s,2],2,lambdaS_diff_mean_only[s,7],2,length=0,lwd=2,col=alpha(spp_cols[s],spp_alpha))
  arrows(lambdaS_diff_mean_only[s,3],2,lambdaS_diff_mean_only[s,6],2,length=0,lwd=4,col=alpha(spp_cols[s],spp_alpha))
  arrows(lambdaS_diff_mean_only[s,4],2,lambdaS_diff_mean_only[s,5],2,length=0,lwd=6,col=alpha(spp_cols[s],spp_alpha))
  points(lambdaS_diff_mean_only[s,1],2,cex=3,pch=0,col=alpha(spp_cols[s],spp_alpha))
  
  arrows(lambdaS_diff_var_only[s,2],3,lambdaS_diff_var_only[s,7],3,length=0,lwd=2,col=alpha(spp_cols[s],spp_alpha))
  arrows(lambdaS_diff_var_only[s,3],3,lambdaS_diff_var_only[s,6],3,length=0,lwd=4,col=alpha(spp_cols[s],spp_alpha))
  arrows(lambdaS_diff_var_only[s,4],3,lambdaS_diff_var_only[s,5],3,length=0,lwd=6,col=alpha(spp_cols[s],spp_alpha))
  points(lambdaS_diff_var_only[s,1],3,cex=3,pch=2,col=alpha(spp_cols[s],spp_alpha))
  
}


win.graph()
par(mar=c(5,5,4,1),xpd=T)
plot(rep(0,8),1:8,type="n",xlim=c(-1.2,1.2),axes=F,ylab=" ",xlab=expression(paste("Endophyte effect on ",lambda,"S")),cex.lab=1.4)
axis(side=1)
axis(side=2,at=(1:8)+0.5,labels=spp_names[1:8],las=1,cex.axis=1.6)
arrows(0,1,0,9,lty=2,col="gray",length=0)
arrows(lambdaS_diff[,2],1:8,lambdaS_diff[,7],1:8,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(lambdaS_diff[,3],1:8,lambdaS_diff[,6],1:8,length=0,lwd=4,col=alpha(spp_cols,spp_alpha))
arrows(lambdaS_diff[,4],1:8,lambdaS_diff[,5],1:8,length=0,lwd=6,col=alpha(spp_cols,spp_alpha))
points(lambdaS_diff[,1],1:8,cex=3,pch=16,col=alpha(spp_cols,spp_alpha))

arrows(lambdaS_diff_mean_only[,2],(1:8)+0.33,lambdaS_diff_mean_only[,7],(1:8)+0.33,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(lambdaS_diff_mean_only[,3],(1:8)+0.33,lambdaS_diff_mean_only[,6],(1:8)+0.33,length=0,lwd=4,col=alpha(spp_cols,spp_alpha))
arrows(lambdaS_diff_mean_only[,4],(1:8)+0.33,lambdaS_diff_mean_only[,5],(1:8)+0.33,length=0,lwd=6,col=alpha(spp_cols,spp_alpha))
points(lambdaS_diff_mean_only[,1],(1:8)+0.33,cex=3,pch=15,col=alpha(spp_cols,spp_alpha))

arrows(lambdaS_diff_var_only[,2],(1:8)+0.66,lambdaS_diff_var_only[,7],(1:8)+0.66,length=0,lwd=2,col=alpha(spp_cols,spp_alpha))
arrows(lambdaS_diff_var_only[,3],(1:8)+0.66,lambdaS_diff_var_only[,6],(1:8)+0.66,length=0,lwd=4,col=alpha(spp_cols,spp_alpha))
arrows(lambdaS_diff_var_only[,4],(1:8)+0.66,lambdaS_diff_var_only[,5],(1:8)+0.66,length=0,lwd=6,col=alpha(spp_cols,spp_alpha))
points(lambdaS_diff_var_only[,1],(1:8)+0.66,cex=3,pch=17,col=alpha(spp_cols,spp_alpha))

