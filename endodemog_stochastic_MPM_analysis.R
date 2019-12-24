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

