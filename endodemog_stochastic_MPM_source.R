library(tidyverse)
library(rstan)
rstan_options( auto_write = TRUE )
options( mc.cores = parallel::detectCores() )
library(bayesplot)
library(scales)
library(popbio)
library(lme4)

# misc functions -----------------------------------------------------------
invlogit<-function(x){exp(x)/(1+exp(x))}

quote_bare <- function( ... ){
  substitute( alist(...) ) %>% 
    eval( ) %>% 
    sapply( deparse )
}

# Parameter assembly function ---------------------------------------------
make_params <- function(species,endo_mean,endo_var,draw,original=0,rfx=F,year=NULL,max_size,
                        surv_par,grow_par,flow_par,fert_par,spike_par,seed_par,recruit_par){
  
  if(rfx==F){rfx_surv <- rfx_grow <- rfx_flow <- rfx_fert <- rfx_spike <- 0}
  if(rfx==T){
    ## timing and survival and growth (size_t / y_t1) is meant to line up with reproduction (size_t1 / y_t1)
    rfx_surv <- surv_par$tau_year[draw,species,(endo_var+1),(year+1)]; 
    rfx_grow <- grow_par$tau_year[draw,species,(endo_var+1),(year+1)];
    rfx_flow <- flow_par$tau_year[draw,species,(endo_var+1),year];
    #rfx_fert <- fert_par$tau_year[draw,species,year]; #no endo effects on variance here
    rfx_spike <- spike_par$tau_year[draw,year]; #no endo effects or species differences
  }
  
  params <- c()
  #survival
  params$surv_int <- surv_par$beta0[draw,species] + 
    endo_mean * surv_par$betaendo[draw,species] + 
    original * surv_par$betaorigin[draw,species] + rfx_surv
  params$surv_slope <- surv_par$betasize[draw,species]
  #growth
  params$grow_int <- grow_par$beta0[draw,species] + 
    endo_mean * grow_par$betaendo[draw,species] + 
    original * grow_par$betaorigin[draw,species] + rfx_grow
  params$grow_slope <- grow_par$betasize[draw,species]  
  params$grow_phi <- grow_par$phi[draw,species] 
  #flowering
  params$flow_int <- flow_par$beta0[draw,species] + 
    endo_mean * flow_par$betaendo[draw,species] + 
    original * flow_par$betaorigin[draw,species] + rfx_flow
  params$flow_slope <- flow_par$betasize[draw,species]  
  #fertility
  #params$fert_int <- fert_par$beta0[draw,species] + 
  #  endo_mean * fert_par$betaendo[draw,species] + 
  #  original * fert_par$betaorigin[draw,species] + rfx_fert
  #params$fert_slope <- fert_par$betasize[draw,species]  
  ## hacking my way through fertility for now with a glmm with spp intercet and other params shared
  params$fert_int <- fert_par[species]
  params$fert_slope <- fert_par[8] + original * fert_par[9]
  #spikelets
  params$spike_int <- spike_par$beta0[draw,species]  + 
    original * spike_par$betaorigin[draw,species] + rfx_spike
  params$spike_slope <- spike_par$betasize[draw,species]  
  #seeds per spikelet
  params$seeds_per_spike <- seed_par$mean_seeds[species] #no posterior sampling here
  #recruits per seed
  params$recruits_per_seed <- recruit_par$mean_rec[species] #no posterior sampling here
  #tack on max size
  params$max_size <- max_size$max_size[species]
  
  return(params)
}

# Vital rate functions ----------------------------------------------------
sx<-function(x,params){
  invlogit(params$surv_int + params$surv_slope*log(x))
}

gxy <- function(x,y,params){
  grow_mean <- params$grow_int + params$grow_slope*log(x)
  grow<-dnbinom(x=y,mu=exp(grow_mean),size=exp(params$grow_phi),log=F)
  truncLower<-dnbinom(x=0,mu=exp(grow_mean),size=exp(params$grow_phi),log=F)
  truncUpper<-sum(dnbinom(x=params$max_size:10000,mu=exp(grow_mean),size=exp(params$grow_phi),log=F))
  return(grow/(1-(truncLower+truncUpper)))
}

pxy<-function(x,y,params){
  sx(x,params) * gxy(x,y,params)
}

fx<-function(x, params){
    flw <- invlogit(params$flow_int + params$flow_slope*log(x))
    fert <- exp(params$fert_int + params$fert_slope*log(x))
    spike <- exp(params$spike_int + params$spike_slope*log(x))
  seedlings <- flw * fert * spike * params$seeds_per_spike * params$recruits_per_seed
  return(seedlings)
}

# Bigmatrix function ------------------------------------------------------
# note from Tom: we should think about adding a reproductive delay
bigmatrix<-function(params){   
  matdim<-params$max_size
  y <- 1:matdim 
  Fmat <- matrix(0,matdim,matdim)
  Fmat[1,]<-fx(x = y, params)
  Tmat <-matrix(0,matdim,matdim)
  Tmat<-t(outer(y,y,pxy,params))
  MPMmat<-Tmat + Fmat
  return(list(MPMmat = MPMmat, Tmat = Tmat, Fmat = Fmat))
}

# lambdaS function##########################################################
lambdaSim<-function(mat_list, ## a list of transition matrices, each corresponding to a study year
                    max_yrs=500 ## how many years the simulation runs (arbitrarily large)
){
  ## grab the dimension of the projection matrix
  matdim<-dim(mat_list[[1]])[1]
  ## grab the number of study years / matrices we have available
  n_years <- length(mat_list)
  ## vector that will hold year-by-year growth rates
  rtracker <- rep(0,max_yrs)
  ## initial vector of population structure -- note that this sums to one, which will be convenient
  n0 <- rep(1/matdim,matdim)
  for(t in 1:max_yrs){ #Start loop
    ## for each year, randomly sample one of the matrices
    A_t <- mat_list[[sample.int(n=n_years,size=1)]]
    ## project the population one step forward
    n0 <- A_t %*% n0
    ## total population size after one year of growth
    N  <- sum(n0)
    ## calculate r as log(N_t+1 / N_t), note that here N_t=1
    rtracker[t]<-log(N)
    ## rescale population vector to sum to one, so the same trick works again next time step
    n0 <-n0/N
  }
  #discard initial values (to get rid of transient)
  burnin    <- round(max_yrs*0.1)
  #Finish and return
  log_lambdaS <- mean(rtracker[-c(1:burnin)])
  lambdaS<-exp(log_lambdaS)
  return(list(log_lambdaS=log_lambdaS,lambdaS=lambdaS))
}