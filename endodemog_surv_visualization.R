## Title: Grass endophyte vital rate variance with a bayesian framework
## Purpose: Visualizes output of endophyte effect on year variance from 
## model outputs for Survival models 
## Authors: Joshua and Tom
#############################################################
library(tidyverse)
library(rstan)
library(StanHeaders)
library(shinystan)
library(bayesplot)
library(devtools)
library(reshape2)
library(ggplot2)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

#  source in raw data from endodemog_surv_mixed_stan.R 

## Read in the survival model output for all species
survPOAL <- readRDS(file = "endodemog_surv_full_POAL.rds")
survPOSY <- readRDS(file = "endodemog_surv_full_POSY.rds")
survLOAR <- readRDS(file = "endodemog_surv_full_LOAR.rds")
survELVI <- readRDS(file = "endodemog_surv_full_ELVI.rds")
survELRI <- readRDS(file = "endodemog_surv_full_ELRI.rds")
survFESU <- readRDS(file = "endodemog_surv_full_FESU.rds")
survAGPE <- readRDS(file = "endodemog_surv_full_AGPE.rds")



# save posteriors within dataframe
params = c("alpha", "beta", "tau_year", "sigma_0")

post_survPOAL <- as.data.frame(survPOAL, pars = params)
post_survPOSY <- as.data.frame(survPOSY, pars = params)
post_survLOAR <- as.data.frame(survLOAR, pars = params)
post_survELVI <- as.data.frame(survELVI, pars = params)
post_survELRI <- as.data.frame(survELRI, pars = params)
post_survFESU <- as.data.frame(survFESU, pars = params)
post_survAGPE <- as.data.frame(survAGPE, pars = params)


yearcolors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
                "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a",
                "#ffff99")
colors2 <- c("#e31a1c","#ff7f00")
# POAL survival

# actual data points for probability of survival
surv_bin0 <- POAL_data1 %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, minus = 0, plus = 1)) %>% 
  filter(endo == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T))

surv_bin1 <- POAL_data1 %>% 
  select(surv_t1, logsize_t, year_t, endo) %>%
  mutate(endo = recode(endo, minus = 0, plus = 1)) %>% 
  filter(endo == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T))


# fitted model line
nvalues <- length(POAL_size_dat$logsize_t)
xdummy <- seq(min(POAL_size_dat$logsize_t), max(POAL_size_dat$logsize_t), length.out = nvalues)

# E- 
ydummy_eminus <- as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*0 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                    + mean(post_survPOAL$'beta[4]')*xdummy*0
                                    + mean(mean(post_survPOAL$'tau_year[1,1]'),mean(post_survPOAL$'tau_year[1,2]'),mean(post_survPOAL$'tau_year[1,3]'),mean(post_survPOAL$'tau_year[1,4]'),mean(post_survPOAL$'tau_year[1,5]'),mean(post_survPOAL$'tau_year[1,6]'),mean(post_survPOAL$'tau_year[1,7]'),mean(post_survPOAL$'tau_year[1,8]'),mean(post_survPOAL$'tau_year[1,9]'),mean(post_survPOAL$'tau_year[1,10]'),mean(post_survPOAL$'tau_year[1,11]'))))
ydummy_eminus_y1 <- as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*0 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                             + mean(post_survPOAL$'beta[4]')*xdummy*0
                             + mean(post_survPOAL$'tau_year[1,1]')))
ydummy_eminus_y2 <- as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*0 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                             + mean(post_survPOAL$'beta[4]')*xdummy*0
                             + mean(post_survPOAL$'tau_year[1,2]')))
ydummy_eminus_y3 <- as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*0 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                             + mean(post_survPOAL$'beta[4]')*xdummy*0
                             + mean(post_survPOAL$'tau_year[1,3]')))
ydummy_eminus_y4 <- as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*0 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                             + mean(post_survPOAL$'beta[4]')*xdummy*0
                             + mean(post_survPOAL$'tau_year[1,4]')))
ydummy_eminus_y5 <- as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*0 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                             + mean(post_survPOAL$'beta[4]')*xdummy*0
                             + mean(post_survPOAL$'tau_year[1,5]')))
ydummy_eminus_y6 <- as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*0 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                             + mean(post_survPOAL$'beta[4]')*xdummy*0
                             + mean(post_survPOAL$'tau_year[1,6]')))
ydummy_eminus_y7 <- as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*0 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                             + mean(post_survPOAL$'beta[4]')*xdummy*0
                             + mean(post_survPOAL$'tau_year[1,7]')))
ydummy_eminus_y8 <- as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*0 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                             + mean(post_survPOAL$'beta[4]')*xdummy*0
                             + mean(post_survPOAL$'tau_year[1,8]')))
ydummy_eminus_y9 <- as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*0 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                             + mean(post_survPOAL$'beta[4]')*xdummy*0
                             + mean(post_survPOAL$'tau_year[1,9]')))
ydummy_eminus_y10 <- as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*0 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                             + mean(post_survPOAL$'beta[4]')*xdummy*0
                             + mean(post_survPOAL$'tau_year[1,10]')))
ydummy_eminus_y11 <- as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*0 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                             + mean(post_survPOAL$'beta[4]')*xdummy*0
                             + mean(post_survPOAL$'tau_year[1,11]')))
POALfitsminus0 <- as.data.frame(cbind(xdummy, ydummy_eminus))
POALfitsminus1 <- as.data.frame(cbind(xdummy,ydummy_eminus_y1,ydummy_eminus_y2, ydummy_eminus_y3, ydummy_eminus_y4, ydummy_eminus_y5, ydummy_eminus_y6, ydummy_eminus_y7, ydummy_eminus_y8, ydummy_eminus_y9, ydummy_eminus_y10, ydummy_eminus_y11))
POALfitsminus <- melt(POALfitsminus1, id = "xdummy",
                      measure.vars = c("ydummy_eminus_y1","ydummy_eminus_y2", 
                                       "ydummy_eminus_y3", "ydummy_eminus_y4", 
                                       "ydummy_eminus_y5", "ydummy_eminus_y6", 
                                       "ydummy_eminus_y7", "ydummy_eminus_y8", 
                                       "ydummy_eminus_y9", "ydummy_eminus_y10",
                                       "ydummy_eminus_y11")) %>% 
  mutate(Year = recode(variable, "ydummy_eminus_y1" = '1',"ydummy_eminus_y2" = "2", 
                                  "ydummy_eminus_y3" = "3", "ydummy_eminus_y4" = "4", 
                                  "ydummy_eminus_y5" = "5", "ydummy_eminus_y6" = "6", 
                                  "ydummy_eminus_y7" = "7", "ydummy_eminus_y8" = "8", 
                                  "ydummy_eminus_y9" = "9", "ydummy_eminus_y10" = "10",
                                  "ydummy_eminus_y11" = "11") )
ggplot(data = POALfitsminus1)+
  geom_line(aes(x = xdummy, y = ydummy_eminus_y1))+
  geom_line(aes(x = xdummy, y = ydummy_eminus_y2))+
  geom_line(aes(x = xdummy, y = ydummy_eminus_y3), col = "red") +
  geom_line(aes(x = xdummy, y = ydummy_eminus_y4))+
  geom_line(aes(x = xdummy, y = ydummy_eminus_y5))+
  geom_line(aes(x = xdummy, y = ydummy_eminus_y6))+
  geom_line(aes(x = xdummy, y = ydummy_eminus_y7))+
  geom_line(aes(x = xdummy, y = ydummy_eminus_y8))+
  geom_line(aes(x = xdummy, y = ydummy_eminus_y9))+
  geom_line(aes(x = xdummy, y = ydummy_eminus_y10))+
  geom_line(aes(x = xdummy, y = ydummy_eminus_y11))
  

# with data points
ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .7) +
  geom_point(data = surv_bin0, aes(x = mean_size, y = mean_surv, color = Year)) +
  labs(title = "POAL E- Survival Probability", x = "log(size_t)", y = "Prob. of Survival") +
  scale_color_manual(values=yearcolors)
# without datapoints
ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .7) +
  labs(title = "POAL E- Survival Probability", x = "log(size_t)", y = "Prob. of Survival") +
  scale_color_manual(values=yearcolors)

# endo effect on mean
ggplot(data = POALfitsminus0) +
  geom_line(aes(x = xdummy, y = ydummy_eminus)) +
  labs(title = "POAL E- Survival Probability", x = "log(size_t)", y = "Prob. of Survival") +
  scale_color_manual(values=colors2)

# E+
ydummy_eplus <- as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*1 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                   + mean(post_survPOAL$'beta[4]')*xdummy*1
                                   + mean(mean(post_survPOAL$'tau_year[2,1]'),mean(post_survPOAL$'tau_year[2,2]'),mean(post_survPOAL$'tau_year[2,3]'),mean(post_survPOAL$'tau_year[2,4]'),mean(post_survPOAL$'tau_year[2,5]'),mean(post_survPOAL$'tau_year[2,6]'),mean(post_survPOAL$'tau_year[2,7]'),mean(post_survPOAL$'tau_year[2,8]'),mean(post_survPOAL$'tau_year[2,9]'),mean(post_survPOAL$'tau_year[2,10]'),mean(post_survPOAL$'tau_year[2,11]'))))
ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*1 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                             + mean(post_survPOAL$'beta[4]')*xdummy*1
                             + mean(post_survPOAL$'tau_year[2,1]')))
ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*1 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                             + mean(post_survPOAL$'beta[4]')*xdummy*1
                             + mean(post_survPOAL$'tau_year[2,2')))
ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*1 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                             + mean(post_survPOAL$'beta[4]')*xdummy*1
                             + mean(post_survPOAL$'tau_year[2,3]')))
ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*1 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                             + mean(post_survPOAL$'beta[4]')*xdummy*1
                             + mean(post_survPOAL$'tau_year[2,4]')))
ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*1 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                             + mean(post_survPOAL$'beta[4]')*xdummy*1
                             + mean(post_survPOAL$'tau_year[2,5]')))
ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*1 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                             + mean(post_survPOAL$'beta[4]')*xdummy*1
                             + mean(post_survPOAL$'tau_year[2,6]')))
ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*1 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                             + mean(post_survPOAL$'beta[4]')*xdummy*1
                             + mean(post_survPOAL$'tau_year[2,7]')))
ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*1 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                             + mean(post_survPOAL$'beta[4]')*xdummy*1
                             + mean(post_survPOAL$'tau_year[2,8]')))
ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*1 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                             + mean(post_survPOAL$'beta[4]')*xdummy*1
                             + mean(post_survPOAL$'tau_year[2,9]')))
ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*1 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                             + mean(post_survPOAL$'beta[4]')*xdummy*1
                             + mean(post_survPOAL$'tau_year[2,10]')))
ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_survPOAL$alpha) + mean(post_survPOAL$'beta[1]')*xdummy + mean(post_survPOAL$'beta[2]')*1 + mean(post_survPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                             + mean(post_survPOAL$'beta[4]')*xdummy*1
                             + mean(post_survPOAL$'tau_year[2,11]')))
POAL_fitsplus0 <- as.data.frame(cbind(xdummy, ydummy_eplus))
POALfitsplus <- as.data.frame(cbind(xdummy,ydummy_eplus_y1,ydummy_eplus_y2, ydummy_eplus_y3, ydummy_eplus_y4, ydummy_eplus_y5, ydummy_eplus_y6, ydummy_eplus_y7, ydummy_eplus_y8, ydummy_eplus_y9, ydummy_eplus_y10, ydummy_eplus_y11))
POALfitsplus <- melt(POALfitsplus, id = "xdummy",
                      measure.vars = c("ydummy_eplus_y1","ydummy_eplus_y2", 
                                       "ydummy_eplus_y3", "ydummy_eplus_y4", 
                                       "ydummy_eplus_y5", "ydummy_eplus_y6", 
                                       "ydummy_eplus_y7", "ydummy_eplus_y8", 
                                       "ydummy_eplus_y9", "ydummy_eplus_y10",
                                       "ydummy_eplus_y11")) %>% 
  mutate(Year = recode(variable, "ydummy_eplus_y1" = '1',"ydummy_eplus_y2" = "2", 
                       "ydummy_eplus_y3" = "3", "ydummy_eplus_y4" = "4", 
                       "ydummy_eplus_y5" = "5", "ydummy_eplus_y6" = "6", 
                       "ydummy_eplus_y7" = "7", "ydummy_eplus_y8" = "8", 
                       "ydummy_eplus_y9" = "9", "ydummy_eplus_y10" = "10",
                       "ydummy_eplus_y11" = "11") )

ggplot(data = POALfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  geom_point(data = surv_bin1, aes(x = mean_size, y = mean_surv, color = Year), lwd = 2) +
  ggtitle("POAL E+ Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=yearcolors)
ggplot(data = POALfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("POAL E+ Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=yearcolors)



# effect on mean


POALfits <- as.data.frame(cbind(xdummy, ydummy_eminus, ydummy_eplus))
POALfits <- melt(POALfits, id = "xdummy", measure.vars = c("ydummy_eplus", "ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ydummy_eminus" = "E-", "ydummy_eplus" = "E+"))

surv_bin <- POAL_data1 %>% 
  select(surv_t1, logsize_t, endo) %>%
  mutate(endo = recode(endo, minus = 'E-', plus = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T))
# without point
ggplot(data = POALfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("POAL Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

# with points
ggplot(data = POALfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = surv_bin, aes(x = mean_size, y = mean_surv, col = endo), lwd = 3) +
  ggtitle("POAL Survival Probability") + xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

ggplot(data = POALfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = surv_bin, aes(x = mean_size, y = mean_surv, col = endo), lwd = 3) +
  xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)


# effect on variance as density plot

stan_dens(survPOAL, pars = c("sigma_0[1]", "sigma_0[2]"))

mcmc_areas(post_survPOAL, pars = c("sigma_0[1]", "sigma_0[2]"),
           prob = 0.8, prob_outer = .99 )


ggplot(data = post_survPOAL,) +
  geom_density( aes(`sigma_0[1]`), col = "#ff7f00", lwd = .8) +
  geom_density(aes(`sigma_0[2]`), col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")

ggplot(data = post_survPOAL,) +
  geom_histogram( aes(`sigma_0[1]`),fill = "#ff7f00", alpha = .8, bins = 100) +
  geom_histogram(aes(`sigma_0[2]`),fill = "#e31a1c", alpha = .8, bins = 100) +
  labs(x = expression(σ[year]), y = "Posterior Count")


ggplot(data = post_survPOAL,) +
  stat_density( aes(`sigma_0[1]`), geom = "line", col = "#ff7f00", lwd = .8) +
  stat_density(aes(`sigma_0[2]`), geom = "line", col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")
ggplot(data = post_survPOAL,) +
  geom_line( aes(`sigma_0[1]`), stat = "density", col = "#ff7f00", lwd = .8) +
  geom_line(aes(`sigma_0[2]`), stat = "density", col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")


# POSY#######################

# POSY survival
nvalues <- length(POSY_size_dat$logsize_t)
xdummy <- seq(min(POSY_size_dat$logsize_t), max(POSY_size_dat$logsize_t), length.out = nvalues)

# E- 
ydummy_eminus <- as.vector(invlogit(mean(post_survPOSY$alpha) + mean(post_survPOSY$'beta[1]')*xdummy + mean(post_survPOSY$'beta[2]')*0 + mean(post_survPOSY$'beta[3]')*mean(POSY_origin_dat$origin1)
                                    + mean(post_survPOSY$'beta[4]')*xdummy*0
                                    + mean(mean(post_survPOSY$'tau_year[1,1]'),mean(post_survPOSY$'tau_year[1,2]'),mean(post_survPOSY$'tau_year[1,3]'),mean(post_survPOSY$'tau_year[1,4]'),mean(post_survPOSY$'tau_year[1,5]'),mean(post_survPOSY$'tau_year[1,6]'),mean(post_survPOSY$'tau_year[1,7]'),mean(post_survPOSY$'tau_year[1,8]'),mean(post_survPOSY$'tau_year[1,9]'),mean(post_survPOSY$'tau_year[1,10]'),mean(post_survPOSY$'tau_year[1,11]'))))

# E+
ydummy_eplus <- as.vector(invlogit(mean(post_survPOSY$alpha) + mean(post_survPOSY$'beta[1]')*xdummy + mean(post_survPOSY$'beta[2]')*1 + mean(post_survPOSY$'beta[3]')*mean(POSY_origin_dat$origin1)
                                   + mean(post_survPOSY$'beta[4]')*xdummy*1
                                   + mean(mean(post_survPOSY$'tau_year[2,1]'),mean(post_survPOSY$'tau_year[2,2]'),mean(post_survPOSY$'tau_year[2,3]'),mean(post_survPOSY$'tau_year[2,4]'),mean(post_survPOSY$'tau_year[2,5]'),mean(post_survPOSY$'tau_year[2,6]'),mean(post_survPOSY$'tau_year[2,7]'),mean(post_survPOSY$'tau_year[2,8]'),mean(post_survPOSY$'tau_year[2,9]'),mean(post_survPOSY$'tau_year[2,10]'),mean(post_survPOSY$'tau_year[2,11]'))))
POSYfits <- as.data.frame(cbind(xdummy, ydummy_eminus, ydummy_eplus))
POSYfits <- melt(POSYfits, id = "xdummy", measure.vars = c("ydummy_eplus", "ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ydummy_eminus" = "E-", "ydummy_eplus" = "E+"))

surv_bin <- POSY_data1 %>% 
  select(surv_t1, logsize_t, endo) %>% 
  mutate(endo = recode(endo, "minus" = 'E-', "plus" = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T))
# with points
ggplot(data = POSYfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = surv_bin, aes(x = mean_size, y = mean_surv, col = endo), lwd = 3) +
  xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

ggplot(data = post_survPOSY,) +
  geom_density( aes(`sigma_0[1]`), col = "#ff7f00", lwd = .8) +
  geom_density(aes(`sigma_0[2]`), col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")



########################

# LOAR survival
nvalues <- length(LOAR_size_dat$logsize_t)
xdummy <- seq(min(LOAR_size_dat$logsize_t), max(LOAR_size_dat$logsize_t), length.out = nvalues)

# E- 
ydummy_eminus <- as.vector(invlogit(mean(post_survLOAR$alpha) + mean(post_survLOAR$'beta[1]')*xdummy + mean(post_survLOAR$'beta[2]')*0 + mean(post_survLOAR$'beta[3]')*mean(LOAR_origin_dat$origin)
                                    + mean(post_survLOAR$'beta[4]')*xdummy*0
                                    + mean(mean(post_survLOAR$'tau_year[1,1]'),mean(post_survLOAR$'tau_year[1,2]'),mean(post_survLOAR$'tau_year[1,3]'),mean(post_survLOAR$'tau_year[1,4]'),mean(post_survLOAR$'tau_year[1,5]'),mean(post_survLOAR$'tau_year[1,6]'),mean(post_survLOAR$'tau_year[1,7]'),mean(post_survLOAR$'tau_year[1,8]'),mean(post_survLOAR$'tau_year[1,9]'),mean(post_survLOAR$'tau_year[1,10]'),mean(post_survLOAR$'tau_year[1,11]'))))

# E+
ydummy_eplus <- as.vector(invlogit(mean(post_survLOAR$alpha) + mean(post_survLOAR$'beta[1]')*xdummy + mean(post_survLOAR$'beta[2]')*1 + mean(post_survLOAR$'beta[3]')*mean(LOAR_origin_dat$origin)
                                   + mean(post_survLOAR$'beta[4]')*xdummy*1
                                   + mean(mean(post_survLOAR$'tau_year[2,1]'),mean(post_survLOAR$'tau_year[2,2]'),mean(post_survLOAR$'tau_year[2,3]'),mean(post_survLOAR$'tau_year[2,4]'),mean(post_survLOAR$'tau_year[2,5]'),mean(post_survLOAR$'tau_year[2,6]'),mean(post_survLOAR$'tau_year[2,7]'),mean(post_survLOAR$'tau_year[2,8]'),mean(post_survLOAR$'tau_year[2,9]'),mean(post_survLOAR$'tau_year[2,10]'),mean(post_survLOAR$'tau_year[2,11]'))))
LOARfits <- as.data.frame(cbind(xdummy, ydummy_eminus, ydummy_eplus))
LOARfits <- melt(LOARfits, id = "xdummy", measure.vars = c("ydummy_eplus", "ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ydummy_eminus" = "E-", "ydummy_eplus" = "E+"))

surv_bin <- LOAR_data1 %>% 
  select(surv_t1, logsize_t, endo) %>% 
  mutate(endo = recode(endo, "0" = 'E-', "1" = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T))
# with points
ggplot(data = LOARfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = surv_bin, aes(x = mean_size, y = mean_surv, col = endo), lwd = 3) +
  xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

ggplot(data = post_survLOAR,) +
  geom_density( aes(`sigma_0[1]`), col = "#ff7f00", lwd = .8) +
  geom_density(aes(`sigma_0[2]`), col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")



# ELVI survival
nvalues <- length(ELVI_size_dat$logsize_t)
xdummy <- seq(min(ELVI_size_dat$logsize_t), max(ELVI_size_dat$logsize_t), length.out = nvalues)

# E- 
ydummy_eminus <- as.vector(invlogit(mean(post_survELVI$alpha) + mean(post_survELVI$'beta[1]')*xdummy + mean(post_survELVI$'beta[2]')*0 + mean(post_survELVI$'beta[3]')*mean(ELVI_origin_dat$origin1)
                                    + mean(post_survELVI$'beta[4]')*xdummy*0
                                    + mean(mean(post_survELVI$'tau_year[1,1]'),mean(post_survELVI$'tau_year[1,2]'),mean(post_survELVI$'tau_year[1,3]'),mean(post_survELVI$'tau_year[1,4]'),mean(post_survELVI$'tau_year[1,5]'),mean(post_survELVI$'tau_year[1,6]'),mean(post_survELVI$'tau_year[1,7]'),mean(post_survELVI$'tau_year[1,8]'),mean(post_survELVI$'tau_year[1,9]'),mean(post_survELVI$'tau_year[1,10]'),mean(post_survELVI$'tau_year[1,11]'))))

# E+
ydummy_eplus <- as.vector(invlogit(mean(post_survELVI$alpha) + mean(post_survELVI$'beta[1]')*xdummy + mean(post_survELVI$'beta[2]')*1 + mean(post_survELVI$'beta[3]')*mean(ELVI_origin_dat$origin1)
                                   + mean(post_survELVI$'beta[4]')*xdummy*1
                                   + mean(mean(post_survELVI$'tau_year[2,1]'),mean(post_survELVI$'tau_year[2,2]'),mean(post_survELVI$'tau_year[2,3]'),mean(post_survELVI$'tau_year[2,4]'),mean(post_survELVI$'tau_year[2,5]'),mean(post_survELVI$'tau_year[2,6]'),mean(post_survELVI$'tau_year[2,7]'),mean(post_survELVI$'tau_year[2,8]'),mean(post_survELVI$'tau_year[2,9]'),mean(post_survELVI$'tau_year[2,10]'),mean(post_survELVI$'tau_year[2,11]'))))
ELVIfits <- as.data.frame(cbind(xdummy, ydummy_eminus, ydummy_eplus))
ELVIfits <- melt(ELVIfits, id = "xdummy", measure.vars = c("ydummy_eplus", "ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ydummy_eminus" = "E-", "ydummy_eplus" = "E+"))

surv_bin <- ELVI_data1 %>% 
  select(surv_t1, logsize_t, endo) %>% 
  mutate(endo = recode(endo, "minus" = 'E-', "plus" = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T))
# with points
ggplot(data = ELVIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = surv_bin, aes(x = mean_size, y = mean_surv, col = endo), lwd = 3) +
  xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

ggplot(data = post_survELVI,) +
  geom_density( aes(`sigma_0[1]`), col = "#ff7f00", lwd = .8) +
  geom_density(aes(`sigma_0[2]`), col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")




# ELRI survival
nvalues <- length(ELRI_size_dat$logsize_t)
xdummy <- seq(min(ELRI_size_dat$logsize_t), max(ELRI_size_dat$logsize_t), length.out = nvalues)

# E- 
ydummy_eminus <- as.vector(invlogit(mean(post_survELRI$alpha) + mean(post_survELRI$'beta[1]')*xdummy + mean(post_survELRI$'beta[2]')*0 + mean(post_survELRI$'beta[3]')*mean(ELRI_origin_dat$origin1)
                                    + mean(post_survELRI$'beta[4]')*xdummy*0
                                    + mean(mean(post_survELRI$'tau_year[1,1]'),mean(post_survELRI$'tau_year[1,2]'),mean(post_survELRI$'tau_year[1,3]'),mean(post_survELRI$'tau_year[1,4]'),mean(post_survELRI$'tau_year[1,5]'),mean(post_survELRI$'tau_year[1,6]'),mean(post_survELRI$'tau_year[1,7]'),mean(post_survELRI$'tau_year[1,8]'),mean(post_survELRI$'tau_year[1,9]'),mean(post_survELRI$'tau_year[1,10]'),mean(post_survELRI$'tau_year[1,11]'))))

# E+
ydummy_eplus <- as.vector(invlogit(mean(post_survELRI$alpha) + mean(post_survELRI$'beta[1]')*xdummy + mean(post_survELRI$'beta[2]')*1 + mean(post_survELRI$'beta[3]')*mean(ELRI_origin_dat$origin1)
                                   + mean(post_survELRI$'beta[4]')*xdummy*1
                                   + mean(mean(post_survELRI$'tau_year[2,1]'),mean(post_survELRI$'tau_year[2,2]'),mean(post_survELRI$'tau_year[2,3]'),mean(post_survELRI$'tau_year[2,4]'),mean(post_survELRI$'tau_year[2,5]'),mean(post_survELRI$'tau_year[2,6]'),mean(post_survELRI$'tau_year[2,7]'),mean(post_survELRI$'tau_year[2,8]'),mean(post_survELRI$'tau_year[2,9]'),mean(post_survELRI$'tau_year[2,10]'),mean(post_survELRI$'tau_year[2,11]'))))
ELRIfits <- as.data.frame(cbind(xdummy, ydummy_eminus, ydummy_eplus))
ELRIfits <- melt(ELRIfits, id = "xdummy", measure.vars = c("ydummy_eplus", "ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ydummy_eminus" = "E-", "ydummy_eplus" = "E+"))

surv_bin <- ELRI_data1 %>% 
  select(surv_t1, logsize_t, endo) %>% 
  mutate(endo = recode(endo, "minus" = 'E-', "plus" = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T))
# with points
ggplot(data = ELRIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = surv_bin, aes(x = mean_size, y = mean_surv, col = endo), lwd = 3) +
  xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

ggplot(data = post_survELRI,) +
  geom_density( aes(`sigma_0[1]`), col = "#ff7f00", lwd = .8) +
  geom_density(aes(`sigma_0[2]`), col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")



# FESU survival
nvalues <- length(FESU_size_dat$logsize_t)
xdummy <- seq(min(FESU_size_dat$logsize_t), max(FESU_size_dat$logsize_t), length.out = nvalues)

# E- 
ydummy_eminus <- as.vector(invlogit(mean(post_survFESU$alpha) + mean(post_survFESU$'beta[1]')*xdummy + mean(post_survFESU$'beta[2]')*0 + mean(post_survFESU$'beta[3]')*mean(FESU_origin_dat$origin1)
                                    + mean(post_survFESU$'beta[4]')*xdummy*0
                                    + mean(mean(post_survFESU$'tau_year[1,1]'),mean(post_survFESU$'tau_year[1,2]'),mean(post_survFESU$'tau_year[1,3]'),mean(post_survFESU$'tau_year[1,4]'),mean(post_survFESU$'tau_year[1,5]'),mean(post_survFESU$'tau_year[1,6]'),mean(post_survFESU$'tau_year[1,7]'),mean(post_survFESU$'tau_year[1,8]'),mean(post_survFESU$'tau_year[1,9]'),mean(post_survFESU$'tau_year[1,10]'),mean(post_survFESU$'tau_year[1,11]'))))

# E+
ydummy_eplus <- as.vector(invlogit(mean(post_survFESU$alpha) + mean(post_survFESU$'beta[1]')*xdummy + mean(post_survFESU$'beta[2]')*1 + mean(post_survFESU$'beta[3]')*mean(FESU_origin_dat$origin1)
                                   + mean(post_survFESU$'beta[4]')*xdummy*1
                                   + mean(mean(post_survFESU$'tau_year[2,1]'),mean(post_survFESU$'tau_year[2,2]'),mean(post_survFESU$'tau_year[2,3]'),mean(post_survFESU$'tau_year[2,4]'),mean(post_survFESU$'tau_year[2,5]'),mean(post_survFESU$'tau_year[2,6]'),mean(post_survFESU$'tau_year[2,7]'),mean(post_survFESU$'tau_year[2,8]'),mean(post_survFESU$'tau_year[2,9]'),mean(post_survFESU$'tau_year[2,10]'),mean(post_survFESU$'tau_year[2,11]'))))
FESUfits <- as.data.frame(cbind(xdummy, ydummy_eminus, ydummy_eplus))
FESUfits <- melt(FESUfits, id = "xdummy", measure.vars = c("ydummy_eplus", "ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ydummy_eminus" = "E-", "ydummy_eplus" = "E+"))

surv_bin <- FESU_data1 %>% 
  select(surv_t1, logsize_t, endo) %>% 
  mutate(endo = recode(endo, "minus" = 'E-', "plus" = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T))
# with points
ggplot(data = FESUfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = surv_bin, aes(x = mean_size, y = mean_surv, col = endo), lwd = 3) +
  xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

ggplot(data = post_survFESU,) +
  geom_density( aes(`sigma_0[1]`), col = "#ff7f00", lwd = .8) +
  geom_density(aes(`sigma_0[2]`), col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")




# AGPE survival
nvalues <- length(AGPE_size_dat$logsize_t)
xdummy <- seq(min(AGPE_size_dat$logsize_t), max(AGPE_size_dat$logsize_t), length.out = nvalues)

# E- 
ydummy_eminus <- as.vector(invlogit(mean(post_survAGPE$alpha) + mean(post_survAGPE$'beta[1]')*xdummy + mean(post_survAGPE$'beta[2]')*0 + mean(post_survAGPE$'beta[3]')*mean(AGPE_origin_dat$origin1)
                                    + mean(post_survAGPE$'beta[4]')*xdummy*0
                                    + mean(mean(post_survAGPE$'tau_year[1,1]'),mean(post_survAGPE$'tau_year[1,2]'),mean(post_survAGPE$'tau_year[1,3]'),mean(post_survAGPE$'tau_year[1,4]'),mean(post_survAGPE$'tau_year[1,5]'),mean(post_survAGPE$'tau_year[1,6]'),mean(post_survAGPE$'tau_year[1,7]'),mean(post_survAGPE$'tau_year[1,8]'),mean(post_survAGPE$'tau_year[1,9]'),mean(post_survAGPE$'tau_year[1,10]'),mean(post_survAGPE$'tau_year[1,11]'))))

# E+
ydummy_eplus <- as.vector(invlogit(mean(post_survAGPE$alpha) + mean(post_survAGPE$'beta[1]')*xdummy + mean(post_survAGPE$'beta[2]')*1 + mean(post_survAGPE$'beta[3]')*mean(AGPE_origin_dat$origin1)
                                   + mean(post_survAGPE$'beta[4]')*xdummy*1
                                   + mean(mean(post_survAGPE$'tau_year[2,1]'),mean(post_survAGPE$'tau_year[2,2]'),mean(post_survAGPE$'tau_year[2,3]'),mean(post_survAGPE$'tau_year[2,4]'),mean(post_survAGPE$'tau_year[2,5]'),mean(post_survAGPE$'tau_year[2,6]'),mean(post_survAGPE$'tau_year[2,7]'),mean(post_survAGPE$'tau_year[2,8]'),mean(post_survAGPE$'tau_year[2,9]'),mean(post_survAGPE$'tau_year[2,10]'),mean(post_survAGPE$'tau_year[2,11]'))))
AGPEfits <- as.data.frame(cbind(xdummy, ydummy_eminus, ydummy_eplus))
AGPEfits <- melt(AGPEfits, id = "xdummy", measure.vars = c("ydummy_eplus", "ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ydummy_eminus" = "E-", "ydummy_eplus" = "E+"))

surv_bin <- AGPE_data1 %>% 
  select(surv_t1, logsize_t, endo) %>% 
  mutate(endo = recode(endo, "minus" = "E-", "0" = 'E-', "plus" = "E+", "1" = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T))
# with points
ggplot(data = AGPEfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = surv_bin, aes(x = mean_size, y = mean_surv, col = endo), lwd = 3) +
  xlab("log(size_t)") + ylab("Prob. of Survival") +   
  scale_color_manual(values=colors2)

ggplot(data = post_survAGPE,) +
  geom_density( aes(`sigma_0[1]`), col = "#ff7f00", lwd = .8) +
  geom_density(aes(`sigma_0[2]`), col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")

