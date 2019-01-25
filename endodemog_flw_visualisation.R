## Title: Grass endophyte vital rate variance with a bayesian framework
## Purpose: Visualizes output of endophyte effect on year variance from 
## model outputs for Flowering models 
## Authors: Joshua and Tom
#############################################################
library(tidyverse)
library(rstan)
library(StanHeaders)
library(shinystan)
library(bayesplot)
library(devtools)
library(reshape2)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

#  source in raw data from endodemog_flw_mixed_stan.R

## Read in the flower model output for all species
flwPOAL <- readRDS(file = "endodemog_flw_full_POAL.rds")
flwPOSY <- readRDS(file = "endodemog_flw_full_POSY.rds")
flwLOAR <- readRDS(file = "endodemog_flw_full_LOAR.rds")
flwELVI <- readRDS(file = "endodemog_flw_full_ELVI.rds")
flwELRI <- readRDS(file = "endodemog_flw_full_ELRI.rds")
flwFESU <- readRDS(file = "endodemog_flw_full_FESU.rds")
flwAGPE <- readRDS(file = "endodemog_flw_full_AGPE.rds")



params = c("alpha", "beta", "tau_year", "sigma_0")
# save posteriors within dataframe
post_flwPOAL <- as.data.frame(flwPOAL, pars = params)
post_flwPOSY <- as.data.frame(flwPOSY, pars = params)
post_flwLOAR <- as.data.frame(flwLOAR, pars = params)
post_flwELVI <- as.data.frame(flwELVI, pars = params)
post_flwELRI <- as.data.frame(flwELRI, pars = params)
post_flwFESU <- as.data.frame(flwFESU, pars = params)
post_flwAGPE <- as.data.frame(flwAGPE, pars = params)


yearcolors <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99",
                "#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#6a3d9a",
                "#ffff99")
colors2 <- c("#e31a1c","#ff7f00")

# POAL flowering
nvalues <- length(POAL_size_dat$logsize_t)
xdummy <- seq(min(POAL_size_dat$logsize_t), max(POAL_size_dat$logsize_t), length.out = nvalues)

# actual data points for probability of flowering
flwbin0 <- POAL_data1 %>% 
  select(seed_t1, logsize_t, year_t, endo) %>%
  mutate(flw_t1 = as.integer(case_when(seed_t1 == 0 ~ 0, seed_t1 > 0 ~ 1, is.na(seed_t1) ~ 0))) %>% 
  mutate(endo = recode(endo, minus = 0, plus = 1)) %>% 
  filter(endo == 0) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(flw_t1,na.rm=T))

flwbin1 <- POAL_data1 %>% 
  select(seed_t1, logsize_t, year_t, endo) %>%
  mutate(flw_t1 = as.integer(case_when(seed_t1 == 0 ~ 0, seed_t1 > 0 ~ 1, is.na(seed_t1) ~ 0))) %>% 
  mutate(endo = recode(endo, minus = 0, plus = 1)) %>% 
  filter(endo == 1) %>% 
  mutate(Year = as.factor(as.integer(as.factor(year_t)))) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 15)) %>% 
  group_by(size_bin, Year) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(flw_t1,na.rm=T))



# E- 
ydummy_eminus <- as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*0 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                    + mean(post_flwPOAL$'beta[4]')*xdummy*0
                                    + mean(mean(post_flwPOAL$'tau_year[1,1]'),mean(post_flwPOAL$'tau_year[1,2]'),mean(post_flwPOAL$'tau_year[1,3]'),mean(post_flwPOAL$'tau_year[1,4]'),mean(post_flwPOAL$'tau_year[1,5]'),mean(post_flwPOAL$'tau_year[1,6]'),mean(post_flwPOAL$'tau_year[1,7]'),mean(post_flwPOAL$'tau_year[1,8]'),mean(post_flwPOAL$'tau_year[1,9]'),mean(post_flwPOAL$'tau_year[1,10]'),mean(post_flwPOAL$'tau_year[1,11]'))))
ydummy_eminus_y1 <- as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*0 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*0
                                       + mean(post_flwPOAL$'tau_year[1,1]')))
ydummy_eminus_y2 <- as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*0 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*0
                                       + mean(post_flwPOAL$'tau_year[1,2]')))
ydummy_eminus_y3 <- as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*0 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*0
                                       + mean(post_flwPOAL$'tau_year[1,3]')))
ydummy_eminus_y4 <- as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*0 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*0
                                       + mean(post_flwPOAL$'tau_year[1,4]')))
ydummy_eminus_y5 <- as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*0 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*0
                                       + mean(post_flwPOAL$'tau_year[1,5]')))
ydummy_eminus_y6 <- as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*0 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*0
                                       + mean(post_flwPOAL$'tau_year[1,6]')))
ydummy_eminus_y7 <- as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*0 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*0
                                       + mean(post_flwPOAL$'tau_year[1,7]')))
ydummy_eminus_y8 <- as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*0 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*0
                                       + mean(post_flwPOAL$'tau_year[1,8]')))
ydummy_eminus_y9 <- as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*0 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*0
                                       + mean(post_flwPOAL$'tau_year[1,9]')))
ydummy_eminus_y10 <- as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*0 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                        + mean(post_flwPOAL$'beta[4]')*xdummy*0
                                        + mean(post_flwPOAL$'tau_year[1,10]')))
ydummy_eminus_y11 <- as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*0 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                        + mean(post_flwPOAL$'beta[4]')*xdummy*0
                                        + mean(post_flwPOAL$'tau_year[1,11]')))
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
  geom_point(data = flwbin0, aes(x = mean_size, y = mean_surv, color = Year)) +
  labs(title = "POAL E- Flowering Probability", x = "log(size_t)", y = "Prob. of Flowering") +
  scale_color_manual(values=yearcolors)
# without datapoints
ggplot(data = POALfitsminus) +
  geom_line(aes(x = xdummy, y = value, color = Year), lwd = .7) +
  labs(title = "POAL E- Flowering Probability", x = "log(size_t)", y = "Prob. of Flowering") +
  scale_color_manual(values=yearcolors)

# endo effect on mean
ggplot(data = POALfitsminus0) +
  geom_line(aes(x = xdummy, y = ydummy_eminus)) +
  labs(title = "POAL E- Flowering Probability", x = "log(size_t)", y = "Prob. of Flowering") +
  scale_color_manual(values=colors2)

# E+
ydummy_eplus <- as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*1 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1)
                                   + mean(post_flwPOAL$'beta[4]')*xdummy*1
                                   + mean(mean(post_flwPOAL$'tau_year[2,1]'),mean(post_flwPOAL$'tau_year[2,2]'),mean(post_flwPOAL$'tau_year[2,3]'),mean(post_flwPOAL$'tau_year[2,4]'),mean(post_flwPOAL$'tau_year[2,5]'),mean(post_flwPOAL$'tau_year[2,6]'),mean(post_flwPOAL$'tau_year[2,7]'),mean(post_flwPOAL$'tau_year[2,8]'),mean(post_flwPOAL$'tau_year[2,9]'),mean(post_flwPOAL$'tau_year[2,10]'),mean(post_flwPOAL$'tau_year[2,11]'))))
ydummy_eplus_y1 <-  as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*1 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*1
                                       + mean(post_flwPOAL$'tau_year[2,1]')))
ydummy_eplus_y2 <-  as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*1 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*1
                                       + mean(post_flwPOAL$'tau_year[2,2')))
ydummy_eplus_y3 <-  as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*1 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*1
                                       + mean(post_flwPOAL$'tau_year[2,3]')))
ydummy_eplus_y4 <-  as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*1 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*1
                                       + mean(post_flwPOAL$'tau_year[2,4]')))
ydummy_eplus_y5 <-  as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*1 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*1
                                       + mean(post_flwPOAL$'tau_year[2,5]')))
ydummy_eplus_y6 <-  as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*1 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*1
                                       + mean(post_flwPOAL$'tau_year[2,6]')))
ydummy_eplus_y7 <-  as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*1 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*1
                                       + mean(post_flwPOAL$'tau_year[2,7]')))
ydummy_eplus_y8 <-  as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*1 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*1
                                       + mean(post_flwPOAL$'tau_year[2,8]')))
ydummy_eplus_y9 <-  as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*1 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                                       + mean(post_flwPOAL$'beta[4]')*xdummy*1
                                       + mean(post_flwPOAL$'tau_year[2,9]')))
ydummy_eplus_y10 <-  as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*1 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                                        + mean(post_flwPOAL$'beta[4]')*xdummy*1
                                        + mean(post_flwPOAL$'tau_year[2,10]')))
ydummy_eplus_y11 <-  as.vector(invlogit(mean(post_flwPOAL$alpha) + mean(post_flwPOAL$'beta[1]')*xdummy + mean(post_flwPOAL$'beta[2]')*1 + mean(post_flwPOAL$'beta[3]')*mean(POAL_origin_dat$origin1) 
                                        + mean(post_flwPOAL$'beta[4]')*xdummy*1
                                        + mean(post_flwPOAL$'tau_year[2,11]')))
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
  geom_point(data = flwbin1, aes(x = mean_size, y = mean_surv, color = Year), lwd = 2) +
  ggtitle("POAL E+ Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +   
  scale_color_manual(values=yearcolors)
ggplot(data = POALfitsplus) +
  geom_line(aes(x = xdummy, y = value, color = Year),  lwd = .8) +
  ggtitle("POAL E+ Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +   
  scale_color_manual(values=yearcolors)



# effect on mean


POALfits <- as.data.frame(cbind(xdummy, ydummy_eminus, ydummy_eplus))
POALfits <- melt(POALfits, id = "xdummy", measure.vars = c("ydummy_eplus", "ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ydummy_eminus" = "E-", "ydummy_eplus" = "E+"))

flwbin <- POAL_data1 %>% 
  select(seed_t1, logsize_t, endo) %>%
  mutate(flw_t1 = as.integer(case_when(seed_t1 == 0 ~ 0, seed_t1 > 0 ~ 1, is.na(seed_t1) ~0))) %>%  
  mutate(endo1 = recode(endo, minus = 'E-', plus = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo1) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_flw = mean(flw_t1,na.rm=T))
# without point
ggplot(data = POALfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  ggtitle("POAL Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +   
  scale_color_manual(values=colors2)
ggplot(data = flwbin)+
  geom_point(aes(x= mean_size, y = mean_flw, col = endo1))

# with points
ggplot(data = POALfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = flwbin, aes(x = mean_size, y = mean_flw, col = endo1), lwd = 3) +
  ggtitle("POAL Flowering Probability") + xlab("log(size_t)") + ylab("Prob. of Flowering") +   
  scale_color_manual(values=colors2)

ggplot(data = POALfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = flwbin, aes(x = mean_size, y = mean_flw, col = endo1), lwd = 3) +
  xlab("log(size_t)") + ylab("Prob. of Flowering") +   
  scale_color_manual(values=colors2)


# effect on variance as density plot

stan_dens(flwPOAL, pars = c("sigma_0[1]", "sigma_0[2]"))

mcmc_areas(post_flwPOAL, pars = c("sigma_0[1]", "sigma_0[2]"),
           prob = 0.8, prob_outer = .99 )


ggplot(data = post_flwPOAL,) +
  geom_density( aes(`sigma_0[1]`), col = "#ff7f00", lwd = .8) +
  geom_density(aes(`sigma_0[2]`), col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")

# POSY#######################

# POSY Flowering
nvalues <- length(POSY_size_dat$logsize_t)
xdummy <- seq(min(POSY_size_dat$logsize_t), max(POSY_size_dat$logsize_t), length.out = nvalues)

# E- 
ydummy_eminus <- as.vector(invlogit(mean(post_flwPOSY$alpha) + mean(post_flwPOSY$'beta[1]')*xdummy + mean(post_flwPOSY$'beta[2]')*0 + mean(post_flwPOSY$'beta[3]')*mean(POSY_origin_dat$origin1)
                                    + mean(post_flwPOSY$'beta[4]')*xdummy*0
                                    + mean(mean(post_flwPOSY$'tau_year[1,1]'),mean(post_flwPOSY$'tau_year[1,2]'),mean(post_flwPOSY$'tau_year[1,3]'),mean(post_flwPOSY$'tau_year[1,4]'),mean(post_flwPOSY$'tau_year[1,5]'),mean(post_flwPOSY$'tau_year[1,6]'),mean(post_flwPOSY$'tau_year[1,7]'),mean(post_flwPOSY$'tau_year[1,8]'),mean(post_flwPOSY$'tau_year[1,9]'),mean(post_flwPOSY$'tau_year[1,10]'),mean(post_flwPOSY$'tau_year[1,11]'))))

# E+
ydummy_eplus <- as.vector(invlogit(mean(post_flwPOSY$alpha) + mean(post_flwPOSY$'beta[1]')*xdummy + mean(post_flwPOSY$'beta[2]')*1+ mean(post_flwPOSY$'beta[3]')*mean(POSY_origin_dat$origin1)
                                   + mean(post_flwPOSY$'beta[4]')*xdummy*1
                                   + mean(mean(post_flwPOSY$'tau_year[2,1]'),mean(post_flwPOSY$'tau_year[2,2]'),mean(post_flwPOSY$'tau_year[2,3]'),mean(post_flwPOSY$'tau_year[2,4]'),mean(post_flwPOSY$'tau_year[2,5]'),mean(post_flwPOSY$'tau_year[2,6]'),mean(post_flwPOSY$'tau_year[2,7]'),mean(post_flwPOSY$'tau_year[2,8]'),mean(post_flwPOSY$'tau_year[2,9]'),mean(post_flwPOSY$'tau_year[2,10]'),mean(post_flwPOSY$'tau_year[2,11]'))))
POSYfits <- as.data.frame(cbind(xdummy, ydummy_eminus, ydummy_eplus))
POSYfits <- melt(POSYfits, id = "xdummy", measure.vars = c("ydummy_eplus", "ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ydummy_eminus" = "E-", "ydummy_eplus" = "E+"))

flwbin <- POSY_data1 %>% 
  select(seed_t1, logsize_t, endo) %>% 
  mutate(flw_t1 = as.integer(case_when(seed_t1 == 0 ~ 0, seed_t1 > 0 ~ 1, is.na(seed_t1) ~ 0))) %>% 
  mutate(endo = recode(endo, "minus" = 'E-', "plus" = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(flw_t1,na.rm=T))
# with points
ggplot(data = POSYfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = flwbin, aes(x = mean_size, y = mean_surv, col = endo), lwd = 3) +
  xlab("log(size_t)") + ylab("Prob. of Flowering") +   
  scale_color_manual(values=colors2)

ggplot(data = post_flwPOSY,) +
  geom_density( aes(`sigma_0[1]`), col = "#ff7f00", lwd = .8) +
  geom_density(aes(`sigma_0[2]`), col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")



########################

# LOAR Flowering
nvalues <- length(LOAR_size_dat$logsize_t)
xdummy <- seq(min(LOAR_size_dat$logsize_t), max(LOAR_size_dat$logsize_t), length.out = nvalues)

# E- 
ydummy_eminus <- as.vector(invlogit(mean(post_flwLOAR$alpha) + mean(post_flwLOAR$'beta[1]')*xdummy + mean(post_flwLOAR$'beta[2]')*0 + mean(post_flwLOAR$'beta[3]')*mean(LOAR_origin_dat$origin)
                                    + mean(post_flwLOAR$'beta[4]')*xdummy*0
                                    + mean(mean(post_flwLOAR$'tau_year[1,1]'),mean(post_flwLOAR$'tau_year[1,2]'),mean(post_flwLOAR$'tau_year[1,3]'),mean(post_flwLOAR$'tau_year[1,4]'),mean(post_flwLOAR$'tau_year[1,5]'),mean(post_flwLOAR$'tau_year[1,6]'),mean(post_flwLOAR$'tau_year[1,7]'),mean(post_flwLOAR$'tau_year[1,8]'),mean(post_flwLOAR$'tau_year[1,9]'),mean(post_flwLOAR$'tau_year[1,10]'),mean(post_flwLOAR$'tau_year[1,11]'))))

# E+
ydummy_eplus <- as.vector(invlogit(mean(post_flwLOAR$alpha) + mean(post_flwLOAR$'beta[1]')*xdummy + mean(post_flwLOAR$'beta[2]')*1 + mean(post_flwLOAR$'beta[3]')*mean(LOAR_origin_dat$origin)
                                   + mean(post_flwLOAR$'beta[4]')*xdummy*1
                                   + mean(mean(post_flwLOAR$'tau_year[2,1]'),mean(post_flwLOAR$'tau_year[2,2]'),mean(post_flwLOAR$'tau_year[2,3]'),mean(post_flwLOAR$'tau_year[2,4]'),mean(post_flwLOAR$'tau_year[2,5]'),mean(post_flwLOAR$'tau_year[2,6]'),mean(post_flwLOAR$'tau_year[2,7]'),mean(post_flwLOAR$'tau_year[2,8]'),mean(post_flwLOAR$'tau_year[2,9]'),mean(post_flwLOAR$'tau_year[2,10]'),mean(post_flwLOAR$'tau_year[2,11]'))))
LOARfits <- as.data.frame(cbind(xdummy, ydummy_eminus, ydummy_eplus))
LOARfits <- melt(LOARfits, id = "xdummy", measure.vars = c("ydummy_eplus", "ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ydummy_eminus" = "E-", "ydummy_eplus" = "E+"))

flwbin <- LOAR_data1 %>% 
  select(seed_t1, logsize_t, endo) %>%
  mutate(flw_t1 = as.integer(case_when(seed_t1 == 0 ~ 0, seed_t1 > 0 ~ 1, is.na(seed_t1) ~ 0))) %>% 
  mutate(endo = recode(endo, "0" = 'E-', "1" = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(flw_t1,na.rm=T))
# with points
ggplot(data = LOARfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = flwbin, aes(x = mean_size, y = mean_surv, col = endo), lwd = 3) +
  xlab("log(size_t)") + ylab("Prob. of Flowering") +   
  scale_color_manual(values=colors2)

ggplot(data = post_flwLOAR,) +
  geom_density( aes(`sigma_0[1]`), col = "#ff7f00", lwd = .8) +
  geom_density(aes(`sigma_0[2]`), col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")



# ELVI Flowering
nvalues <- length(ELVI_size_dat$logsize_t)
xdummy <- seq(min(ELVI_size_dat$logsize_t), max(ELVI_size_dat$logsize_t), length.out = nvalues)

# E- 
ydummy_eminus <- as.vector(invlogit(mean(post_flwELVI$alpha) + mean(post_flwELVI$'beta[1]')*xdummy + mean(post_flwELVI$'beta[2]')*0 + mean(post_flwELVI$'beta[3]')*mean(ELVI_origin_dat$origin1)
                                    + mean(post_flwELVI$'beta[4]')*xdummy*0
                                    + mean(mean(post_flwELVI$'tau_year[1,1]'),mean(post_flwELVI$'tau_year[1,2]'),mean(post_flwELVI$'tau_year[1,3]'),mean(post_flwELVI$'tau_year[1,4]'),mean(post_flwELVI$'tau_year[1,5]'),mean(post_flwELVI$'tau_year[1,6]'),mean(post_flwELVI$'tau_year[1,7]'),mean(post_flwELVI$'tau_year[1,8]'),mean(post_flwELVI$'tau_year[1,9]'),mean(post_flwELVI$'tau_year[1,10]'),mean(post_flwELVI$'tau_year[1,11]'))))

# E+
ydummy_eplus <- as.vector(invlogit(mean(post_flwELVI$alpha) + mean(post_flwELVI$'beta[1]')*xdummy + mean(post_flwELVI$'beta[2]')*1 + mean(post_flwELVI$'beta[3]')*mean(ELVI_origin_dat$origin1)
                                   + mean(post_flwELVI$'beta[4]')*xdummy*1
                                   + mean(mean(post_flwELVI$'tau_year[2,1]'),mean(post_flwELVI$'tau_year[2,2]'),mean(post_flwELVI$'tau_year[2,3]'),mean(post_flwELVI$'tau_year[2,4]'),mean(post_flwELVI$'tau_year[2,5]'),mean(post_flwELVI$'tau_year[2,6]'),mean(post_flwELVI$'tau_year[2,7]'),mean(post_flwELVI$'tau_year[2,8]'),mean(post_flwELVI$'tau_year[2,9]'),mean(post_flwELVI$'tau_year[2,10]'),mean(post_flwELVI$'tau_year[2,11]'))))
ELVIfits <- as.data.frame(cbind(xdummy, ydummy_eminus, ydummy_eplus))
ELVIfits <- melt(ELVIfits, id = "xdummy", measure.vars = c("ydummy_eplus", "ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ydummy_eminus" = "E-", "ydummy_eplus" = "E+"))

flwbin <- ELVI_data1 %>% 
  select(seed_t1, logsize_t, endo) %>%
  mutate(flw_t1 = as.integer(case_when(seed_t1 == 0 ~ 0, seed_t1 > 0 ~ 1, is.na(seed_t1) ~ 0))) %>% 
  mutate(endo = recode(endo, "minus" = 'E-', "plus" = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(flw_t1,na.rm=T))
# with points
ggplot(data = ELVIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = flwbin, aes(x = mean_size, y = mean_surv, col = endo), lwd = 3) +
  xlab("log(size_t)") + ylab("Prob. of Flowering") +   
  scale_color_manual(values=colors2)

ggplot(data = post_flwELVI,) +
  geom_density( aes(`sigma_0[1]`), col = "#ff7f00", lwd = .8) +
  geom_density(aes(`sigma_0[2]`), col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")




# ELRI Flowering
nvalues <- length(ELRI_size_dat$logsize_t)
xdummy <- seq(min(ELRI_size_dat$logsize_t), max(ELRI_size_dat$logsize_t), length.out = nvalues)

# E- 
ydummy_eminus <- as.vector(invlogit(mean(post_flwELRI$alpha) + mean(post_flwELRI$'beta[1]')*xdummy + mean(post_flwELRI$'beta[2]')*0 + mean(post_flwELRI$'beta[3]')*mean(ELRI_origin_dat$origin1)
                                    + mean(post_flwELRI$'beta[4]')*xdummy*0
                                    + mean(mean(post_flwELRI$'tau_year[1,1]'),mean(post_flwELRI$'tau_year[1,2]'),mean(post_flwELRI$'tau_year[1,3]'),mean(post_flwELRI$'tau_year[1,4]'),mean(post_flwELRI$'tau_year[1,5]'),mean(post_flwELRI$'tau_year[1,6]'),mean(post_flwELRI$'tau_year[1,7]'),mean(post_flwELRI$'tau_year[1,8]'),mean(post_flwELRI$'tau_year[1,9]'),mean(post_flwELRI$'tau_year[1,10]'),mean(post_flwELRI$'tau_year[1,11]'))))

# E+
ydummy_eplus <- as.vector(invlogit(mean(post_flwELRI$alpha) + mean(post_flwELRI$'beta[1]')*xdummy + mean(post_flwELRI$'beta[2]')*1 + mean(post_flwELRI$'beta[3]')*mean(ELRI_origin_dat$origin1)
                                   + mean(post_flwELRI$'beta[4]')*xdummy*1
                                   + mean(mean(post_flwELRI$'tau_year[2,1]'),mean(post_flwELRI$'tau_year[2,2]'),mean(post_flwELRI$'tau_year[2,3]'),mean(post_flwELRI$'tau_year[2,4]'),mean(post_flwELRI$'tau_year[2,5]'),mean(post_flwELRI$'tau_year[2,6]'),mean(post_flwELRI$'tau_year[2,7]'),mean(post_flwELRI$'tau_year[2,8]'),mean(post_flwELRI$'tau_year[2,9]'),mean(post_flwELRI$'tau_year[2,10]'),mean(post_flwELRI$'tau_year[2,11]'))))
ELRIfits <- as.data.frame(cbind(xdummy, ydummy_eminus, ydummy_eplus))
ELRIfits <- melt(ELRIfits, id = "xdummy", measure.vars = c("ydummy_eplus", "ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ydummy_eminus" = "E-", "ydummy_eplus" = "E+"))

flwbin <- ELRI_data1 %>% 
  select(seed_t1, logsize_t, endo) %>%
  mutate(flw_t1 = as.integer(case_when(seed_t1 == 0 ~ 0, seed_t1 > 0 ~ 1, is.na(seed_t1) ~ 0))) %>% 
  mutate(endo = recode(endo, "minus" = 'E-', "plus" = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(flw_t1,na.rm=T))
# with points
ggplot(data = ELRIfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = flwbin, aes(x = mean_size, y = mean_surv, col = endo), lwd = 3) +
  xlab("log(size_t)") + ylab("Prob. of Flowering") +   
  scale_color_manual(values=colors2)

ggplot(data = post_flwELRI,) +
  geom_density( aes(`sigma_0[1]`), col = "#ff7f00", lwd = .8) +
  geom_density(aes(`sigma_0[2]`), col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")


# FESU Flowering
nvalues <- length(FESU_size_dat$logsize_t)
xdummy <- seq(min(FESU_size_dat$logsize_t), max(FESU_size_dat$logsize_t), length.out = nvalues)

# E- 
ydummy_eminus <- as.vector(invlogit(mean(post_flwFESU$alpha) + mean(post_flwFESU$'beta[1]')*xdummy + mean(post_flwFESU$'beta[2]')*0 + mean(post_flwFESU$'beta[3]')*mean(FESU_origin_dat$origin1)
                                    + mean(post_flwFESU$'beta[4]')*xdummy*0
                                    + mean(mean(post_flwFESU$'tau_year[1,1]'),mean(post_flwFESU$'tau_year[1,2]'),mean(post_flwFESU$'tau_year[1,3]'),mean(post_flwFESU$'tau_year[1,4]'),mean(post_flwFESU$'tau_year[1,5]'),mean(post_flwFESU$'tau_year[1,6]'),mean(post_flwFESU$'tau_year[1,7]'),mean(post_flwFESU$'tau_year[1,8]'),mean(post_flwFESU$'tau_year[1,9]'),mean(post_flwFESU$'tau_year[1,10]'),mean(post_flwFESU$'tau_year[1,11]'))))

# E+
ydummy_eplus <- as.vector(invlogit(mean(post_flwFESU$alpha) + mean(post_flwFESU$'beta[1]')*xdummy + mean(post_flwFESU$'beta[2]')*1 + mean(post_flwFESU$'beta[3]')*mean(FESU_origin_dat$origin1)
                                   + mean(post_flwFESU$'beta[4]')*xdummy*1
                                   + mean(mean(post_flwFESU$'tau_year[2,1]'),mean(post_flwFESU$'tau_year[2,2]'),mean(post_flwFESU$'tau_year[2,3]'),mean(post_flwFESU$'tau_year[2,4]'),mean(post_flwFESU$'tau_year[2,5]'),mean(post_flwFESU$'tau_year[2,6]'),mean(post_flwFESU$'tau_year[2,7]'),mean(post_flwFESU$'tau_year[2,8]'),mean(post_flwFESU$'tau_year[2,9]'),mean(post_flwFESU$'tau_year[2,10]'),mean(post_flwFESU$'tau_year[2,11]'))))
FESUfits <- as.data.frame(cbind(xdummy, ydummy_eminus, ydummy_eplus))
FESUfits <- melt(FESUfits, id = "xdummy", measure.vars = c("ydummy_eplus", "ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ydummy_eminus" = "E-", "ydummy_eplus" = "E+"))

flwbin <- FESU_data1 %>% 
  select(seed_t1, logsize_t, endo) %>%
  mutate(flw_t1 = as.integer(case_when(seed_t1 == 0 ~ 0, seed_t1 > 0 ~ 1, is.na(seed_t1) ~ 0))) %>% 
  mutate(endo = recode(endo, "minus" = 'E-', "plus" = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(flw_t1,na.rm=T))
# with points
ggplot(data = FESUfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = flwbin, aes(x = mean_size, y = mean_surv, col = endo), lwd = 3) +
  xlab("log(size_t)") + ylab("Prob. of Flowering") +   
  scale_color_manual(values=colors2)

ggplot(data = post_flwFESU,) +
  geom_density( aes(`sigma_0[1]`), col = "#ff7f00", lwd = .8) +
  geom_density(aes(`sigma_0[2]`), col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")




# AGPE Flowering
nvalues <- length(AGPE_size_dat$logsize_t)
xdummy <- seq(min(AGPE_size_dat$logsize_t), max(AGPE_size_dat$logsize_t), length.out = nvalues)

# E- 
ydummy_eminus <- as.vector(invlogit(mean(post_flwAGPE$alpha) + mean(post_flwAGPE$'beta[1]')*xdummy + mean(post_flwAGPE$'beta[2]')*0 + mean(post_flwAGPE$'beta[3]')*mean(AGPE_origin_dat$origin1)
                                    + mean(post_flwAGPE$'beta[4]')*xdummy*0
                                    + mean(mean(post_flwAGPE$'tau_year[1,1]'),mean(post_flwAGPE$'tau_year[1,2]'),mean(post_flwAGPE$'tau_year[1,3]'),mean(post_flwAGPE$'tau_year[1,4]'),mean(post_flwAGPE$'tau_year[1,5]'),mean(post_flwAGPE$'tau_year[1,6]'),mean(post_flwAGPE$'tau_year[1,7]'),mean(post_flwAGPE$'tau_year[1,8]'),mean(post_flwAGPE$'tau_year[1,9]'),mean(post_flwAGPE$'tau_year[1,10]'),mean(post_flwAGPE$'tau_year[1,11]'))))

# E+
ydummy_eplus <- as.vector(invlogit(mean(post_flwAGPE$alpha) + mean(post_flwAGPE$'beta[1]')*xdummy + mean(post_flwAGPE$'beta[2]')*1 + mean(post_flwAGPE$'beta[3]')*mean(AGPE_origin_dat$origin1)
                                   + mean(post_flwAGPE$'beta[4]')*xdummy*1
                                   + mean(mean(post_flwAGPE$'tau_year[2,1]'),mean(post_flwAGPE$'tau_year[2,2]'),mean(post_flwAGPE$'tau_year[2,3]'),mean(post_flwAGPE$'tau_year[2,4]'),mean(post_flwAGPE$'tau_year[2,5]'),mean(post_flwAGPE$'tau_year[2,6]'),mean(post_flwAGPE$'tau_year[2,7]'),mean(post_flwAGPE$'tau_year[2,8]'),mean(post_flwAGPE$'tau_year[2,9]'),mean(post_flwAGPE$'tau_year[2,10]'),mean(post_flwAGPE$'tau_year[2,11]'))))
AGPEfits <- as.data.frame(cbind(xdummy, ydummy_eminus, ydummy_eplus))
AGPEfits <- melt(AGPEfits, id = "xdummy", measure.vars = c("ydummy_eplus", "ydummy_eminus")) %>% 
  mutate(Endo = recode(variable, "ydummy_eminus" = "E-", "ydummy_eplus" = "E+"))

flwbin <- AGPE_data1 %>% 
  select(seed_t1, logsize_t, endo) %>% 
  mutate(flw_t1 = as.integer(case_when(seed_t1 == 0 ~ 0, seed_t1 > 0 ~ 1, is.na(seed_t1) ~ 0))) %>% 
  mutate(endo = recode(endo, "minus" = "E-", "0" = 'E-', "plus" = "E+", "1" = 'E+')) %>% 
  mutate(size_bin = cut_interval(logsize_t, n = 12)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(flw_t1,na.rm=T))
# with points
ggplot(data = AGPEfits) +
  geom_line(aes(x = xdummy, y = value, col = Endo),lwd = .8) +
  geom_point(data = flwbin, aes(x = mean_size, y = mean_surv, col = endo), lwd = 3) +
  xlab("log(size_t)") + ylab("Prob. of Flowering") +   
  scale_color_manual(values=colors2)

ggplot(data = post_flwAGPE,) +
  geom_density( aes(`sigma_0[1]`), col = "#ff7f00", lwd = .8) +
  geom_density(aes(`sigma_0[2]`), col = "#e31a1c", lwd = .8) +
  labs(x = expression(σ[year]), y = "Posterior Density")



