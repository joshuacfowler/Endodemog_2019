library(tidyverse)
library(rstan)
library(StanHeaders)
library(shinystan)
library(bayesplot)
library(devtools)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(scales)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

## Load Josh's model fits
survPOSY <- readRDS(file = "D:/Dropbox/EndodemogData/Model_Runs/endodemog_surv_POSY.rds")
survPOSY@stanmodel
params = c("beta", "tau_year", "sigma_e")
post_survPOSY <- as.data.frame(survPOSY, pars = params)


## Load in full data frame
LTREB_endodemog <- 
  read.csv(file = "D:/Dropbox/EndodemogData/Fulldataplusmetadata/endo_demog_long.csv")
## Josh's clean-up code
LTREB_data <- LTREB_endodemog %>% 
  mutate(size_t, logsize_t = log(size_t)) %>% 
  mutate(size_t1, logsize_t1 = log(size_t1)) %>%  
  mutate(surv_t1 = as.integer(recode(surv_t1, "0" = 0, "1" =1, "2" = 1, "4" = 1))) %>% 
  mutate(endo_01 = as.integer(case_when(endo == "0" | endo == "minus" ~ 0,
                                        endo == "1"| endo =="plus" ~ 1))) %>% 
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
# NA's in survival come from mostly 2017 recruits.
LTREB_data1 <- LTREB_data %>%
  filter(!is.na(surv_t1)) %>% 
  filter(!is.na(logsize_t)) %>% 
  filter(logsize_t >= 0) %>% 
  filter(!is.na(endo_01))

POSY_surv <- LTREB_data1 %>% 
  filter(species == "POSY") %>% 
  select(logsize_t, endo, year_t, surv_t1) %>% 
  mutate(size_bin = as.integer(as.factor(cut_interval(logsize_t, n = 10))))%>%
  mutate(endo = recode(endo, minus = 0, plus = 1)) %>% 
  group_by(size_bin, endo) %>% 
  summarise(mean_size = mean((logsize_t),na.rm=T),
            mean_surv = mean(surv_t1,na.rm=T),
            surv_n = n())

x.dummy<-seq(min(POSY_surv$mean_size),max(POSY_surv$mean_size),0.1)
post_survPOSY_samples <- post_survPOSY[sample.int(n = nrow(post_survPOSY),100),]

plot(POSY_surv$mean_size,POSY_surv$mean_surv,type="n",ylim=c(0,1))
for(i in 1:nrow(post_survPOSY_samples)){
  lines(x.dummy,
        invlogit((post_survPOSY_samples[i,]$`beta[1]`) + (post_survPOSY_samples[i,]$`beta[2]`) * x.dummy), col=alpha("black",0.8), lwd=3, lty=2)
  lines(x.dummy,
        invlogit((post_survPOSY_samples[i,]$`beta[1]`) + (post_survPOSY_samples[i,]$`beta[2]`) * x.dummy + 
                   (post_survPOSY_samples[i,]$`beta[3]`) + (post_survPOSY_samples[i,]$`beta[5]`) * x.dummy), col=alpha("black",0.2), lwd=3, lty=1)
}

lines(x.dummy,
      invlogit(mean(post_survPOSY$`beta[1]`) + mean(post_survPOSY$`beta[2]`) * x.dummy), lwd=1, lty=2)
lines(x.dummy,
      invlogit(mean(post_survPOSY$`beta[1]`) + mean(post_survPOSY$`beta[2]`) * x.dummy + 
                 mean(post_survPOSY$`beta[3]`) + mean(post_survPOSY$`beta[5]`) * x.dummy),lwd=1, lty=1)


points(POSY_surv$mean_size[POSY_surv$endo==0],POSY_surv$mean_surv[POSY_surv$endo==0],
       cex=(log(POSY_surv$surv_n)/(max(log(POSY_surv$surv_n))))*2)
points(POSY_surv$mean_size[POSY_surv$endo==1],POSY_surv$mean_surv[POSY_surv$endo==1],pch=16)
