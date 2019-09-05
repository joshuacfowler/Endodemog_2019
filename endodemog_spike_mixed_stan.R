## Title: Grass endophyte population model with a bayesian framework
## Purpose: Creates spikelet production kernel written in STAN, 
## and does visualisation of posterior predictive checks
## Authors: Joshua and Tom
#############################################################

library(tidyverse)
library(rstan)
library(StanHeaders)
library(shinystan)
library(bayesplot)
library(devtools)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }


#############################################################################################
####### Data manipulation to prepare data as lists for Stan models------------------
#############################################################################################

# spikelet data lists are generated in the endodemog_data_processing.R file, 
# within the section titled "Preparing datalists for Spikelet Kernel"
source("endodemog_data_processing.R")
#########################################################################################################
# Stan model for spikelet production ------------------------------
#########################################################################################################
## run this code to optimize computer system settings for MCMC
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## MCMC settings
ni <-10000
nb <- 5000
nc <- 3

