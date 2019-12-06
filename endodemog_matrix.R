## Title: Grass endophyte population model with a bayesian framework
## Purpose: Assemble matrix model from vital rate estimates
## Authors: Joshua and Tom
## Updated: 11/19/2019
#############################################################
library(tidyverse)
library(rstan)
library(StanHeaders)
library(popbio)
library(bbmle)
library(SPEI)

invlogit<-function(x){exp(x)/(1+exp(x))}
logit = function(x) { log(x/(1-x)) }

#############################################################################################
####### Data manipulation to prepare data as lists for Stan models------------------
#############################################################################################
# filepaths
joshpath <- "/Users/joshuacfowler/Dropbox/EndodemogData/"
tompath <- "C:/Users/tm9/Dropbox/EndodemogData/"
# Growth data lists are generated in the endodemog_data_processing.R file
# within the section titled "Preparing datalists for Growth Kernel"
source("endodemog_data_processing.R")
# endo_demog <- read.csv("~/Dropbox/EndodemogData/fulldataplusmetadata/endo_demog_long.csv")
# loar <- endo_demog %>% 
#   filter(species == "LOAR") %>% 
#   mutate(log_size_t = log(size_t),
#          log_size_t1 = log(size_t1),
#          flow_t = as.integer(seed_t > 0))
# 
# # DATA ISSUES:
# # I noticed that there are about 10 plants with a size of 0 tillers. They are:
# filter(loar,size_t==0)
# # I also noticed that there are a few instances where survival was > 1:
# filter(loar,surv_t1 > 1)
# # We'll need to cut these for the analysis.
# loar <- loar %>% filter(size_t > 0,
#                         surv_t1 <= 1)

# 
# #  Read in the model outputs ----------------------------------------------
source("endodemog_model_outputs.R")

# survLOAR <- read_rds(path = paste(joshpath, "Model_Runs/endodemog_surv_LOAR.rds", sep = ""))
# 
# growLOAR <- read_rds(path = paste(joshpath, "Model_Runs/endodemog_grow_LOAR.rds", sep = ""))
# flwLOAR <- read_rds(path =  paste(joshpath, "Model_Runs/endodemog_flw_LOAR.rds", sep = ""))
# 
# fertLOAR <- read_rds(path =  paste(joshpath, "Model_Runs/endodemog_fert_LOAR.rds", sep = ""))
# 
# spikeLOAR <- read_rds(path = paste(joshpath, "Model_Runs/endodemog_spike_LOAR.rds", sep = ""))
# seedLOAR <- read_rds(path = paste(joshpath, "Model_Runs/endodemog_seed_mean_LOAR.rds", sep = ""))
# s_to_sLOAR <- read_rds(path = paste(joshpath, "Model_Runs/endodemog_s_to_s_LOAR.rds", sep = ""))

#############################################################################################
# Assembling matrix model -------------------------------------------------
#############################################################################################

# collect all the parameters into one vector

# This is the linear predictor for surv, grow, flw, fert, spike
# mu[n] = beta[1] + beta[2]*logsize_t[n] + beta[3]*endo_01[n] +beta[4]*origin_01[n]
# + tau_year[endo_index[n],year_t[n]]
# + tau_plot[plot[n]];

getparams <- function(surv, grow, flw, fert, spike, seed, s_to_s, data){
params <- c()
params[1] <- lapply(rstan::extract(surv, pars = "beta[1]"), FUN = mean); names(params)[1] <- "surv_beta1"
params[2] <- lapply(rstan::extract(surv, pars = "beta[2]"), FUN = mean); names(params)[2] <- "surv_beta2"
params[3] <- lapply(rstan::extract(surv, pars = "beta[3]"), FUN = mean); names(params)[3] <- "surv_beta3"
params[4] <- lapply(rstan::extract(surv, pars = "beta[4]"), FUN = mean); names(params)[4] <- "surv_beta4"
# collect the growth parameters
params[5] <- lapply(rstan::extract(grow, pars = "beta[1]"), FUN = mean); names(params)[5] <- "grow_beta1"
params[6] <- lapply(rstan::extract(grow, pars = "beta[2]"), FUN = mean); names(params)[6] <- "grow_beta2"
params[7] <- lapply(rstan::extract(grow, pars = "beta[3]"), FUN = mean); names(params)[7] <- "grow_beta3"
params[8] <- lapply(rstan::extract(grow, pars = "beta[4]"), FUN = mean); names(params)[8] <- "grow_beta4"
params[9] <- lapply(rstan::extract(grow, pars = "phi"), FUN = mean); names(params)[9] <- "grow_phi"
# collect the flowering parameters
params[10] <- lapply(rstan::extract(flw, pars = "beta[1]"), FUN = mean); names(params)[10] <- "flw_beta1"
params[11] <- lapply(rstan::extract(flw, pars = "beta[2]"), FUN = mean); names(params)[11] <- "flw_beta2"
params[12] <- lapply(rstan::extract(flw, pars = "beta[3]"), FUN = mean); names(params)[12] <- "flw_beta3"
params[13] <- lapply(rstan::extract(flw, pars = "beta[4]"), FUN = mean); names(params)[13] <- "flw_beta4"
# collect the fertility parameters
params[14] <- lapply(rstan::extract(fert, pars = "beta[1]"), FUN = mean); names(params)[14] <- "fert_beta1"
params[15] <- lapply(rstan::extract(fert, pars = "beta[2]"), FUN = mean); names(params)[15] <- "fert_beta2"
params[16] <- lapply(rstan::extract(fert, pars = "beta[3]"), FUN = mean); names(params)[16] <- "fert_beta3"
params[17] <- lapply(rstan::extract(fert, pars = "beta[4]"), FUN = mean); names(params)[17] <- "fert_beta4"
params[18] <- lapply(rstan::extract(fert, pars = "phi"), FUN = mean); names(params)[18] <- "fert_phi"
# collect the spikelets/infl parameters
params[19] <- lapply(rstan::extract(spike, pars = "beta[1]"), FUN = mean); names(params)[19] <- "spike_beta1"
params[20] <- lapply(rstan::extract(spike, pars = "beta[2]"), FUN = mean); names(params)[20] <- "spike_beta2"
params[21] <- lapply(rstan::extract(spike, pars = "beta[3]"), FUN = mean); names(params)[21] <- "spike_beta3"
params[22] <- lapply(rstan::extract(spike, pars = "beta[4]"), FUN = mean); names(params)[22] <- "spike_beta4"
params[23] <- lapply(rstan::extract(spike, pars = "phi"), FUN = mean); names(params)[23] <- "spike_phi"
# collect the seed/spikelet parameters
params[24] <- lapply(rstan::extract(seed, pars = "mu_seed"), FUN = mean); names(params)[24] <- "mu_seed"
# collect the recruitment parameters
params[25] <- lapply(rstan::extract(s_to_s, pars = "beta[1]"), FUN = mean); names(params)[25] <- "s_to_s_beta1"
params[26] <- lapply(rstan::extract(s_to_s, pars = "beta[2]"), FUN = mean); names(params)[26] <- "s_to_s_beta2"
# min and max size
params[27] <- 1; names(params)[27]<-"min_size"
params[28] <- quantile(data$size_t1,0.95,na.rm=T); names(params)[28]<-"max_size"
# note that I define max size as the 95TH pctile of observed sized. The very max sizes observed often have very poor 
# replication, and I find that these few observations (and the corresponding vital rate predictions) can have a strong
# influence on the results. So this approach is conservative, but you can always experiment with this.
# collect the endo specific year parameters
# surv
params[29] <- lapply(rstan::extract(surv, pars = "tau_year[1,1]"), FUN = mean); names(params)[29] <- "surv_eminus_y1"
params[30] <- lapply(rstan::extract(surv, pars = "tau_year[1,2]"), FUN = mean); names(params)[30] <- "surv_eminus_y2"
params[31] <- lapply(rstan::extract(surv, pars = "tau_year[1,3]"), FUN = mean); names(params)[31] <- "surv_eminus_y3"
params[32] <- lapply(rstan::extract(surv, pars = "tau_year[1,4]"), FUN = mean); names(params)[32] <- "surv_eminus_y4"
params[33] <- lapply(rstan::extract(surv, pars = "tau_year[1,5]"), FUN = mean); names(params)[33] <- "surv_eminus_y5"
params[34] <- lapply(rstan::extract(surv, pars = "tau_year[1,6]"), FUN = mean); names(params)[34] <- "surv_eminus_y6"
params[35] <- lapply(rstan::extract(surv, pars = "tau_year[1,7]"), FUN = mean); names(params)[35] <- "surv_eminus_y7"
params[36] <- lapply(rstan::extract(surv, pars = "tau_year[1,8]"), FUN = mean); names(params)[36] <- "surv_eminus_y8"
params[37] <- lapply(rstan::extract(surv, pars = "tau_year[1,9]"), FUN = mean); names(params)[37] <- "surv_eminus_y9"
params[38] <- lapply(rstan::extract(surv, pars = "tau_year[1,10]"), FUN = mean); names(params)[38] <- "surv_eminus_y10"
params[39] <- lapply(rstan::extract(surv, pars = "tau_year[1,11]"), FUN = mean); names(params)[39] <- "surv_eminus_y11"

params[40] <- lapply(rstan::extract(surv, pars = "tau_year[2,1]"), FUN = mean); names(params)[40] <- "surv_eplus_y1"
params[41] <- lapply(rstan::extract(surv, pars = "tau_year[2,2]"), FUN = mean); names(params)[41] <- "surv_eplus_y2"
params[42] <- lapply(rstan::extract(surv, pars = "tau_year[2,3]"), FUN = mean); names(params)[42] <- "surv_eplus_y3"
params[43] <- lapply(rstan::extract(surv, pars = "tau_year[2,4]"), FUN = mean); names(params)[43] <- "surv_eplus_y4"
params[44] <- lapply(rstan::extract(surv, pars = "tau_year[2,5]"), FUN = mean); names(params)[44] <- "surv_eplus_y5"
params[45] <- lapply(rstan::extract(surv, pars = "tau_year[2,6]"), FUN = mean); names(params)[45] <- "surv_eplus_y6"
params[46] <- lapply(rstan::extract(surv, pars = "tau_year[2,7]"), FUN = mean); names(params)[46] <- "surv_eplus_y7"
params[47] <- lapply(rstan::extract(surv, pars = "tau_year[2,8]"), FUN = mean); names(params)[47] <- "surv_eplus_y8"
params[48] <- lapply(rstan::extract(surv, pars = "tau_year[2,9]"), FUN = mean); names(params)[48] <- "surv_eplus_y9"
params[49] <- lapply(rstan::extract(surv, pars = "tau_year[2,10]"), FUN = mean); names(params)[49] <- "surv_eplus_y10"
params[50] <- lapply(rstan::extract(surv, pars = "tau_year[2,11]"), FUN = mean); names(params)[50] <- "surv_eplus_y11"

# grow
params[51] <- lapply(rstan::extract(grow, pars = "tau_year[1,1]"), FUN = mean); names(params)[51] <- "grow_eminus_y1"
params[52] <- lapply(rstan::extract(grow, pars = "tau_year[1,2]"), FUN = mean); names(params)[52] <- "grow_eminus_y2"
params[53] <- lapply(rstan::extract(grow, pars = "tau_year[1,3]"), FUN = mean); names(params)[53] <- "grow_eminus_y3"
params[54] <- lapply(rstan::extract(grow, pars = "tau_year[1,4]"), FUN = mean); names(params)[54] <- "grow_eminus_y4"
params[55] <- lapply(rstan::extract(grow, pars = "tau_year[1,5]"), FUN = mean); names(params)[55] <- "grow_eminus_y5"
params[56] <- lapply(rstan::extract(grow, pars = "tau_year[1,6]"), FUN = mean); names(params)[56] <- "grow_eminus_y6"
params[57] <- lapply(rstan::extract(grow, pars = "tau_year[1,7]"), FUN = mean); names(params)[57] <- "grow_eminus_y7"
params[58] <- lapply(rstan::extract(grow, pars = "tau_year[1,8]"), FUN = mean); names(params)[58] <- "grow_eminus_y8"
params[59] <- lapply(rstan::extract(grow, pars = "tau_year[1,9]"), FUN = mean); names(params)[59] <- "grow_eminus_y9"
params[60] <- lapply(rstan::extract(grow, pars = "tau_year[1,10]"), FUN = mean); names(params)[60] <- "grow_eminus_y10"
params[61] <- lapply(rstan::extract(grow, pars = "tau_year[1,11]"), FUN = mean); names(params)[61] <- "grow_eminus_y11"

params[62] <- lapply(rstan::extract(grow, pars = "tau_year[2,1]"), FUN = mean); names(params)[62] <- "grow_eplus_y1"
params[63] <- lapply(rstan::extract(grow, pars = "tau_year[2,2]"), FUN = mean); names(params)[63] <- "grow_eplus_y2"
params[64] <- lapply(rstan::extract(grow, pars = "tau_year[2,3]"), FUN = mean); names(params)[64] <- "grow_eplus_y3"
params[65] <- lapply(rstan::extract(grow, pars = "tau_year[2,4]"), FUN = mean); names(params)[65] <- "grow_eplus_y4"
params[66] <- lapply(rstan::extract(grow, pars = "tau_year[2,5]"), FUN = mean); names(params)[66] <- "grow_eplus_y5"
params[67] <- lapply(rstan::extract(grow, pars = "tau_year[2,6]"), FUN = mean); names(params)[67] <- "grow_eplus_y6"
params[68] <- lapply(rstan::extract(grow, pars = "tau_year[2,7]"), FUN = mean); names(params)[68] <- "grow_eplus_y7"
params[69] <- lapply(rstan::extract(grow, pars = "tau_year[2,8]"), FUN = mean); names(params)[69] <- "grow_eplus_y8"
params[70] <- lapply(rstan::extract(grow, pars = "tau_year[2,9]"), FUN = mean); names(params)[70] <- "grow_eplus_y9"
params[71] <- lapply(rstan::extract(grow, pars = "tau_year[2,10]"), FUN = mean); names(params)[71] <- "grow_eplus_y10"
params[72] <- lapply(rstan::extract(grow, pars = "tau_year[2,11]"), FUN = mean); names(params)[72] <- "grow_eplus_y11"

# flw
params[73] <- lapply(rstan::extract(flw, pars = "tau_year[1,1]"), FUN = mean); names(params)[73] <- "flw_eminus_y1"
params[74] <- lapply(rstan::extract(flw, pars = "tau_year[1,2]"), FUN = mean); names(params)[74] <- "flw_eminus_y2"
params[75] <- lapply(rstan::extract(flw, pars = "tau_year[1,3]"), FUN = mean); names(params)[75] <- "flw_eminus_y3"
params[76] <- lapply(rstan::extract(flw, pars = "tau_year[1,4]"), FUN = mean); names(params)[76] <- "flw_eminus_y4"
params[77] <- lapply(rstan::extract(flw, pars = "tau_year[1,5]"), FUN = mean); names(params)[77] <- "flw_eminus_y5"
params[78] <- lapply(rstan::extract(flw, pars = "tau_year[1,6]"), FUN = mean); names(params)[78] <- "flw_eminus_y6"
params[79] <- lapply(rstan::extract(flw, pars = "tau_year[1,7]"), FUN = mean); names(params)[79] <- "flw_eminus_y7"
params[80] <- lapply(rstan::extract(flw, pars = "tau_year[1,8]"), FUN = mean); names(params)[80] <- "flw_eminus_y8"
params[81] <- lapply(rstan::extract(flw, pars = "tau_year[1,9]"), FUN = mean); names(params)[81] <- "flw_eminus_y9"
params[82] <- lapply(rstan::extract(flw, pars = "tau_year[1,10]"), FUN = mean); names(params)[82] <- "flw_eminus_y10"
params[83] <- lapply(rstan::extract(flw, pars = "tau_year[1,11]"), FUN = mean); names(params)[83] <- "flw_eminus_y11"

params[84] <- lapply(rstan::extract(flw, pars = "tau_year[2,1]"), FUN = mean); names(params)[84] <- "flw_eplus_y1"
params[85] <- lapply(rstan::extract(flw, pars = "tau_year[2,2]"), FUN = mean); names(params)[85] <- "flw_eplus_y2"
params[86] <- lapply(rstan::extract(flw, pars = "tau_year[2,3]"), FUN = mean); names(params)[86] <- "flw_eplus_y3"
params[87] <- lapply(rstan::extract(flw, pars = "tau_year[2,4]"), FUN = mean); names(params)[87] <- "flw_eplus_y4"
params[88] <- lapply(rstan::extract(flw, pars = "tau_year[2,5]"), FUN = mean); names(params)[88] <- "flw_eplus_y5"
params[89] <- lapply(rstan::extract(flw, pars = "tau_year[2,6]"), FUN = mean); names(params)[89] <- "flw_eplus_y6"
params[90] <- lapply(rstan::extract(flw, pars = "tau_year[2,7]"), FUN = mean); names(params)[90] <- "flw_eplus_y7"
params[91] <- lapply(rstan::extract(flw, pars = "tau_year[2,8]"), FUN = mean); names(params)[91] <- "flw_eplus_y8"
params[92] <- lapply(rstan::extract(flw, pars = "tau_year[2,9]"), FUN = mean); names(params)[92] <- "flw_eplus_y9"
params[93] <- lapply(rstan::extract(flw, pars = "tau_year[2,10]"), FUN = mean); names(params)[93] <- "flw_eplus_y10"
params[94] <- lapply(rstan::extract(flw, pars = "tau_year[2,11]"), FUN = mean); names(params)[94] <- "flw_eplus_y11"

# fert
params[95] <- lapply(rstan::extract(fert, pars = "tau_year[1,1]"), FUN = mean); names(params)[95] <- "fert_eminus_y1"
params[96] <- lapply(rstan::extract(fert, pars = "tau_year[1,2]"), FUN = mean); names(params)[96] <- "fert_eminus_y2"
params[97] <- lapply(rstan::extract(fert, pars = "tau_year[1,3]"), FUN = mean); names(params)[97] <- "fert_eminus_y3"
params[98] <- lapply(rstan::extract(fert, pars = "tau_year[1,4]"), FUN = mean); names(params)[98] <- "fert_eminus_y4"
params[99] <- lapply(rstan::extract(fert, pars = "tau_year[1,5]"), FUN = mean); names(params)[99] <- "fert_eminus_y5"
params[100] <- lapply(rstan::extract(fert, pars = "tau_year[1,6]"), FUN = mean); names(params)[100] <- "fert_eminus_y6"
params[101] <- lapply(rstan::extract(fert, pars = "tau_year[1,7]"), FUN = mean); names(params)[101] <- "fert_eminus_y7"
params[102] <- lapply(rstan::extract(fert, pars = "tau_year[1,8]"), FUN = mean); names(params)[102] <- "fert_eminus_y8"
params[103] <- lapply(rstan::extract(fert, pars = "tau_year[1,9]"), FUN = mean); names(params)[103] <- "fert_eminus_y9"
params[104] <- lapply(rstan::extract(fert, pars = "tau_year[1,10]"), FUN = mean); names(params)[104] <- "fert_eminus_y10"
params[105] <- lapply(rstan::extract(fert, pars = "tau_year[1,11]"), FUN = mean); names(params)[105] <- "fert_eminus_y11"

params[106] <- lapply(rstan::extract(fert, pars = "tau_year[2,1]"), FUN = mean); names(params)[106] <- "fert_eplus_y1"
params[107] <- lapply(rstan::extract(fert, pars = "tau_year[2,2]"), FUN = mean); names(params)[107] <- "fert_eplus_y2"
params[108] <- lapply(rstan::extract(fert, pars = "tau_year[2,3]"), FUN = mean); names(params)[108] <- "fert_eplus_y3"
params[109] <- lapply(rstan::extract(fert, pars = "tau_year[2,4]"), FUN = mean); names(params)[109] <- "fert_eplus_y4"
params[110] <- lapply(rstan::extract(fert, pars = "tau_year[2,5]"), FUN = mean); names(params)[110] <- "fert_eplus_y5"
params[111] <- lapply(rstan::extract(fert, pars = "tau_year[2,6]"), FUN = mean); names(params)[111] <- "fert_eplus_y6"
params[112] <- lapply(rstan::extract(fert, pars = "tau_year[2,7]"), FUN = mean); names(params)[112] <- "fert_eplus_y7"
params[113] <- lapply(rstan::extract(fert, pars = "tau_year[2,8]"), FUN = mean); names(params)[113] <- "fert_eplus_y8"
params[114] <- lapply(rstan::extract(fert, pars = "tau_year[2,9]"), FUN = mean); names(params)[114] <- "fert_eplus_y9"
params[115] <- lapply(rstan::extract(fert, pars = "tau_year[2,10]"), FUN = mean); names(params)[115] <- "fert_eplus_y10"
params[116] <- lapply(rstan::extract(fert, pars = "tau_year[2,11]"), FUN = mean); names(params)[116] <- "fert_eplus_y11"

# spike
params[117] <- lapply(rstan::extract(spike, pars = "tau_year[1,1]"), FUN = mean); names(params)[117] <- "spike_eminus_y1"
params[118] <- lapply(rstan::extract(spike, pars = "tau_year[1,2]"), FUN = mean); names(params)[118] <- "spike_eminus_y2"
params[119] <- lapply(rstan::extract(spike, pars = "tau_year[1,3]"), FUN = mean); names(params)[119] <- "spike_eminus_y3"
params[120] <- lapply(rstan::extract(spike, pars = "tau_year[1,4]"), FUN = mean); names(params)[120] <- "spike_eminus_y4"
params[121] <- lapply(rstan::extract(spike, pars = "tau_year[1,5]"), FUN = mean); names(params)[121] <- "spike_eminus_y5"
params[122] <- lapply(rstan::extract(spike, pars = "tau_year[1,6]"), FUN = mean); names(params)[122] <- "spike_eminus_y6"
params[123] <- lapply(rstan::extract(spike, pars = "tau_year[1,7]"), FUN = mean); names(params)[123] <- "spike_eminus_y7"
params[124] <- lapply(rstan::extract(spike, pars = "tau_year[1,8]"), FUN = mean); names(params)[124] <- "spike_eminus_y8"
params[125] <- lapply(rstan::extract(spike, pars = "tau_year[1,9]"), FUN = mean); names(params)[125] <- "spike_eminus_y9"
params[126] <- lapply(rstan::extract(spike, pars = "tau_year[1,10]"), FUN = mean); names(params)[126] <- "spike_eminus_y10"
params[127] <- lapply(rstan::extract(spike, pars = "tau_year[1,11]"), FUN = mean); names(params)[127] <- "spike_eminus_y11"

params[128] <- lapply(rstan::extract(spike, pars = "tau_year[2,1]"), FUN = mean); names(params)[128] <- "spike_eplus_y1"
params[129] <- lapply(rstan::extract(spike, pars = "tau_year[2,2]"), FUN = mean); names(params)[129] <- "spike_eplus_y2"
params[130] <- lapply(rstan::extract(spike, pars = "tau_year[2,3]"), FUN = mean); names(params)[130] <- "spike_eplus_y3"
params[131] <- lapply(rstan::extract(spike, pars = "tau_year[2,4]"), FUN = mean); names(params)[131] <- "spike_eplus_y4"
params[132] <- lapply(rstan::extract(spike, pars = "tau_year[2,5]"), FUN = mean); names(params)[132] <- "spike_eplus_y5"
params[133] <- lapply(rstan::extract(spike, pars = "tau_year[2,6]"), FUN = mean); names(params)[133] <- "spike_eplus_y6"
params[134] <- lapply(rstan::extract(spike, pars = "tau_year[2,7]"), FUN = mean); names(params)[134] <- "spike_eplus_y7"
params[135] <- lapply(rstan::extract(spike, pars = "tau_year[2,8]"), FUN = mean); names(params)[135] <- "spike_eplus_y8"
params[136] <- lapply(rstan::extract(spike, pars = "tau_year[2,9]"), FUN = mean); names(params)[136] <- "spike_eplus_y9"
params[137] <- lapply(rstan::extract(spike, pars = "tau_year[2,10]"), FUN = mean); names(params)[137] <- "spike_eplus_y10"
params[138] <- lapply(rstan::extract(spike, pars = "tau_year[2,11]"), FUN = mean); names(params)[138] <- "spike_eplus_y11"

# s_to_s
params[139] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,1]"), FUN = mean); names(params)[139] <- "s_to_s_eminus_y1"
params[140] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,2]"), FUN = mean); names(params)[140] <- "s_to_s_eminus_y2"
params[141] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,3]"), FUN = mean); names(params)[141] <- "s_to_s_eminus_y3"
params[142] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,4]"), FUN = mean); names(params)[142] <- "s_to_s_eminus_y4"
params[143] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,5]"), FUN = mean); names(params)[143] <- "s_to_s_eminus_y5"
params[144] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,6]"), FUN = mean); names(params)[144] <- "s_to_s_eminus_y6"
params[145] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,7]"), FUN = mean); names(params)[145] <- "s_to_s_eminus_y7"
params[146] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,8]"), FUN = mean); names(params)[146] <- "s_to_s_eminus_y8"
params[147] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,9]"), FUN = mean); names(params)[147] <- "s_to_s_eminus_y9"
params[148] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,10]"), FUN = mean); names(params)[148] <- "s_to_s_eminus_y10"
params[149] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,11]"), FUN = mean); names(params)[149] <- "s_to_s_eminus_y11"

params[150] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,1]"), FUN = mean); names(params)[150] <- "s_to_s_eplus_y1"
params[151] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,2]"), FUN = mean); names(params)[151] <- "s_to_s_eplus_y2"
params[152] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,3]"), FUN = mean); names(params)[152] <- "s_to_s_eplus_y3"
params[153] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,4]"), FUN = mean); names(params)[153] <- "s_to_s_eplus_y4"
params[154] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,5]"), FUN = mean); names(params)[154] <- "s_to_s_eplus_y5"
params[155] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,6]"), FUN = mean); names(params)[155] <- "s_to_s_eplus_y6"
params[156] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,7]"), FUN = mean); names(params)[156] <- "s_to_s_eplus_y7"
params[157] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,8]"), FUN = mean); names(params)[157] <- "s_to_s_eplus_y8"
params[158] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,9]"), FUN = mean); names(params)[158] <- "s_to_s_eplus_y9"
params[159] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,10]"), FUN = mean); names(params)[159] <- "s_to_s_eplus_y10"
params[160] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,11]"), FUN = mean); names(params)[160] <- "s_to_s_eplus_y11"

params <- unlist(params)
return(params)
}

agpe_params <- getparams(surv = survAGPE, grow = growAGPE, flw = flwAGPE, fert = fertAGPE, spike = spikeAGPE, seed = seedAGPE, s_to_s = s_to_sAGPE, data = AGPE_surv_data)
elri_params <- getparams(surv = survELRI, grow = growELRI, flw = flwELRI, fert = fertELRI, spike = spikeELRI, seed = seedELRI, s_to_s = s_to_sELRI, data = ELRI_surv_data)
elvi_params <- getparams(surv = survELVI, grow = growELVI, flw = flwELVI, fert = fertELVI, spike = spikeELVI, seed = seedELVI, s_to_s = s_to_sELVI, data = ELVI_surv_data)
fesu_params <- getparams(surv = survFESU, grow = growFESU, flw = flwFESU, fert = fertFESU, spike = spikeFESU, seed = seedFESU, s_to_s = s_to_sFESU, data = FESU_surv_data)
loar_params <- getparams(surv = survLOAR, grow = growLOAR, flw = flwLOAR, fert = fertLOAR, spike = spikeLOAR, seed = seedLOAR, s_to_s = s_to_sLOAR, data = LOAR_surv_data)
poal_params <- getparams(surv = survPOAL, grow = growPOAL, flw = flwPOAL, fert = fertPOAL, spike = spikePOAL, seed = seedPOAL, s_to_s = s_to_sPOAL, data = POAL_surv_data)
posy_params <- getparams(surv = survPOSY, grow = growPOSY, flw = flwPOSY, fert = fertPOSY, spike = spikePOSY, seed = seedPOSY, s_to_s = s_to_sPOSY, data = POSY_surv_data)


# define functions that will be used to populate projection matrix
#SURVIVAL AT SIZE X.
# currently this is fitting as if the E- is the intercept
sx<-function(x,params, endo = FALSE, recruit = FALSE, year = FALSE){
  # Eminus Original Plants
  if(endo == FALSE & recruit == FALSE & year == FALSE) surv <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x))
  if(endo == FALSE & recruit == FALSE & year != FALSE) surv <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params[paste0("surv_eminus_y",year)])
  # Eminus recruit Plants
  if(endo == FALSE & recruit == TRUE & year == FALSE) surv <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta4"])
  if(endo == FALSE & recruit == TRUE & year != FALSE) surv <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta4"] + params[paste0("surv_eminus_y",year)])
  # Eplus priginal Plants
  if(endo == TRUE & recruit == FALSE & year == FALSE) surv <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"])
  if(endo == TRUE & recruit == FALSE & year != FALSE) surv <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params[paste0("surv_eplus_y",year)])
  # Eplus recruit Plants
  if(endo == TRUE & recruit == TRUE & year == FALSE) surv <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_beta4"])
  if(endo == TRUE & recruit == TRUE & year != FALSE) surv <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_beta4"] + params[paste0("surv_eplus_y",year)])
return(surv)
}
sx(2, loar_params, endo = FALSE, recruit = FALSE, year = 2)
sx<-function(x,params){
  # eminus original plants
  eminus_surv <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x))
  eminus_surv_y1 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_eminus_y1"])
  eminus_surv_y2 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_eminus_y2"])
  eminus_surv_y3 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_eminus_y3"])
  eminus_surv_y4 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_eminus_y4"])
  eminus_surv_y5 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_eminus_y5"])
  eminus_surv_y6 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_eminus_y6"])
  eminus_surv_y7 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_eminus_y7"])
  eminus_surv_y8 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_eminus_y8"])
  eminus_surv_y9 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_eminus_y9"])
  eminus_surv_y10 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_eminus_y10"])
  eminus_surv_y11 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_eminus_y11"])
  # eplus original plants
  eplus_surv <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"])
  eplus_surv_y1 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_eplus_y1"])
  eplus_surv_y2 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_eplus_y2"])
  eplus_surv_y3 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_eplus_y3"])
  eplus_surv_y4 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_eplus_y4"])
  eplus_surv_y5 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_eplus_y5"])
  eplus_surv_y6 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_eplus_y6"])
  eplus_surv_y7 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_eplus_y7"])
  eplus_surv_y8 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_eplus_y8"])
  eplus_surv_y9 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_eplus_y9"])
  eplus_surv_y10 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_eplus_y10"])
  eplus_surv_y11 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_eplus_y11"])
  # eminus recruit plants
  eminus_surv_rec <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta4"])
  eminus_surv_rec_y1 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta4"] + params["surv_eminus_y1"])
  eminus_surv_rec_y2 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta4"] + params["surv_eminus_y2"])
  eminus_surv_rec_y3 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta4"] + params["surv_eminus_y3"])
  eminus_surv_rec_y4 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta4"] + params["surv_eminus_y4"])
  eminus_surv_rec_y5 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta4"] + params["surv_eminus_y5"])
  eminus_surv_rec_y6 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta4"] + params["surv_eminus_y6"])
  eminus_surv_rec_y7 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta4"] + params["surv_eminus_y7"])
  eminus_surv_rec_y8 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta4"] + params["surv_eminus_y8"])
  eminus_surv_rec_y9 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta4"] + params["surv_eminus_y9"])
  eminus_surv_rec_y10 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta4"] + params["surv_eminus_y10"])
  eminus_surv_rec_y11 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta4"] + params["surv_eminus_y11"])
  # eplus recruit plants
  eplus_surv_rec <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_beta4"])
  eplus_surv_rec_y1 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_beta4"] + params["surv_eplus_y1"])
  eplus_surv_rec_y2 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_beta4"] + params["surv_eplus_y2"])
  eplus_surv_rec_y3 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_beta4"] + params["surv_eplus_y3"])
  eplus_surv_rec_y4 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_beta4"] + params["surv_eplus_y4"])
  eplus_surv_rec_y5 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_beta4"] + params["surv_eplus_y5"])
  eplus_surv_rec_y6 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_beta4"] + params["surv_eplus_y6"])
  eplus_surv_rec_y7 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_beta4"] + params["surv_eplus_y7"])
  eplus_surv_rec_y8 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_beta4"] + params["surv_eplus_y8"])
  eplus_surv_rec_y9 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_beta4"] + params["surv_eplus_y9"])
  eplus_surv_rec_y10 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_beta4"] + params["surv_eplus_y10"])
  eplus_surv_rec_y11 <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"] + params["surv_beta4"] + params["surv_eplus_y11"])
  
  return(list(eminus_surv=eminus_surv,eminus_surv_y1=eminus_surv_y1,eminus_surv_y2=eminus_surv_y2,eminus_surv_y3=eminus_surv_y3,eminus_surv_y4=eminus_surv_y4,eminus_surv_y5=eminus_surv_y5,eminus_surv_y6=eminus_surv_y6,eminus_surv_y7=eminus_surv_y7,eminus_surv_y8=eminus_surv_y8,eminus_surv_y9=eminus_surv_y9,eminus_surv_y10=eminus_surv_y10,eminus_surv_y11=eminus_surv_y11,
              eplus_surv=eplus_surv,eplus_surv_y1=eplus_surv_y1,eplus_surv_y2=eplus_surv_y2,eplus_surv_y3=eplus_surv_y3,eplus_surv_y4=eplus_surv_y4,eplus_surv_y5=eplus_surv_y5,eplus_surv_y6=eplus_surv_y6,eplus_surv_y7=eplus_surv_y7,eplus_surv_y8=eplus_surv_y8,eplus_surv_y9=eplus_surv_y9,eplus_surv_y10=eplus_surv_y10,eplus_surv_y11=eplus_surv_y11,
              eminus_surv_rec=eminus_surv_rec,eminus_surv_rec_y1=eminus_surv_rec_y1,eminus_surv_rec_y2=eminus_surv_rec_y2,eminus_surv_rec_y3=eminus_surv_rec_y3,eminus_surv_rec_y4=eminus_surv_rec_y4,eminus_surv_rec_y5=eminus_surv_rec_y5,eminus_surv_rec_y6=eminus_surv_rec_y6,eminus_surv_rec_y7=eminus_surv_rec_y7,eminus_surv_rec_y8=eminus_surv_rec_y8,eminus_surv_rec_y9=eminus_surv_rec_y9,eminus_surv_rec_y10=eminus_surv_rec_y10,eminus_surv_rec_y11=eminus_surv_rec_y11,
              eplus_surv_rec=eplus_surv_rec,eplus_surv_rec_y1=eplus_surv_rec_y1,eplus_surv_rec_y2=eplus_surv_rec_y2,eplus_surv_rec_y3=eplus_surv_rec_y3,eplus_surv_rec_y4=eplus_surv_rec_y4,eplus_surv_rec_y5=eplus_surv_rec_y5,eplus_surv_rec_y6=eplus_surv_rec_y6,eplus_surv_rec_y7=eplus_surv_rec_y7,eplus_surv_rec_y8=eplus_surv_rec_y8,eplus_surv_rec_y9=eplus_surv_rec_y9,eplus_surv_rec_y10=eplus_surv_rec_y10,eplus_surv_rec_y11=eplus_surv_rec_y11))
}
# sx(x = 10, params = poal_params) # we can test out our functions for different sizes of x
gxy <- function(x,y,params,endo = FALSE, recruit = FALSE, year = FALSE){
  # eminus original plants
  if(endo == FALSE & recruit == FALSE & year == FALSE) grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x)
  if(endo == FALSE & recruit == FALSE & year != FALSE) grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params[paste0("grow_eminus_y",year)]
  # eminus recruit plants
  if(endo == FALSE & recruit == TRUE & year == FALSE) grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta4"]
  if(endo == FALSE & recruit == TRUE & year != FALSE) grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta4"] + params[paste0("grow_eminus_y",year)]
  # eplus original plants
  if(endo == TRUE & recruit == FALSE & year == FALSE) grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"]
  if(endo == TRUE & recruit == FALSE & year != FALSE) grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x)  + params["grow_beta3"] + params[paste0("grow_eplus_y",year)]
  # eplus recruit plants
  if(endo == TRUE & recruit == TRUE & year == FALSE) grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_beta4"]
  if(endo == TRUE & recruit == TRUE & year != FALSE) grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x)  + params["grow_beta3"] + params["grow_beta4"] + params[paste0("grow_eplus_y",year)]
  
  grow <- dnbinom(x=y, mu = exp(grow.mean), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(grow.mean), size = params["grow_phi"])))
  return(grow)
}
# gxy(1,3,loar_params,endo = FALSE, recruit = TRUE, year = FALSE)
gxy <- function(x,y,params){
  # eminus original plants
  eminus_grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x)
  eminus_pr_grow <- dnbinom(x=y, mu = exp(eminus_grow.mean), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean), size = params["grow_phi"])))
  eminus_grow.mean_y1 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_eminus_y1"]
  eminus_pr_grow_y1 <- dnbinom(x=y, mu = exp(eminus_grow.mean_y1), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_y1), size = params["grow_phi"])))
  eminus_grow.mean_y2 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_eminus_y2"]
  eminus_pr_grow_y2 <- dnbinom(x=y, mu = exp(eminus_grow.mean_y2), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_y2), size = params["grow_phi"])))
  eminus_grow.mean_y3 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_eminus_y3"]
  eminus_pr_grow_y3 <- dnbinom(x=y, mu = exp(eminus_grow.mean_y3), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_y3), size = params["grow_phi"])))
  eminus_grow.mean_y4 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_eminus_y4"]
  eminus_pr_grow_y4 <- dnbinom(x=y, mu = exp(eminus_grow.mean_y4), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_y4), size = params["grow_phi"])))
  eminus_grow.mean_y5 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_eminus_y5"]
  eminus_pr_grow_y5 <- dnbinom(x=y, mu = exp(eminus_grow.mean_y5), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_y5), size = params["grow_phi"])))
  eminus_grow.mean_y6 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_eminus_y6"]
  eminus_pr_grow_y6 <- dnbinom(x=y, mu = exp(eminus_grow.mean_y6), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_y6), size = params["grow_phi"])))
  eminus_grow.mean_y7 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_eminus_y7"]
  eminus_pr_grow_y7 <- dnbinom(x=y, mu = exp(eminus_grow.mean_y7), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_y7), size = params["grow_phi"])))
  eminus_grow.mean_y8 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_eminus_y8"]
  eminus_pr_grow_y8 <- dnbinom(x=y, mu = exp(eminus_grow.mean_y8), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_y8), size = params["grow_phi"])))
  eminus_grow.mean_y9 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_eminus_y9"]
  eminus_pr_grow_y9 <- dnbinom(x=y, mu = exp(eminus_grow.mean_y9), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_y9), size = params["grow_phi"])))
  eminus_grow.mean_y10 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_eminus_y10"]
  eminus_pr_grow_y10 <- dnbinom(x=y, mu = exp(eminus_grow.mean_y10), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_y10), size = params["grow_phi"])))
  eminus_grow.mean_y11 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_eminus_y11"]
  eminus_pr_grow_y11 <- dnbinom(x=y, mu = exp(eminus_grow.mean_y11), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_y11), size = params["grow_phi"])))
  # eplus original plants
  eplus_grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"]
  eplus_pr_grow <- dnbinom(x=y, mu = exp(eplus_grow.mean), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean), size = params["grow_phi"])))
  eplus_grow.mean_y1 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_eplus_y1"]
  eplus_pr_grow_y1 <- dnbinom(x=y, mu = exp(eplus_grow.mean_y1), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_y1), size = params["grow_phi"])))
  eplus_grow.mean_y2 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_eplus_y2"]
  eplus_pr_grow_y2 <- dnbinom(x=y, mu = exp(eplus_grow.mean_y2), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_y2), size = params["grow_phi"])))
  eplus_grow.mean_y3 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_eplus_y3"]
  eplus_pr_grow_y3 <- dnbinom(x=y, mu = exp(eplus_grow.mean_y3), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_y3), size = params["grow_phi"])))
  eplus_grow.mean_y4 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_eplus_y4"]
  eplus_pr_grow_y4 <- dnbinom(x=y, mu = exp(eplus_grow.mean_y4), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_y4), size = params["grow_phi"])))
  eplus_grow.mean_y5 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_eplus_y5"]
  eplus_pr_grow_y5 <- dnbinom(x=y, mu = exp(eplus_grow.mean_y5), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_y5), size = params["grow_phi"])))
  eplus_grow.mean_y6 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_eplus_y6"]
  eplus_pr_grow_y6 <- dnbinom(x=y, mu = exp(eplus_grow.mean_y6), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_y6), size = params["grow_phi"])))
  eplus_grow.mean_y7 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_eplus_y7"]
  eplus_pr_grow_y7 <- dnbinom(x=y, mu = exp(eplus_grow.mean_y7), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_y7), size = params["grow_phi"])))
  eplus_grow.mean_y8 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_eplus_y8"]
  eplus_pr_grow_y8 <- dnbinom(x=y, mu = exp(eplus_grow.mean_y8), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_y8), size = params["grow_phi"])))
  eplus_grow.mean_y9 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_eplus_y9"]
  eplus_pr_grow_y9 <- dnbinom(x=y, mu = exp(eplus_grow.mean_y9), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_y9), size = params["grow_phi"])))
  eplus_grow.mean_y10 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_eplus_y10"]
  eplus_pr_grow_y10 <- dnbinom(x=y, mu = exp(eplus_grow.mean_y10), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_y10), size = params["grow_phi"])))
  eplus_grow.mean_y11 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_eplus_y11"]
  eplus_pr_grow_y11 <- dnbinom(x=y, mu = exp(eplus_grow.mean_y11), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_y11), size = params["grow_phi"])))
  # eminus recruit plants
  eminus_grow.mean_rec <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta4"]
  eminus_pr_grow_rec <- dnbinom(x=y, mu = exp(eminus_grow.mean_rec), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_rec), size = params["grow_phi"])))
  eminus_grow.mean_rec_y1 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta4"] + params["grow_eminus_y1"]
  eminus_pr_grow_rec_y1 <- dnbinom(x=y, mu = exp(eminus_grow.mean_rec_y1), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_rec_y1), size = params["grow_phi"])))
  eminus_grow.mean_rec_y2 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta4"] + params["grow_eminus_y2"]
  eminus_pr_grow_rec_y2 <- dnbinom(x=y, mu = exp(eminus_grow.mean_rec_y2), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_rec_y2), size = params["grow_phi"])))
  eminus_grow.mean_rec_y3 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta4"] + params["grow_eminus_y3"]
  eminus_pr_grow_rec_y3 <- dnbinom(x=y, mu = exp(eminus_grow.mean_rec_y3), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_rec_y3), size = params["grow_phi"])))
  eminus_grow.mean_rec_y4 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta4"] + params["grow_eminus_y4"]
  eminus_pr_grow_rec_y4 <- dnbinom(x=y, mu = exp(eminus_grow.mean_rec_y4), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_rec_y4), size = params["grow_phi"])))
  eminus_grow.mean_rec_y5 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta4"] + params["grow_eminus_y5"]
  eminus_pr_grow_rec_y5 <- dnbinom(x=y, mu = exp(eminus_grow.mean_rec_y5), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_rec_y5), size = params["grow_phi"])))
  eminus_grow.mean_rec_y6 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta4"] + params["grow_eminus_y6"]
  eminus_pr_grow_rec_y6 <- dnbinom(x=y, mu = exp(eminus_grow.mean_rec_y6), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_rec_y6), size = params["grow_phi"])))
  eminus_grow.mean_rec_y7 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta4"] + params["grow_eminus_y7"]
  eminus_pr_grow_rec_y7 <- dnbinom(x=y, mu = exp(eminus_grow.mean_rec_y7), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_rec_y7), size = params["grow_phi"])))
  eminus_grow.mean_rec_y8 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta4"] + params["grow_eminus_y8"]
  eminus_pr_grow_rec_y8 <- dnbinom(x=y, mu = exp(eminus_grow.mean_rec_y8), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_rec_y8), size = params["grow_phi"])))
  eminus_grow.mean_rec_y9 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta4"] + params["grow_eminus_y9"]
  eminus_pr_grow_rec_y9 <- dnbinom(x=y, mu = exp(eminus_grow.mean_rec_y9), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_rec_y9), size = params["grow_phi"])))
  eminus_grow.mean_rec_y10 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta4"] + params["grow_eminus_y10"]
  eminus_pr_grow_rec_y10 <- dnbinom(x=y, mu = exp(eminus_grow.mean_rec_y10), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_rec_y10), size = params["grow_phi"])))
  eminus_grow.mean_rec_y11 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta4"] + params["grow_eminus_y11"]
  eminus_pr_grow_rec_y11 <- dnbinom(x=y, mu = exp(eminus_grow.mean_rec_y11), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean_rec_y11), size = params["grow_phi"])))
  # eplus recruit plants
  eplus_grow.mean_rec <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_beta4"]
  eplus_pr_grow_rec <- dnbinom(x=y, mu = exp(eplus_grow.mean_rec), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_rec), size = params["grow_phi"])))
  eplus_grow.mean_rec_y1 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_beta4"] + params["grow_eplus_y1"]
  eplus_pr_grow_rec_y1 <- dnbinom(x=y, mu = exp(eplus_grow.mean_rec_y1), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_rec_y1), size = params["grow_phi"])))
  eplus_grow.mean_rec_y2 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_beta4"] + params["grow_eplus_y2"]
  eplus_pr_grow_rec_y2 <- dnbinom(x=y, mu = exp(eplus_grow.mean_rec_y2), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_rec_y2), size = params["grow_phi"])))
  eplus_grow.mean_rec_y3 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_beta4"] + params["grow_eplus_y3"]
  eplus_pr_grow_rec_y3 <- dnbinom(x=y, mu = exp(eplus_grow.mean_rec_y3), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_rec_y3), size = params["grow_phi"])))
  eplus_grow.mean_rec_y4 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_beta4"] + params["grow_eplus_y4"]
  eplus_pr_grow_rec_y4 <- dnbinom(x=y, mu = exp(eplus_grow.mean_rec_y4), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_rec_y4), size = params["grow_phi"])))
  eplus_grow.mean_rec_y5 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_beta4"] + params["grow_eplus_y5"]
  eplus_pr_grow_rec_y5 <- dnbinom(x=y, mu = exp(eplus_grow.mean_rec_y5), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_rec_y5), size = params["grow_phi"])))
  eplus_grow.mean_rec_y6 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_beta4"] + params["grow_eplus_y6"]
  eplus_pr_grow_rec_y6 <- dnbinom(x=y, mu = exp(eplus_grow.mean_rec_y6), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_rec_y6), size = params["grow_phi"])))
  eplus_grow.mean_rec_y7 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_beta4"] + params["grow_eplus_y7"]
  eplus_pr_grow_rec_y7 <- dnbinom(x=y, mu = exp(eplus_grow.mean_rec_y7), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_rec_y7), size = params["grow_phi"])))
  eplus_grow.mean_rec_y8 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_beta4"] + params["grow_eplus_y8"]
  eplus_pr_grow_rec_y8 <- dnbinom(x=y, mu = exp(eplus_grow.mean_rec_y8), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_rec_y8), size = params["grow_phi"])))
  eplus_grow.mean_rec_y9 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_beta4"] + params["grow_eplus_y9"]
  eplus_pr_grow_rec_y9 <- dnbinom(x=y, mu = exp(eplus_grow.mean_rec_y9), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_rec_y9), size = params["grow_phi"])))
  eplus_grow.mean_rec_y10 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_beta4"] + params["grow_eplus_y10"]
  eplus_pr_grow_rec_y10 <- dnbinom(x=y, mu = exp(eplus_grow.mean_rec_y10), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_rec_y10), size = params["grow_phi"])))
  eplus_grow.mean_rec_y11 <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"] + params["grow_beta4"] + params["grow_eplus_y11"]
  eplus_pr_grow_rec_y11 <- dnbinom(x=y, mu = exp(eplus_grow.mean_rec_y11), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean_rec_y11), size = params["grow_phi"])))

  return(list(eminus_pr_grow=eminus_pr_grow,eminus_pr_grow_y1=eminus_pr_grow_y1,eminus_pr_grow_y2=eminus_pr_grow_y2,eminus_pr_grow_y3=eminus_pr_grow_y3,eminus_pr_grow_y4=eminus_pr_grow_y4,eminus_pr_grow_y5=eminus_pr_grow_y5,eminus_pr_grow_y6=eminus_pr_grow_y6,eminus_pr_grow_y7=eminus_pr_grow_y7,eminus_pr_grow_y8=eminus_pr_grow_y8,eminus_pr_grow_y9=eminus_pr_grow_y9,eminus_pr_grow_y10=eminus_pr_grow_y10,eminus_pr_grow_y11=eminus_pr_grow_y11, 
            eplus_pr_grow=eplus_pr_grow,eplus_pr_grow_y1=eplus_pr_grow_y1,eplus_pr_grow_y2=eplus_pr_grow_y2,eplus_pr_grow_y3=eplus_pr_grow_y3,eplus_pr_grow_y4=eplus_pr_grow_y4,eplus_pr_grow_y5=eplus_pr_grow_y5,eplus_pr_grow_y6=eplus_pr_grow_y6,eplus_pr_grow_y7=eplus_pr_grow_y7,eplus_pr_grow_y8=eplus_pr_grow_y8,eplus_pr_grow_y9=eplus_pr_grow_y9,eplus_pr_grow_y10=eplus_pr_grow_y10,eplus_pr_grow_y11=eplus_pr_grow_y11, 
            eminus_pr_grow_rec=eminus_pr_grow_rec, eminus_pr_grow_rec_y1=eminus_pr_grow_rec_y1,eminus_pr_grow_rec_y2=eminus_pr_grow_rec_y2,eminus_pr_grow_rec_y3=eminus_pr_grow_rec_y3,eminus_pr_grow_rec_y4=eminus_pr_grow_rec_y4,eminus_pr_grow_rec_y5=eminus_pr_grow_rec_y5,eminus_pr_grow_rec_y6=eminus_pr_grow_rec_y6,eminus_pr_grow_rec_y7=eminus_pr_grow_rec_y7,eminus_pr_grow_rec_y8=eminus_pr_grow_rec_y8,eminus_pr_grow_rec_y9=eminus_pr_grow_rec_y9,eminus_pr_grow_rec_y10=eminus_pr_grow_rec_y10,eminus_pr_grow_rec_y11=eminus_pr_grow_rec_y11, 
            eplus_pr_grow_rec=eplus_pr_grow_rec,eplus_pr_grow_rec_y1=eplus_pr_grow_rec_y1,eplus_pr_grow_rec_y2=eplus_pr_grow_rec_y2,eplus_pr_grow_rec_y3=eplus_pr_grow_rec_y3,eplus_pr_grow_rec_y4=eplus_pr_grow_rec_y4,eplus_pr_grow_rec_y5=eplus_pr_grow_rec_y5,eplus_pr_grow_rec_y6=eplus_pr_grow_rec_y6,eplus_pr_grow_rec_y7=eplus_pr_grow_rec_y7,eplus_pr_grow_rec_y8=eplus_pr_grow_rec_y8,eplus_pr_grow_rec_y9=eplus_pr_grow_rec_y9,eplus_pr_grow_rec_y10=eplus_pr_grow_rec_y10,eplus_pr_grow_rec_y11=eplus_pr_grow_rec_y11))
}
# gxy(x = 10, y = 11, params = loar_params)
# I'll truncate this later...
# prob=dnbinom(1:n_post_draws, mu = exp(mu[i,j]), size = phi[i]/(1-dnbinom(0, mu = exp(mu[i,j]), size = phi[i]))))


#SURVIVAL*GROWTH
# This function generates the transition matrix by multiplying growth and survival.
pxy<-function(x,y,params,endo = FALSE, recruit = FALSE, year = FALSE){
  sx(x,params, endo, recruit, year) * gxy(x,y,params, endo, recruit, year)
}
# pxy(x = 1,y = 4, params = loar_params, endo = TRUE, recruit = FALSE, year = FALSE)

pxy<-function(x,y,params,output = ""){
  # eminus original plants
  pxy_eminus <- sx(x,params)$eminus_surv * gxy(x,y,params)$eminus_pr_grow
  pxy_eminus_y1 <- sx(x,params)$eminus_surv_y1 * gxy(x,y,params)$eminus_pr_grow_y1
  pxy_eminus_y2 <- sx(x,params)$eminus_surv_y1 * gxy(x,y,params)$eminus_pr_grow_y2
  pxy_eminus_y3 <- sx(x,params)$eminus_surv_y1 * gxy(x,y,params)$eminus_pr_grow_y3
  pxy_eminus_y4 <- sx(x,params)$eminus_surv_y1 * gxy(x,y,params)$eminus_pr_grow_y4
  pxy_eminus_y5 <- sx(x,params)$eminus_surv_y1 * gxy(x,y,params)$eminus_pr_grow_y5
  pxy_eminus_y6 <- sx(x,params)$eminus_surv_y1 * gxy(x,y,params)$eminus_pr_grow_y6
  pxy_eminus_y7 <- sx(x,params)$eminus_surv_y1 * gxy(x,y,params)$eminus_pr_grow_y7
  pxy_eminus_y8 <- sx(x,params)$eminus_surv_y1 * gxy(x,y,params)$eminus_pr_grow_y8
  pxy_eminus_y9 <- sx(x,params)$eminus_surv_y1 * gxy(x,y,params)$eminus_pr_grow_y9
  pxy_eminus_y10 <- sx(x,params)$eminus_surv_y1 * gxy(x,y,params)$eminus_pr_grow_y10
  pxy_eminus_y11 <- sx(x,params)$eminus_surv_y1 * gxy(x,y,params)$eminus_pr_grow_y11
  
  # eplus original plants
  pxy_eplus <- sx(x,params)$eplus_surv * gxy(x,y,params)$eplus_pr_grow
  pxy_eplus_y1 <- sx(x,params)$eplus_surv_y1 * gxy(x,y,params)$eplus_pr_grow_y1
  pxy_eplus_y2 <- sx(x,params)$eplus_surv_y1 * gxy(x,y,params)$eplus_pr_grow_y2
  pxy_eplus_y3 <- sx(x,params)$eplus_surv_y1 * gxy(x,y,params)$eplus_pr_grow_y3
  pxy_eplus_y4 <- sx(x,params)$eplus_surv_y1 * gxy(x,y,params)$eplus_pr_grow_y4
  pxy_eplus_y5 <- sx(x,params)$eplus_surv_y1 * gxy(x,y,params)$eplus_pr_grow_y5
  pxy_eplus_y6 <- sx(x,params)$eplus_surv_y1 * gxy(x,y,params)$eplus_pr_grow_y6
  pxy_eplus_y7 <- sx(x,params)$eplus_surv_y1 * gxy(x,y,params)$eplus_pr_grow_y7
  pxy_eplus_y8 <- sx(x,params)$eplus_surv_y1 * gxy(x,y,params)$eplus_pr_grow_y8
  pxy_eplus_y9 <- sx(x,params)$eplus_surv_y1 * gxy(x,y,params)$eplus_pr_grow_y9
  pxy_eplus_y10 <- sx(x,params)$eplus_surv_y1 * gxy(x,y,params)$eplus_pr_grow_y10
  pxy_eplus_y11 <- sx(x,params)$eplus_surv_y1 * gxy(x,y,params)$eplus_pr_grow_y11
  # eminus recruit plants
  pxy_eminus_rec <- sx(x,params)$eminus_surv_rec * gxy(x,y,params)$eminus_pr_grow_rec
  pxy_eminus_rec_y1 <- sx(x,params)$eminus_surv_rec_y1 * gxy(x,y,params)$eminus_pr_grow_rec_y1
  pxy_eminus_rec_y2 <- sx(x,params)$eminus_surv_rec_y2 * gxy(x,y,params)$eminus_pr_grow_rec_y2
  pxy_eminus_rec_y3 <- sx(x,params)$eminus_surv_rec_y3 * gxy(x,y,params)$eminus_pr_grow_rec_y3
  pxy_eminus_rec_y4 <- sx(x,params)$eminus_surv_rec_y4 * gxy(x,y,params)$eminus_pr_grow_rec_y4
  pxy_eminus_rec_y5 <- sx(x,params)$eminus_surv_rec_y5 * gxy(x,y,params)$eminus_pr_grow_rec_y5
  pxy_eminus_rec_y6 <- sx(x,params)$eminus_surv_rec_y6 * gxy(x,y,params)$eminus_pr_grow_rec_y6
  pxy_eminus_rec_y7 <- sx(x,params)$eminus_surv_rec_y7 * gxy(x,y,params)$eminus_pr_grow_rec_y7
  pxy_eminus_rec_y8 <- sx(x,params)$eminus_surv_rec_y8 * gxy(x,y,params)$eminus_pr_grow_rec_y8
  pxy_eminus_rec_y9 <- sx(x,params)$eminus_surv_rec_y9 * gxy(x,y,params)$eminus_pr_grow_rec_y9
  pxy_eminus_rec_y10 <- sx(x,params)$eminus_surv_rec_y10 * gxy(x,y,params)$eminus_pr_grow_rec_y10
  pxy_eminus_rec_y11 <- sx(x,params)$eminus_surv_rec_y11 * gxy(x,y,params)$eminus_pr_grow_rec_y11
  # eplus recruit plants
  pxy_eplus_rec <- sx(x,params)$eplus_surv_rec * gxy(x,y,params)$eplus_pr_grow_rec
  pxy_eplus_rec_y1 <- sx(x,params)$eplus_surv_rec_y1 * gxy(x,y,params)$eplus_pr_grow_rec_y1
  pxy_eplus_rec_y2 <- sx(x,params)$eplus_surv_rec_y2 * gxy(x,y,params)$eplus_pr_grow_rec_y2
  pxy_eplus_rec_y3 <- sx(x,params)$eplus_surv_rec_y3 * gxy(x,y,params)$eplus_pr_grow_rec_y3
  pxy_eplus_rec_y4 <- sx(x,params)$eplus_surv_rec_y4 * gxy(x,y,params)$eplus_pr_grow_rec_y4
  pxy_eplus_rec_y5 <- sx(x,params)$eplus_surv_rec_y5 * gxy(x,y,params)$eplus_pr_grow_rec_y5
  pxy_eplus_rec_y6 <- sx(x,params)$eplus_surv_rec_y6 * gxy(x,y,params)$eplus_pr_grow_rec_y6
  pxy_eplus_rec_y7 <- sx(x,params)$eplus_surv_rec_y7 * gxy(x,y,params)$eplus_pr_grow_rec_y7
  pxy_eplus_rec_y8 <- sx(x,params)$eplus_surv_rec_y8 * gxy(x,y,params)$eplus_pr_grow_rec_y8
  pxy_eplus_rec_y9 <- sx(x,params)$eplus_surv_rec_y9 * gxy(x,y,params)$eplus_pr_grow_rec_y9
  pxy_eplus_rec_y10 <- sx(x,params)$eplus_surv_rec_y10 * gxy(x,y,params)$eplus_pr_grow_rec_y10
  pxy_eplus_rec_y11 <- sx(x,params)$eplus_surv_rec_y11 * gxy(x,y,params)$eplus_pr_grow_rec_y11
  
  result <- list(pxy_eminus=pxy_eminus,pxy_eminus_y1=pxy_eminus_y1,pxy_eminus_y2=pxy_eminus_y2,pxy_eminus_y3=pxy_eminus_y3,pxy_eminus_y4=pxy_eminus_y4,pxy_eminus_y5=pxy_eminus_y5,pxy_eminus_y6=pxy_eminus_y6,pxy_eminus_y7=pxy_eminus_y7,pxy_eminus_y8=pxy_eminus_y8,pxy_eminus_y9=pxy_eminus_y9,pxy_eminus_y10=pxy_eminus_y10,pxy_eminus_y11=pxy_eminus_y11,
                 pxy_eplus=pxy_eplus,pxy_eplus_y1=pxy_eplus_y1,pxy_eplus_y2=pxy_eplus_y2,pxy_eplus_y3=pxy_eplus_y3,pxy_eplus_y4=pxy_eplus_y4,pxy_eplus_y5=pxy_eplus_y5,pxy_eplus_y6=pxy_eplus_y6,pxy_eplus_y7=pxy_eplus_y7,pxy_eplus_y8=pxy_eplus_y8,pxy_eplus_y9=pxy_eplus_y9,pxy_eplus_y10=pxy_eplus_y10,pxy_eplus_y11=pxy_eplus_y11,
                 pxy_eminus_rec=pxy_eminus_rec,pxy_eminus_rec_y1=pxy_eminus_rec_y1,pxy_eminus_rec_y2=pxy_eminus_rec_y2,pxy_eminus_rec_y3=pxy_eminus_rec_y3,pxy_eminus_rec_y4=pxy_eminus_rec_y4,pxy_eminus_rec_y5=pxy_eminus_rec_y5,pxy_eminus_rec_y6=pxy_eminus_rec_y6,pxy_eminus_rec_y7=pxy_eminus_rec_y7,pxy_eminus_rec_y8=pxy_eminus_rec_y8,pxy_eminus_rec_y9=pxy_eminus_rec_y9,pxy_eminus_rec_y10=pxy_eminus_rec_y10,pxy_eminus_rec_y11=pxy_eminus_rec_y11,
                 pxy_eplus_rec=pxy_eplus_rec,pxy_eplus_rec_y1=pxy_eplus_rec_y1,pxy_eplus_rec_y2=pxy_eplus_rec_y2,pxy_eplus_rec_y3=pxy_eplus_rec_y3,pxy_eplus_rec_y4=pxy_eplus_rec_y4,pxy_eplus_rec_y5=pxy_eplus_rec_y5,pxy_eplus_rec_y6=pxy_eplus_rec_y6,pxy_eplus_rec_y7=pxy_eplus_rec_y7,pxy_eplus_rec_y8=pxy_eplus_rec_y8,pxy_eplus_rec_y9=pxy_eplus_rec_y9,pxy_eplus_rec_y10=pxy_eplus_rec_y10,pxy_eplus_rec_y11=pxy_eplus_rec_y11)
         if(output == "pxy_eminus") return(result$pxy_eminus)
         if(output == "pxy_eminus_y1") return(result$pxy_eminus_y1)
         if(output == "pxy_eminus_y2") return(result$pxy_eminus_y2)
         if(output == "pxy_eminus_y3") return(result$pxy_eminus_y3)
         if(output == "pxy_eminus_y4") return(result$pxy_eminus_y4)
         if(output == "pxy_eminus_y5") return(result$pxy_eminus_y5)
         if(output == "pxy_eminus_y6") return(result$pxy_eminus_y6)
         if(output == "pxy_eminus_y7") return(result$pxy_eminus_y7)
         if(output == "pxy_eminus_y8") return(result$pxy_eminus_y8)
         if(output == "pxy_eminus_y9") return(result$pxy_eminus_y9)
         if(output == "pxy_eminus_y10") return(result$pxy_eminus_y10)
         if(output == "pxy_eminus_y11") return(result$pxy_eminus_y11)
  
         if(output == "pxy_eplus") return(result$pxy_eplus) 
         if(output == "pxy_eplus_y1") return(result$pxy_eplus_y1)
         if(output == "pxy_eplus_y2") return(result$pxy_eplus_y2)
         if(output == "pxy_eplus_y3") return(result$pxy_eplus_y3)
         if(output == "pxy_eplus_y4") return(result$pxy_eplus_y4)
         if(output == "pxy_eplus_y5") return(result$pxy_eplus_y5)
         if(output == "pxy_eplus_y6") return(result$pxy_eplus_y6)
         if(output == "pxy_eplus_y7") return(result$pxy_eplus_y7)
         if(output == "pxy_eplus_y8") return(result$pxy_eplus_y8)
         if(output == "pxy_eplus_y9") return(result$pxy_eplus_y9)
         if(output == "pxy_eplus_y10") return(result$pxy_eplus_y10)
         if(output == "pxy_eplus_y11") return(result$pxy_eplus_y11)
  
         if(output == "pxy_eminus_rec") return(result$pxy_eminus_rec)
         if(output == "pxy_eminus_rec_y1") return(result$pxy_eminus_rec_y1)
         if(output == "pxy_eminus_rec_y2") return(result$pxy_eminus_rec_y2)
         if(output == "pxy_eminus_rec_y3") return(result$pxy_eminus_rec_y3)
         if(output == "pxy_eminus_rec_y4") return(result$pxy_eminus_rec_y4)
         if(output == "pxy_eminus_rec_y5") return(result$pxy_eminus_rec_y5)
         if(output == "pxy_eminus_rec_y6") return(result$pxy_eminus_rec_y6)
         if(output == "pxy_eminus_rec_y7") return(result$pxy_eminus_rec_y7)
         if(output == "pxy_eminus_rec_y8") return(result$pxy_eminus_rec_y8)
         if(output == "pxy_eminus_rec_y9") return(result$pxy_eminus_rec_y9)
         if(output == "pxy_eminus_rec_y10") return(result$pxy_eminus_rec_y10)
         if(output == "pxy_eminus_rec_y11") return(result$pxy_eminus_rec_y11)
  
         if(output == "pxy_eplus_rec") return(result$pxy_eplus_rec) 
         if(output == "pxy_eplus_rec_y1") return(result$pxy_eplus_rec_y1)
         if(output == "pxy_eplus_rec_y2") return(result$pxy_eplus_rec_y2)
         if(output == "pxy_eplus_rec_y3") return(result$pxy_eplus_rec_y3)
         if(output == "pxy_eplus_rec_y4") return(result$pxy_eplus_rec_y4)
         if(output == "pxy_eplus_rec_y5") return(result$pxy_eplus_rec_y5)
         if(output == "pxy_eplus_rec_y6") return(result$pxy_eplus_rec_y6)
         if(output == "pxy_eplus_rec_y7") return(result$pxy_eplus_rec_y7)
         if(output == "pxy_eplus_rec_y8") return(result$pxy_eplus_rec_y8)
         if(output == "pxy_eplus_rec_y9") return(result$pxy_eplus_rec_y9)
         if(output == "pxy_eplus_rec_y10") return(result$pxy_eplus_rec_y10)
         if(output == "pxy_eplus_rec_y11") return(result$pxy_eplus_rec_y11)
  else return(result)
}
pxy(x = 10,y = 1,params = loar_params, output = "pxy_eplus_rec_y1")


#FERTILITY--returns number of seedlings, which we will assume (for now) to be 1-tiller, produced by size X
fx<-function(x, params, endo = FALSE, recruit = FALSE, year = FALSE){
# flowering
  # Eminus original plant
if(endo == FALSE & recruit == FALSE & year == FALSE){
  flw <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x))
  fert <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x))
  spike <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x))
  recruitment <- invlogit(params["s_to_s_beta1"])
}
if(endo == FALSE & recruit == FALSE & year != FALSE) {
  flw <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params[paste0("flw_eminus_y",year)])
  fert <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params[paste0("fert_eminus_y", year)])
  spike <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params[paste0("spike_eminus_y", year)])
  recruitment <- invlogit(params["s_to_s_beta1"] + params[paste0("s_to_s_eminus_y", year)])
}
  # Eminus recruit plant
if(endo == FALSE & recruit == TRUE & year == FALSE) {
  flw <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta4"])
  fert <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta4"])
  spike <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta4"])
  recruitment <- invlogit(params["s_to_s_beta1"])
}
if(endo == FALSE & recruit == TRUE & year != FALSE) {
  flw <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta4"] + params[paste0("flw_eminus_y",year)])
  fert <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta4"] + params[paste0("fert_eminus_y",year)])
  spike <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta4"] + params[paste0("spike_eminus_y",year)])
  recruitment <- invlogit(params["s_to_s_beta1"] + params[paste0("s_to_s_eminus_y", year)])
}
  # Eplus original plant
if(endo == TRUE & recruit == FALSE & year == FALSE) {
  flw <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"])
  fert <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"])
  spike <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"])
  recruitment <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"])
}
if(endo == TRUE & recruit == FALSE & year != FALSE) {
  flw <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params[paste0("flw_eplus_y",year)])
  fert <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params[paste0("fert_eplus_y",year)])
  spike <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params[paste0("spike_eplus_y",year)])
  recruitment <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"]+ params[paste0("s_to_s_eplus_y", year)])
}
  # EPlus recruit plant
if(endo == TRUE & recruit == TRUE & year == FALSE) {
  flw <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_beta4"])
  fert <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_beta4"])
  spike <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_beta4"])
  recruitment <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"])
}
if(endo == TRUE & recruit == TRUE & year != FALSE) {
  flw <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_beta4"] + params[paste0("flw_eplus_y",year)])
  fert <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_beta4"] + params[paste0("fert_eplus_y",year)])
  spike <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_beta4"] + params[paste0("spike_eplus_y",year)])
  recruitment <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"]+ params[paste0("s_to_s_eplus_y", year)])
}
# mean seed per spikelet
seed.mean <- params["mu_seed"]

seedlings <- flw * fert * spike * seed.mean * recruitment
return(seedlings)
}
fx(x=3, params = loar_params, endo = TRUE, year = 1)


fx<-function(x, params){
  # flowering
    # eminus original plants
    eminus_flw <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x)) 
    eminus_flw_y1 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_eminus_y1"]) 
    eminus_flw_y2 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_eminus_y2"]) 
    eminus_flw_y3 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_eminus_y3"]) 
    eminus_flw_y4 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_eminus_y4"]) 
    eminus_flw_y5 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_eminus_y5"]) 
    eminus_flw_y6 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_eminus_y6"]) 
    eminus_flw_y7 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_eminus_y7"]) 
    eminus_flw_y8 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_eminus_y8"]) 
    eminus_flw_y9 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_eminus_y9"]) 
    eminus_flw_y10 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_eminus_y10"]) 
    eminus_flw_y11 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_eminus_y11"]) 
    # eplus original plants
    eplus_flw <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"])
    eplus_flw_y1 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_eplus_y1"])
    eplus_flw_y2 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_eplus_y2"])
    eplus_flw_y3 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_eplus_y3"])
    eplus_flw_y4 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_eplus_y4"])
    eplus_flw_y5 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_eplus_y5"])
    eplus_flw_y6 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_eplus_y6"])
    eplus_flw_y7 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_eplus_y7"])
    eplus_flw_y8 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_eplus_y8"])
    eplus_flw_y9 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_eplus_y9"])
    eplus_flw_y10 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_eplus_y10"])
    eplus_flw_y11 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_eplus_y11"])
    # eminus recruit plants
    eminus_flw_rec <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta4"]) 
    eminus_flw_rec_y1 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta4"] + params["flw_eminus_y1"]) 
    eminus_flw_rec_y2 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta4"] + params["flw_eminus_y2"]) 
    eminus_flw_rec_y3 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta4"] + params["flw_eminus_y3"]) 
    eminus_flw_rec_y4 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta4"] + params["flw_eminus_y4"]) 
    eminus_flw_rec_y5 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta4"] + params["flw_eminus_y5"]) 
    eminus_flw_rec_y6 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta4"] + params["flw_eminus_y6"]) 
    eminus_flw_rec_y7 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta4"] + params["flw_eminus_y7"]) 
    eminus_flw_rec_y8 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta4"] + params["flw_eminus_y8"]) 
    eminus_flw_rec_y9 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta4"] + params["flw_eminus_y9"]) 
    eminus_flw_rec_y10 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta4"] + params["flw_eminus_y10"]) 
    eminus_flw_rec_y11 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta4"] + params["flw_eminus_y11"]) 
    # eplus recruit plants
    eplus_flw_rec <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_beta4"])
    eplus_flw_rec_y1 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_beta4"] + params["flw_eplus_y1"])
    eplus_flw_rec_y2 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_beta4"] + params["flw_eplus_y2"])
    eplus_flw_rec_y3 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_beta4"] + params["flw_eplus_y3"])
    eplus_flw_rec_y4 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_beta4"] + params["flw_eplus_y4"])
    eplus_flw_rec_y5 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_beta4"] + params["flw_eplus_y5"])
    eplus_flw_rec_y6 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_beta4"] + params["flw_eplus_y6"])
    eplus_flw_rec_y7 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_beta4"] + params["flw_eplus_y7"])
    eplus_flw_rec_y8 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_beta4"] + params["flw_eplus_y8"])
    eplus_flw_rec_y9 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_beta4"] + params["flw_eplus_y9"])
    eplus_flw_rec_y10 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_beta4"] + params["flw_eplus_y10"])
    eplus_flw_rec_y11 <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"] + params["flw_beta4"] + params["flw_eplus_y11"])

  # flw tiller no.
    # eminus original plants
    eminus_fert <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x))
    eminus_fert_y1 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_eminus_y1"])
    eminus_fert_y2 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_eminus_y2"])
    eminus_fert_y3 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_eminus_y3"])
    eminus_fert_y4 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_eminus_y4"])
    eminus_fert_y5 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_eminus_y5"])
    eminus_fert_y6 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_eminus_y6"])
    eminus_fert_y7 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_eminus_y7"])
    eminus_fert_y8 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_eminus_y8"])
    eminus_fert_y9 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_eminus_y9"])
    eminus_fert_y10 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_eminus_y10"])
    eminus_fert_y11 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_eminus_y11"])
    # eplus original plants
    eplus_fert <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"])
    eplus_fert_y1 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_eplus_y1"])
    eplus_fert_y2 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_eplus_y2"])
    eplus_fert_y3 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_eplus_y3"])
    eplus_fert_y4 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_eplus_y4"])
    eplus_fert_y5 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_eplus_y5"])
    eplus_fert_y6 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_eplus_y6"])
    eplus_fert_y7 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_eplus_y7"])
    eplus_fert_y8 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_eplus_y8"])
    eplus_fert_y9 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_eplus_y9"])
    eplus_fert_y10 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_eplus_y10"])
    eplus_fert_y11 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_eplus_y11"])
    # eminus recruit plants
    eminus_fert_rec <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x)+ params["fert_beta4"])
    eminus_fert_rec_y1 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x)+ params["fert_beta4"] + params["fert_eminus_y1"])
    eminus_fert_rec_y2 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x)+ params["fert_beta4"] + params["fert_eminus_y2"])
    eminus_fert_rec_y3 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x)+ params["fert_beta4"] + params["fert_eminus_y3"])
    eminus_fert_rec_y4 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x)+ params["fert_beta4"] + params["fert_eminus_y4"])
    eminus_fert_rec_y5 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x)+ params["fert_beta4"] + params["fert_eminus_y5"])
    eminus_fert_rec_y6 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x)+ params["fert_beta4"] + params["fert_eminus_y6"])
    eminus_fert_rec_y7 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x)+ params["fert_beta4"] + params["fert_eminus_y7"])
    eminus_fert_rec_y8 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x)+ params["fert_beta4"] + params["fert_eminus_y8"])
    eminus_fert_rec_y9 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x)+ params["fert_beta4"] + params["fert_eminus_y9"])
    eminus_fert_rec_y10 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x)+ params["fert_beta4"] + params["fert_eminus_y10"])
    eminus_fert_rec_y11 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x)+ params["fert_beta4"] + params["fert_eminus_y11"])
    # eplus recruit plants
    eplus_fert_rec <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_beta4"])
    eplus_fert_rec_y1 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_beta4"] + params["fert_eplus_y1"])
    eplus_fert_rec_y2 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_beta4"] + params["fert_eplus_y2"])
    eplus_fert_rec_y3 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_beta4"] + params["fert_eplus_y3"])
    eplus_fert_rec_y4 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_beta4"] + params["fert_eplus_y4"])
    eplus_fert_rec_y5 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_beta4"] + params["fert_eplus_y5"])
    eplus_fert_rec_y6 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_beta4"] + params["fert_eplus_y6"])
    eplus_fert_rec_y7 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_beta4"] + params["fert_eplus_y7"])
    eplus_fert_rec_y8 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_beta4"] + params["fert_eplus_y8"])
    eplus_fert_rec_y9 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_beta4"] + params["fert_eplus_y9"])
    eplus_fert_rec_y10 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_beta4"] + params["fert_eplus_y10"])
    eplus_fert_rec_y11 <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"] + params["fert_beta4"] + params["fert_eplus_y11"])
    
  # spikelet per infl
    # eminus original plants
    eminus_spike <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x))
    eminus_spike_y1 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_eminus_y1"])
    eminus_spike_y2 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_eminus_y2"])
    eminus_spike_y3 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_eminus_y3"])
    eminus_spike_y4 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_eminus_y4"])
    eminus_spike_y5 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_eminus_y5"])
    eminus_spike_y6 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_eminus_y6"])
    eminus_spike_y7 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_eminus_y7"])
    eminus_spike_y8 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_eminus_y8"])
    eminus_spike_y9 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_eminus_y9"])
    eminus_spike_y10 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_eminus_y10"])
    eminus_spike_y11 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_eminus_y11"])
    # eplus original plants
    eplus_spike <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"])
    eplus_spike_y1 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_eplus_y1"])
    eplus_spike_y2 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_eplus_y2"])
    eplus_spike_y3 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_eplus_y3"])
    eplus_spike_y4 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_eplus_y4"])
    eplus_spike_y5 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_eplus_y5"])
    eplus_spike_y6 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_eplus_y6"])
    eplus_spike_y7 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_eplus_y7"])
    eplus_spike_y8 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_eplus_y8"])
    eplus_spike_y9 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_eplus_y9"])
    eplus_spike_y10 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_eplus_y10"])
    eplus_spike_y11 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_eplus_y11"])
    # eminus recruit plants
    eminus_spike_rec <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta4"])
    eminus_spike_rec_y1 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta4"] + params["spike_eminus_y1"])
    eminus_spike_rec_y2 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta4"] + params["spike_eminus_y2"])
    eminus_spike_rec_y3 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta4"] + params["spike_eminus_y3"])
    eminus_spike_rec_y4 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta4"] + params["spike_eminus_y4"])
    eminus_spike_rec_y5 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta4"] + params["spike_eminus_y5"])
    eminus_spike_rec_y6 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta4"] + params["spike_eminus_y6"])
    eminus_spike_rec_y7 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta4"] + params["spike_eminus_y7"])
    eminus_spike_rec_y8 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta4"] + params["spike_eminus_y8"])
    eminus_spike_rec_y9 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta4"] + params["spike_eminus_y9"])
    eminus_spike_rec_y10 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta4"] + params["spike_eminus_y10"])
    eminus_spike_rec_y11 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta4"] + params["spike_eminus_y11"])
    # eplus recruit plants
    eplus_spike_rec <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_beta4"])
    eplus_spike_rec_y1 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_beta4"] + params["spike_eplus_y1"])
    eplus_spike_rec_y2 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_beta4"] + params["spike_eplus_y2"])
    eplus_spike_rec_y3 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_beta4"] + params["spike_eplus_y3"])
    eplus_spike_rec_y4 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_beta4"] + params["spike_eplus_y4"])
    eplus_spike_rec_y5 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_beta4"] + params["spike_eplus_y5"])
    eplus_spike_rec_y6 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_beta4"] + params["spike_eplus_y6"])
    eplus_spike_rec_y7 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_beta4"] + params["spike_eplus_y7"])
    eplus_spike_rec_y8 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_beta4"] + params["spike_eplus_y8"])
    eplus_spike_rec_y9 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_beta4"] + params["spike_eplus_y9"])
    eplus_spike_rec_y10 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_beta4"] + params["spike_eplus_y10"])
    eplus_spike_rec_y11 <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"] + params["spike_beta4"] + params["spike_eplus_y11"])
    
  # recruitment probability
    # eminus
    eminus_recruitment <- invlogit(params["s_to_s_beta1"])
    eminus_recruitment_y1 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_eminus_y1"])
    eminus_recruitment_y2 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_eminus_y2"])
    eminus_recruitment_y3 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_eminus_y3"])
    eminus_recruitment_y4 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_eminus_y4"])
    eminus_recruitment_y5 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_eminus_y5"])
    eminus_recruitment_y6 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_eminus_y6"])
    eminus_recruitment_y7 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_eminus_y7"])
    eminus_recruitment_y8 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_eminus_y8"])
    eminus_recruitment_y9 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_eminus_y9"])
    eminus_recruitment_y10 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_eminus_y10"])
    eminus_recruitment_y11 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_eminus_y11"])

    # eplus
    eplus_recruitment <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"])
    eplus_recruitment_y1 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"] + params["s_to_s_eplus_y1"])
    eplus_recruitment_y2 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"] + params["s_to_s_eplus_y2"])
    eplus_recruitment_y3 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"] + params["s_to_s_eplus_y3"])
    eplus_recruitment_y4 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"] + params["s_to_s_eplus_y4"])
    eplus_recruitment_y5 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"] + params["s_to_s_eplus_y5"])
    eplus_recruitment_y6 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"] + params["s_to_s_eplus_y6"])
    eplus_recruitment_y7 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"] + params["s_to_s_eplus_y7"])
    eplus_recruitment_y8 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"] + params["s_to_s_eplus_y8"])
    eplus_recruitment_y9 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"] + params["s_to_s_eplus_y9"])
    eplus_recruitment_y10 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"] + params["s_to_s_eplus_y10"])
    eplus_recruitment_y11 <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"] + params["s_to_s_eplus_y11"])

  # mean seed per spikelet
    seed.mean <- params["mu_seed"]
  
  eminus_seedlings <- eminus_flw * eminus_fert * eminus_spike * seed.mean * eminus_recruitment
  eminus_seedlings_y1 <- eminus_flw_y1 * eminus_fert_y1 * eminus_spike_y1 * seed.mean * eminus_recruitment_y1
  eminus_seedlings_y2 <- eminus_flw_y2 * eminus_fert_y2 * eminus_spike_y2 * seed.mean * eminus_recruitment_y2
  eminus_seedlings_y3 <- eminus_flw_y3 * eminus_fert_y3 * eminus_spike_y3 * seed.mean * eminus_recruitment_y3
  eminus_seedlings_y4 <- eminus_flw_y4 * eminus_fert_y4 * eminus_spike_y4 * seed.mean * eminus_recruitment_y4
  eminus_seedlings_y5 <- eminus_flw_y5 * eminus_fert_y5 * eminus_spike_y5 * seed.mean * eminus_recruitment_y5
  eminus_seedlings_y6 <- eminus_flw_y6 * eminus_fert_y6 * eminus_spike_y6 * seed.mean * eminus_recruitment_y6
  eminus_seedlings_y7 <- eminus_flw_y7 * eminus_fert_y7 * eminus_spike_y7 * seed.mean * eminus_recruitment_y7
  eminus_seedlings_y8 <- eminus_flw_y8 * eminus_fert_y8 * eminus_spike_y8 * seed.mean * eminus_recruitment_y8
  eminus_seedlings_y9 <- eminus_flw_y9 * eminus_fert_y9 * eminus_spike_y9 * seed.mean * eminus_recruitment_y9
  eminus_seedlings_y10 <- eminus_flw_y10 * eminus_fert_y10 * eminus_spike_y10 * seed.mean * eminus_recruitment_y10
  eminus_seedlings_y11 <- eminus_flw_y11 * eminus_fert_y11 * eminus_spike_y11 * seed.mean * eminus_recruitment_y11
  
  eplus_seedlings <- eplus_flw * eplus_fert * eplus_spike * seed.mean * eplus_recruitment
  eplus_seedlings_y1 <- eplus_flw_y1 * eplus_fert_y1 * eplus_spike_y1 * seed.mean * eplus_recruitment_y1
  eplus_seedlings_y2 <- eplus_flw_y2 * eplus_fert_y2 * eplus_spike_y2 * seed.mean * eplus_recruitment_y2
  eplus_seedlings_y3 <- eplus_flw_y3 * eplus_fert_y3 * eplus_spike_y3 * seed.mean * eplus_recruitment_y3
  eplus_seedlings_y4 <- eplus_flw_y4 * eplus_fert_y4 * eplus_spike_y4 * seed.mean * eplus_recruitment_y4
  eplus_seedlings_y5 <- eplus_flw_y5 * eplus_fert_y5 * eplus_spike_y5 * seed.mean * eplus_recruitment_y5
  eplus_seedlings_y6 <- eplus_flw_y6 * eplus_fert_y6 * eplus_spike_y6 * seed.mean * eplus_recruitment_y6
  eplus_seedlings_y7 <- eplus_flw_y7 * eplus_fert_y7 * eplus_spike_y7 * seed.mean * eplus_recruitment_y7
  eplus_seedlings_y8 <- eplus_flw_y8 * eplus_fert_y8 * eplus_spike_y8 * seed.mean * eplus_recruitment_y8
  eplus_seedlings_y9 <- eplus_flw_y9 * eplus_fert_y9 * eplus_spike_y9 * seed.mean * eplus_recruitment_y9
  eplus_seedlings_y10 <- eplus_flw_y10 * eplus_fert_y10 * eplus_spike_y10 * seed.mean * eplus_recruitment_y10
  eplus_seedlings_y11 <- eplus_flw_y11 * eplus_fert_y11 * eplus_spike_y11 * seed.mean * eplus_recruitment_y11
  
  eminus_seedlings_rec <- eminus_flw_rec * eminus_fert_rec * eminus_spike_rec * seed.mean * eminus_recruitment
  eminus_seedlings_rec_y1 <- eminus_flw_rec_y1 * eminus_fert_rec_y1 * eminus_spike_rec_y1 * seed.mean * eminus_recruitment_y1
  eminus_seedlings_rec_y2 <- eminus_flw_rec_y2 * eminus_fert_rec_y2 * eminus_spike_rec_y2 * seed.mean * eminus_recruitment_y2
  eminus_seedlings_rec_y3 <- eminus_flw_rec_y3 * eminus_fert_rec_y3 * eminus_spike_rec_y3 * seed.mean * eminus_recruitment_y3
  eminus_seedlings_rec_y4 <- eminus_flw_rec_y4 * eminus_fert_rec_y4 * eminus_spike_rec_y4 * seed.mean * eminus_recruitment_y4
  eminus_seedlings_rec_y5 <- eminus_flw_rec_y5 * eminus_fert_rec_y5 * eminus_spike_rec_y5 * seed.mean * eminus_recruitment_y5
  eminus_seedlings_rec_y6 <- eminus_flw_rec_y6 * eminus_fert_rec_y6 * eminus_spike_rec_y6 * seed.mean * eminus_recruitment_y6
  eminus_seedlings_rec_y7 <- eminus_flw_rec_y7 * eminus_fert_rec_y7 * eminus_spike_rec_y7 * seed.mean * eminus_recruitment_y7
  eminus_seedlings_rec_y8 <- eminus_flw_rec_y8 * eminus_fert_rec_y8 * eminus_spike_rec_y8 * seed.mean * eminus_recruitment_y8
  eminus_seedlings_rec_y9 <- eminus_flw_rec_y9 * eminus_fert_rec_y9 * eminus_spike_rec_y9 * seed.mean * eminus_recruitment_y9
  eminus_seedlings_rec_y10 <- eminus_flw_rec_y10 * eminus_fert_rec_y10 * eminus_spike_rec_y10 * seed.mean * eminus_recruitment_y10
  eminus_seedlings_rec_y11 <- eminus_flw_rec_y11 * eminus_fert_rec_y11 * eminus_spike_rec_y11 * seed.mean * eminus_recruitment_y11
  
  eplus_seedlings_rec <- eplus_flw_rec * eplus_fert_rec * eplus_spike_rec * seed.mean * eplus_recruitment
  eplus_seedlings_rec_y1 <- eplus_flw_rec_y1 * eplus_fert_rec_y1 * eplus_spike_rec_y1 * seed.mean * eplus_recruitment_y1
  eplus_seedlings_rec_y2 <- eplus_flw_rec_y2 * eplus_fert_rec_y2 * eplus_spike_rec_y2 * seed.mean * eplus_recruitment_y2
  eplus_seedlings_rec_y3 <- eplus_flw_rec_y3 * eplus_fert_rec_y3 * eplus_spike_rec_y3 * seed.mean * eplus_recruitment_y3
  eplus_seedlings_rec_y4 <- eplus_flw_rec_y4 * eplus_fert_rec_y4 * eplus_spike_rec_y4 * seed.mean * eplus_recruitment_y4
  eplus_seedlings_rec_y5 <- eplus_flw_rec_y5 * eplus_fert_rec_y5 * eplus_spike_rec_y5 * seed.mean * eplus_recruitment_y5
  eplus_seedlings_rec_y6 <- eplus_flw_rec_y6 * eplus_fert_rec_y6 * eplus_spike_rec_y6 * seed.mean * eplus_recruitment_y6
  eplus_seedlings_rec_y7 <- eplus_flw_rec_y7 * eplus_fert_rec_y7 * eplus_spike_rec_y7 * seed.mean * eplus_recruitment_y7
  eplus_seedlings_rec_y8 <- eplus_flw_rec_y8 * eplus_fert_rec_y8 * eplus_spike_rec_y8 * seed.mean * eplus_recruitment_y8
  eplus_seedlings_rec_y9 <- eplus_flw_rec_y9 * eplus_fert_rec_y9 * eplus_spike_rec_y9 * seed.mean * eplus_recruitment_y9
  eplus_seedlings_rec_y10 <- eplus_flw_rec_y10 * eplus_fert_rec_y10 * eplus_spike_rec_y10 * seed.mean * eplus_recruitment_y10
  eplus_seedlings_rec_y11 <- eplus_flw_rec_y11 * eplus_fert_rec_y11 * eplus_spike_rec_y11 * seed.mean * eplus_recruitment_y11
  
  
  return(list(eminus_seedlings=eminus_seedlings,eminus_seedlings_y1=eminus_seedlings_y1,eminus_seedlings_y2=eminus_seedlings_y2,eminus_seedlings_y3=eminus_seedlings_y3,eminus_seedlings_y4=eminus_seedlings_y4,eminus_seedlings_y5=eminus_seedlings_y5,eminus_seedlings_y6=eminus_seedlings_y6,eminus_seedlings_y7=eminus_seedlings_y7,eminus_seedlings_y8=eminus_seedlings_y8,eminus_seedlings_y9=eminus_seedlings_y9,eminus_seedlings_y10=eminus_seedlings_y10,eminus_seedlings_y11=eminus_seedlings_y11,
              eplus_seedlings=eplus_seedlings,eplus_seedlings_y1=eplus_seedlings_y1,eplus_seedlings_y2=eplus_seedlings_y2,eplus_seedlings_y3=eplus_seedlings_y3,eplus_seedlings_y4=eplus_seedlings_y4,eplus_seedlings_y5=eplus_seedlings_y5,eplus_seedlings_y6=eplus_seedlings_y6,eplus_seedlings_y7=eplus_seedlings_y7,eplus_seedlings_y8=eplus_seedlings_y8,eplus_seedlings_y9=eplus_seedlings_y9,eplus_seedlings_y10=eplus_seedlings_y10,eplus_seedlings_y11=eplus_seedlings_y11,
              eminus_seedlings_rec=eminus_seedlings_rec,eminus_seedlings_rec_y1=eminus_seedlings_rec_y1,eminus_seedlings_rec_y2=eminus_seedlings_rec_y2,eminus_seedlings_rec_y3=eminus_seedlings_rec_y3,eminus_seedlings_rec_y4=eminus_seedlings_rec_y4,eminus_seedlings_rec_y5=eminus_seedlings_rec_y5,eminus_seedlings_rec_y6=eminus_seedlings_rec_y6,eminus_seedlings_rec_y7=eminus_seedlings_rec_y7,eminus_seedlings_rec_y8=eminus_seedlings_rec_y8,eminus_seedlings_rec_y9=eminus_seedlings_rec_y9,eminus_seedlings_rec_y10=eminus_seedlings_rec_y10,eminus_seedlings_rec_y11=eminus_seedlings_rec_y11,
              eplus_seedlings_rec=eplus_seedlings_rec,eplus_seedlings_rec_y1=eplus_seedlings_rec_y1,eplus_seedlings_rec_y2=eplus_seedlings_rec_y2,eplus_seedlings_rec_y3=eplus_seedlings_rec_y3,eplus_seedlings_rec_y4=eplus_seedlings_rec_y4,eplus_seedlings_rec_y5=eplus_seedlings_rec_y5,eplus_seedlings_rec_y6=eplus_seedlings_rec_y6,eplus_seedlings_rec_y7=eplus_seedlings_rec_y7,eplus_seedlings_rec_y8=eplus_seedlings_rec_y8,eplus_seedlings_rec_y9=eplus_seedlings_rec_y9,eplus_seedlings_rec_y10=eplus_seedlings_rec_y10,eplus_seedlings_rec_y11=eplus_seedlings_rec_y11))
}
# fx(x=1, params = loar_params)

## note from Tom: we should think about adding a reproductive delay
# finally, here is the function that takes the parameter vector and assembles the matrix model from all of the pieces
bigmatrix<-function(params, endo  = FALSE, recruit = FALSE, year = FALSE){   
  matdim<-params["max_size"]## matrix dimension
  y <- 1:params["max_size"]## size (tiller number) associated with each class
  # Fertility matrix
  Fmat <- matrix(0,matdim,matdim)
  # all seedlings get dumped into top row (1-tiller)
  Fmat[1,]<-fx(x = y, params, endo, recruit, year)
  
  # Growth/survival transition matrix
  Tmat <-matrix(0,matdim,matdim)
  # Filling the transition matrix 
  Tmat<-t(outer(y,y,pxy,params, endo, recruit, year))
  
  
  # # Put it all together
  # #sum the Tmat & Fmat to get the whole matrix

  MPMmat<-Tmat + Fmat
  return(list(MPMmat = MPMmat, Tmat = Tmat, Fmat = Fmat))
}
lambda(bigmatrix(loar_params,endo = FALSE)$MPMmat)

bigmatrix<-function(params){   
  
  matdim<-params["max_size"]## matrix dimension
  y <- 1:params["max_size"]## size (tiller number) associated with each class
  # Fertility matrix
  Fmat_eminus<-matrix(0,matdim,matdim)
  Fmat_eminus_y1<-matrix(0,matdim,matdim)
  Fmat_eminus_y2<-matrix(0,matdim,matdim)
  Fmat_eminus_y3<-matrix(0,matdim,matdim)
  Fmat_eminus_y4<-matrix(0,matdim,matdim)
  Fmat_eminus_y5<-matrix(0,matdim,matdim)
  Fmat_eminus_y6<-matrix(0,matdim,matdim)
  Fmat_eminus_y7<-matrix(0,matdim,matdim)
  Fmat_eminus_y8<-matrix(0,matdim,matdim)
  Fmat_eminus_y9<-matrix(0,matdim,matdim)
  Fmat_eminus_y10<-matrix(0,matdim,matdim)
  Fmat_eminus_y11<-matrix(0,matdim,matdim)
  
  Fmat_eplus<-matrix(0,matdim,matdim)
  Fmat_eplus_y1<-matrix(0,matdim,matdim)
  Fmat_eplus_y2<-matrix(0,matdim,matdim)
  Fmat_eplus_y3<-matrix(0,matdim,matdim)
  Fmat_eplus_y4<-matrix(0,matdim,matdim)
  Fmat_eplus_y5<-matrix(0,matdim,matdim)
  Fmat_eplus_y6<-matrix(0,matdim,matdim)
  Fmat_eplus_y7<-matrix(0,matdim,matdim)
  Fmat_eplus_y8<-matrix(0,matdim,matdim)
  Fmat_eplus_y9<-matrix(0,matdim,matdim)
  Fmat_eplus_y10<-matrix(0,matdim,matdim)
  Fmat_eplus_y11<-matrix(0,matdim,matdim)
  
  Fmat_eminus_rec<-matrix(0,matdim,matdim)
  Fmat_eminus_rec_y1<-matrix(0,matdim,matdim)
  Fmat_eminus_rec_y2<-matrix(0,matdim,matdim)
  Fmat_eminus_rec_y3<-matrix(0,matdim,matdim)
  Fmat_eminus_rec_y4<-matrix(0,matdim,matdim)
  Fmat_eminus_rec_y5<-matrix(0,matdim,matdim)
  Fmat_eminus_rec_y6<-matrix(0,matdim,matdim)
  Fmat_eminus_rec_y7<-matrix(0,matdim,matdim)
  Fmat_eminus_rec_y8<-matrix(0,matdim,matdim)
  Fmat_eminus_rec_y9<-matrix(0,matdim,matdim)
  Fmat_eminus_rec_y10<-matrix(0,matdim,matdim)
  Fmat_eminus_rec_y11<-matrix(0,matdim,matdim)
  
  Fmat_eplus_rec<-matrix(0,matdim,matdim)
  Fmat_eplus_rec_y1<-matrix(0,matdim,matdim)
  Fmat_eplus_rec_y2<-matrix(0,matdim,matdim)
  Fmat_eplus_rec_y3<-matrix(0,matdim,matdim)
  Fmat_eplus_rec_y4<-matrix(0,matdim,matdim)
  Fmat_eplus_rec_y5<-matrix(0,matdim,matdim)
  Fmat_eplus_rec_y6<-matrix(0,matdim,matdim)
  Fmat_eplus_rec_y7<-matrix(0,matdim,matdim)
  Fmat_eplus_rec_y8<-matrix(0,matdim,matdim)
  Fmat_eplus_rec_y9<-matrix(0,matdim,matdim)
  Fmat_eplus_rec_y10<-matrix(0,matdim,matdim)
  Fmat_eplus_rec_y11<-matrix(0,matdim,matdim)
  
  # all seedlings get dumped into top row (1-tiller)
  Fmat_eminus[1,]<-fx(x = y, params=params)$eminus_seedlings 
  Fmat_eminus_y1[1,]<-fx(x = y, params=params)$eminus_seedlings_y1 
  Fmat_eminus_y2[1,]<-fx(x = y, params=params)$eminus_seedlings_y2 
  Fmat_eminus_y3[1,]<-fx(x = y, params=params)$eminus_seedlings_y3 
  Fmat_eminus_y4[1,]<-fx(x = y, params=params)$eminus_seedlings_y4 
  Fmat_eminus_y5[1,]<-fx(x = y, params=params)$eminus_seedlings_y5 
  Fmat_eminus_y6[1,]<-fx(x = y, params=params)$eminus_seedlings_y6 
  Fmat_eminus_y7[1,]<-fx(x = y, params=params)$eminus_seedlings_y7 
  Fmat_eminus_y8[1,]<-fx(x = y, params=params)$eminus_seedlings_y8 
  Fmat_eminus_y9[1,]<-fx(x = y, params=params)$eminus_seedlings_y9 
  Fmat_eminus_y10[1,]<-fx(x = y, params=params)$eminus_seedlings_y10 
  Fmat_eminus_y11[1,]<-fx(x = y, params=params)$eminus_seedlings_y11 
  
  Fmat_eplus[1,]<-fx(x = y, params=params)$eplus_seedlings 
  Fmat_eplus_y1[1,]<-fx(x = y, params=params)$eplus_seedlings_y1 
  Fmat_eplus_y2[1,]<-fx(x = y, params=params)$eplus_seedlings_y2 
  Fmat_eplus_y3[1,]<-fx(x = y, params=params)$eplus_seedlings_y3 
  Fmat_eplus_y4[1,]<-fx(x = y, params=params)$eplus_seedlings_y4 
  Fmat_eplus_y5[1,]<-fx(x = y, params=params)$eplus_seedlings_y5 
  Fmat_eplus_y6[1,]<-fx(x = y, params=params)$eplus_seedlings_y6 
  Fmat_eplus_y7[1,]<-fx(x = y, params=params)$eplus_seedlings_y7 
  Fmat_eplus_y8[1,]<-fx(x = y, params=params)$eplus_seedlings_y8 
  Fmat_eplus_y9[1,]<-fx(x = y, params=params)$eplus_seedlings_y9 
  Fmat_eplus_y10[1,]<-fx(x = y, params=params)$eplus_seedlings_y10 
  Fmat_eplus_y11[1,]<-fx(x = y, params=params)$eplus_seedlings_y11 
  
  Fmat_eminus_rec[1,]<-fx(x = y, params=params)$eminus_seedlings_rec 
  Fmat_eminus_rec_y1[1,]<-fx(x = y, params=params)$eminus_seedlings_rec_y1
  Fmat_eminus_rec_y2[1,]<-fx(x = y, params=params)$eminus_seedlings_rec_y2 
  Fmat_eminus_rec_y3[1,]<-fx(x = y, params=params)$eminus_seedlings_rec_y3
  Fmat_eminus_rec_y4[1,]<-fx(x = y, params=params)$eminus_seedlings_rec_y4
  Fmat_eminus_rec_y5[1,]<-fx(x = y, params=params)$eminus_seedlings_rec_y5
  Fmat_eminus_rec_y6[1,]<-fx(x = y, params=params)$eminus_seedlings_rec_y6
  Fmat_eminus_rec_y7[1,]<-fx(x = y, params=params)$eminus_seedlings_rec_y7
  Fmat_eminus_rec_y8[1,]<-fx(x = y, params=params)$eminus_seedlings_rec_y8
  Fmat_eminus_rec_y9[1,]<-fx(x = y, params=params)$eminus_seedlings_rec_y9
  Fmat_eminus_rec_y10[1,]<-fx(x = y, params=params)$eminus_seedlings_rec_y10 
  Fmat_eminus_rec_y11[1,]<-fx(x = y, params=params)$eminus_seedlings_rec_y11
  
  Fmat_eplus_rec[1,]<-fx(x = y, params=params)$eplus_seedlings_rec 
  Fmat_eplus_rec_y1[1,]<-fx(x = y, params=params)$eplus_seedlings_rec_y1
  Fmat_eplus_rec_y2[1,]<-fx(x = y, params=params)$eplus_seedlings_rec_y2 
  Fmat_eplus_rec_y3[1,]<-fx(x = y, params=params)$eplus_seedlings_rec_y3
  Fmat_eplus_rec_y4[1,]<-fx(x = y, params=params)$eplus_seedlings_rec_y4
  Fmat_eplus_rec_y5[1,]<-fx(x = y, params=params)$eplus_seedlings_rec_y5
  Fmat_eplus_rec_y6[1,]<-fx(x = y, params=params)$eplus_seedlings_rec_y6
  Fmat_eplus_rec_y7[1,]<-fx(x = y, params=params)$eplus_seedlings_rec_y7
  Fmat_eplus_rec_y8[1,]<-fx(x = y, params=params)$eplus_seedlings_rec_y8
  Fmat_eplus_rec_y9[1,]<-fx(x = y, params=params)$eplus_seedlings_rec_y9
  Fmat_eplus_rec_y10[1,]<-fx(x = y, params=params)$eplus_seedlings_rec_y10 
  Fmat_eplus_rec_y11[1,]<-fx(x = y, params=params)$eplus_seedlings_rec_y11
  
  # Growth/survival transition matrix
  Tmat_eminus<-matrix(0,matdim,matdim)
  Tmat_eminus_y1<-matrix(0,matdim,matdim)
  Tmat_eminus_y2<-matrix(0,matdim,matdim)
  Tmat_eminus_y3<-matrix(0,matdim,matdim)
  Tmat_eminus_y4<-matrix(0,matdim,matdim)
  Tmat_eminus_y5<-matrix(0,matdim,matdim)
  Tmat_eminus_y6<-matrix(0,matdim,matdim)
  Tmat_eminus_y7<-matrix(0,matdim,matdim)
  Tmat_eminus_y8<-matrix(0,matdim,matdim)
  Tmat_eminus_y9<-matrix(0,matdim,matdim)
  Tmat_eminus_y10<-matrix(0,matdim,matdim)
  Tmat_eminus_y11<-matrix(0,matdim,matdim)
  
  Tmat_eplus<-matrix(0,matdim,matdim)
  Tmat_eplus_y1<-matrix(0,matdim,matdim)
  Tmat_eplus_y2<-matrix(0,matdim,matdim)
  Tmat_eplus_y3<-matrix(0,matdim,matdim)
  Tmat_eplus_y4<-matrix(0,matdim,matdim)
  Tmat_eplus_y5<-matrix(0,matdim,matdim)
  Tmat_eplus_y6<-matrix(0,matdim,matdim)
  Tmat_eplus_y7<-matrix(0,matdim,matdim)
  Tmat_eplus_y8<-matrix(0,matdim,matdim)
  Tmat_eplus_y9<-matrix(0,matdim,matdim)
  Tmat_eplus_y10<-matrix(0,matdim,matdim)
  Tmat_eplus_y11<-matrix(0,matdim,matdim)
  
  Tmat_eminus_rec<-matrix(0,matdim,matdim)
  Tmat_eminus_rec_y1<-matrix(0,matdim,matdim)
  Tmat_eminus_rec_y2<-matrix(0,matdim,matdim)
  Tmat_eminus_rec_y3<-matrix(0,matdim,matdim)
  Tmat_eminus_rec_y4<-matrix(0,matdim,matdim)
  Tmat_eminus_rec_y5<-matrix(0,matdim,matdim)
  Tmat_eminus_rec_y6<-matrix(0,matdim,matdim)
  Tmat_eminus_rec_y7<-matrix(0,matdim,matdim)
  Tmat_eminus_rec_y8<-matrix(0,matdim,matdim)
  Tmat_eminus_rec_y9<-matrix(0,matdim,matdim)
  Tmat_eminus_rec_y10<-matrix(0,matdim,matdim)
  Tmat_eminus_rec_y11<-matrix(0,matdim,matdim)
  
  Tmat_eplus_rec<-matrix(0,matdim,matdim)
  Tmat_eplus_rec_y1<-matrix(0,matdim,matdim)
  Tmat_eplus_rec_y2<-matrix(0,matdim,matdim)
  Tmat_eplus_rec_y3<-matrix(0,matdim,matdim)
  Tmat_eplus_rec_y4<-matrix(0,matdim,matdim)
  Tmat_eplus_rec_y5<-matrix(0,matdim,matdim)
  Tmat_eplus_rec_y6<-matrix(0,matdim,matdim)
  Tmat_eplus_rec_y7<-matrix(0,matdim,matdim)
  Tmat_eplus_rec_y8<-matrix(0,matdim,matdim)
  Tmat_eplus_rec_y9<-matrix(0,matdim,matdim)
  Tmat_eplus_rec_y10<-matrix(0,matdim,matdim)
  Tmat_eplus_rec_y11<-matrix(0,matdim,matdim)
  
  Tmat_eminus<-t(outer(y,y,pxy,params=params,output="pxy_eminus"))
  Tmat_eminus_y1<-t(outer(y,y,pxy,params=params,output="pxy_eminus_y1"))
  Tmat_eminus_y2<-t(outer(y,y,pxy,params=params,output="pxy_eminus_y2"))
  Tmat_eminus_y3<-t(outer(y,y,pxy,params=params,output="pxy_eminus_y3"))
  Tmat_eminus_y4<-t(outer(y,y,pxy,params=params,output="pxy_eminus_y4"))
  Tmat_eminus_y5<-t(outer(y,y,pxy,params=params,output="pxy_eminus_y5"))
  Tmat_eminus_y6<-t(outer(y,y,pxy,params=params,output="pxy_eminus_y6"))
  Tmat_eminus_y7<-t(outer(y,y,pxy,params=params,output="pxy_eminus_y7"))
  Tmat_eminus_y8<-t(outer(y,y,pxy,params=params,output="pxy_eminus_y8"))
  Tmat_eminus_y9<-t(outer(y,y,pxy,params=params,output="pxy_eminus_y9"))
  Tmat_eminus_y10<-t(outer(y,y,pxy,params=params,output="pxy_eminus_y10"))
  Tmat_eminus_y11<-t(outer(y,y,pxy,params=params,output="pxy_eminus_y11"))
  
  Tmat_eplus<-t(outer(y,y,pxy,params=params,output="pxy_eplus"))
  Tmat_eplus_y1<-t(outer(y,y,pxy,params=params,output="pxy_eplus_y1"))
  Tmat_eplus_y2<-t(outer(y,y,pxy,params=params,output="pxy_eplus_y2"))
  Tmat_eplus_y3<-t(outer(y,y,pxy,params=params,output="pxy_eplus_y3"))
  Tmat_eplus_y4<-t(outer(y,y,pxy,params=params,output="pxy_eplus_y4"))
  Tmat_eplus_y5<-t(outer(y,y,pxy,params=params,output="pxy_eplus_y5"))
  Tmat_eplus_y6<-t(outer(y,y,pxy,params=params,output="pxy_eplus_y6"))
  Tmat_eplus_y7<-t(outer(y,y,pxy,params=params,output="pxy_eplus_y7"))
  Tmat_eplus_y8<-t(outer(y,y,pxy,params=params,output="pxy_eplus_y8"))
  Tmat_eplus_y9<-t(outer(y,y,pxy,params=params,output="pxy_eplus_y9"))
  Tmat_eplus_y10<-t(outer(y,y,pxy,params=params,output="pxy_eplus_y10"))
  Tmat_eplus_y11<-t(outer(y,y,pxy,params=params,output="pxy_eplus_y11"))
  
  Tmat_eminus_rec<-t(outer(y,y,pxy,params=params,output="pxy_eminus_rec"))
  Tmat_eminus_rec_y1<-t(outer(y,y,pxy,params=params,output="pxy_eminus_rec_y1"))
  Tmat_eminus_rec_y2<-t(outer(y,y,pxy,params=params,output="pxy_eminus_rec_y2"))
  Tmat_eminus_rec_y3<-t(outer(y,y,pxy,params=params,output="pxy_eminus_rec_y3"))
  Tmat_eminus_rec_y4<-t(outer(y,y,pxy,params=params,output="pxy_eminus_rec_y4"))
  Tmat_eminus_rec_y5<-t(outer(y,y,pxy,params=params,output="pxy_eminus_rec_y5"))
  Tmat_eminus_rec_y6<-t(outer(y,y,pxy,params=params,output="pxy_eminus_rec_y6"))
  Tmat_eminus_rec_y7<-t(outer(y,y,pxy,params=params,output="pxy_eminus_rec_y7"))
  Tmat_eminus_rec_y8<-t(outer(y,y,pxy,params=params,output="pxy_eminus_rec_y8"))
  Tmat_eminus_rec_y9<-t(outer(y,y,pxy,params=params,output="pxy_eminus_rec_y9"))
  Tmat_eminus_rec_y10<-t(outer(y,y,pxy,params=params,output="pxy_eminus_rec_y10"))
  Tmat_eminus_rec_y11<-t(outer(y,y,pxy,params=params,output="pxy_eminus_rec_y11"))
  
  Tmat_eplus_rec<-t(outer(y,y,pxy,params=params,output="pxy_eplus_rec"))
  Tmat_eplus_rec_y1<-t(outer(y,y,pxy,params=params,output="pxy_eplus_rec_y1"))
  Tmat_eplus_rec_y2<-t(outer(y,y,pxy,params=params,output="pxy_eplus_rec_y2"))
  Tmat_eplus_rec_y3<-t(outer(y,y,pxy,params=params,output="pxy_eplus_rec_y3"))
  Tmat_eplus_rec_y4<-t(outer(y,y,pxy,params=params,output="pxy_eplus_rec_y4"))
  Tmat_eplus_rec_y5<-t(outer(y,y,pxy,params=params,output="pxy_eplus_rec_y5"))
  Tmat_eplus_rec_y6<-t(outer(y,y,pxy,params=params,output="pxy_eplus_rec_y6"))
  Tmat_eplus_rec_y7<-t(outer(y,y,pxy,params=params,output="pxy_eplus_rec_y7"))
  Tmat_eplus_rec_y8<-t(outer(y,y,pxy,params=params,output="pxy_eplus_rec_y8"))
  Tmat_eplus_rec_y9<-t(outer(y,y,pxy,params=params,output="pxy_eplus_rec_y9"))
  Tmat_eplus_rec_y10<-t(outer(y,y,pxy,params=params,output="pxy_eplus_rec_y10"))
  Tmat_eplus_rec_y11<-t(outer(y,y,pxy,params=params,output="pxy_eplus_rec_y11"))
  
  # Put it all together
  #sum the Tmat & Fmat to get the whole matrix
  # eminus original
  MPMmat_eminus<-Tmat_eminus+Fmat_eminus 
  MPMmat_eminus_y1<-Tmat_eminus_y1+Fmat_eminus_y1 
  MPMmat_eminus_y2<-Tmat_eminus_y2+Fmat_eminus_y2
  MPMmat_eminus_y3<-Tmat_eminus_y3+Fmat_eminus_y3 
  MPMmat_eminus_y4<-Tmat_eminus_y4+Fmat_eminus_y4 
  MPMmat_eminus_y5<-Tmat_eminus_y5+Fmat_eminus_y5 
  MPMmat_eminus_y6<-Tmat_eminus_y6+Fmat_eminus_y6 
  MPMmat_eminus_y7<-Tmat_eminus_y7+Fmat_eminus_y7 
  MPMmat_eminus_y8<-Tmat_eminus_y8+Fmat_eminus_y8 
  MPMmat_eminus_y9<-Tmat_eminus_y9+Fmat_eminus_y9 
  MPMmat_eminus_y10<-Tmat_eminus_y10+Fmat_eminus_y10 
  MPMmat_eminus_y11<-Tmat_eminus_y11+Fmat_eminus_y11
  
  # eplus original
  MPMmat_eplus<-Tmat_eplus+Fmat_eplus 
  MPMmat_eplus_y1<-Tmat_eplus_y1+Fmat_eplus_y1 
  MPMmat_eplus_y2<-Tmat_eplus_y2+Fmat_eplus_y2
  MPMmat_eplus_y3<-Tmat_eplus_y3+Fmat_eplus_y3 
  MPMmat_eplus_y4<-Tmat_eplus_y4+Fmat_eplus_y4 
  MPMmat_eplus_y5<-Tmat_eplus_y5+Fmat_eplus_y5 
  MPMmat_eplus_y6<-Tmat_eplus_y6+Fmat_eplus_y6 
  MPMmat_eplus_y7<-Tmat_eplus_y7+Fmat_eplus_y7 
  MPMmat_eplus_y8<-Tmat_eplus_y8+Fmat_eplus_y8 
  MPMmat_eplus_y9<-Tmat_eplus_y9+Fmat_eplus_y9 
  MPMmat_eplus_y10<-Tmat_eplus_y10+Fmat_eplus_y10 
  MPMmat_eplus_y11<-Tmat_eplus_y11+Fmat_eplus_y11
  
  # eminus recruit
  MPMmat_eminus_rec<-Tmat_eminus_rec+Fmat_eminus_rec 
  MPMmat_eminus_rec_y1<-Tmat_eminus_rec_y1+Fmat_eminus_rec_y1 
  MPMmat_eminus_rec_y2<-Tmat_eminus_rec_y2+Fmat_eminus_rec_y2
  MPMmat_eminus_rec_y3<-Tmat_eminus_rec_y3+Fmat_eminus_rec_y3 
  MPMmat_eminus_rec_y4<-Tmat_eminus_rec_y4+Fmat_eminus_rec_y4 
  MPMmat_eminus_rec_y5<-Tmat_eminus_rec_y5+Fmat_eminus_rec_y5 
  MPMmat_eminus_rec_y6<-Tmat_eminus_rec_y6+Fmat_eminus_rec_y6 
  MPMmat_eminus_rec_y7<-Tmat_eminus_rec_y7+Fmat_eminus_rec_y7 
  MPMmat_eminus_rec_y8<-Tmat_eminus_rec_y8+Fmat_eminus_rec_y8 
  MPMmat_eminus_rec_y9<-Tmat_eminus_rec_y9+Fmat_eminus_rec_y9 
  MPMmat_eminus_rec_y10<-Tmat_eminus_rec_y10+Fmat_eminus_rec_y10 
  MPMmat_eminus_rec_y11<-Tmat_eminus_rec_y11+Fmat_eminus_rec_y11
  
  # eplus recruit
  MPMmat_eplus_rec<-Tmat_eplus+Fmat_eplus_rec 
  MPMmat_eplus_rec_y1<-Tmat_eplus_rec_y1+Fmat_eplus_rec_y1 
  MPMmat_eplus_rec_y2<-Tmat_eplus_rec_y2+Fmat_eplus_rec_y2
  MPMmat_eplus_rec_y3<-Tmat_eplus_rec_y3+Fmat_eplus_rec_y3 
  MPMmat_eplus_rec_y4<-Tmat_eplus_rec_y4+Fmat_eplus_rec_y4 
  MPMmat_eplus_rec_y5<-Tmat_eplus_rec_y5+Fmat_eplus_rec_y5 
  MPMmat_eplus_rec_y6<-Tmat_eplus_rec_y6+Fmat_eplus_rec_y6 
  MPMmat_eplus_rec_y7<-Tmat_eplus_rec_y7+Fmat_eplus_rec_y7 
  MPMmat_eplus_rec_y8<-Tmat_eplus_rec_y8+Fmat_eplus_rec_y8 
  MPMmat_eplus_rec_y9<-Tmat_eplus_rec_y9+Fmat_eplus_rec_y9 
  MPMmat_eplus_rec_y10<-Tmat_eplus_rec_y10+Fmat_eplus_rec_y10 
  MPMmat_eplus_rec_y11<-Tmat_eplus_rec_y11+Fmat_eplus_rec_y11
  
  
  return(list(MPMmat_eminus=MPMmat_eminus,MPMmat_eminus_y1=MPMmat_eminus_y1,MPMmat_eminus_y2=MPMmat_eminus_y2,MPMmat_eminus_y3=MPMmat_eminus_y3,MPMmat_eminus_y4=MPMmat_eminus_y4,MPMmat_eminus_y5=MPMmat_eminus_y5,MPMmat_eminus_y6=MPMmat_eminus_y6,MPMmat_eminus_y7=MPMmat_eminus_y7,MPMmat_eminus_y8=MPMmat_eminus_y8,MPMmat_eminus_y9=MPMmat_eminus_y9,MPMmat_eminus_y10=MPMmat_eminus_y10,MPMmat_eminus_y11=MPMmat_eminus_y11,
              MPMmat_eplus=MPMmat_eplus,MPMmat_eplus_y1=MPMmat_eplus_y1,MPMmat_eplus_y2=MPMmat_eplus_y2,MPMmat_eplus_y3=MPMmat_eplus_y3,MPMmat_eplus_y4=MPMmat_eplus_y4,MPMmat_eplus_y5=MPMmat_eplus_y5,MPMmat_eplus_y6=MPMmat_eplus_y6,MPMmat_eplus_y7=MPMmat_eplus_y7,MPMmat_eplus_y8=MPMmat_eplus_y8,MPMmat_eplus_y9=MPMmat_eplus_y9,MPMmat_eplus_y10=MPMmat_eplus_y10,MPMmat_eplus_y11=MPMmat_eplus_y11,
              MPMmat_eminus_rec=MPMmat_eminus_rec,MPMmat_eminus_rec_y1=MPMmat_eminus_rec_y1,MPMmat_eminus_rec_y2=MPMmat_eminus_rec_y2,MPMmat_eminus_rec_y3=MPMmat_eminus_rec_y3,MPMmat_eminus_rec_y4=MPMmat_eminus_rec_y4,MPMmat_eminus_rec_y5=MPMmat_eminus_rec_y5,MPMmat_eminus_rec_y6=MPMmat_eminus_rec_y6,MPMmat_eminus_rec_y7=MPMmat_eminus_rec_y7,MPMmat_eminus_rec_y8=MPMmat_eminus_rec_y8,MPMmat_eminus_rec_y9=MPMmat_eminus_rec_y9,MPMmat_eminus_rec_y10=MPMmat_eminus_rec_y10,MPMmat_eminus_rec_y11=MPMmat_eminus_rec_y11,
              MPMmat_eplus_rec=MPMmat_eplus_rec,MPMmat_eplus_rec_y1=MPMmat_eplus_rec_y1,MPMmat_eplus_rec_y2=MPMmat_eplus_rec_y2,MPMmat_eplus_rec_y3=MPMmat_eplus_rec_y3,MPMmat_eplus_rec_y4=MPMmat_eplus_rec_y4,MPMmat_eplus_rec_y5=MPMmat_eplus_rec_y5,MPMmat_eplus_rec_y6=MPMmat_eplus_rec_y6,MPMmat_eplus_rec_y7=MPMmat_eplus_rec_y7,MPMmat_eplus_rec_y8=MPMmat_eplus_rec_y8,MPMmat_eplus_rec_y9=MPMmat_eplus_rec_y9,MPMmat_eplus_rec_y10=MPMmat_eplus_rec_y10,MPMmat_eplus_rec_y11=MPMmat_eplus_rec_y11,
              Fmat_eminus=Fmat_eminus,Fmat_eminus_y1=Fmat_eminus_y1,Fmat_eminus_y2=Fmat_eminus_y2,Fmat_eminus_y3=Fmat_eminus_y3,Fmat_eminus_y4=Fmat_eminus_y4,Fmat_eminus_y5=Fmat_eminus_y5,Fmat_eminus_y6=Fmat_eminus_y6,Fmat_eminus_y7=Fmat_eminus_y7,Fmat_eminus_y8=Fmat_eminus_y8,Fmat_eminus_y9=Fmat_eminus_y9,Fmat_eminus_y10=Fmat_eminus_y10,Fmat_eminus_y11=Fmat_eminus_y11,
              Fmat_eplus=Fmat_eplus,Fmat_eplus_y1=Fmat_eplus_y1,Fmat_eplus_y2=Fmat_eplus_y2,Fmat_eplus_y3=Fmat_eplus_y3,Fmat_eplus_y4=Fmat_eplus_y4,Fmat_eplus_y5=Fmat_eplus_y5,Fmat_eplus_y6=Fmat_eplus_y6,Fmat_eplus_y7=Fmat_eplus_y7,Fmat_eplus_y8=Fmat_eplus_y8,Fmat_eplus_y9=Fmat_eplus_y9,Fmat_eplus_y10=Fmat_eplus_y10,Fmat_eplus_y11=Fmat_eplus_y11,
              Fmat_eminus_rec=Fmat_eminus_rec,Fmat_eminus_rec_y1=Fmat_eminus_rec_y1,Fmat_eminus_rec_y2=Fmat_eminus_rec_y2,Fmat_eminus_rec_y3=Fmat_eminus_rec_y3,Fmat_eminus_rec_y4=Fmat_eminus_rec_y4,Fmat_eminus_rec_y5=Fmat_eminus_rec_y5,Fmat_eminus_rec_y6=Fmat_eminus_rec_y6,Fmat_eminus_rec_y7=Fmat_eminus_rec_y7,Fmat_eminus_rec_y8=Fmat_eminus_rec_y8,Fmat_eminus_rec_y9=Fmat_eminus_rec_y9,Fmat_eminus_rec_y10=Fmat_eminus_rec_y10,Fmat_eminus_rec_y11=Fmat_eminus_rec_y11,
              Fmat_eplus_rec=Fmat_eplus_rec,Fmat_eplus_rec_y1=Fmat_eplus_rec_y1,Fmat_eplus_rec_y2=Fmat_eplus_rec_y2,Fmat_eplus_rec_y3=Fmat_eplus_rec_y3,Fmat_eplus_rec_y4=Fmat_eplus_rec_y4,Fmat_eplus_rec_y5=Fmat_eplus_rec_y5,Fmat_eplus_rec_y6=Fmat_eplus_rec_y6,Fmat_eplus_rec_y7=Fmat_eplus_rec_y7,Fmat_eplus_rec_y8=Fmat_eplus_rec_y8,Fmat_eplus_rec_y9=Fmat_eplus_rec_y9,Fmat_eplus_rec_y10=Fmat_eplus_rec_y10,Fmat_eplus_rec_y11=Fmat_eplus_rec_y11,
              Tmat_eminus=Tmat_eminus,Tmat_eminus_y1=Tmat_eminus_y1,Tmat_eminus_y2=Tmat_eminus_y2,Tmat_eminus_y3=Tmat_eminus_y3,Tmat_eminus_y4=Tmat_eminus_y4,Tmat_eminus_y5=Tmat_eminus_y5,Tmat_eminus_y6=Tmat_eminus_y6,Tmat_eminus_y7=Tmat_eminus_y7,Tmat_eminus_y8=Tmat_eminus_y8,Tmat_eminus_y9=Tmat_eminus_y9,Tmat_eminus_y10=Tmat_eminus_y10,Tmat_eminus_y11=Tmat_eminus_y11,
              Tmat_eplus=Tmat_eplus,Tmat_eplus_y1=Tmat_eplus_y1,Tmat_eplus_y2=Tmat_eplus_y2,Tmat_eplus_y3=Tmat_eplus_y3,Tmat_eplus_y4=Tmat_eplus_y4,Tmat_eplus_y5=Tmat_eplus_y5,Tmat_eplus_y6=Tmat_eplus_y6,Tmat_eplus_y7=Tmat_eplus_y7,Tmat_eplus_y8=Tmat_eplus_y8,Tmat_eplus_y9=Tmat_eplus_y9,Tmat_eplus_y10=Tmat_eplus_y10,Tmat_eplus_y11=Tmat_eplus_y11,
              Tmat_eminus_rec=Tmat_eminus_rec,Tmat_eminus_rec_y1=Tmat_eminus_rec_y1,Tmat_eminus_rec_y2=Tmat_eminus_rec_y2,Tmat_eminus_rec_y3=Tmat_eminus_rec_y3,Tmat_eminus_rec_y4=Tmat_eminus_rec_y4,Tmat_eminus_rec_y5=Tmat_eminus_rec_y5,Tmat_eminus_rec_y6=Tmat_eminus_rec_y6,Tmat_eminus_rec_y7=Tmat_eminus_rec_y7,Tmat_eminus_rec_y8=Tmat_eminus_rec_y8,Tmat_eminus_rec_y9=Tmat_eminus_rec_y9,Tmat_eminus_rec_y10=Tmat_eminus_rec_y10,Tmat_eminus_rec_y11=Tmat_eminus_rec_y11,
              Tmat_eplus_rec=Tmat_eplus_rec,Tmat_eplus_rec_y1=Tmat_eplus_rec_y1,Tmat_eplus_rec_y2=Tmat_eplus_rec_y2,Tmat_eplus_rec_y3=Tmat_eplus_rec_y3,Tmat_eplus_rec_y4=Tmat_eplus_rec_y4,Tmat_eplus_rec_y5=Tmat_eplus_rec_y5,Tmat_eplus_rec_y6=Tmat_eplus_rec_y6,Tmat_eplus_rec_y7=Tmat_eplus_rec_y7,Tmat_eplus_rec_y8=Tmat_eplus_rec_y8,Tmat_eplus_rec_y9=Tmat_eplus_rec_y9,Tmat_eplus_rec_y10=Tmat_eplus_rec_y10,Tmat_eplus_rec_y11=Tmat_eplus_rec_y11))
}

## population growth rate (eigenalaysis of the projection matrix)
# matrix <- bigmatrix(loar_params)
# calculating the matrices
agpe_mat <- bigmatrix(agpe_params)
elri_mat <- bigmatrix(elri_params)
elvi_mat <- bigmatrix(elvi_params)
fesu_mat <- bigmatrix(fesu_params)
loar_mat <- bigmatrix(loar_params)
poal_mat <- bigmatrix(poal_params)
posy_mat <- bigmatrix(posy_params)


# avg lambdas
lambda(agpe_mat$MPMmat_eminus)
lambda(agpe_mat$MPMmat_eplus)
lambda(fesu_mat$MPMmat_eminus)
lambda(fesu_mat$MPMmat_eplus)
lambda(loar_mat$MPMmat_eminus)
lambda(loar_mat$MPMmat_eplus)
lambda(poal_mat$MPMmat_eminus)
lambda(poal_mat$MPMmat_eplus)
lambda(posy_mat$MPMmat_eminus)
lambda(posy_mat$MPMmat_eplus)

# saving the yearly growth rates to vectors 
yearly <- function(mat){
  yearly_eminus <- c()
  yearly_eplus <- c()
  yearly_eminus_rec <- c()
  yearly_eplus_rec <- c()
  year_t1 <- c()
  # eminus original
  yearly_eminus[1] <- lambda(mat$MPMmat_eminus_y1)
  yearly_eminus[2] <- lambda(mat$MPMmat_eminus_y2)
  yearly_eminus[3] <- lambda(mat$MPMmat_eminus_y3)
  yearly_eminus[4] <- lambda(mat$MPMmat_eminus_y4)
  yearly_eminus[5] <- lambda(mat$MPMmat_eminus_y5)
  yearly_eminus[6] <- lambda(mat$MPMmat_eminus_y6)
  yearly_eminus[7] <- lambda(mat$MPMmat_eminus_y7)
  yearly_eminus[8] <- lambda(mat$MPMmat_eminus_y8)
  yearly_eminus[9] <- lambda(mat$MPMmat_eminus_y9)
  yearly_eminus[10] <- lambda(mat$MPMmat_eminus_y10)
  yearly_eminus[11] <- lambda(mat$MPMmat_eminus_y11)
  # eplus original
  yearly_eplus[1] <- lambda(mat$MPMmat_eplus_y1)
  yearly_eplus[2] <- lambda(mat$MPMmat_eplus_y2)
  yearly_eplus[3] <- lambda(mat$MPMmat_eplus_y3)
  yearly_eplus[4] <- lambda(mat$MPMmat_eplus_y4)
  yearly_eplus[5] <- lambda(mat$MPMmat_eplus_y5)
  yearly_eplus[6] <- lambda(mat$MPMmat_eplus_y6)
  yearly_eplus[7] <- lambda(mat$MPMmat_eplus_y7)
  yearly_eplus[8] <- lambda(mat$MPMmat_eplus_y8)
  yearly_eplus[9] <- lambda(mat$MPMmat_eplus_y9)
  yearly_eplus[10] <- lambda(mat$MPMmat_eplus_y10)
  yearly_eplus[11] <- lambda(mat$MPMmat_eplus_y11)
  # eminus recruit
  yearly_eminus_rec[1] <- lambda(mat$MPMmat_eminus_rec_y1)
  yearly_eminus_rec[2] <- lambda(mat$MPMmat_eminus_rec_y2)
  yearly_eminus_rec[3] <- lambda(mat$MPMmat_eminus_rec_y3)
  yearly_eminus_rec[4] <- lambda(mat$MPMmat_eminus_rec_y4)
  yearly_eminus_rec[5] <- lambda(mat$MPMmat_eminus_rec_y5)
  yearly_eminus_rec[6] <- lambda(mat$MPMmat_eminus_rec_y6)
  yearly_eminus_rec[7] <- lambda(mat$MPMmat_eminus_rec_y7)
  yearly_eminus_rec[8] <- lambda(mat$MPMmat_eminus_rec_y8)
  yearly_eminus_rec[9] <- lambda(mat$MPMmat_eminus_rec_y9)
  yearly_eminus_rec[10] <- lambda(mat$MPMmat_eminus_rec_y10)
  yearly_eminus_rec[11] <- lambda(mat$MPMmat_eminus_rec_y11)
  # eplus recruit
  yearly_eplus_rec[1] <- lambda(mat$MPMmat_eplus_rec_y1)
  yearly_eplus_rec[2] <- lambda(mat$MPMmat_eplus_rec_y2)
  yearly_eplus_rec[3] <- lambda(mat$MPMmat_eplus_rec_y3)
  yearly_eplus_rec[4] <- lambda(mat$MPMmat_eplus_rec_y4)
  yearly_eplus_rec[5] <- lambda(mat$MPMmat_eplus_rec_y5)
  yearly_eplus_rec[6] <- lambda(mat$MPMmat_eplus_rec_y6)
  yearly_eplus_rec[7] <- lambda(mat$MPMmat_eplus_rec_y7)
  yearly_eplus_rec[8] <- lambda(mat$MPMmat_eplus_rec_y8)
  yearly_eplus_rec[9] <- lambda(mat$MPMmat_eplus_rec_y9)
  yearly_eplus_rec[10] <- lambda(mat$MPMmat_eplus_rec_y10)
  yearly_eplus_rec[11] <- lambda(mat$MPMmat_eplus_rec_y11)
  # year date
  year_t1[1] <- 2008
  year_t1[2] <- 2009
  year_t1[3] <- 2010
  year_t1[4] <- 2011
  year_t1[5] <- 2012
  year_t1[6] <- 2013
  year_t1[7] <- 2014
  year_t1[8] <- 2015
  year_t1[9] <- 2016
  year_t1[10] <- 2017
  year_t1[11] <- 2018

  lambdas <- as_tibble(cbind(yearly_eminus, yearly_eplus, yearly_eminus_rec, yearly_eplus_rec, year_t1))
  
  return(lambdas)
}
# yearly(agpe_mat)


# some starter histograms
# eminus = orange, eplus = purple
ggplot(data = yearly(agpe_mat))+
  geom_histogram(aes(yearly_eminus_rec), bins = 8, fill = "#ff7f00", alpha = .6) +
  geom_histogram(aes(yearly_eplus_rec), bins = 8, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") +labs(title = "AGPE") + theme_classic(base_size = 20)

ggplot(data = yearly(elri_mat))+
  geom_histogram(aes(yearly_eminus_rec), bins = 8, fill = "#ff7f00", alpha = .6) +
  geom_histogram(aes(yearly_eplus_rec), bins = 8, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") +labs(title = "ELRI") + theme_classic(base_size = 20)

ggplot(data = yearly(elvi_mat))+
  geom_histogram(aes(yearly_eminus_rec), bins = 8, fill = "#ff7f00", alpha = .6) +
  geom_histogram(aes(yearly_eplus_rec), bins = 8, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") +labs(title = "ELVI") + theme_classic(base_size = 20)

ggplot(data = yearly(fesu_mat))+
  geom_histogram(aes(yearly_eminus_rec), bins = 8, fill = "#ff7f00", alpha = .6) +
  geom_histogram(aes(yearly_eplus_rec), bins = 8, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") +labs(title = "FESU") + theme_classic(base_size = 20)

ggplot(data = yearly(loar_mat))+
  geom_histogram(aes(yearly_eminus_rec), bins = 8, fill = "#ff7f00", alpha = .6) +
  geom_histogram(aes(yearly_eplus_rec), bins = 8, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") +labs(title = "LOAR") + theme_classic(base_size = 20)

ggplot(data = yearly(poal_mat))+
  geom_histogram(aes(yearly_eminus_rec), bins = 8, fill = "#ff7f00", alpha = .6) +
  geom_histogram(aes(yearly_eplus_rec), bins = 8, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") +labs(title = "POAL") + theme_classic(base_size = 20)

ggplot(data = yearly(posy_mat))+
  geom_histogram(aes(yearly_eminus_rec), bins = 8, fill = "#ff7f00", alpha = .6) +
  geom_histogram(aes(yearly_eplus_rec), bins = 8, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") +labs(title = "POSY") + theme_classic(base_size = 20) 

# option of a density plot
ggplot(data = yearly(posy_mat))+
  geom_density(aes(yearly_eminus_rec), bw = .1, fill = "#ff7f00", alpha = .6) +
  geom_density(aes(yearly_eplus_rec), bw = .1, fill = "#6a3d9a", alpha = .6) + 
  labs(x = "population growth rate") + labs(title = "POSY") +theme_classic(base_size = 20)

# some starter time series
ggplot(data = yearly(agpe_mat))+
  geom_point(aes(x = year_t1, y = yearly_eminus_rec),size = 4, col = "#ff7f00") +
  geom_point(aes(x = year_t1, y = yearly_eplus_rec),size = 4, col = "#6a3d9a") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(y = "population growth rate", x = "year_t1") + ggtitle("AGPE") + theme_classic(base_size = 20) + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))
ggplot(data = yearly(elri_mat))+
  geom_point(aes(x = year_t1, y = yearly_eminus_rec),size = 4, col = "#ff7f00") +
  geom_point(aes(x = year_t1, y = yearly_eplus_rec),size = 4, col = "#6a3d9a") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(y = "population growth rate", x = "year_t1") + ggtitle("ELRI") + theme_classic(base_size = 20) + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))
ggplot(data = yearly(elvi_mat))+
  geom_point(aes(x = year_t1, y = yearly_eminus_rec),size = 4, col = "#ff7f00") +
  geom_point(aes(x = year_t1, y = yearly_eplus_rec),size = 4, col = "#6a3d9a") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(y = "population growth rate", x = "year_t1") + ggtitle("ELVI") + theme_classic(base_size = 20) + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))

ggplot(data = yearly(fesu_mat))+
  geom_point(aes(x = year_t1, y = yearly_eminus_rec),size = 4, col = "#ff7f00") +
  geom_point(aes(x = year_t1, y = yearly_eplus_rec),size = 4, col = "#6a3d9a") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(y = "population growth rate", x = "year_t1") + ggtitle("FESU") + theme_classic(base_size = 20) + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))

ggplot(data = yearly(loar_mat))+
  geom_point(aes(x = year_t1, y = yearly_eminus_rec),size = 4, col = "#ff7f00") +
  geom_point(aes(x = year_t1, y = yearly_eplus_rec),size = 4, col = "#6a3d9a") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(y = "population growth rate", x = "year_t1") + ggtitle("LOAR") + theme_classic(base_size = 20) + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))

ggplot(data = yearly(poal_mat))+
  geom_point(aes(x = year_t1, y = yearly_eminus_rec),size = 4, col = "#ff7f00") +
  geom_point(aes(x = year_t1, y = yearly_eplus_rec),size = 4, col = "#6a3d9a") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(y = "population growth rate", x = "year_t1") + ggtitle("POAL") + theme_classic(base_size = 20) + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))

ggplot(data = yearly(posy_mat))+
  geom_point(aes(x = year_t1, y = yearly_eminus_rec),size = 4, col = "#ff7f00") +
  geom_point(aes(x = year_t1, y = yearly_eplus_rec),size = 4, col = "#6a3d9a") +
  geom_hline(yintercept = 1, linetype = "dotted") +
  labs(y = "population growth rate", x = "year_t1") + ggtitle("POSY") + theme_classic(base_size = 20) + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))



# I'm gonna read in env. data here
# data downloaded from PRISM, daily ppt, tmin, tmean, tmax for GPS point 39.235900000000,-86.218100000000 from Jan 1, 2006 to Jan 1, 2019
climate <- read_csv(file = "~/Dropbox/EndodemogData/PRISMClimateData_BrownCo.csv") %>% 
  mutate(year = year(Date), month = month(Date), day = day(Date)) %>% 
  rename(ppt = `ppt (mm)`, tmean = `tmean (degrees C)`) %>% 
  mutate(site_lat = 39.235900000000, site_long = -86.218100000000)
  

AGPE_climate <- climate %>% 
  mutate(census_month = 9, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>% 
  filter(climate_year != 2006) %>% 
  group_by(climate_year) %>% 
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean))
AGPE <- yearly(agpe_mat) %>% 
  mutate(species = "AGPE") %>% 
  merge(AGPE_climate, by.x = c("year_t1"), by.y = c("climate_year"))
ELRI_climate <- climate %>%
  mutate(census_month = 7, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>%
  filter(climate_year != 2006) %>%
  group_by(climate_year) %>%
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean))
ELRI <- yearly(elri_mat) %>%
  mutate(species = "ELRI") %>%
  merge(ELRI_climate, by.x = c("year_t1"), by.y = c("climate_year"))
ELVI_climate <- climate %>%
  mutate(census_month = 7, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>%
  filter(climate_year != 2006) %>%
  group_by(climate_year) %>%
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean))
ELVI <- yearly(elvi_mat) %>%
  mutate(species = "ELVI") %>%
  merge(ELVI_climate, by.x = c("year_t1"), by.y = c("climate_year"))
FESU_climate <- climate %>% 
  mutate(census_month = 6, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>% 
  filter(climate_year != 2006) %>% 
  group_by(climate_year) %>% 
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean))
FESU <- yearly(fesu_mat) %>% 
  mutate(species = "FESU") %>% 
  merge(FESU_climate, by.x = c("year_t1"), by.y = c("climate_year"))
LOAR_climate <- climate %>% 
  mutate(census_month = 7, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>% 
  filter(climate_year != 2006) %>% 
  group_by(climate_year) %>% 
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean))
LOAR <- yearly(loar_mat) %>% 
  mutate(species = "LOAR") %>% 
  merge(LOAR_climate, by.x = c("year_t1"), by.y = c("climate_year"))
POAL_climate <- climate %>% 
  mutate(census_month = 5, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>% 
  filter(climate_year != 2006) %>% 
  group_by(climate_year) %>% 
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean))
POAL <- yearly(poal_mat) %>% 
  mutate(species = "POAL") %>% 
  merge(POAL_climate, by.x = c("year_t1"), by.y = c("climate_year"))
POSY_climate <- climate %>% 
  mutate(census_month = 5, climate_year = as.numeric(ifelse(month >=census_month, year+1, year))) %>% 
  filter(climate_year != 2006) %>% 
  group_by(climate_year) %>% 
  summarize('Cumulative PPT (mm)' = sum(ppt),
            'Mean Temp. (C˚)' = mean(tmean))
POSY <- yearly(posy_mat) %>% 
  mutate(species = "POSY") %>% 
  merge(POSY_climate, by.x = c("year_t1"), by.y = c("climate_year"))

endo_climate <-  rbind(AGPE,ELRI,ELVI,FESU,LOAR,POAL,POSY) %>% 
  melt(id.vars = c("species", "year_t1", 'Cumulative PPT (mm)', 'Mean Temp. (C˚)')) %>% 
  filter(variable != "yearly_eminus", variable != "yearly_eplus")


ggplot(data = endo_climate) +
  geom_smooth(aes(x = `Cumulative PPT (mm)`, y = value),color = "grey", method = "glm", se = FALSE) +
  geom_point(aes(x = `Cumulative PPT (mm)`, y = value, color = variable, shape = species)) +
  facet_grid(rows = vars(variable), cols = vars(species)) +
  theme_classic() + scale_color_manual(values = c("#ff7f00", "#6a3d9a")) +scale_shape_manual(values=c(0:6))+ labs(y = "Population Growth") + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))


ggplot(data = endo_climate) +
  geom_smooth(aes(x = `Cumulative PPT (mm)`, y = value), color = "grey", method = "glm", se = FALSE) +
  geom_point(aes(x = `Cumulative PPT (mm)`, y = value, color = variable, shape = species)) +
  facet_grid(cols = vars(species)) +
  theme_classic() + scale_color_manual(values = c("#ff7f00", "#6a3d9a"))+ scale_shape_manual(values=c(0:6))+ labs(y = "Population Growth") + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))


ggplot(data = endo_climate) +
  geom_smooth(aes(x = `Mean Temp. (C˚)`, y = value),color = "grey", method = "glm", se = FALSE) +
  geom_point(aes(x = `Mean Temp. (C˚)`, y = value, color = variable, shape = species)) +
  facet_grid(rows = vars(variable), cols = vars(species)) +
  theme_classic() + scale_color_manual(values = c("#ff7f00", "#6a3d9a")) +scale_shape_manual(values=c(0:6))+ labs(y = "Population Growth") + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))


ggplot(data = endo_climate) +
  geom_smooth(aes(x = `Mean Temp. (C˚)`, y = value), color = "grey", method = "glm", se = FALSE) +
  geom_point(aes(x = `Mean Temp. (C˚)`, y = value, color = variable, shape = species)) +
  facet_grid(cols = vars(species)) +
  theme_classic() + scale_color_manual(values = c("#ff7f00", "#6a3d9a"))+ scale_shape_manual(values=c(0:6))+ labs(y = "Population Growth") + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))

endo_climate %>% 
  filter(species == "AGPE") %>% 
  ggplot()+
  geom_point(aes(x = year_t1, y = `Cumulative PPT (mm)`), color = "black") +
  theme_classic() + ylim(0,1400) + scale_color_manual(values = c("#ff7f00", "#6a3d9a"))+ scale_shape_manual(values=c(0:6))+ labs(y = "Cumulative Precipitation (mm)") + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))

endo_climate %>% 
  filter(species == "AGPE") %>% 
  ggplot()+
  geom_point(aes(x = year_t1, y = `Mean Temp. (C˚)`), color = "black") +
  theme_classic() + ylim(0,15) + scale_color_manual(values = c("#ff7f00", "#6a3d9a"))+ scale_shape_manual(values=c(0:6))+ labs(y = "Mean Temp. (C˚)") + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))




ggplot(data = endo_climate) +
  geom_point(aes(x = year_t1, y = value, color = variable)) +
  facet_grid(rows = vars(species)) +
  theme_classic() + scale_color_manual(values = c("#ff7f00","#6a3d9a"))+ labs(y = "Population Growth") + scale_x_continuous(breaks = c(2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019))

# Calculate effect size and sd for lambdas

endo <- endo_climate %>% 
  mutate(variable = recode(variable, yearly_eplus_rec = "E+", yearly_eminus_rec = "E-")) %>% 
  group_by(species,variable) %>% 
  summarize(mean_l = mean(value),
            sd_l = sd(value),
            lowCI = mean_l - 2*sd_l,
            highCI = mean_l + 2*sd_l)

CatPlot <- ggplot(data = endo, aes(x = variable,
                                   y = mean_l,
                                   ymin = lowCI, 
                                   ymax = highCI,
                                   color = variable)) +
  geom_pointrange(size = 1.0)  +
  facet_grid(rows = vars(species)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  xlab("Endophyte") + ylab("Mean Pop. Growth Rate (± 2 sd)") +
  coord_flip() +
  theme_classic(base_size = 15) + scale_color_manual(values = c("#ff7f00","#6a3d9a"))+theme(axis.title = element_text(size = 30), axis.text = element_text(size = 20)) 








# some mean and sd of the lambda values,
mean(yearly(agpe_mat)$yearly_eminus)
mean(yearly(agpe_mat)$yearly_eplus)
sd(yearly(agpe_mat)$yearly_eminus)
sd(yearly(agpe_mat)$yearly_eplus)
mean(yearly(agpe_mat)$yearly_eminus_rec)
mean(yearly(agpe_mat)$yearly_eplus_rec)
sd(yearly(agpe_mat)$yearly_eminus_rec)
sd(yearly(agpe_mat)$yearly_eplus_rec)











image(bigmatrix(agpe_params)$Fmat_eplus_y1)
image(bigmatrix(agpe_params)$Fmat_eminus)
image(bigmatrix(agpe_params)$Tmat_eplus)
image(bigmatrix(agpe_params)$Tmat_eminus)
image(bigmatrix(agpe_params)$MPMmat_eplus)
image(bigmatrix(agpe_params)$MPMmat_eminus)
lambda(bigmatrix(agpe_params)$MPMmat_eplus_rec)
lambda(bigmatrix(agpe_params)$MPMmat_eminus)

image(bigmatrix(fesu_params)$Fmat_eplus)
image(bigmatrix(fesu_params)$Fmat_eminus)
image(bigmatrix(fesu_params)$Tmat_eplus)
image(bigmatrix(fesu_params)$Tmat_eminus)
image(bigmatrix(fesu_params)$MPMmat_eplus)
image(bigmatrix(fesu_params)$MPMmat_eminus)
lambda(bigmatrix(fesu_params)$MPMmat_eplus)
lambda(bigmatrix(fesu_params)$MPMmat_eminus)

image(bigmatrix(loar_params)$Fmat_eplus)
image(bigmatrix(loar_params)$Fmat_eminus)
image(bigmatrix(loar_params)$Tmat_eplus)
image(bigmatrix(loar_params)$Tmat_eminus)
image(bigmatrix(loar_params)$MPMmat_eplus)
image(bigmatrix(loar_params)$MPMmat_eminus)
lambda(bigmatrix(loar_params)$MPMmat_eplus)
lambda(bigmatrix(loar_params)$MPMmat_eminus)

image(bigmatrix(poal_params)$Fmat_eplus)
image(bigmatrix(poal_params)$Fmat_eminus)
image(bigmatrix(poal_params)$Tmat_eplus)
image(bigmatrix(poal_params)$Tmat_eminus)
image(bigmatrix(poal_params)$MPMmat_eplus)
image(bigmatrix(poal_params)$MPMmat_eminus)
lambda(bigmatrix(poal_params)$MPMmat_eplus)
lambda(bigmatrix(poal_params)$MPMmat_eminus)

image(bigmatrix(posy_params)$Fmat_eplus)
image(bigmatrix(posy_params)$Fmat_eminus)
image(bigmatrix(posy_params)$Tmat_eplus)
image(bigmatrix(posy_params)$Tmat_eminus)
image(bigmatrix(posy_params)$MPMmat_eplus)
image(bigmatrix(posy_params)$MPMmat_eminus)
lambda(bigmatrix(posy_params)$MPMmat_eplus)
lambda(bigmatrix(posy_params)$MPMmat_eminus)




