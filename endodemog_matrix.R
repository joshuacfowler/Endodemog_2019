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

# spike
params[95] <- lapply(rstan::extract(spike, pars = "tau_year[1,1]"), FUN = mean); names(params)[95] <- "spike_eminus_y1"
params[96] <- lapply(rstan::extract(spike, pars = "tau_year[1,2]"), FUN = mean); names(params)[96] <- "spike_eminus_y2"
params[97] <- lapply(rstan::extract(spike, pars = "tau_year[1,3]"), FUN = mean); names(params)[97] <- "spike_eminus_y3"
params[98] <- lapply(rstan::extract(spike, pars = "tau_year[1,4]"), FUN = mean); names(params)[98] <- "spike_eminus_y4"
params[99] <- lapply(rstan::extract(spike, pars = "tau_year[1,5]"), FUN = mean); names(params)[99] <- "spike_eminus_y5"
params[100] <- lapply(rstan::extract(spike, pars = "tau_year[1,6]"), FUN = mean); names(params)[100] <- "spike_eminus_y6"
params[101] <- lapply(rstan::extract(spike, pars = "tau_year[1,7]"), FUN = mean); names(params)[101] <- "spike_eminus_y7"
params[102] <- lapply(rstan::extract(spike, pars = "tau_year[1,8]"), FUN = mean); names(params)[102] <- "spike_eminus_y8"
params[103] <- lapply(rstan::extract(spike, pars = "tau_year[1,9]"), FUN = mean); names(params)[103] <- "spike_eminus_y9"
params[104] <- lapply(rstan::extract(spike, pars = "tau_year[1,10]"), FUN = mean); names(params)[104] <- "spike_eminus_y10"
params[105] <- lapply(rstan::extract(spike, pars = "tau_year[1,11]"), FUN = mean); names(params)[105] <- "spike_eminus_y11"

params[106] <- lapply(rstan::extract(spike, pars = "tau_year[2,1]"), FUN = mean); names(params)[106] <- "spike_eplus_y1"
params[107] <- lapply(rstan::extract(spike, pars = "tau_year[2,2]"), FUN = mean); names(params)[107] <- "spike_eplus_y2"
params[108] <- lapply(rstan::extract(spike, pars = "tau_year[2,3]"), FUN = mean); names(params)[108] <- "spike_eplus_y3"
params[109] <- lapply(rstan::extract(spike, pars = "tau_year[2,4]"), FUN = mean); names(params)[109] <- "spike_eplus_y4"
params[110] <- lapply(rstan::extract(spike, pars = "tau_year[2,5]"), FUN = mean); names(params)[110] <- "spike_eplus_y5"
params[111] <- lapply(rstan::extract(spike, pars = "tau_year[2,6]"), FUN = mean); names(params)[111] <- "spike_eplus_y6"
params[112] <- lapply(rstan::extract(spike, pars = "tau_year[2,7]"), FUN = mean); names(params)[112] <- "spike_eplus_y7"
params[113] <- lapply(rstan::extract(spike, pars = "tau_year[2,8]"), FUN = mean); names(params)[113] <- "spike_eplus_y8"
params[114] <- lapply(rstan::extract(spike, pars = "tau_year[2,9]"), FUN = mean); names(params)[114] <- "spike_eplus_y9"
params[115] <- lapply(rstan::extract(spike, pars = "tau_year[2,10]"), FUN = mean); names(params)[115] <- "spike_eplus_y10"
params[116] <- lapply(rstan::extract(spike, pars = "tau_year[2,11]"), FUN = mean); names(params)[116] <- "spike_eplus_y11"

# s_to_s
params[117] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,1]"), FUN = mean); names(params)[117] <- "s_to_s_eminus_y1"
params[118] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,2]"), FUN = mean); names(params)[118] <- "s_to_s_eminus_y2"
params[119] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,3]"), FUN = mean); names(params)[119] <- "s_to_s_eminus_y3"
params[120] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,4]"), FUN = mean); names(params)[120] <- "s_to_s_eminus_y4"
params[121] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,5]"), FUN = mean); names(params)[121] <- "s_to_s_eminus_y5"
params[122] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,6]"), FUN = mean); names(params)[122] <- "s_to_s_eminus_y6"
params[123] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,7]"), FUN = mean); names(params)[123] <- "s_to_s_eminus_y7"
params[124] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,8]"), FUN = mean); names(params)[124] <- "s_to_s_eminus_y8"
params[125] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,9]"), FUN = mean); names(params)[125] <- "s_to_s_eminus_y9"
params[126] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,10]"), FUN = mean); names(params)[126] <- "s_to_s_eminus_y10"
params[127] <- lapply(rstan::extract(s_to_s, pars = "tau_year[1,11]"), FUN = mean); names(params)[127] <- "s_to_s_eminus_y11"

params[128] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,1]"), FUN = mean); names(params)[128] <- "s_to_s_eplus_y1"
params[129] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,2]"), FUN = mean); names(params)[129] <- "s_to_s_eplus_y2"
params[130] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,3]"), FUN = mean); names(params)[130] <- "s_to_s_eplus_y3"
params[131] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,4]"), FUN = mean); names(params)[131] <- "s_to_s_eplus_y4"
params[132] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,5]"), FUN = mean); names(params)[132] <- "s_to_s_eplus_y5"
params[133] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,6]"), FUN = mean); names(params)[133] <- "s_to_s_eplus_y6"
params[134] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,7]"), FUN = mean); names(params)[134] <- "s_to_s_eplus_y7"
params[135] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,8]"), FUN = mean); names(params)[135] <- "s_to_s_eplus_y8"
params[136] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,9]"), FUN = mean); names(params)[136] <- "s_to_s_eplus_y9"
params[137] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,10]"), FUN = mean); names(params)[137] <- "s_to_s_eplus_y10"
params[138] <- lapply(rstan::extract(s_to_s, pars = "tau_year[2,11]"), FUN = mean); names(params)[138] <- "s_to_s_eplus_y11"

params <- unlist(params)
return(params)
}

agpe_params <- getparams(surv = survAGPE, grow = growAGPE, flw = flwAGPE, fert = fertAGPE, spike = spikeAGPE, seed = seedAGPE, s_to_s = s_to_sAGPE, data = AGPE_surv_data)
fesu_params <- getparams(surv = survFESU, grow = growFESU, flw = flwFESU, fert = fertFESU, spike = spikeFESU, seed = seedFESU, s_to_s = s_to_sFESU, data = FESU_surv_data)
loar_params <- getparams(surv = survLOAR, grow = growLOAR, flw = flwLOAR, fert = fertLOAR, spike = spikeLOAR, seed = seedLOAR, s_to_s = s_to_sLOAR, data = LOAR_surv_data)
poal_params <- getparams(surv = survPOAL, grow = growPOAL, flw = flwPOAL, fert = fertPOAL, spike = spikePOAL, seed = seedPOAL, s_to_s = s_to_sPOAL, data = POAL_surv_data)
posy_params <- getparams(surv = survPOSY, grow = growPOSY, flw = flwPOSY, fert = fertPOSY, spike = spikePOSY, seed = seedPOSY, s_to_s = s_to_sPOSY, data = POSY_surv_data)


# define functions that will be used to populate projection matrix
#SURVIVAL AT SIZE X.
# currently this is fitting as if the E- is the intercept

sx<-function(x,params){
  eminus_surv <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x))
  eplus_surv <- invlogit(params["surv_beta1"] + params["surv_beta2"]*log(x) + params["surv_beta3"])
  return(list(eminus_surv = eminus_surv, eplus_surv = eplus_surv))
}
# sx(x = 10, params = poal_params) # we can test out our functions for different sizes of x
gxy <- function(x,y,params){
  eminus_grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x)
  eminus_pr_grow <- dnbinom(x=y, mu = exp(eminus_grow.mean), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eminus_grow.mean), size = params["grow_phi"])))
  eplus_grow.mean <- params["grow_beta1"] + params["grow_beta2"]*log(x) + params["grow_beta3"]
  eplus_pr_grow <- dnbinom(x=y, mu = exp(eplus_grow.mean), size = params["grow_phi"]/(1-dnbinom(0, mu = exp(eplus_grow.mean), size = params["grow_phi"])))
  
return(list(eminus_pr_grow  = eminus_pr_grow, eplus_pr_grow = eplus_pr_grow))
}
# gxy(x = 10, y = 11, params = loar_params)$eplus_pr_grow
# I'll truncate this later...
# prob=dnbinom(1:n_post_draws, mu = exp(mu[i,j]), size = phi[i]/(1-dnbinom(0, mu = exp(mu[i,j]), size = phi[i]))))


#SURVIVAL*GROWTH
pxy_eminus<-function(x,y,params){
  sx(x,params)$eminus_surv * gxy(x,y,params)$eminus_pr_grow
}
pxy_eplus<-function(x,y,params){
  sx(x,params)$eplus_surv * gxy(x,y,params)$eplus_pr_grow
}
# pxy_eminus(x = 2, y = 10, params = loar_params)
# pxy_eplus(x = 2, y = 10, params = loar_params)

#FERTILITY--returns number of seedlings, which we will assume (for now) to be 1-tiller, produced by size X
fx<-function(x, params){
  eminus_p_flw <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x)) 
  eminus_fert.mean <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x))
  eminus_spike.mean <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x))
  eminus_seed.mean <- params["mu_seed"]
  eminus_p_rec <- invlogit(params["s_to_s_beta1"])
  
  eplus_p_flw <- invlogit(params["flw_beta1"] + params["flw_beta2"]*log(x) + params["flw_beta3"])
  eplus_fert.mean <- exp(params["fert_beta1"] + params["fert_beta2"]*log(x) + params["fert_beta3"])
  eplus_spike.mean <- exp(params["spike_beta1"] + params["spike_beta2"]*log(x) + params["spike_beta3"])
  eplus_seed.mean <- params["mu_seed"]
  eplus_p_rec <- invlogit(params["s_to_s_beta1"] + params["s_to_s_beta2"])
  
  eminus_seedlings <- eminus_p_flw * eminus_fert.mean * eminus_spike.mean * eminus_seed.mean * eminus_p_rec
  eplus_seedlings <- eplus_p_flw * eplus_fert.mean * eplus_spike.mean * eplus_seed.mean * eplus_p_rec
  return(list(eminus_seedlings = eminus_seedlings, eplus_seedlings = eplus_seedlings))
}

## note from Tom: we should think about adding a reproductive delay
# finally, here is the function that takes the parameter vector and assembles the matrix model from all of the pieces
bigmatrix<-function(params){   
  
  matdim<-params["max_size"]## matrix dimension
  y <- 1:params["max_size"]## size (tiller number) associated with each class
  # Fertility matrix
  Fmat_eminus<-matrix(0,matdim,matdim)
  Fmat_eplus<-matrix(0,matdim,matdim)
  
  # all seedlings get dumped into top row (1-tiller)
  Fmat_eminus[1,]<-fx(x = y, params=params)$eminus_seedlings 
  Fmat_eplus[1,]<-fx(x = y, params=params)$eplus_seedlings 
  
  # Growth/survival transition matrix
  Tmat_eminus<-matrix(0,matdim,matdim)
  Tmat_eminus<-t(outer(y,y,pxy_eminus,params=params))
  
  Tmat_eplus<-matrix(0,matdim,matdim)
  Tmat_eplus<-t(outer(y,y,pxy_eplus,params=params))
  
  # Put it all together
  MPMmat_eminus<-Tmat_eminus+Fmat_eminus #sum the Tmat & Fmat to get the whole matrix
  MPMmat_eplus<-Tmat_eplus+Fmat_eplus #sum the Tmat & Fmat to get the whole matrix
  
  return(list(MPMmat_eminus=MPMmat_eminus,MPMmat_eplus=MPMmat_eplus,Fmat_eminus=Fmat_eminus,Fmat_eplus=Fmat_eplus,Tmat_eminus=Tmat_eminus,Tmat_eplus=Tmat_eplus))
}

## population growth rate (eigenalaysis of the projection matrix)
# matrix <- bigmatrix(loar_params)
image(bigmatrix(agpe_params)$Fmat_eplus)
image(bigmatrix(agpe_params)$Fmat_eminus)
image(bigmatrix(agpe_params)$Tmat_eplus)
image(bigmatrix(agpe_params)$Tmat_eminus)
image(bigmatrix(agpe_params)$MPMmat_eplus)
image(bigmatrix(agpe_params)$MPMmat_eminus)
lambda(bigmatrix(agpe_params)$MPMmat_eplus)
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




