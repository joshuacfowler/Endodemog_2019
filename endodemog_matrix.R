## Title: Grass endophyte population model with a bayesian framework
## Purpose: Assemble matrix model from vital rate estimates
## Authors: Joshua and Tom
## Updated: 11/19/2019
#############################################################


# Assembling matrix model -------------------------------------------------

# LOAR
# collect all the parameters into one vector
loar_params <- c()
loar_params[1] <- fixef(surv_mod)[1]; names(loar_params)[1]<-"surv_intercept"
loar_params[2] <- fixef(surv_mod)[2]; names(loar_params)[2]<-"surv_slope"
loar_params[3] <- fixef(grow_mod)[1]; names(loar_params)[3]<-"grow_intercept"
loar_params[4] <- fixef(grow_mod)[2]; names(loar_params)[4]<-"grow_slope"
loar_params[5] <- fixef(flow_mod)[1]; names(loar_params)[5]<-"flow_intercept"
loar_params[6] <- fixef(flow_mod)[2]; names(loar_params)[6]<-"flow_slope"
loar_params[7] <- fixef(seed_mod)[1]; names(loar_params)[7]<-"seed_intercept"
loar_params[8] <- fixef(seed_mod)[2]; names(loar_params)[8]<-"seed_slope"
loar_params[9] <- prob_recruit; names(loar_params)[9]<-"prob_recruit"
loar_params[10] <- 1; names(loar_params)[10]<-"min_size"
loar_params[11] <- quantile(loar$size_t1,0.95,na.rm=T); names(loar_params)[11]<-"max_size"
# note that I define max size as the 95TH pctile of observed sized. The very max sizes observed often have very poor 
# replication, and I find that these few observations (and the corresponding vital rate predictions) can have a strong
# influence on the results. So this approach is conservative, but you can always experiment with this.

# define functions that will be used to populate projection matrix
#SURVIVAL AT SIZE X.
sx<-function(x,params){
  ilogit(params["surv_intercept"] + params["surv_slope"]*log(x))
}
