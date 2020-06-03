## purpose calculate stochastic growth rates, estimate endophyte effects,
## decompose endophyte effects on means versus variances.

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

  rtracker  <- rtracker[-c(1:burnin)]
  #Finish and return
  log_lambdaS <- mean(rtracker)
  lambdaS<-exp(log_lambdaS)
  return(list(log_lambdaS=log_lambdaS,lambdaS=lambdaS))
}

## toy example to show that this works
toy_matrices <- list(matrix(c(1,0.5,2,0.1),2,2),matrix(c(0.8,0.7,5,0.3),2,2),matrix(c(0.4,0.6,0.4,0.9),2,2))
lambdaSim(toy_matrices)


# Analysis ----------------------------------------------------------------
## For each species, we want four lists of matrices:
## 1. E-plus matrices with endo effects on mean (apply beta3) and variance (use E+ random year effects)
## 2. E-plus matrices with endo effects on mean (apply beta3) but not variance (use E- random year effects)
## 3. E-plus matrices with endo effects variance (use E+ random year effects) but not mean (do not apply beta3)
## 4. E-minus matrices with no endo effects (do not apply beta3 and use E- random year effects)

## to start, you can just use the posterior mean parameter estimates (don't worry about posterior sampling)
## so for AGPE as an example:
AGPE_lambdaS <- c()

AGPE_Eplus_mean_var <- list() ## assemble this from parameters
AGPE_lambdaS[1]<-lambdaSim(AGPE_Eplus_mean_var)$log_lambdaS

AGPE_Eplus_mean_only <- list() ## assemble this from parameters
AGPE_lambdaS[2]<-lambdaSim(AGPE_Eplus_mean_only)$log_lambdaS

AGPE_Eplus_var_only <- list() ## assemble this from parameters
AGPE_lambdaS[3]<-lambdaSim(AGPE_Eplus_var_only)$log_lambdaS

AGPE_Eminus <- list() ## assemble this from parameters
AGPE_Eplus_mean_var <- list(bigmatrix(agpe_1m_1v_rec_params, year = 1)$MPMmat,
                            bigmatrix(agpe_1m_1v_rec_params, year = 2)$MPMmat,
                            bigmatrix(agpe_1m_1v_rec_params, year = 3)$MPMmat,
                            bigmatrix(agpe_1m_1v_rec_params, year = 4)$MPMmat,
                            bigmatrix(agpe_1m_1v_rec_params, year = 5)$MPMmat,
                            bigmatrix(agpe_1m_1v_rec_params, year = 6)$MPMmat,
                            bigmatrix(agpe_1m_1v_rec_params, year = 7)$MPMmat,
                            bigmatrix(agpe_1m_1v_rec_params, year = 8)$MPMmat,
                            bigmatrix(agpe_1m_1v_rec_params, year = 9)$MPMmat,
                            bigmatrix(agpe_1m_1v_rec_params, year = 10)$MPMmat,
                            bigmatrix(agpe_1m_1v_rec_params, year = 11)$MPMmat) ## assemble this from parameters
AGPE_lambdaS[1]<-lambdaSim(AGPE_Eplus_mean_var)$log_lambdaS

AGPE_Eplus_mean_only <- list(bigmatrix(agpe_1m_0v_rec_params, year = 1)$MPMmat,
                             bigmatrix(agpe_1m_0v_rec_params, year = 2)$MPMmat,
                             bigmatrix(agpe_1m_0v_rec_params, year = 3)$MPMmat,
                             bigmatrix(agpe_1m_0v_rec_params, year = 4)$MPMmat,
                             bigmatrix(agpe_1m_0v_rec_params, year = 5)$MPMmat,
                             bigmatrix(agpe_1m_0v_rec_params, year = 6)$MPMmat,
                             bigmatrix(agpe_1m_0v_rec_params, year = 7)$MPMmat,
                             bigmatrix(agpe_1m_0v_rec_params, year = 8)$MPMmat,
                             bigmatrix(agpe_1m_0v_rec_params, year = 9)$MPMmat,
                             bigmatrix(agpe_1m_0v_rec_params, year = 10)$MPMmat,
                             bigmatrix(agpe_1m_0v_rec_params, year = 11)$MPMmat) ## assemble this from parameters
AGPE_lambdaS[2]<-lambdaSim(AGPE_Eplus_mean_only)$log_lambdaS

AGPE_Eplus_var_only <- list(bigmatrix(agpe_0m_1v_rec_params, year = 1)$MPMmat,
                            bigmatrix(agpe_0m_1v_rec_params, year = 2)$MPMmat,
                            bigmatrix(agpe_0m_1v_rec_params, year = 3)$MPMmat,
                            bigmatrix(agpe_0m_1v_rec_params, year = 4)$MPMmat,
                            bigmatrix(agpe_0m_1v_rec_params, year = 5)$MPMmat,
                            bigmatrix(agpe_0m_1v_rec_params, year = 6)$MPMmat,
                            bigmatrix(agpe_0m_1v_rec_params, year = 7)$MPMmat,
                            bigmatrix(agpe_0m_1v_rec_params, year = 8)$MPMmat,
                            bigmatrix(agpe_0m_1v_rec_params, year = 9)$MPMmat,
                            bigmatrix(agpe_0m_1v_rec_params, year = 10)$MPMmat,
                            bigmatrix(agpe_0m_1v_rec_params, year = 11)$MPMmat) ## assemble this from parameters
AGPE_lambdaS[3]<-lambdaSim(AGPE_Eplus_var_only)$log_lambdaS

AGPE_Eminus <- list(bigmatrix(agpe_0m_0v_rec_params, year = 1)$MPMmat,
                    bigmatrix(agpe_0m_0v_rec_params, year = 2)$MPMmat,
                    bigmatrix(agpe_0m_0v_rec_params, year = 3)$MPMmat,
                    bigmatrix(agpe_0m_0v_rec_params, year = 4)$MPMmat,
                    bigmatrix(agpe_0m_0v_rec_params, year = 5)$MPMmat,
                    bigmatrix(agpe_0m_0v_rec_params, year = 6)$MPMmat,
                    bigmatrix(agpe_0m_0v_rec_params, year = 7)$MPMmat,
                    bigmatrix(agpe_0m_0v_rec_params, year = 8)$MPMmat,
                    bigmatrix(agpe_0m_0v_rec_params, year = 9)$MPMmat,
                    bigmatrix(agpe_0m_0v_rec_params, year = 10)$MPMmat,
                    bigmatrix(agpe_0m_0v_rec_params, year = 11)$MPMmat) ## assemble this from parameters
AGPE_lambdaS[4]<-lambdaSim(AGPE_Eminus)$log_lambdaS

## compute E+ / E- difference for different types of E+ effects (both, mean, variance)
AGPE_lambdaS_effect <- AGPE_lambdaS[1:3]-AGPE_lambdaS[4]

## visualize results to see how much the endo effect is due to mean vs variance
barplot(AGPE_lambdaS_effect,col=rainbow(3),
        names.arg = c("Total effect","Mean effect","Variance effect"))

## so for FESU as an example:
FESU_lambdaS <- c()

FESU_Eplus_mean_var <- list(bigmatrix(fesu_1m_1v_rec_params, year = 1)$MPMmat,
                            bigmatrix(fesu_1m_1v_rec_params, year = 2)$MPMmat,
                            bigmatrix(fesu_1m_1v_rec_params, year = 3)$MPMmat,
                            bigmatrix(fesu_1m_1v_rec_params, year = 4)$MPMmat,
                            bigmatrix(fesu_1m_1v_rec_params, year = 5)$MPMmat,
                            bigmatrix(fesu_1m_1v_rec_params, year = 6)$MPMmat,
                            bigmatrix(fesu_1m_1v_rec_params, year = 7)$MPMmat,
                            bigmatrix(fesu_1m_1v_rec_params, year = 8)$MPMmat,
                            bigmatrix(fesu_1m_1v_rec_params, year = 9)$MPMmat,
                            bigmatrix(fesu_1m_1v_rec_params, year = 10)$MPMmat,
                            bigmatrix(fesu_1m_1v_rec_params, year = 11)$MPMmat) ## assemble this from parameters
FESU_lambdaS[1]<-lambdaSim(FESU_Eplus_mean_var)$log_lambdaS

FESU_Eplus_mean_only <- list(bigmatrix(fesu_1m_0v_rec_params, year = 1)$MPMmat,
                             bigmatrix(fesu_1m_0v_rec_params, year = 2)$MPMmat,
                             bigmatrix(fesu_1m_0v_rec_params, year = 3)$MPMmat,
                             bigmatrix(fesu_1m_0v_rec_params, year = 4)$MPMmat,
                             bigmatrix(fesu_1m_0v_rec_params, year = 5)$MPMmat,
                             bigmatrix(fesu_1m_0v_rec_params, year = 6)$MPMmat,
                             bigmatrix(fesu_1m_0v_rec_params, year = 7)$MPMmat,
                             bigmatrix(fesu_1m_0v_rec_params, year = 8)$MPMmat,
                             bigmatrix(fesu_1m_0v_rec_params, year = 9)$MPMmat,
                             bigmatrix(fesu_1m_0v_rec_params, year = 10)$MPMmat,
                             bigmatrix(fesu_1m_0v_rec_params, year = 11)$MPMmat) ## assemble this from parameters
FESU_lambdaS[2]<-lambdaSim(FESU_Eplus_mean_only)$log_lambdaS

FESU_Eplus_var_only <- list(bigmatrix(fesu_0m_1v_rec_params, year = 1)$MPMmat,
                            bigmatrix(fesu_0m_1v_rec_params, year = 2)$MPMmat,
                            bigmatrix(fesu_0m_1v_rec_params, year = 3)$MPMmat,
                            bigmatrix(fesu_0m_1v_rec_params, year = 4)$MPMmat,
                            bigmatrix(fesu_0m_1v_rec_params, year = 5)$MPMmat,
                            bigmatrix(fesu_0m_1v_rec_params, year = 6)$MPMmat,
                            bigmatrix(fesu_0m_1v_rec_params, year = 7)$MPMmat,
                            bigmatrix(fesu_0m_1v_rec_params, year = 8)$MPMmat,
                            bigmatrix(fesu_0m_1v_rec_params, year = 9)$MPMmat,
                            bigmatrix(fesu_0m_1v_rec_params, year = 10)$MPMmat,
                            bigmatrix(fesu_0m_1v_rec_params, year = 11)$MPMmat) ## assemble this from parameters
FESU_lambdaS[3]<-lambdaSim(FESU_Eplus_var_only)$log_lambdaS

FESU_Eminus <- list(bigmatrix(fesu_0m_0v_rec_params, year = 1)$MPMmat,
                    bigmatrix(fesu_0m_0v_rec_params, year = 2)$MPMmat,
                    bigmatrix(fesu_0m_0v_rec_params, year = 3)$MPMmat,
                    bigmatrix(fesu_0m_0v_rec_params, year = 4)$MPMmat,
                    bigmatrix(fesu_0m_0v_rec_params, year = 5)$MPMmat,
                    bigmatrix(fesu_0m_0v_rec_params, year = 6)$MPMmat,
                    bigmatrix(fesu_0m_0v_rec_params, year = 7)$MPMmat,
                    bigmatrix(fesu_0m_0v_rec_params, year = 8)$MPMmat,
                    bigmatrix(fesu_0m_0v_rec_params, year = 9)$MPMmat,
                    bigmatrix(fesu_0m_0v_rec_params, year = 10)$MPMmat,
                    bigmatrix(fesu_0m_0v_rec_params, year = 11)$MPMmat) ## assemble this from parameters
FESU_lambdaS[4]<-lambdaSim(FESU_Eminus)$log_lambdaS

## compute E+ / E- difference for different types of E+ effects (both, mean, variance)
FESU_lambdaS_effect <- FESU_lambdaS[1:3]-FESU_lambdaS[4]

## visualize results to see how much the endo effe
barplot(FESU_lambdaS_effect,col=rainbow(3),
        names.arg = c("Total effect","Mean effect","Variance effect"))





## so for LOAR as an example:
LOAR_lambdaS <- c()

LOAR_Eplus_mean_var <- list(bigmatrix(loar_1m_1v_rec_params, year = 1)$MPMmat,
                            bigmatrix(loar_1m_1v_rec_params, year = 2)$MPMmat,
                            bigmatrix(loar_1m_1v_rec_params, year = 3)$MPMmat,
                            bigmatrix(loar_1m_1v_rec_params, year = 4)$MPMmat,
                            bigmatrix(loar_1m_1v_rec_params, year = 5)$MPMmat,
                            bigmatrix(loar_1m_1v_rec_params, year = 6)$MPMmat,
                            bigmatrix(loar_1m_1v_rec_params, year = 7)$MPMmat,
                            bigmatrix(loar_1m_1v_rec_params, year = 8)$MPMmat,
                            bigmatrix(loar_1m_1v_rec_params, year = 9)$MPMmat,
                            bigmatrix(loar_1m_1v_rec_params, year = 10)$MPMmat,
                            bigmatrix(loar_1m_1v_rec_params, year = 11)$MPMmat) ## assemble this from parameters
LOAR_lambdaS[1]<-lambdaSim(LOAR_Eplus_mean_var)$log_lambdaS

LOAR_Eplus_mean_only <- list(bigmatrix(loar_1m_0v_rec_params, year = 1)$MPMmat,
                             bigmatrix(loar_1m_0v_rec_params, year = 2)$MPMmat,
                             bigmatrix(loar_1m_0v_rec_params, year = 3)$MPMmat,
                             bigmatrix(loar_1m_0v_rec_params, year = 4)$MPMmat,
                             bigmatrix(loar_1m_0v_rec_params, year = 5)$MPMmat,
                             bigmatrix(loar_1m_0v_rec_params, year = 6)$MPMmat,
                             bigmatrix(loar_1m_0v_rec_params, year = 7)$MPMmat,
                             bigmatrix(loar_1m_0v_rec_params, year = 8)$MPMmat,
                             bigmatrix(loar_1m_0v_rec_params, year = 9)$MPMmat,
                             bigmatrix(loar_1m_0v_rec_params, year = 10)$MPMmat,
                             bigmatrix(loar_1m_0v_rec_params, year = 11)$MPMmat) ## assemble this from parameters
LOAR_lambdaS[2]<-lambdaSim(LOAR_Eplus_mean_only)$log_lambdaS

LOAR_Eplus_var_only <- list(bigmatrix(loar_0m_1v_rec_params, year = 1)$MPMmat,
                            bigmatrix(loar_0m_1v_rec_params, year = 2)$MPMmat,
                            bigmatrix(loar_0m_1v_rec_params, year = 3)$MPMmat,
                            bigmatrix(loar_0m_1v_rec_params, year = 4)$MPMmat,
                            bigmatrix(loar_0m_1v_rec_params, year = 5)$MPMmat,
                            bigmatrix(loar_0m_1v_rec_params, year = 6)$MPMmat,
                            bigmatrix(loar_0m_1v_rec_params, year = 7)$MPMmat,
                            bigmatrix(loar_0m_1v_rec_params, year = 8)$MPMmat,
                            bigmatrix(loar_0m_1v_rec_params, year = 9)$MPMmat,
                            bigmatrix(loar_0m_1v_rec_params, year = 10)$MPMmat,
                            bigmatrix(loar_0m_1v_rec_params, year = 11)$MPMmat) ## assemble this from parameters
LOAR_lambdaS[3]<-lambdaSim(LOAR_Eplus_var_only)$log_lambdaS

LOAR_Eminus <- list(bigmatrix(loar_0m_0v_rec_params, year = 1)$MPMmat,
                    bigmatrix(loar_0m_0v_rec_params, year = 2)$MPMmat,
                    bigmatrix(loar_0m_0v_rec_params, year = 3)$MPMmat,
                    bigmatrix(loar_0m_0v_rec_params, year = 4)$MPMmat,
                    bigmatrix(loar_0m_0v_rec_params, year = 5)$MPMmat,
                    bigmatrix(loar_0m_0v_rec_params, year = 6)$MPMmat,
                    bigmatrix(loar_0m_0v_rec_params, year = 7)$MPMmat,
                    bigmatrix(loar_0m_0v_rec_params, year = 8)$MPMmat,
                    bigmatrix(loar_0m_0v_rec_params, year = 9)$MPMmat,
                    bigmatrix(loar_0m_0v_rec_params, year = 10)$MPMmat,
                    bigmatrix(loar_0m_0v_rec_params, year = 11)$MPMmat) ## assemble this from parameters
LOAR_lambdaS[4]<-lambdaSim(LOAR_Eminus)$log_lambdaS

## compute E+ / E- difference for different types of E+ effects (both, mean, variance)
LOAR_lambdaS_effect <- LOAR_lambdaS[1:3]-LOAR_lambdaS[4]

## visualize results to see how much the endo effe
barplot(LOAR_lambdaS_effect,col=rainbow(3),
        names.arg = c("Total effect","Mean effect","Variance effect"))
