
    data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> K;                       // number of predictors
    
    int<lower=0> nyear;                       // number of years (used as index)
    int<lower=0> year_t[N];                      // year of observation
    int<lower=0> nEndo;                       // number of endo treatments
    
    int<lower=0, upper=1> surv_t1[N];      // plant survival at time t+1 and target variable (response)
    vector<lower=-1>[N] logsize_t;                  // log of plant size at time t (predictor)
    int<lower=0, upper=1> endo[N];            // endophyte status 
    int<lower=1, upper=2> endo_index[N];       // index for endophyte effect
    int<lower=0, upper=1> origin[N];            // origin status
    }
    
    parameters {
    real alpha;                  // intercept
    vector[K] beta;              // predictor parameters
    matrix[nEndo, nyear] tau_year;      // random year effect
      
    real<lower=0> sigma_0[nEndo];        //year variance intercept
    }
    

    model {
    vector[N] mu;
   
       // Linear Predictor
    for(n in 1:N){
       mu[n] = alpha + beta[1]*logsize_t[n] + beta[2]*endo[n] + beta[3]*origin[n] 
       + beta[4]*logsize_t[n]*endo[n] +
       + tau_year[endo_index[n], year_t[n]];
    }
    // Priors
    alpha ~ normal(0,1e6);      // prior for fixed intercept
    beta ~ normal(0,1e6);      // prior for predictor intercepts
    for(n in 1:nyear){
    for(e in 1:nEndo){
    tau_year[e,n] ~ normal(0,sigma_0[e]);   // prior for year random effects
    }}
    // Likelihood
      surv_t1 ~ bernoulli_logit(mu);
    }
    
   generated quantities{
    int yrep[N];
    vector[N] mu;
    
    // for posterior predictive check
    for(n in 1:N){
      mu[n] = alpha + beta[1]*logsize_t[n] + beta[2]*endo[n] + beta[3]*origin[n] 
      +beta[4]*logsize_t[n]*endo[n]
      + tau_year[endo_index[n], year_t[n]];
      
      yrep[n] = bernoulli_logit_rng(mu[n]);
    }
    
    }
  
    
