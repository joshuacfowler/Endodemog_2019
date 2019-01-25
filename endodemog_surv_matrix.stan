
    data { 
    int<lower=0> N;                       // number of observations
    int<lower=0> K;                       // number of predictors
    
    int<lower=0> nyear;                       // number of years (used as index)
    int<lower=0> year_t[N];                      // year of observation
    int<lower=0> nEndo;                       // number of endo treatments
    int<lower=1, upper=2> endo_index[N];       // index for endophyte effect

    int<lower=0, upper=1> surv_t1[N];      // plant survival at time t+1 and target variable (response)
    matrix[N,K] Xs;                  //  predictor matrix - surv_t1~logsize_t+endo+origin+logsize_t*endo
    }
    
    parameters {
    vector[K] beta;              // predictor parameters
    matrix[nEndo, nyear] tau_year;      // random year effect
      
    real<lower=0> sigma_0[nEndo];        //year variance intercept
    }
    

    model {
    vector[N] mu;
   
       // Linear Predictor
    for(n in 1:N){
       mu[n] = Xs[n]*beta
       + tau_year[endo_index[n], year_t[n]];
    }
    // Priors
    beta ~ normal(0,1e6);      // prior for predictor intercepts
    for(n in 1:nyear){
    for(e in 1:nEndo){
    tau_year[e,n] ~ normal(0,sigma_0[e]);   // prior for year random effects
    }}
    // Likelihood
      surv_t1 ~ bernoulli_logit(mu);
    }
    
   //generated quantities{
    //int yrep[N];
    //vector[N] mu;
    
    // for posterior predictive check
    //for(n in 1:N){
     // mu[n] = Xs[n]*beta
    //   + tau_year[endo_index[n], year_t[n]];
      
    //  yrep[n] = bernoulli_logit_rng(mu[n]);
    //}
    
   // }
  
    
