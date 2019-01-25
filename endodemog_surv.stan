
    data { 
    int<lower=0> N;                   // number of observations
    int surv_t1[N];  // plant survival at time t+1 and target variable
    vector[N] logsize_t;        // log of plant size at time t
    }
    
    parameters {
    real alpha; //fixed intercept
    real beta; // fixed slope 
    }

    model {
    vector[N] mu;
  
  //Priors
    //alpha ~normal(0,10);
    //beta ~ normal(0,10);
    mu = alpha + beta*logsize_t;  // linear predictor
    surv_t1 ~ bernoulli_logit(mu); // likelihood
    }
    
