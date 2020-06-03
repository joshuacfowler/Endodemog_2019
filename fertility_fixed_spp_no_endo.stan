
data {
    int<lower=0> nYear;                       // number of years
    int<lower=0> nPlot;                         // number of plots
    int<lower=0> nSpp;                        // number of host species
    int<lower=0> N;                       // number of observations for surv model
    int<lower=0, upper=nYear> year_t[N];         // year of observation for surv model
    int<lower=0> plot[N];                   // plot of observation for surv model
    int<lower=0, upper=nSpp> spp[N];         // year of observation for surv model
    int<lower=0> y[N];      // plant survival or flowering at time t+1
    vector<lower=0>[N] logsize_t;             // plant size at time t for surv model
    int<lower=0,upper=1> origin_01[N];          // plant origin status for surv model
}

parameters {
    vector[nSpp] beta0; 
    vector[nSpp] betasize;                     
    vector[nSpp] betaorigin;  

    vector[nYear] tau_year;        // random plot effect
    real<lower=0> sigma_year;          // plot variance effect
    vector[nPlot] tau_plot;        // random plot effect
    real<lower=0> sigma_plot;          // plot variance effect
    
    //neg bin overdispersion
    vector[nSpp] phi; 
}

transformed parameters {
    real lambda[N]; 
    real od[N];

    // surv Linear Predictor
    for(n in 1:N){
    lambda[n] = beta0[spp[n]] + 
    betasize[spp[n]]*logsize_t[n] + 
    betaorigin[spp[n]]*origin_01[n] + 
    tau_year[year_t[n]] + tau_plot[plot[n]];
    
    od[n] = exp(phi[spp[n]]);
    }
}

model {
    // priors

    sigma_plot ~ inv_gamma(0.001, 0.001);
    for(i in 1:nPlot){
      tau_plot[i] ~ normal(0,sigma_plot);
    }
    sigma_year ~ inv_gamma(0.001, 0.001);
    for(i in 1:nYear){
      tau_year[i] ~ normal(0,sigma_year);
    }
    
    for(s in 1:nSpp){
      beta0[s] ~ normal(0,100); 
      betasize[s] ~ normal(0,100); 
      betaorigin[s] ~ normal(0,100); 
      phi[s] ~ normal(0,100);
    }

    for(n in 1:N){
      y[n] ~ neg_binomial_2_log(lambda[n],od[n]);
      target += -log1m(neg_binomial_2_log_lpmf(0 | lambda[n], od[n])); 
     }
}

