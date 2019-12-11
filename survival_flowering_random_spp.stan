//right now this model has random species effects on the coefficients but pooled year and plot rfx 
//and no endo effects on variance

data {
    // indices
    int<lower=0> nYear;                       // number of years
    int<lower=0> nPlot;                         // number of plots
    int<lower=0> nSpp;                        // number of host species

    // surv data
    int<lower=0> N;                       // number of observations for surv model
    int<lower=0, upper=nYear> year_t[N];         // year of observation for surv model
    int<lower=0> plot[N];                   // plot of observation for surv model
    int<lower=0, upper=nSpp> spp[N];         // year of observation for surv model
    
    int<lower=0, upper=1> y[N];      // plant survival or flowering at time t+1
    vector<lower=0>[N] logsize_t;             // plant size at time t for surv model
    int<lower=0,upper=1> endo_01[N];            // plant endophyte status for surv model
    int<lower=0,upper=1> origin_01[N];          // plant origin status for surv model
}

parameters {
    // surv params
    vector[nSpp] beta0; //spp-specific                     // predictor parameters as grand means and spp rfx
    real mu_beta0; //spp mean
    real<lower=0> sigma_beta0; //spp variance
    vector[nSpp] betasize;                     
    real mu_betasize; 
    real<lower=0> sigma_betasize;
    vector[nSpp] betaendo;                     
    real mu_betaendo;                     
    real<lower=0> sigma_betaendo;
    vector[nSpp] betaorigin;                     
    real mu_betaorigin;                     
    real<lower=0> sigma_betaorigin;
    
    vector[nYear] tau_year;      // random year effect
    real<lower=0> sigma_year;        //year variance by endophyte effect
    vector[nPlot] tau_plot;        // random plot effect
    real<lower=0> sigma_plot;          // plot variance effect
}

transformed parameters {
    real mu[N];                           // surv Linear Predictor

       for(n in 1:N){
    mu[n] = beta0[spp[n]] + betasize[spp[n]]*logsize_t[n] + betaendo[spp[n]]*endo_01[n] + betaorigin[spp[n]]*origin_01[n]
    + tau_year[year_t[n]] 
    + tau_plot[plot[n]];
       }
    }

model {
    // surv priors
    mu_beta0 ~ normal(0,100);      // prior for predictor intercepts
    mu_betasize ~ normal(0,100);      // prior for predictor intercepts
    mu_betaendo ~ normal(0,100);      // prior for predictor intercepts
    mu_betaorigin ~ normal(0,100);      // prior for predictor intercepts
    //sigma priors
    sigma_beta0 ~ inv_gamma(0.001, 0.001);
    sigma_betasize ~ inv_gamma(0.001, 0.001);
    sigma_betaendo ~ inv_gamma(0.001, 0.001);
    sigma_betaorigin ~ inv_gamma(0.001, 0.001);
    sigma_year ~ inv_gamma(0.001, 0.001);
    sigma_plot ~ inv_gamma(0.001, 0.001);
    //sample random effects
    for(i in 1:nYear){
      tau_year[i] ~ normal(0,sigma_year);
    }
    for(i in 1:nPlot){
      tau_plot[i] ~ normal(0,sigma_plot);
    }
    for(i in 1:nSpp){
      beta0[i] ~ normal(0,sigma_beta0);
      betasize[i] ~ normal(0,sigma_betasize);
      betaendo[i] ~ normal(0,sigma_betaendo);
      betaorigin[i] ~ normal(0,sigma_betaorigin);
    }

    
    y ~ bernoulli_logit(mu);
}

