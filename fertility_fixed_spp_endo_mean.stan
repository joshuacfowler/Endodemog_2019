//right now this model has random species effects on the coefficients but pooled year and plot rfx 
//and no endo effects on variance

data {
    // indices
    int<lower=0> nYear;                       // number of years
    int<lower=0> nPlot;                         // number of plots
    int<lower=0> nSpp;                        // number of host species
    int<lower=0> nEndo;                        // number of endophyte levels
    // surv data
    int<lower=0> N;                       // number of observations for surv model
    int<lower=0, upper=nYear> year_t[N];         // year of observation for surv model
    int<lower=0> plot[N];                   // plot of observation for surv model
    int<lower=0, upper=nSpp> spp[N];         // year of observation for surv model
    int<lower=0> y[N];      // plant survival or flowering at time t+1
    vector<lower=0>[N] logsize_t;             // plant size at time t for surv model
    int<lower=0,upper=1> endo_01[N];            // plant endophyte status for surv model
    int<lower=0,upper=1> origin_01[N];          // plant origin status for surv model
}

parameters {
    // surv params
    vector[nSpp] beta0; //spp-specific      // predictor parameters as grand means and spp rfx
    //real mu_beta0; //spp mean
    //real<lower=0> sigma_beta0; //spp variance
    vector[nSpp] betasize;                     
    //real mu_betasize; 
    //real<lower=0> sigma_betasize;
    vector[nSpp] betaendo;                     
    //real mu_betaendo;                     
    //real<lower=0> sigma_betaendo;
    vector[nSpp] betaorigin;  
    
    real tau_year[nSpp,nYear];      // random year effect, unique to species and endo
    //real mu_sigma0;
    //real<lower=0> sigma_mu_sigma0;
    vector[nSpp] sigma0;
    //real mu_sigmaendo;
    //vector[nSpp] sigmaendo;
    //real<lower=0> sigma_mu_sigmaendo;
    
    vector[nPlot] tau_plot;        // random plot effect
    real<lower=0> sigma_plot;          // plot variance effect
    
    //neg bin overdispersion
    //real mu_phi;
    vector[nSpp] phi; 
    //real<lower=0> sigma_phi;
}

transformed parameters {
    real lambda[N]; 
    real od[N];
    real sigma_year[nSpp];

    // surv Linear Predictor
    for(n in 1:N){
    lambda[n] = beta0[spp[n]] + betasize[spp[n]]*logsize_t[n] + betaendo[spp[n]]*endo_01[n] +
    betaorigin[spp[n]]*origin_01[n]
    + tau_year[spp[n],year_t[n]] + tau_plot[plot[n]];
    
    od[n] = exp(phi[spp[n]]);
    }
    
    // endo effect on variance
    for(s in 1:nSpp){
      //for(d in 1:nEndo){
        sigma_year[s] = exp(sigma0[s]);
      //}
    }
}

model {
    // priors
    // this coefficient is shared across spp
         // prior for predictor intercepts
    //the are the grand means across spp
    //mu_beta0 ~ normal(0,100);      // prior for predictor intercepts
    //mu_betasize ~ normal(0,100);      // prior for predictor intercepts
    //mu_betaendo ~ normal(0,100);      // prior for predictor intercepts
    //mu_sigma0 ~ normal(0,100);
    //mu_sigmaendo ~ normal(0,100);
    //mu_phi ~ normal(0,100);
    //these are all spp variances
    //sigma_beta0 ~ inv_gamma(0.001, 0.001);
    //sigma_betasize ~ inv_gamma(0.001, 0.001);
    //sigma_betaendo ~ inv_gamma(0.001, 0.001);
    //sigma_mu_sigma0 ~ inv_gamma(0.001, 0.001);
    //sigma_mu_sigmaendo ~ inv_gamma(0.001, 0.001);
    //sigma_phi ~ inv_gamma(0.001, 0.001);
    //this is plot variance
    sigma_plot ~ inv_gamma(0.001, 0.001);
    for(i in 1:nPlot){
      tau_plot[i] ~ normal(0,sigma_plot);
    }
    for(s in 1:nSpp){
      beta0[s] ~ normal(0,100); 
      betasize[s] ~ normal(0,100); 
      betaendo[s] ~ normal(0,100); 
      betaorigin[s] ~ normal(0,100); 
      sigma0[s] ~ normal(0,100); 
      //sigmaendo[s] ~ normal(0,100); 
      phi[s] ~ normal(0,100);
      //for(d in 1:nEndo){
        for(t in 1:nYear){
          tau_year[s,t] ~ normal(0,sigma_year[s]);
        }
      //}
    }

    for(n in 1:N){
      y[n] ~ neg_binomial_2_log(lambda[n],od[n]);
      target += -log1m(neg_binomial_2_log_lpmf(0 | lambda[n], od[n])); 
     }
}

