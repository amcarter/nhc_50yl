data {
  int<lower=0> N; // number of observations
  vector<upper=0>[N] R_obs; // ER
  vector<lower=0>[N] P_obs; // GPP
  vector<lower=0>[N] tempK; // water temperature in K
  vector<lower=0>[N] Q; // discharge
  vector<lower=0>[N] litter; // daily carbon input
  real<lower=0> C0; // initial carbon storage
}

parameters {
  //real<lower=0> R20;
  real<lower=0,upper=1> beta;
  real<lower=0> sigma_proc;
  real<lower=0> sigma_obs;
  vector<lower=0>[N] C; // organic carbon
}

transformed parameters{
  
}

model {
  // Process model
  real k_b;     // Boltzman's constant in eV/K
  real AR;      // Fraction of GPP that is respired by autotrophic respiration 
  real E;       // activation energy
  real K_20;
  vector[N] R;  // actual respiration (unobserved)
  vector[N] K;  // decay rate of organic carbon
  
  k_b = 8.6173 * 10^-5;
  AR = 0.5;
  E = 0.63;
  K_20 = 0.01;
  
  //initialize
  C[1] ~ lognormal(log(C0), sigma_proc);
  R[1] = R_obs[1];
  
  
  // iterate through timesteps:
  for(i in 2:N){
    C[i] ~ lognormal(log((C[i-1] + litter[i] + R[i-1] + P_obs[i-1])*(1 - beta*Q[i])), sigma_proc);
    K[i] = K_20 * exp(-E/k_b * (1/tempK[i] - 0.003412969));
    R[i] = -AR * P_obs[i] - K[i] * C[i];
  }
  
  // Observation Model:
  R_obs ~ normal(R, sigma_obs);
  
  // Priors:
  beta ~ uniform(0, 1);
  sigma_obs ~ normal(0.08, 0.02);
  sigma_proc ~ normal(0,0.05);
  C ~ lognormal(4.5, 0.75);
}
