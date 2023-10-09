data {
  int<lower=0> N; // number of observations
  vector<lower=0>[N] R_obs; // ER
  vector<lower=0>[N] P_obs; // GPP
  vector<lower=0>[N] tempK; // water temperature in K
  vector<lower=0>[N] C;
}

parameters {
  real<lower=0> b1;
  real<lower=0> E; // activation energy
  real<lower=0> sigma_obs;
}

model {
  // Process model
  // initialize latent states:
  real k_b; // Boltzman's constant in eV/K
  real AR; // Fraction of GPP that is respired by autotrophic respiration 
  vector[N] R; // actual respiration (unobserved)
  
  k_b = 8.6173 * 10^-5;
  AR = 0.5;
  R[1] = R_obs[1];
  
  // iterate through timesteps:
  for(i in 2:N){
    R[i] = 10^8 *b1 * C[i] * exp(-E/(k_b * tempK[i])) + AR * P_obs[i];
  }
  
  // Observation Model:
  R_obs ~ normal(R, sigma_obs);
  
  // Priors:
  b1 ~ normal(0, 1);
  //b2 ~ normal(0, 1);
  E ~ normal(0.63, 0.1);
  //sigma_proc ~ normal(0,1);
  sigma_obs ~ normal(0,1);
}
