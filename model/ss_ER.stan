data {
  int<lower=0> N; // number of observations
  vector<lower=0>[N] R_obs; // ER
  vector<lower=0>[N] P_obs; // GPP
  vector<lower=0>[N] tempK; // water temperature in K
  vector<lower=0>[N] Q; // discharge
  vector<lower=0>[N] litter; // daily carbon input
  real<lower=0> C0; // initial carbon storage
}

parameters {
  real<lower=0> b1;
  real<lower=0,upper=1> b2;
  real<lower=0> E; // activation energy
  vector<lower=0>[N] C; // organic carbon
  real<lower=0> sigma_proc;
  real<lower=0> sigma_obs;
}

model {
  // Process model
  // initialize latent states:
  real k_b; // Boltzman's constant
  real AR; // Fraction of autotrophic respiration 
  vector[N] R; // underlying respiration
  
  k_b = 8.6173 * 10^-5;
  AR = 0.5;
  C[1] ~ normal(C0, sigma_proc);
  R[1] ~ normal(R_obs[1], sigma_obs);
  
  // iterate through timesteps:
  for(i in 2:N){
    C[i] ~ normal((C[i-1] + litter[i] - R[i-1] + P_obs[i-1])*(1 - b2*Q[i]), sigma_proc);
    R[i] = 10^8 *b1 * C[i] * exp(-E/(k_b * tempK[i])) + AR * P_obs[i];
  }
  
  // Observation Model:
  R_obs ~ normal(R, sigma_obs);
  
  // Priors:
  b1 ~ normal(0, 1);
  b2 ~ normal(0, 1);
  E ~ normal(0.65, 0.1);
  sigma_proc ~ normal(0,1);
  sigma_obs ~ normal(0,1);
}
