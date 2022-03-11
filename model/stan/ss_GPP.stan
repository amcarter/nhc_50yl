// Stan model to predict GPP as an autoregressive function of light and temp

data {
  int<lower=0> N;
  vector[N] P_obs;
  vector[N] light;
  vector[N] tempK;
  vector[N] Q;
}

parameters {
  real<lower=0,upper=1> phi;
  real<lower=0> b_light;
  real<lower=0,upper=1> b_q;
  real<lower=0> sigma_obs;
  real<lower=0> sigma_proc;
  vector<lower=0>[N] P;
}

model {

  P[1] ~ normal(P_obs[1], sigma_proc);
  
  print(" P:", P[1], "light", light[1], "b_l", b_light, " phi ", phi);

  for(i in 2:N){
    P[i] ~ normal((phi * P[i-1] + b_light * light[i])*(1-b_q * Q[i])* exp(0.003413 - 1/tempK[i]), sigma_proc);

  }
  
  P_obs ~ normal(P, sigma_obs);
  
  
  phi ~ beta(2.975,1.275); //mean = 0.7, sd = 0.8
  b_q ~ uniform(0,1);
  b_light ~ lognormal(0,1);
  sigma_proc ~ normal(0, 0.05);
  sigma_obs ~ normal(0.08, 0.02);
}

