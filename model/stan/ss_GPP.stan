// Stan model to predict GPP as an autoregressive function of light and temp

data {
  int<lower=0> N;
  vector[N] P_obs;
  vector[N] light;
  vector[N] tempK;
  vector[N] Q; // standardized discharge
}

parameters {
  real<lower=0,upper=1> beta;
  vector<lower=0>[N] B; // biomass
  real r; // growth rate
  real lambda; //r/K - growth rate/carying capacity
  real<lower=0> P_20;
  
  real<lower=0> sigma_obs;
  real<lower=0> sigma_proc;
  
}

transformed parameters{
  real pred_P [N];
  real E;
  real k_b;
 
  k_b = 8.6173 * 10^-5;
  E = 0.35;
 
  
  //real<lower=0> P [N];

  for(i in 1:N){
    pred_P[i] = light[i]*exp(B[i]) * P_20 * exp(E/k_b * (0.003413 - 1/tempK[i]));
  }

}
model {

  B[1] ~ normal(log(P_obs[1]/light[1]), 1);
  
  for(j in 2:N){
    B[j] ~ normal((B[j-1] + r + lambda*exp(B[j-1])) * (1-beta*Q[j]), sigma_proc);

  }
  
  P_obs ~ normal(pred_P, sigma_obs);
  
  beta ~ uniform(0,1);
  sigma_proc ~ normal(0, 0.05);
  sigma_obs ~ normal(0.08, 0.02);
  
  r ~ normal(0,1);
  lambda ~ normal(0,1)T[,0];
  P_20 ~ normal(1,0.5)T[,0];
}

