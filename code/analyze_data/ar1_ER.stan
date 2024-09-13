//
//
/*----------------------- Data --------------------------*/
  data {
      int<lower = 1> N_site;        // number of sites in the model
      int<lower = 1> N_obs;         // Number of observed data points
      int<lower = 0> N_mis;         // Number of missing observations
      int<lower = 1, upper = N_obs + N_mis> ii_obs[N_obs]; // indices of observations
      int<lower = 1, upper = N_obs + N_mis> ii_mis[N_mis]; // indices of missing obs
      vector[N_obs] R_obs;           // observed GPP
      vector[N_obs + N_mis] log_GPP;   // covariate
      vector[N_obs + N_mis] temperature;     // covariate
      vector[N_obs + N_mis] log_Q;        // covariate
      int<lower = 1, upper = N_site> site[N_obs + N_mis];    // the site for each observation coded as integers
  }

/*----------------------- Transformed data---------------*/
  transformed data{
      int<lower = 0> N = N_obs + N_mis;
  }

/*----------------------- Parameters --------------------*/
  parameters {
      vector[N_mis] R_mis;          // estimates of missing observations
      vector[N_site] beta_0;        // site specific intercepts
      real<lower = 0> site_sigma;   // the variance on the site specific betas
      vector[5] beta;               // global intercept and covariate coefficients
      real<lower = 0, upper=1> phi; // Auto-regressive parameter
      real<lower=0> sigma;          // Standard deviation of the random innovations
  }

/*----------------------- Transformed parameters---------*/
  transformed parameters {
      vector[N] y;
      vector[N] y2;
      vector[N] mu;

      y[ii_obs] = R_obs;
      y[ii_mis] = R_mis;

      mu = beta_0[site] + beta[2]*log_GPP + beta[3]*temperature + beta[4]*log_Q + beta[5]*(temperature.*log_Q);

      y2 = y - mu;
  }

/*----------------------- Model -------------------------*/
  model {
      // Prior distributions
      phi ~ beta(1,1);
      beta ~ normal(0, 5);
      site_sigma ~ normal(0, 1);
      beta_0 ~ normal(beta[1], site_sigma);
      sigma ~ normal(0, 1);

      // likelihood
      for(i in 2:N){

          y2[i] ~ normal(phi * y2[i-1], sigma);

      }


  }
