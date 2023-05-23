functions{
    void add_iter();
    void reset_iter();
    int get_iter();
}

data {
  int<lower=1> S; // number of supplied samples
  int<lower=1> M; // number of samples for integration
  int<lower=1> G; // number of latent classes
  int<lower=1> N;  // number of observations
  vector[N] y;  // response vector

  int<lower=1> I; //number of strata
  int<lower=1, upper=I> ii[N]; //stratum of each observation
  int<lower=1, upper=N+1> offsets[I+1]; //strata observation start indices

  int<lower=1> B_mu;  // number of fixed effects
  matrix[N, B_mu] x_mu;  // population-level design matrix

  // parameters
  vector<lower=0, upper=1>[G] theta[S];  // mixing proportions
  vector[G] intercept_mu[S];  // to identify mixtures
  matrix[B_mu, G] beta_mu[S];  // population-level effects
  real<lower=0> sdz_mu[S];
  real<lower=0> sigma[S];
}

transformed data {
  real log_M = log(M);
  reset_iter();
}

model {
  // empty
}

generated quantities {
  int s = get_iter();
  vector[I] log_lik;
  
  {
    vector[G] log_theta_s = log(theta[s]);
    vector[G] intercept_s = intercept_mu[s];
    matrix[B_mu, G] beta_mu_s = beta_mu[s];
    real sdz_mu_s = sdz_mu[s];
    real sigma_s = sigma[s];
    vector[M] log_obs_m;
    vector[G] log_lik_g;
    real r_mu[M] = normal_rng(rep_vector(0, M), sdz_mu_s);

    for (i in 1:I) {
      for (g in 1:G) {
        for (m in 1:M) {
          log_obs_m[m] = normal_id_glm_lpdf(y[offsets[i]:(offsets[i + 1] - 1)] | x_mu[offsets[i]:(offsets[i + 1] - 1), ], intercept_s[g] + r_mu[m], beta_mu_s[, g], sigma_s);
        }
        
        // compute cluster-specific log liks (mean of log_obs)
        log_lik_g[g] = log_sum_exp(log_obs_m) - log_M;
      }
      
      // compute overall mixture likelihood for the trajectory
      log_lik[i] = log_sum_exp(log_lik_g + log_theta_s);
    }
  }
  
  add_iter();
}
