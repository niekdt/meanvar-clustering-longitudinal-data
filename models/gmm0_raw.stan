// random-intercept GMM
// https://mc-stan.org/docs/2_21/reference-manual/ordered-vector.html

data {
  int<lower=1> G; // number of latent classes
  int<lower=1> N;  // number of observations
  vector[N] y;  // response vector
  int<lower=1> B_mu;  // number of fixed effects
  matrix[N, B_mu] x_mu;  // population-level design matrix
  int<lower=1> I; //number of strata
  int<lower=1> ii[N]; //stratum of each observation
  int<lower=1, upper=N+1> offsets[I+1]; //strata observation start indices
  
  vector[G] prior_theta;  // prior concentration
  real prior_interceptMu_mean[G];
  real prior_interceptMu_sigma[G];
  real prior_sdzMu_nu;
  real prior_sdzMu_sigma;
  real prior_sigma_nu;
  real prior_sigma_sigma;
  
  int use_priors;
}

parameters {
  simplex[G] theta;  // mixing proportions
  ordered[G] intercept_mu;
  matrix[B_mu, G] beta_mu;  // population-level effects
  real<lower=0> sdz_mu;
  vector[I] r_mu[G];  // random intercepts
  real<lower=0> sigma;  // residual SD
}

model {
  vector[G] log_theta = log(theta);
  vector[G] ps;
  
  for (g in 1:G) {
    r_mu[g] ~ normal(intercept_mu[g], sdz_mu);
  }
  
  for (i in 1:I) {
    for (g in 1:G) {
      ps[g] = normal_id_glm_lpdf(y[offsets[i]:(offsets[i+1] - 1)] | x_mu[offsets[i]:(offsets[i+1] - 1),], r_mu[g][i], beta_mu[,g], sigma);
    }
    target += log_sum_exp(log_theta + ps);
  }
  
  if (use_priors) {
    theta ~ dirichlet(prior_theta);
    intercept_mu ~ normal(prior_interceptMu_mean, prior_interceptMu_sigma);
    to_vector(beta_mu) ~ std_normal();
    sdz_mu ~ student_t(prior_sdzMu_nu, 0, prior_sdzMu_sigma);
    sigma ~ student_t(prior_sigma_nu, 0, prior_sigma_sigma);
  }
}

generated quantities {
  matrix[G, I] pp;
  matrix[N, G] mus;
  
  for (g in 1:G) {
    mus[, g] = r_mu[g][ii] + x_mu * beta_mu[,g];
  }

  {
    vector[G] log_theta = log(theta);
    vector[G] ps;
    vector[G] pst;
    
    for (i in 1:I) {
      for (g in 1:G) {
        ps[g] = normal_id_glm_lpdf(y[offsets[i]:(offsets[i+1] - 1)] | x_mu[offsets[i]:(offsets[i+1] - 1),], r_mu[g][i], beta_mu[,g], sigma);
      }
      pst = ps + log_theta;
      pp[, i] = exp(pst - log_sum_exp(pst));
    }
  }
}
