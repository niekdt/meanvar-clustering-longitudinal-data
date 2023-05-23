// Group-based trajectory model with group-specific sigma
// Niek Den Teuling 2020
data {
  int<lower=1> G; // number of latent classes
  int<lower=1> N;  // number of observations
  vector[N] y;  // response vector
  int<lower=1> I; //number of strata
  int<lower=1, upper=N+1> offsets[I+1]; //strata observation start indices
  
  int<lower=1> B_mu;  // number of fixed effects
  matrix[N, B_mu] x_mu;  // population-level design matrix
  
  vector[G] prior_theta;  // prior concentration
  real prior_interceptMu_mean[G];
  real prior_interceptMu_sigma[G];
  real prior_sigma_nu;
  real prior_sigma_sigma;
  
  int use_priors;
  int use_sigmas; // whether to use group-specific sigma (or shared)
}

parameters {
  simplex[G] theta;  // mixing proportions
  ordered[G] intercept_mu;  // to identify mixtures
  matrix[B_mu, G] beta_mu;  // population-level effects
  real<lower=0> sigma[use_sigmas ? G : 1];  // residual SD
}

model {
  vector[G] log_theta = log(theta);
  vector[G] ps;
  for (i in 1:I) {
    for (g in 1:G) {
      ps[g] = normal_id_glm_lpdf(y[offsets[i]:(offsets[i+1] - 1)] | x_mu[offsets[i]:(offsets[i+1] - 1),], intercept_mu[g], beta_mu[,g], sigma[use_sigmas ? g : 1]);
    }
    target += log_sum_exp(log_theta + ps);
  }

  // priors including all constants
  if (use_priors) {
    theta ~ dirichlet(prior_theta);
    sigma ~ student_t(prior_sigma_nu, 0, prior_sigma_sigma);
    to_vector(beta_mu) ~ std_normal();
    intercept_mu ~ normal(prior_interceptMu_mean, prior_interceptMu_sigma);
  }
}

generated quantities {
  vector[I] log_lik;
  matrix[G, I] pp;
  
  {
    vector[G] log_theta = log(theta);
    vector[G] ps;
    for (i in 1:I) {
      ps = log_theta;
      for (g in 1:G) {
        ps[g] += normal_id_glm_lpdf(y[offsets[i]:(offsets[i+1] - 1)] | x_mu[offsets[i]:(offsets[i+1] - 1),], intercept_mu[g], beta_mu[,g], sigma[use_sigmas ? g : 1]);
      }
      log_lik[i] = log_sum_exp(ps);
      pp[,i] = exp(ps - log_lik[i]);
    }
  }
}
