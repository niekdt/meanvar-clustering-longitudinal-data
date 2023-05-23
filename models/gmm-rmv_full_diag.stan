// https://mc-stan.org/users/documentation/case-studies/identifying_mixture_models.html
// https://mc-stan.org/docs/2_18/stan-users-guide/vectorizing-mixtures.html
// Latent-class mixed model with random intercept and shared variance
// https://mc-stan.org/docs/2_19/functions-reference/composed-functions.html
// https://mc-stan.org/misc/warnings.html

// GMM with QR decomposition, diagonal vcov, random variance, and mean-variance relation
// Niek Den Teuling 2022

data {
  int<lower=1> G; // number of latent classes
  int<lower=1> N;  // number of observations
  vector[N] y;  // response vector
  
  int<lower=1> I; //number of strata
  int<lower=1> ii[N]; //stratum of each observation
  int<lower=1, upper=N+1> offsets[I+1]; //strata observation start indices
  
  int<lower=1> B_mu;  // number of fixed effects
  matrix[N, B_mu] x_mu;  // population-level design matrix
  
  int use_priors;
  int use_sdzs;
  int use_sigmas;
  
  // priors
  vector[G] prior_theta;
  real prior_cfc;
  real prior_sdzMu_nu;
  real prior_sdzMu_sigma;
  real prior_sdzSigma_nu;
  real prior_sdzSigma_sigma;
  real prior_interceptSigma_mean;
  real prior_interceptSigma_sigma;
  real prior_cv_nu;
  real prior_cv_sigma;
}

transformed data {
  real ybar = mean(y);
  vector[N] yc = y - ybar;
  vector[B_mu] xbar_mu;
  matrix[N, B_mu] xc_mu;
  matrix[N, B_mu] Q;
  matrix[B_mu, B_mu] R;
  matrix[B_mu, B_mu] R_inv;

  print("Centering...");
  for (i in 1:B_mu) {
    xbar_mu[i] = mean(x_mu[,i]);
    xc_mu[,i] = x_mu[,i] - xbar_mu[i];
  }
  print("xbar_mu:", xbar_mu);
  
  print("Computing QR decomposition with auto-scaling...");
  Q = qr_thin_Q(xc_mu) * sqrt(N - 1);
  R = qr_thin_R(xc_mu) / sqrt(N - 1);
  R_inv = inverse(R);
}

parameters {
  simplex[G] theta;  // mixing proportions
  ordered[G] alpha;  // to identify mixtures
  matrix[B_mu, G] b_mu;  // population-level effects
  
  vector<lower=0>[use_sdzs ? G : 1] sdz0_mu;
  matrix<lower=0>[use_sdzs ? G : 1, B_mu] sdz_mu;  // group-specific standard deviations
  
  vector[I] r0_mu[G];  // unscaled random effects per group
  matrix[B_mu, I] r_mu[G];  // unscaled random effects per group
  
  vector<lower=0>[I] r0_sigma[G];
  vector<lower=0>[use_sdzs ? G : 1] sdz0_sigma;
  vector[use_sigmas ? G : 1] intercept_sigma;
  vector[G] cv;
}


transformed parameters {
  matrix[N, G] mus;
  matrix<lower=0>[N, G] sigmas;
  
  for (g in 1:G) {
    mus[, g] = r0_mu[g][ii] + rows_dot_product(Q, r_mu[g][, ii]') + ybar;
    sigmas[, g] = r0_sigma[g][ii] .* exp(cv[g] * mus[,g]);
  }
}


model {
  vector[G] log_theta = log(theta);
  vector[G] ps;
  
  for (g in 1:G) {
    r0_mu[g] ~ normal(alpha[g], sdz0_mu[use_sdzs ? g : 1]);    
    r0_sigma[g] ~ lognormal(intercept_sigma[use_sigmas ? g : 1], sdz0_sigma[use_sdzs ? g : 1]);
    for (l in 1:B_mu) {
      r_mu[g][l,] ~ normal(b_mu[l,g], sdz_mu[use_sdzs ? g : 1, l]);
    }
  }
  
  for (i in 1:I) {
    for (g in 1:G) {
      ps[g] = normal_lpdf(y[offsets[i]:(offsets[i+1] - 1)] | mus[offsets[i]:(offsets[i+1] - 1), g], sigmas[offsets[i]:(offsets[i+1] - 1), g]);
    }
    target += log_sum_exp(log_theta + ps);
  }
  
  if (use_priors) {
    theta ~ dirichlet(prior_theta);
    alpha ~ std_normal();
    to_vector(b_mu) ~ std_normal();
    sdz0_mu ~ student_t(prior_sdzMu_nu, 0, prior_sdzMu_sigma);
    to_vector(sdz_mu) ~ student_t(prior_sdzMu_nu, 0, prior_sdzMu_sigma);
    sdz0_sigma ~ student_t(prior_sdzSigma_nu, 0, prior_sdzSigma_sigma);
    intercept_sigma ~ normal(prior_interceptSigma_mean, prior_interceptSigma_sigma);
    cv ~ student_t(prior_cv_nu, 0, prior_cv_sigma);
  }
}

generated quantities {
  vector[G] intercept_mu;
  matrix[B_mu, G] beta_mu;
  vector[use_sigmas ? G : 1] sigma = exp(intercept_sigma);
  matrix[G, I] pp;
  
  for (g in 1:G) {
    beta_mu[,g] = R_inv * b_mu[,g];
    intercept_mu[g] = alpha[g] - dot_product(xbar_mu, beta_mu[,g]) + ybar;
  }

  {
    vector[G] ps;
    vector[G] log_theta = log(theta);

    for (i in 1:I) {
      ps = log_theta;
      for (g in 1:G) {
        ps[g] += normal_lpdf(y[offsets[i]:(offsets[i+1] - 1)] | mus[offsets[i]:(offsets[i+1] - 1), g], sigmas[offsets[i]:(offsets[i+1] - 1), g]);
      }
      pp[,i] = exp(ps - log_sum_exp(ps));
    }
  }
}
