// https://mc-stan.org/users/documentation/case-studies/identifying_mixture_models.html
// https://mc-stan.org/docs/2_18/stan-users-guide/vectorizing-mixtures.html
// Latent-class mixed model with random intercept and shared variance
// https://mc-stan.org/docs/2_19/functions-reference/composed-functions.html
// https://mc-stan.org/misc/warnings.html
// log_sum_exp(a,b) = log(exp(a) + exp(b))

data {
  int<lower=1> G; // number of latent classes
  int<lower=1> N;  // number of observations
  vector[N] y;  // response vector
  
  int<lower=1> I; //number of strata
  int<lower=1> ii[N]; //stratum of each observation
  int<lower=1, upper=N+1> offsets[I+1]; //strata observation start indices
  
  int<lower=1> B_mu;  // number of fixed effects
  matrix[N, B_mu] x_mu;  // population-level design matrix
  
  int<lower=1> L_mu; // number of random effects
  matrix[N, L_mu] u_mu; // random effects design matrix
  
  // priors
  vector[G] prior_theta;
  real prior_cfc;
  real prior_interceptMu_mean[G];
  real prior_interceptMu_sigma[G];
  real prior_sdzMu_nu;
  real prior_sdzMu_sigma;
  real prior_sigma_nu;
  real prior_sigma_sigma;
  
  int use_priors;
  int use_sigmas;
  int use_cor;
}

transformed data {
  matrix[L_mu, L_mu] corDiag = diag_matrix(rep_vector(1, L_mu));
}

parameters {
  simplex[G] theta;  // mixing proportions
  ordered[G] intercept_mu;  // to identify mixtures
  matrix[B_mu, G] beta_mu;  // population-level effects
  matrix<lower=0>[G, L_mu] sdz_mu;  // group-specific standard deviations
  matrix[L_mu, I] z_mu[G];  // unscaled random effects per group
  real<lower=0> sigma[use_sigmas ? G : 1];  // residual SD
  cholesky_factor_corr[use_cor ? L_mu : 0] cfc;
}

transformed parameters {
  //matrix[G, N] mu;
  // for (g in 1:G) {
  //   matrix [L_mu, I] r_mu = diag_pre_multiply(sdz_mu[g], use_cor ? cfc : corDiag) * z_mu[g];
  //   mu[g, ] = (intercept_mu[g] + x_mu * beta_mu[,g] + rows_dot_product(u_mu, r_mu[, ii]'))';
  // }
}

model {
  vector[G] log_theta = log(theta);
  vector[G] ps;
  matrix[G, N] mu;
  
  for (g in 1:G) {
    matrix [L_mu, I] r_mu = diag_pre_multiply(sdz_mu[g], use_cor ? cfc : corDiag) * z_mu[g];
    mu[g, ] = (intercept_mu[g] + x_mu * beta_mu[,g] + rows_dot_product(u_mu, r_mu[, ii]'))';
  }
  
  for (i in 1:I) {
    for (g in 1:G) {
      ps[g] = normal_lpdf(y[offsets[i]:(offsets[i+1] - 1)] | mu[g, offsets[i]:(offsets[i+1] - 1)], sigma[use_sigmas ? g : 1]);
    }
    target += log_sum_exp(log_theta + ps);
  }
  
  for (g in 1:G) {
    to_vector(z_mu[g]) ~ std_normal();
  }
  
  if (use_priors) {
    theta ~ dirichlet(prior_theta);
    to_vector(sdz_mu) ~ student_t(prior_sdzMu_nu, 0, prior_sdzMu_sigma);
    sigma ~ student_t(prior_sigma_nu, 0, prior_sigma_sigma);
    cfc ~ lkj_corr_cholesky(prior_cfc);
    to_vector(beta_mu) ~ std_normal();
    intercept_mu ~ normal(prior_interceptMu_mean, prior_interceptMu_sigma);
  }
}

generated quantities {
  real log_lik;
  real pp;
  real cor;
  /*
  vector[I] log_lik;
  matrix[G, I] pp;
  vector<lower=-1, upper=1>[use_cor ? L_mu - 1 : 0] cor;
  
  {
    vector[G] ps;
    vector[G] log_theta = log(theta);
    for (i in 1:I) {
      ps = log_theta;
      for (g in 1:G) {
        //ps[g] += normal_id_glm_lpdf(y[offsets[i]:(offsets[i+1] - 1)] | Q[offsets[i]:(offsets[i+1] - 1),], alpha_i[g,i], b[,g], sigma[g]);
        ps[g] += normal_lpdf(y[offsets[i]:(offsets[i+1] - 1)] | mu[g, offsets[i]:(offsets[i+1] - 1)], sigma[use_sigmas ? g : 1]);
      }
      log_lik[i] = log_sum_exp(ps);
      pp[,i] = exp(ps - log_lik[i]);
    }
  }
  
  if (use_cor && L_mu > 1) {
    matrix[L_mu, L_mu] cormat = multiply_lower_tri_self_transpose(cfc);
    for (k in 1:L_mu) {
      for (j in 1:(k - 1)) {
        cor[choose(k - 1, 2) + j] = cormat[j, k];
      }
    }
  }
  */
}
