// https://mc-stan.org/users/documentation/case-studies/identifying_mixture_models.html
// https://mc-stan.org/docs/2_18/stan-users-guide/vectorizing-mixtures.html
// Latent-class mixed model with random intercept and shared variance
// https://mc-stan.org/docs/2_19/functions-reference/composed-functions.html
// https://mc-stan.org/misc/warnings.html

// GMM with QR decomposition. Vcov is not yet optimized!

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
  real prior_sdzMu_nu;
  real prior_sdzMu_sigma;
  real prior_sigma_nu;
  real prior_sigma_sigma;
  
  int use_priors;
  int use_sigmas;
  int use_sdzs;
  int use_cor;
}

transformed data {
  real ybar = mean(y);
  vector[N] yc = y - ybar;
  vector[B_mu] xbar_mu;
  matrix[N, B_mu] xc_mu;
  matrix[N, B_mu] Q;
  matrix[B_mu, B_mu] R;
  matrix[B_mu, B_mu] R_inv;
  vector[L_mu] l0 = rep_vector(0, L_mu);
  matrix[L_mu, L_mu] corDiag = diag_matrix(rep_vector(1, L_mu));

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
  matrix<lower=0>[use_sdzs ? G : 1, L_mu] sdz_mu;  // group-specific standard deviations
  matrix[L_mu, I] r_mu[G];  // unscaled random effects per group
  real<lower=0> sigma[use_sigmas ? G : 1];  // residual SD
  cholesky_factor_corr[use_cor ? L_mu : 0] cfc;
}

transformed parameters {
  matrix[N, G] alpha_re;
  for (g in 1:G) {
    alpha_re[, g] = alpha[g] + rows_dot_product(u_mu, r_mu[g][, ii]');
  }
}

model {
  vector[G] log_theta = log(theta);
  vector[G] ps;
  matrix[L_mu, L_mu] cov = use_cor ? cfc : corDiag;
  matrix[L_mu, L_mu] l_mu;
  
  for (g in 1:G) {
    l_mu = diag_pre_multiply(sdz_mu[use_sdzs ? g : 1,], cov);
    for (i in 1:I) {
      r_mu[g][, i] ~ multi_normal_cholesky(l0, l_mu);
    }
  }
  
  for (i in 1:I) {
    for (g in 1:G) {
      ps[g] = normal_id_glm_lpdf(yc[offsets[i]:(offsets[i+1] - 1)] | Q[offsets[i]:(offsets[i+1] - 1),], alpha_re[offsets[i]:(offsets[i+1] - 1), g], b_mu[,g], sigma[use_sigmas ? g : 1]);
    }
    target += log_sum_exp(log_theta + ps);
  }
  
  if (use_priors) {
    theta ~ dirichlet(prior_theta);
    alpha ~ std_normal();
    to_vector(b_mu) ~ std_normal();
    to_vector(sdz_mu) ~ student_t(prior_sdzMu_nu, 0, prior_sdzMu_sigma);
    sigma ~ student_t(prior_sigma_nu, 0, prior_sigma_sigma);
    cfc ~ lkj_corr_cholesky(prior_cfc);
  }
}

generated quantities {
  vector[G] intercept_mu;
  matrix[B_mu, G] beta_mu;
  vector[I] log_lik;
  matrix[G, I] pp;
  // corr_matrix[L_mu] cor = multiply_lower_tri_self_transpose(cfc);
  // cov_matrix[L_mu] cov = quad_form_diag(cor, sdz_mu);
  vector<lower=-1, upper=1>[use_cor ? L_mu - 1 : 0] cor;
  
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
        ps[g] += normal_id_glm_lpdf(yc[offsets[i]:(offsets[i+1] - 1)] | Q[offsets[i]:(offsets[i+1] - 1),], alpha_re[offsets[i]:(offsets[i+1] - 1), g], b_mu[,g], sigma[use_sigmas ? g : 1]);
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
}
