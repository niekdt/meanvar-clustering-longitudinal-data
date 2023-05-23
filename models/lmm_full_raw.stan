// Fast linear mixed effects model with random effects on all fixed effects
// Niek Den Teuling 2019
// https://mc-stan.org/docs/2_19/stan-users-guide/multivariate-hierarchical-priors-section.html
// https://discourse.mc-stan.org/t/fitting-a-large-data-hierarchical-model-can-i-trust-variational-inference-algorithm/1780/8

data {
  int<lower=1> N;  // number of observations
  vector[N] y;  // response variable
  int<lower=1> B;  // number of population-level effects, excluding the intercept
  matrix[N, B] x;  // population-level design matrix, excluding the intercept
  int<lower=1> I;
  int<lower=1, upper=I> ii[N];
  int<lower=1, upper=N+1> offsets[I+1]; //strata observation start indices
  
  int use_priors;
  real prior_sdz_nu;
  real prior_sdz_mu;
  real prior_sdz_sigma;
  real prior_sigma_nu;
  real prior_sigma_mu;
  real prior_sigma_sigma;
  real prior_cfc;
}

transformed data {
  int L = B + 1; // number of random effects
  vector[B] xbar;
  real ybar = mean(y);
  vector[N] yc = y - ybar;
  
  for (i in 1:B) {
    xbar[i] = mean(x[,i]);
  }
}

parameters {
  real intercept;
  vector[B] beta;  // coefficients of Q
  real<lower=0> sigma;  // residual SD
  vector<lower=0,upper=pi()/2>[L] tau_unif;
  cholesky_factor_corr[L] L_Omega;
  matrix[L, I] z;
}

transformed parameters {
  vector<lower=0>[L] tau;     // prior scale
  matrix[L, I] r;
  for (l in 1:L) {
    tau[l] = 2.5 * tan(tau_unif[l]);
  }
  
  for (i in 1:I) {
    r[, i] = diag_pre_multiply(tau, L_Omega) * z[, i];
  }
}

model {
  to_vector(z) ~ std_normal();
  
  for (i in 1:I) {
    y[offsets[i]:(offsets[i+1] - 1)] ~ normal_id_glm(x[offsets[i]:(offsets[i+1] - 1), ], intercept + r[1, i], beta + r[2:L, i], sigma);
  }
  
  if (use_priors) {
    intercept ~ std_normal();
    beta ~ std_normal();
    sigma ~ student_t(prior_sigma_nu, prior_sigma_mu, prior_sigma_sigma);
    L_Omega ~ lkj_corr_cholesky(prior_cfc);
  }
}

generated quantities {
  //real intercept = ybar - dot_product(xbar, beta);
  corr_matrix[L] Omega = L_Omega * L_Omega';
  cov_matrix[L] cov = diag_pre_multiply(tau, L_Omega) * diag_pre_multiply(tau, L_Omega)';
  vector[L] sdz = diagonal(cov);
}
