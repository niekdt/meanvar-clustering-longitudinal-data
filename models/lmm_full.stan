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
  matrix[N, B] xc;
  matrix[N, B] Q;
  matrix[B, B] R;
  matrix[B, B] R_inv;
  
  // centering
  print("Centering...");
  for (i in 1:B) {
    xbar[i] = mean(x[,i]);
    xc[,i] = x[,i] - xbar[i];
  }
  
  // QR decomposition
  print("Computing QR decomposition with auto-scaling...");
  Q = qr_thin_Q(xc) * sqrt(N - 1);
  R = qr_thin_R(xc) / sqrt(N - 1);
  R_inv = inverse(R);
}

parameters {
  vector[B] b;  // coefficients of Q
  real<lower=0> sigma;  // residual SD
  real<lower=0> sdz;  // group-level standard deviations
  vector[I] r; // random intercepts
  cholesky_factor_corr[L] cfc;
  matrix[L, L] cov;
}

model {
  r ~ normal(0, sdz);
  
  yc ~ normal_id_glm(Q, r[ii], b, sigma);
  
  if (use_priors) {
    b ~ std_normal();
    sdz ~ student_t(prior_sdz_nu, prior_sdz_mu, prior_sdz_sigma);
    sigma ~ student_t(prior_sigma_nu, prior_sigma_mu, prior_sigma_sigma);
    cfc ~ lkj_corr_cholesky(prior_cfc);
  }
}

generated quantities {
  vector[B] beta = R_inv * b;
  real intercept = ybar - dot_product(xbar, beta);
  matrix[L, L] cov = diag_pre_multiply(tau, L_Omega) * diag_pre_multiply(tau, L_Omega)';
}
