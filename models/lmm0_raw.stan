// Random intercept model for Gaussian response with identity link

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
}

parameters {
  real intercept;
  vector[B] beta;  // coefficients of Q
  real<lower=0> sigma;  // residual SD
  real<lower=0> sdz;  // group-level standard deviations
  vector[I] r; // random intercepts
}

model {
  r ~ normal(intercept, sdz);
  
  y ~ normal_id_glm(x, r[ii], beta, sigma);
  
  if (use_priors) {
    intercept ~ std_normal();
    beta ~ std_normal();
    sdz ~ student_t(prior_sdz_nu, prior_sdz_mu, prior_sdz_sigma);
    sigma ~ student_t(prior_sigma_nu, prior_sigma_mu, prior_sigma_sigma);
  }
}
