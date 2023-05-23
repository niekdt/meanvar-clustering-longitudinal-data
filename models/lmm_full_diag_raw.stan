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
  //real intercept;
  vector[B] beta;  // fixed effects coefficients
  real<lower=0> sigma;  // residual SD
  real<lower=0> sdz0;  // group-level standard deviations
  vector<lower=0>[B] sdz;
  vector[I] r0; // random intercepts
  matrix[B, I] r;
}

model {
  r0 ~ normal(0, sdz0);
  for (ib in 1:B) {
    r[ib, ] ~ normal(beta[ib], sdz[ib]);
  }
  
  for (i in 1:I) {
    y[offsets[i]:(offsets[i+1] - 1)] ~ normal_id_glm(x[offsets[i]:(offsets[i+1] - 1), ], r0[i], r[, i], sigma);
  }
  
  
  if (use_priors) {
    //intercept ~ std_normal();
    beta ~ std_normal();
    sdz0 ~ student_t(prior_sdz_nu, prior_sdz_mu, prior_sdz_sigma);
    sdz ~ student_t(prior_sdz_nu, prior_sdz_mu, prior_sdz_sigma);
    sigma ~ student_t(prior_sigma_nu, prior_sigma_mu, prior_sigma_sigma);
  }
}

generated quantities {
  real intercept = ybar - dot_product(xbar, beta);
}
