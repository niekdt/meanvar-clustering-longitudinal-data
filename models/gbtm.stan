// Group-based trajectory model with group-specific sigma using QR decomposition
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
  vector[G] prior_interceptMu_mean;
  real prior_interceptMu_sigma[G];
  real prior_sigma_nu;
  real prior_sigma_sigma;
  
  int use_priors;
  int use_sigmas; // whether to use group-specific sigma (or shared)
}

transformed data {
  real ybar = mean(y);
  vector[N] yc = y - ybar;
  vector[B_mu] xbar;
  matrix[B_mu, G] xbarmat;
  matrix[N, B_mu] xc_mu;
  matrix[N, B_mu] Q;
  matrix[B_mu, B_mu] R;
  matrix[B_mu, B_mu] Rinv;

  print("Centering...");
  for (b in 1:B_mu) {
    xbar[b] = mean(x_mu[,b]);
    xc_mu[,b] = x_mu[,b] - xbar[b];
  }
  xbarmat = rep_matrix(xbar, G);
  print("xbar:", xbar);
  
  print("Computing QR decomposition with auto-scaling...");
  Q = qr_thin_Q(xc_mu) * sqrt(N - 1);
  R = qr_thin_R(xc_mu) / sqrt(N - 1);
  Rinv = inverse(R);
}

parameters {
  simplex[G] theta;  // mixing proportions
  ordered[G] alpha;  // to identify mixtures
  //ordered[G] intercept_mu;
  matrix[B_mu, G] b_mu;  // population-level effects
  real<lower=0> sigma[use_sigmas ? G : 1];  // residual SD
}

model {
  vector[G] log_theta = log(theta);
  vector[G] ps;
  for (i in 1:I) {
    for (g in 1:G) {
      ps[g] = normal_id_glm_lpdf(yc[offsets[i]:(offsets[i+1] - 1)] | Q[offsets[i]:(offsets[i+1] - 1),], alpha[g], b_mu[,g], sigma[use_sigmas ? g : 1]);
    }
    target += log_sum_exp(log_theta + ps);
  }

  // priors including all constants
  if (use_priors) {
    theta ~ dirichlet(prior_theta);
    sigma ~ student_t(prior_sigma_nu, 0, prior_sigma_sigma);
    to_vector(b_mu) ~ std_normal();
    alpha ~ std_normal();
  }
}

generated quantities {
  vector[G] intercept_mu;
  matrix[B_mu, G] beta_mu;
  vector[I] log_lik;
  vector[N] log_lik_n;
  matrix[G, I] pp;
  
  for (g in 1:G) {
    beta_mu[,g] = Rinv * b_mu[,g];
    intercept_mu[g] = alpha[g] - dot_product(xbar, beta_mu[,g]) + ybar;
  }
  
  {
    vector[G] log_theta = log(theta);
    vector[G] ps;
    
    for (n in 1:N) {
      for (g in 1:G) {
        ps[g] = normal_lpdf(yc[n] | alpha[g] + Q[n, ] * b_mu[,g], sigma[use_sigmas ? g : 1]);
      }
      log_lik_n[n] = log_sum_exp(log_theta + ps);
    }
    
    for (i in 1:I) {
      ps = log_theta;
      for (g in 1:G) {
        ps[g] += normal_id_glm_lpdf(yc[offsets[i]:(offsets[i+1] - 1)] | Q[offsets[i]:(offsets[i+1] - 1),], alpha[g], b_mu[,g], sigma[use_sigmas ? g : 1]);
      }
      log_lik[i] = log_sum_exp(ps);
      pp[,i] = exp(ps - log_lik[i]);
    }
  }
}
