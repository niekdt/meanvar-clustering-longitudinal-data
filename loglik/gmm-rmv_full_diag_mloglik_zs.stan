functions{
    void add_iter();
    void reset_iter();
    int get_iter();
}

data {
  int<lower=1> S; // number of supplied samples
  int<lower=1> M; // number of samples for integration
  int<lower=1> G; // number of latent classes
  int<lower=1> N;  // number of observations
  vector[N] y;  // response vector

  int<lower=1> I; //number of strata
  int<lower=1, upper=I> ii[N]; //stratum of each observation
  int<lower=1, upper=N+1> offsets[I+1]; //strata observation start indices

  int<lower=1> B_mu;  // number of fixed effects
  matrix[N, B_mu] x_mu;  // population-level design matrix

  matrix[M, 1 + B_mu] z; // standard random draws for integration
  int g_is[I,S]; // class allocations per sample

  int use_sdzs;
  int use_sigmas;

  // parameters
  vector<lower=0, upper=1>[G] theta[S];  // mixing proportions
  vector[G] alpha[S];  // to identify mixtures
  matrix[B_mu, G] b_mu[S];  // population-level effects
  vector<lower=0>[1] sdz0_mu[S];
  matrix<lower=0>[1, B_mu] sdz_mu[S];
  
  vector[1] intercept_sigma[S];
  vector<lower=0>[1] sdz0_sigma[S];
  vector[G] cv[S];
}

transformed data {
  real log_M = log(M);
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
  
  reset_iter();
}

model {
  // empty
}

generated quantities {
  int s = get_iter();
  vector[I] log_lik;
  
  {
    int g_i[I] = g_is[,s];
    vector[G] log_theta_s = log(theta[s]);
    vector[G] alpha_s = alpha[s];
    matrix[B_mu, G] b_mu_s = b_mu[s];
    real sdz0_mu_s = sdz0_mu[s][1];
    vector[B_mu] sdz_mu_s = sdz_mu[s][1, ]';
    real intercept_sigma_s = intercept_sigma[s][1];
    real sdz0_sigma_s = sdz0_sigma[s][1];
    vector[G] cv_s = cv[s];
    vector[M] log_obs_m;
    vector[M] r0_mu = z[, 1] * sdz0_mu_s;
    real r0_sigma[M] = lognormal_rng(rep_vector(intercept_sigma_s, M), sdz0_sigma_s);
    matrix[B_mu, M] r_mu;
    int g;
    
    for (b in 1:B_mu) {
      r_mu[b, ] = z[, b + 1]' * sdz_mu_s[b];
    }

    for (i in 1:I) {
      int N_i = offsets[i + 1] - offsets[i];
      vector[N_i] ys = y[offsets[i]:(offsets[i+1] - 1)];
      vector[N_i] mus;
      vector[N_i] sigmas;
      g = g_i[i];
      
      for (m in 1:M) {
        mus = ybar + alpha_s[g] + r0_mu[m] + Q[offsets[i]:(offsets[i+1] - 1), ] * (b_mu_s[, g] + r_mu[, m]);
        sigmas = r0_sigma[m] * exp(cv_s[g] * mus);
        log_obs_m[m] = normal_lpdf(ys | mus, sigmas);
      }
      
      // compute overall mixture likelihood for the trajectory
      log_lik[i] = log_sum_exp(log_obs_m) - log_M + log_theta_s[g];
    }
  }
  
  add_iter();
}
