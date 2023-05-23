functions{
    void add_iter();
    void reset_iter();
    int get_iter();
}

data {
  int<lower=1> S; // number of supplied samples
  int<lower=1> M; // number of samples for integration
  int<lower=1> N;  // number of observations
  vector[N] y;  // response variable
  int<lower=1> B;  // number of population-level effects, excluding the intercept
  matrix[N, B] x;  // population-level design matrix, excluding the intercept
  int<lower=1> I;
  int<lower=1, upper=I> ii[N];
  int<lower=1, upper=N+1> offsets[I+1]; //strata observation start indices
  real alpha[S];
  vector[B] b[S];
  real<lower=0> sdz0[S];
  vector<lower=0>[B] sdz[S];
  real<lower=0> sigma[S];
}

transformed data {
  real log_M = log(M);
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
  
  reset_iter();
}

model {
  // empty
}

generated quantities {
  int s = get_iter();
  vector[I] log_lik;
  
  {
    real alpha_s = alpha[s];
    vector[B] b_s = b[s];
    real sdz0_s = sdz0[s];
    vector[B] sdz_s = sdz[s];
    real sigma_s = sigma[s];
    vector[M] log_obs_m;
    real r0[M] = normal_rng(rep_vector(0, M), sdz0_s);
    matrix[B, M] r;
    for (ib in 1:B) {
      r[ib, ] = to_vector(normal_rng(rep_vector(0, M), rep_vector(sdz_s[ib], M)))';
    }

    for (i in 1:I) {
      for (m in 1:M) {
        log_obs_m[m] = normal_id_glm_lpdf(yc[offsets[i]:(offsets[i + 1] - 1)] | Q[offsets[i]:(offsets[i + 1] - 1), ], alpha_s + r0[m], b_s + r[, m], sigma_s);
      }
      
      log_lik[i] = log_sum_exp(log_obs_m) - log_M;
    }
  }
  
  add_iter();
}
