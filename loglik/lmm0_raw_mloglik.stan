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
  real intercept[S];
  vector[B] beta[S];
  real<lower=0> sdz[S];
  real<lower=0> sigma[S];
  real llMin;
  vector[M] z; // the standard random draws for approximating the loglik
}

transformed data {
  real log_M = log(M);
  reset_iter();
}

model {
  // empty
}

generated quantities {
  int<lower=1, upper=S> s = get_iter();
  vector[I] log_lik;

  {
    vector[B] beta_s = beta[s];
    real sigma_s = sigma[s];
    vector[M] log_obs_m;
    real log_obs_im;
    vector[M] r = intercept[s] + z * sdz[s];
    //real r[M] = normal_rng(rep_vector(intercept[s], M), sdz[s]);

    for (i in 1:I) {
      for (m in 1:M) {
        log_obs_im = normal_id_glm_lpdf(y[offsets[i]:(offsets[i + 1] - 1)] | x[offsets[i]:(offsets[i + 1] - 1), ], r[m], beta_s, sigma_s);
        log_obs_m[m] = log_obs_im > llMin ? log_obs_im : llMin;
      }
      log_lik[i] = log_sum_exp(log_obs_m) - log_M;
    }
  }
  
  add_iter();
}
