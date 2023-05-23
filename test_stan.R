# Script for quickly testing Stan and whether we can compile models on this system

library(rstan)

message('== Stan self-test ==')

message('Quick compilation test...')
example(stan_model, package = "rstan", run.dontrun = TRUE)

# workaround...
Sys.setenv(PKG_LIBS = Sys.getenv("LOCAL_LIBS"))

example(stan_model, package = "rstan", run.dontrun = TRUE)

message('Trying to compile Stan model...')
cmod2 = stan_model(
  file = file.path(getwd(), 'models', 'gmm.stan'), 
  allow_undefined = TRUE,
  auto_write = TRUE,
  save_dso = TRUE,
  verbose = TRUE
)

message('Trying to compile stan model with external code...')
cmod = stan_model(
  file = file.path(getwd(), 'loglik', 'gmm_full_diag_mloglik_z.stan'), 
  allow_undefined = TRUE,
  auto_write = TRUE,
  save_dso = TRUE,
  verbose = TRUE,
  includes = sprintf('\n#include "%s"\n', file.path(getwd(), 'loglik', 'get_iter.hpp'))
)

message('Looks good!')
