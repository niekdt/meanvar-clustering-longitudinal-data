message('Hello.')
message(sprintf('Initial working directory: %s', getwd()))

if(.Platform$OS.type == 'unix') {
    workerId = as.integer(Sys.getenv('PBS_ARRAYID'))
    message(sprintf('-- Worker %d --', workerId))
}

queue = Sys.getenv('JOBQUEUE')

source('include.R')

message('Stan self-test')
cmod = stan_model(
  file = file.path(getwd(), 'loglik', 'gmm_full_diag_mloglik_z.stan'), 
  allow_undefined = TRUE,
  auto_write = TRUE,
  save_dso = TRUE,
  includes = sprintf('\n#include "%s"\n', file.path(getwd(), 'loglik', 'get_iter.hpp'))
)
message('Stan self-test success')

sim_init()

message(sprintf('Active on job queue "%s"', queue))

simworkr(queue = queue)

message('No jobs left.')
withTimeout(disconnect(), timeout=5)
message('Goodbye.')
