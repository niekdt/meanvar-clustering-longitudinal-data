sim_init()
# deleteExperiment('final')

# Define methods ####
chains = 1
seeds = 1 # 1:10
# GMM
methodGmm2 = generate_method_cases('lcMethodStanGMM', 
  model = 'gmm_full_diag', 
  nClusters = 2, 
  warmup = 2e3,
  samples = 1e3,
  adapt_delta = .7,
  use_sigmas = FALSE, 
  use_sdzs = FALSE, 
  use_cor = FALSE,
  chains = chains,
  seed = seeds
)

# MV
methodMv2 = generate_method_cases('lcMethodStanGMMMeanVar', 
  model = 'gmm-mv_full_diag', 
  nClusters = 2, 
  warmup = 2e3,
  samples = 1e3,
  adapt_delta = .9,
  use_sigmas = FALSE,
  use_sdzs = FALSE, 
  use_cor = FALSE,
  chains = chains,
  seed = seeds
)

# RV
methodRv2 = generate_method_cases('lcMethodStanGMMRV', 
  model = 'gmm-rv_full_diag', 
  nClusters = 2, 
  warmup = 2e3,
  samples = 1e3,
  adapt_delta = .9,
  use_sigmas = FALSE,
  use_sdzs = FALSE, 
  use_cor = FALSE,
  chains = chains,
  seed = seeds
)

# RMV
methodRmv2 = generate_method_cases('lcMethodStanGMMRMV', 
  model = 'gmm-rmv_full_diag', 
  nClusters = 2, 
  warmup = 2e3,
  samples = 1e3,
  adapt_delta = .95,
  use_sigmas = FALSE,
  use_sdzs = FALSE, 
  use_cor = FALSE,
  chains = chains,
  seed = seeds
)

# Define data scenarios ####
dataSeeds = 1:500
dataCases = generate_data_cases(
  numGroups = 2,
  numTraj = 250,
  numObs = 10,
  fixed = c('Partial', 'Equal'),
  random = 'Med2',
  randomSigma = c('None', 'High'),
  sigma = 'Med',
  cv = c('None', 'HighDesc'),
  seed = dataSeeds
) %>%
  .[data.fixed != 'Equal' | data.cv != 'None']

# Create combined cases
# 2 groups
cases2 = list(
  create_cases(methodMv2, dataCases),
  create_cases(methodGmm2, dataCases),
  create_cases(methodRv2, dataCases),
  create_cases(methodRmv2, dataCases)
) %>% 
  rbindlist(fill = TRUE)
print_cases(cases2)
simworkr::submitJobs('finalc', cases2)

# 3 groups ####
cases3 = copy(cases2) %>%
  .[, data.numGroups := 3] %>%
  .[, nClusters := 3]
simworkr::submitJobs('finalc', cases3)
