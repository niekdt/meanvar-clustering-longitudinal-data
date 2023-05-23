sim_init()
# clearJobs()

nRep = 20 # 100
ad = .8
warmup = 1500
samples = 1500

dataCases2mv = generate_data_cases(
  numGroups = 2,
  numTraj = 250,
  numObs = 10,
  fixed = c('Partial', 'Equal'),
  random = c('Low', 'High'),
  randomSigma = c('None', 'High'),
  sigma = c('Low', 'High'),
  cv = c('LowDesc', 'HighDesc'),
  seed = 1:nRep)

dataCases2n = generate_data_cases(
  numGroups = 2,
  numTraj = 250,
  numObs = 10,
  fixed = c('Partial'),
  random = c('Low', 'High'),
  randomSigma = c('None', 'High'),
  sigma = c('Low', 'High'),
  cv = 'None',
  seed = 1:nRep)

dataCases2 = rbind(dataCases2mv, dataCases2n)

# GMM ####
gmmMethodCases2 = generate_method_cases('lcMethodStanGMM', 
  model = 'gmm_full_diag', 
  nClusters = 2, 
  warmup = warmup,
  samples = samples,
  adapt_delta = ad,
  use_sigmas = FALSE, 
  use_sdzs = FALSE, 
  use_cor = FALSE,
  chains = 1)

casesGmm2 = create_cases(gmmMethodCases2, dataCases2) %T>% print_cases()
simworkr::submitJobs('sim_gmm', casesGmm2)

casesGmm3 = copy(casesGmm2) %>%
  .[, data.numGroups := 3] %>%
  .[, nClusters := 3]
simworkr::submitJobs('sim_gmm', casesGmm3)

 # GMM-MV ####
mvMethodCases2 = generate_method_cases('lcMethodStanGMMMeanVar', 
  model = 'gmm-mv_full_diag', 
  nClusters = 2, 
  warmup = warmup,
  samples = samples,
  adapt_delta = ad,
  use_sigmas = FALSE,
  use_sdzs = FALSE, 
  use_cor = FALSE,
  chains = 1)
casesMv2 = create_cases(mvMethodCases2, dataCases2) %T>% print_cases()
simworkr::submitJobs('sim_gmm_mv', casesMv2)

# 3 groups
casesMv3 = copy(casesMv2) %>%
  .[, data.numGroups := 3] %>%
  .[, nClusters := 3]
simworkr::submitJobs('sim_gmm_mv', casesMv3)


# GMM-RV ####
rvMethodCases2 = generate_method_cases('lcMethodStanGMMRV', 
  model = 'gmm-rv_full_diag', 
  nClusters = 2, 
  warmup = warmup,
  samples = samples,
  adapt_delta = ad,
  use_sigmas = FALSE,
  use_sdzs = FALSE, 
  use_cor = FALSE,
  chains = 1)
casesRv2 = create_cases(rvMethodCases2, dataCases2) %T>% print_cases()
simworkr::submitJobs('sim_gmm_rv', casesRv2)

# 3 groups
casesRv3 = copy(casesRv2) %>%
  .[, data.numGroups := 3] %>%
  .[, nClusters := 3]
simworkr::submitJobs('sim_gmm_rv', casesRv3)

# GMM-RMV ####
rmvMethodCases2 = generate_method_cases('lcMethodStanGMMRMV', 
  model = 'gmm-rmv_full_diag', 
  nClusters = 2, 
  warmup = warmup,
  samples = samples,
  adapt_delta = ad,
  use_sigmas = FALSE,
  use_sdzs = FALSE, 
  use_cor = FALSE,
  chains = 1)
casesRmv2 = create_cases(rmvMethodCases2, dataCases2) %T>% print_cases()
simworkr::submitJobs('sim_gmm_rmv', casesRmv2)

# 3 groups
casesRmv3 = copy(casesRmv2) %>%
  .[, data.numGroups := 3] %>%
  .[, nClusters := 3]
simworkr::submitJobs('sim_gmm_rmv', casesRmv3)
