sim_init()

# Methods ####
gmmMethodCases = generate_method_cases(
  'lcMethodStanGMM', 
  model = 'gmm_full_diag', 
  formula = 'Value ~0 + Time',
  nClusters = 1:3, 
  warmup = 4e3,
  samples = 1e3,
  adapt_delta = .7,
  use_sigmas = FALSE,
  use_sdzs = FALSE, 
  use_cor = FALSE,
  chains = 1,
  margLogLikIter = 5e3,
  margLogLikMixType = 'smax.alloc'
) %>% 
  .[seed <= nClusters]

mvMethodCases = generate_method_cases(
  'lcMethodStanGMMMeanVar', 
  model = 'gmm-mv_full_diag', 
  formula = 'Value ~0 + Time',
  nClusters = 1:3, 
  warmup = 4e3,
  samples = 1e3,
  adapt_delta = .7,
  use_sigmas = FALSE,
  use_sdzs = FALSE, 
  use_cor = FALSE,
  chains = 1,
  margLogLikIter = 5e3,
  margLogLikMixType = 'smax.alloc'
) %>% 
  .[seed <= nClusters]

rmvMethodCases = generate_method_cases(
  'lcMethodStanGMMRMV', 
  model = 'gmm-rmv_full_diag', 
  formula = 'Value ~0 + Time',
  nClusters = 1:3, 
  warmup = 4e3,
  samples = 1e3,
  adapt_delta = .7,
  use_sigmas = FALSE,
  use_sdzs = FALSE, 
  use_cor = FALSE,
  chains = 1,
  margLogLikIter = 5e3,
  margLogLikMixType = 'smax.alloc'
) %>% 
  .[seed <= nClusters]

# Data scenarios ####
dataCasesMV = generate_data_cases(
  numGroups = 2,
  numTraj = 250,
  numObs = 20,
  fixed = c('PartialLinear', 'EqualLinear'),
  random = 'HighLinear',
  randomSigma = 'None',
  sigma = 'High',
  cv = 'HighDesc',
  seed = 1:100
)

dataCasesRV = generate_data_cases(
  numGroups = 2,
  numTraj = 250,
  numObs = 20,
  fixed = c('PartialLinear', 'EqualLinear'),
  random = 'HighLinear',
  randomSigma = 'High',
  sigma = 'High',
  cv = 'None',
  seed = 1:100
)

gen_gtsdata(
  fixed = 'EqualLinear', 
  numGroups=2, 
  numTraj=250, 
  numObs=20, 
  random = 'HighLinear', 
  randomSigma = 'High',
  sigma = 'High', 
  cv = 'None', 
  seed = 1
) %T>% 
  plot_gtsdata() %>% 
  plotTrajectories(response = 'Value', cluster = 'Group')


# GMM
gmmCasesMV = create_cases(gmmMethodCases, dataCasesMV); nrow(gmmCasesMV)
simworkr::submitJobs('class_gmm_final2', gmmCasesMV)
gmmCasesRV = create_cases(gmmMethodCases, dataCasesRV); nrow(gmmCasesRV)
simworkr::submitJobs('class_gmm_final2', gmmCasesRV)
# MV-GMM
mvCasesMV = create_cases(mvMethodCases, dataCasesMV); nrow(mvCasesMV)
simworkr::submitJobs('class_gmmv_final2', mvCasesMV)
mvCasesRV = create_cases(mvMethodCases, dataCasesRV); nrow(mvCasesRV)
simworkr::submitJobs('class_gmmv_final2', mvCasesRV)
# RMV-GMM
rmvCasesMV = create_cases(rmvMethodCases, dataCasesMV); nrow(rmvCasesMV)
simworkr::submitJobs('class_gmmrmv_final2', rmvCasesMV)
rmvCasesRV = create_cases(rmvMethodCases, dataCasesRV); nrow(rmvCasesRV)
simworkr::submitJobs('class_gmmrmv_final2', rmvCasesRV)


