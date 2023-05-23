sim_init()
# clearJobs()

# covidExp = 'logcovid0'
covidExp = 'logcovid0-final'
generate_covid_data_cases = function(state, minDate, maxDate, excludeZeroCounties, minNormNew) {
  mc = stackoverflow::match.call.defaults()
  args = as.list(mc)[-1]
  
  dataCases = do.call(expand.grid, args) %>% 
    as.data.table() %>%
    setnames(paste0('data.', names(args)))
  dataCases[]
}

dataCases = generate_covid_data_cases(
  state = 'all',
  minDate = '2020-06-01',
  maxDate = '2020-09-13',
  excludeZeroCounties = TRUE,
  minNormNew = .5
)

# Methods ####
{
  f = 'LogNormNewConfirmed ~ 0 + bs(Time)'
  nClus = 1:4
  ad = .7 #c(.7, .9)
  use_sigmas = c(FALSE)
  warmup = 3e3
  samples = 1e3
  seeds = 1:10
}

# GMM ####
gmmCases = generate_method_cases(
  'lcMethodStanGMM', 
  model = 'gmm_full_diag', 
  nClusters = nClus, 
  warmup = warmup,
  samples = samples,
  adapt_delta = .7,
  use_sigmas = TRUE,
  use_sdzs = FALSE, 
  use_cor = FALSE,
  chains = 1,
  margLogLikIter = 5e3,
  margLogLikMixType = 'smax.alloc',
  margLogLikCache = FALSE,
  maxFits = 1,
  seed = seeds) %>%
  .[, formula := rep(f)] %>% 
  .[nClusters > 1 | seed == 1] %>%
  .[]

allGmmCases = create_cases(
  methodCases = gmmCases,
  dataCases = dataCases, 
  dataFun = 'covidWeekData'
) %T>% print_cases()

simworkr::submitJobs(covidExp, allGmmCases)

# GMM-MV ####
gmmMvCases = generate_method_cases(
  'lcMethodStanGMMMeanVar', 
  model = 'gmm-mv_full_diag', 
  nClusters = nClus, 
  warmup = warmup,
  samples = samples,
  adapt_delta = ad,
  use_sigmas = use_sigmas,
  use_sdzs = FALSE, 
  use_cor = FALSE,
  chains = 1,
  margLogLikIter = 5e3,
  margLogLikMixType = 'smax.alloc',
  margLogLikCache = FALSE,
  maxFits = 1,
  seed = seeds,
) %>%
  .[, formula := rep(f)] %>% 
  .[nClusters > 1 | seed == 1] %>%
  .[]

allGmmMvCases = create_cases(
  methodCases = gmmMvCases,
  dataCases = dataCases, 
  dataFun = 'covidWeekData'
) %T>% print_cases()

simworkr::submitJobs(covidExp, allGmmMvCases)

# GMM-RV ####
# mloglik not implemented yet
# gmmRvCases = generate_method_cases(
#   'lcMethodStanGMMRV', 
#   model = 'gmm-rv_full_diag', 
#   nClusters = nClus, 
#   warmup = warmup,
#   samples = samples,
#   adapt_delta = ad,
#   use_sigmas = use_sigmas,
#   use_sdzs = FALSE, 
#   use_cor = FALSE,
#   chains = 1,
#   margLogLikIter = 5e3,
#   margLogLikMixType = 'smax.alloc',
#   margLogLikCache = FALSE,
#   maxFits = 1,
#   seed = seeds,
# ) %>%
#   .[, formula := rep(f)] %>% 
#   .[nClusters > 1 | seed == 1] %>%
#   .[]
# 
# allGmmRvCases = create_cases(
#   methodCases = gmmRvCases,
#   dataCases = dataCases, 
#   dataFun = 'covidWeekData'
# ) %T>% print_cases()
# 
# simworkr::submitJobs(covidExp, allGmmRvCases)

# GMM-RMV ####
gmmRmvCases = generate_method_cases(
  'lcMethodStanGMMRMV', 
  model = 'gmm-rmv_full_diag', 
  nClusters = nClus, 
  warmup = warmup,
  samples = samples,
  adapt_delta = ad,
  use_sigmas = use_sigmas,
  use_sdzs = FALSE, 
  use_cor = FALSE,
  chains = 1,
  margLogLikIter = 5e3,
  margLogLikMixType = 'smax.alloc',
  margLogLikCache = FALSE,
  maxFits = 1,
  seed = seeds
) %>%
  .[, formula := rep(f)] %>% 
  .[nClusters > 1 | seed == 1] %>%
  .[]

allGmmRmvCases = create_cases(
  methodCases = gmmRmvCases, 
  dataCases = dataCases, 
  dataFun = 'covidWeekData'
) %T>% print_cases()

simworkr::submitJobs(covidExp, allGmmRmvCases)
