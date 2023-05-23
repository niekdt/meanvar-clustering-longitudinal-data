gen_gtsdata = function(numGroups = 2, 
                       fixed = 'Distinct', 
                       random = 'Low',
                       sigma = 'None', 
                       randomSigma = 'None',
                       cv = 'None', 
                       numTraj = 250,
                       numObs = 20,
                       tStart = 0,
                       tEnd = 2,
                       seed = 1) {
  set.seed(seed)
  assert_that(has_name(dataOptions, as.character(numGroups)))
  opts = dataOptions[[as.character(numGroups)]]
  
  props = sqrt(1:numGroups) / sum(sqrt(1:numGroups))
  sizes = round(props * numTraj)
  
  if(sum(sizes) != numTraj) {
    sizes[which.max(sizes)] = sizes[which.max(sizes)] + (numTraj - sum(sizes))
  }
  assert_that(is.count(numTraj),
              is.count(numGroups),
              sum(sizes) == numTraj,
              is.count(numObs),
              is.numeric(tStart),
              is.numeric(tEnd),
              is.character(fixed) || is.factor(fixed),
              is.character(random) || is.factor(random),
              is.character(sigma) || is.factor(sigma),
              is.character(randomSigma) || is.factor(randomSigma),
              is.character(cv) || is.factor(cv),
              has_name(opts$fixed, fixed),
              has_name(opts$random, random),
              has_name(opts$sigma, sigma),
              has_name(opts$randomSigma, randomSigma),
              has_name(opts$cv, cv))
  
  data = gtsdata_create_repmeas(sizes = sizes, numObs = numObs, tStart = tStart, tEnd = tEnd) %>%
    gtsdata_add_fixed(formula = ~ Time + I(Time^2), coefs = opts$fixed[[as.character(fixed)]]) %>%
    gtsdata_add_random(formula = ~ Time + I(Time^2), covMat = opts$random[[as.character(random)]]) %>%
    gtsdata_add_fixed(formula = ~ 1, coefs = opts$sigma[[as.character(sigma)]], what='Sigma') %>%
    gtsdata_add_random_sigma_intercept(sigma = opts$randomSigma[[as.character(randomSigma)]]) %>%
    gtsdata_add_meanVariance(cv = opts$cv[[as.character(cv)]]) %>%
    tsdata_response(family = rNO)
  
  setattr(data, 'args', match.call()[-1] %>% as.list())
  
  return(data)
}
