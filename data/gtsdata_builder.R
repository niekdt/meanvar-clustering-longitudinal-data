library(mvnfast)

#' @param numObs Shared number of observations
#' @param fixed Shared fixed effects
#' @param tStart Shared start time
#' @param tEnd Shared end time
#' @param fe List of group-specific fixed effects
#' @param sigma Scalar (or vector) of (group-specific) sigma.
gtsdata_create_repmeas = function(sizes, 
                           numObs, 
                           fixed = ~ Time,
                           tStart = 0, 
                           tEnd = 1,
                           groupNames = LETTERS[seq_along(sizes)]) {
  assert_that(length(sizes) > 0,
              all(sapply(sizes, is.count)),
              all(sizes > 0),
              is.count(numObs),
              is.number(tStart),
              is.number(tEnd),
              tStart < tEnd,
              is.character(groupNames),
              length(groupNames) == length(sizes),
              uniqueN(groupNames) == length(sizes))
  
  G = length(sizes)

  gengroup = function(size_g, name_g) {
    tsdata_create_repmeas(group=name_g, numTraj=size_g, numObs=numObs, tStart=tStart, tEnd=tEnd)
  }
  
  gtsdata = mapply(gengroup, 
                  size_g = sizes,
                  name_g = groupNames,
                  SIMPLIFY = FALSE) %>% 
    tsdata_merge()
  
  return(gtsdata[])
}

gtsdata_apply = function(gtsdata, FUN, ...) {
  groupDataList = split(gtsdata, by='Group')
  groupTrends = lapply(groupDataList, function(tsdata) {
    trendData = tsdata_trends(tsdata) %>%
      .[Group %in% unique(tsdata$Group)]
  })
  mapply(setattr, groupDataList, 'trend', groupTrends)
  
  G = length(groupDataList)
  args = list(...)
  assert_that(all(sapply(args, length) %in% c(1, G)), msg='some gtsdata input arguments have the wrong length (should be 1 or G))')
  
  newdata = do.call(mapply, c(FUN = FUN, data = list(groupDataList), args, SIMPLIFY = FALSE)) %>%
    tsdata_merge()
  return(newdata)
}

# Mu ####
gtsdata_add_fixed = function(gtsdata, formula, coefs, what='Mu') {
  gtsdata_apply(gtsdata, FUN = function(data, size_g, fe_g) {
      tsdata_add_fixed(data, formula, coefs = fe_g, what = what)
    }, 
    fe_g = coefs
  )
}

#' @title Add random effects to the data
#' @param covMat covariance matrix or list of group-specific covariance matrices
gtsdata_add_random = function(gtsdata, formula, covMat) {
  if(is.matrix(covMat)) {
    covMat = list(covMat)
  }
  
  gtsdata_apply(gtsdata, FUN = function(data, size_g, cov_g) {
    tsdata_add_random(data, formula, cov_g)
  }, 
  size_g = gtsdata[, .(N = uniqueN(Id)), keyby=Group]$N,
  cov_g = covMat
  )
}

# Sigma ####
#' @param sigma standard deviation for the random effects
gtsdata_add_random_sigma_intercept = function(gtsdata, sigma) {
  assert_that(is.numeric(sigma),
              all(sigma >= 0),
              all(gtsdata$Sigma > 0))
  
  gtsdata_apply(gtsdata, FUN = function(data, sigma_g) {
      tsdata_add_random_sigma_intercept(data, sigma = sigma_g)
    }, 
    sigma_g = sigma
  )
}

#' @title Add group-specific heteroskedasticity to the data
#' @param cv Group-specific coefficient of variation, or scalar
gtsdata_add_meanVariance = function(gtsdata, cv) {
  assert_that(is.numeric(cv),
              all(gtsdata$Sigma > 0))

  gtsdata_apply(gtsdata, FUN = function(data, cv_g) {
      data[, Sigma := exp(log(Sigma) + cv_g * Mu)]
    }, 
    cv_g = cv
  )
}


# old
library(mvnfast)
library(gamlss.dist)
library(R.utils)
gen_tsdata_multilevel = function(seed=NULL, 
                                 numTraj=100, 
                                 numObs=10, 
                                 tStart=0,
                                 tEnd=1,
                                 mu=~1 + I(Time^2),
                                 muCoefs=c(1, -1),
                                 muRandom=~1,
                                 muRandomCov=.1^2,
                                 muLink=I,
                                 sigma=~1,
                                 sigmaCoefs=.01,
                                 sigmaRandom=~1,
                                 sigmaRandomCov=.1^2,
                                 sigmaLink=exp,
                                 response=rNO,
                                 name='All',
                                 ... #additional arguments passed to response()
) {
  set.seed(seed)
  muGroupR = rmvn(numTraj, mu=0, sigma=muRandomCov)
  sigmaGroupR = rmvn(numTraj, mu=0, sigma=sigmaRandomCov)
  
  tsdata = tsdata_create_repmeas(group=name, numTraj=numTraj, numObs=numObs, tStart=tStart, tEnd=tEnd) %>%
    tsdata_add_formula_fixed(mu, coefs=muCoefs, what='muRaw') %>%
    tsdata_add_formula_random(muRandom, tsCoefs=muGroupR, what='muRaw') %>%
    tsdata_add_formula_fixed(sigma, coefs=sigmaCoefs, what='sigmaRaw') %>%
    tsdata_add_formula_random(sigmaRandom, tsCoefs=sigmaGroupR, what='sigmaRaw')
  
  tsdata[, mu := muLink(muRaw)]
  tsdata[, sigma := sigmaLink(sigmaRaw)]
  tsdata[, Value := doCall(response, n=.N, args=c(.SD, list(...)))]
  
  return(tsdata)
}
