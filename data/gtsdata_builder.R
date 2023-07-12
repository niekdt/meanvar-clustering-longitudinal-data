library(data.table)
library(assertthat)
library(magrittr)

#' @param numObs Shared number of observations
#' @param fixed Shared fixed effects
#' @param tStart Shared start time
#' @param tEnd Shared end time
#' @param fe List of group-specific fixed effects
#' @param sigma Scalar (or vector) of (group-specific) sigma.
gtsdata_create_repmeas = function(
    sizes, 
    numObs, 
    fixed = ~ Time,
    tStart = 0, 
    tEnd = 1,
    groupNames = LETTERS[seq_along(sizes)]
) {
  assert_that(
    length(sizes) > 0,
    all(sapply(sizes, is.count)),
    all(sizes > 0),
    is.count(numObs),
    is.number(tStart),
    is.number(tEnd),
    tStart < tEnd,
    is.character(groupNames),
    length(groupNames) == length(sizes),
    uniqueN(groupNames) == length(sizes)
  )
  
  G = length(sizes)

  gengroup = function(size_g, name_g) {
    tsdata_create_repmeas(group=name_g, numTraj=size_g, numObs=numObs, tStart=tStart, tEnd=tEnd)
  }
  
  gtsdata = mapply(
    gengroup, 
    size_g = sizes,
    name_g = groupNames,
    SIMPLIFY = FALSE
  ) %>% 
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
  assert_that(
    is.numeric(sigma),
    all(sigma >= 0),
    all(gtsdata$Sigma > 0)
  )
  
  gtsdata_apply(gtsdata, FUN = function(data, sigma_g) {
      tsdata_add_random_sigma_intercept(data, sigma = sigma_g)
    }, 
    sigma_g = sigma
  )
}

#' @title Add group-specific heteroskedasticity to the data
#' @param cv Group-specific coefficient of variation, or scalar
gtsdata_add_meanVariance = function(gtsdata, cv) {
  assert_that(
    is.numeric(cv),
    all(gtsdata$Sigma > 0)
  )

  gtsdata_apply(gtsdata, FUN = function(data, cv_g) {
      data[, Sigma := exp(log(Sigma) + cv_g * Mu)]
    }, 
    cv_g = cv
  )
}
