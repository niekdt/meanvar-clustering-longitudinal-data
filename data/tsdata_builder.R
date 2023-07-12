library(data.table)
library(assertthat)
library(magrittr)

#' @description  Creates dataset comprising repeated measures data (regularly-spaced time series)
#' @param group Name for this group of time series
#' @param numTraj Number of trajectories to generate
#' @param numObs Number of repeated measurements per trajectory
#' @param tStart
#' @param tEnd
#' @param trajLabels character vector of the trajectory names
tsdata_create_repmeas = function(
  group='ts', 
  numTraj, 
  numObs, 
  tStart=0, 
  tEnd=numObs-1, 
  trajLabels=paste0(group, '.', 1:numTraj)
) {
  assert_that(
    is.count(numObs),
    is.count(numTraj)
  )
    
  if(length(trajLabels) == 1 && numTraj > 1) {
    trajLabels = paste0(trajLabels, '.', 1:numTraj)
  }

  if(!is.factor(trajLabels)) {
    trajLabels = factor(trajLabels, levels=unique(trajLabels))
  }

  dt_meas = data.table(
      Group=factor(group),
      Id=rep(trajLabels, each=numObs),
      Time=seq(tStart*1.0, tEnd*1.0, length.out=numObs) %>% rep(numTraj),
      Mu=0,
      Sigma=0,
      Value=0,
      key=c('Group', 'Id', 'Time')
  )

  #group trend
  #TODO include average of covariates
  dt_trend = data.table(
    Group=factor(group),
    Time=seq(tStart*1.0, tEnd*1.0, length.out=numObs),
    Mu=0.0,
    Sigma=0.0,
    Value=0.0
  )
  setattr(dt_meas, 'trend', dt_trend)
  return(dt_meas)
}

#' @description Generate the response from the distributional parameters (mu and sigma) at each observation time.
tsdata_response = function(data, family = gamlss.dist::rNO) {
  data[, Value := family(.N, mu = Mu, sigma = Sigma)][]
}

# Mu ####
#' @description Add population-wide effects to the data
#' @param formula the formula to evaluate on the data, e.g. ~ poly(Time, 2, raw=T)
#' @param coefmat the coefficients for each of the formula terms
tsdata_add_fixed = function(data, formula, coefs, what='Mu') {
  numTraj = uniqueN(data$Id)
  assert_that(
    is(formula, 'formula'), 
    !attr(terms(formula), 'response'),
    is.numeric(coefs),
    !is.matrix(coefs)
  )
  
  X = model.matrix(formula, data=data)
  assert_that(
    length(coefs) == ncol(X), 
    msg='number of coefficients does not match the model'
  )
  
  if(has_name(data, what)) {
    data[, c(what) := get(what) + colSums(coefs * t(X))]
  } else {
    data[, c(what) := colSums(coefs * t(X))]
  }
  
  dt_trend = tsdata_trends(data) %>% copy()
  assert_that(!is.null(dt_trend), msg='missing trend data')
  Xtrend = model.matrix(formula, data=dt_trend)
  if(has_name(dt_trend, what)) {
    dt_trend[, c(what) := get(what) + colSums(coefs * t(Xtrend))]
  } else {
    dt_trend[, c(what) := colSums(coefs * t(Xtrend))]
  }
  setattr(data, 'trend', dt_trend)
  
  return(data[])
}


#' @description Add a patient-specific random intercept
#' @param sigma Standard deviation of the normal distribution
tsdata_add_random_intercept = function(data, sigma, what='Mu') {
  assert_that(has_name(data, what))
  numTraj = uniqueN(data$Id)
  tsIntercepts = rnorm(numTraj, mean = 0, sd = sigma)

  data[, c(what) := get(what) + tsIntercepts[rleidv(Id)]]
  return(data[])
}


#' @description Add trajectory-specific effects to the data
#' @param formula the formula to evaluate on the data, e.g. ~ poly(Time, 2, raw=T)
tsdata_add_random = function(data, formula, covMat, what='Mu') {
  assert_that(
    is(formula, 'formula'), 
    !attr(terms(formula), 'response')
  )
  numTraj = uniqueN(data$Id)
  X = model.matrix(formula, data=data)
  tsCoefs = mvnfast::rmvn(numTraj, mu=rep(0, nrow(covMat)), sigma=covMat)
  
  assert_that(
    is.matrix(tsCoefs),
    has_name(data, what)
  )
  assert_that(nrow(tsCoefs) == numTraj, msg='number of coefficient rows does not match the number of trajectories')
  assert_that(ncol(tsCoefs) == ncol(X), msg='number of coefficients does not match the model')
  
  rowcoefmat = tsCoefs[rep(seq_len(numTraj), data[, .N, by=Id]$N),]
  stopifnot(nrow(rowcoefmat) == nrow(data))
  
  data[, c(what) := get(what) + rowSums(rowcoefmat * X)]
  return(data[])
}

# Sigma ####
#' @title Add random (id) variance intercept to the data
#' @details Currently only supports random sigma intercept
#' @param sigmaCovMat The shared variance-covariance matrix for the sigma random effects (1 x 1)
tsdata_add_random_sigma_intercept = function(data, sigma) {
  assert_that(
    is.number(sigma),
    sigma >= 0,
    all(data$Sigma > 0)
  )
  
  numTraj = uniqueN(data$Id)
  tsIntercepts = rnorm(numTraj, mean = 0, sd = sigma)

  newdata = copy(data) %>%  
    .[, SigmaInt := tsIntercepts[rleidv(Id)]] %>%
    .[, Sigma := exp(log(Sigma) + SigmaInt)]
  return(newdata[])
}

#' @title Add group-specific heteroskedasticity to the data
#' @param cv Coefficient of variation
tsdata_add_meanVariance = function(data, cv) {
  assert_that(
    is.numeric(cv),
    all(data$Sigma > 0)
  )
  
 data[, Sigma := exp(log(Sigma) + cv * Mu)]
 return(data[])
}

# Noise ####
tsdata_add_normnoise_cv = function(data, cv=1) {
  stopifnot(cv >= 0)
  data[, Value := Value + rnorm(.N, mean=0, sd=data$Value * cv)]
}

tsdata_add_noise = function(data, rfun=rnorm, ...) {
  data[, Value := Value + rfun(.N, ...)]
}

tsdata_add_lognoise = function(data, meanlog, sdlog, center=FALSE) {
  stopifnot(sdlog >= 0)
  data[, Value := Value + rlnorm(.N, meanlog=meanlog, sdlog=sdlog)]
  if(center) {
      data[, Value := Value - qlnorm(.5, meanlog=meanlog, sdlog=sdlog)]
  }
}

tsdata_add_noise_arima = function(data, ar=NULL, ma=NULL, sd=1, order=NULL) {
  data[, Value := Value + stats::arima.sim(n=.N, list(ar=ar, ma=ma, sd=sd)) %>% as.numeric(), by=Id]
}
