library(latrend)
library(data.table)
library(magrittr)
library(assertthat)
library(rstan)
library(HDInterval)

setClass('lcModelStan', representation=list(permutations='list'), contains='lcModel')

# dimension order: samples, row, columns, chain
as.array.lcModelStan = function(
  x, 
  par, 
  samples = seq(nWarmup(x) + 1, nSamples(x)), 
  permuteK = length(x@permutations) > 0, 
  chains = seq_len(nChains(x)),
  ...
) {
  assert_that(
    is.lcModel(x),
    is.scalar(par),
    is.character(par),
    all(samples %in% seq_len(nSamples(x))),
    all(chains %in% seq_len(nChains(x)))
  )
  assert_that(
    has_name(x, par),
    msg = sprintf('model does not contain samples for parameter "%s"', par)
  )
  stanmod = x@model
  nSamples = nSamples(x)
  nChains = nChains(x)
  pardim = stanmod@par_dims[[par]]
  sampleParNames = names(stanmod@sim$samples[[1]])
  assert_that(all(chains %in% seq_len(nChains)), msg='invalid chain indices in argument "chains"')

  if(length(pardim) == 0) {
    # scalar. cannot permute
    sampleMat = do.call(cbind, lapply(stanmod@sim$samples, '[[', par))
    dimnames(sampleMat) = list(Sample = NULL, Chain = NULL)
    
    assert_that(
      nrow(sampleMat) == nSamples,
      ncol(sampleMat) == nChains
    )
    return (sampleMat[samples, chains, drop = FALSE])
  }
  else if(length(pardim) == 1) {
    if(pardim == 0) {
      # empty par
      return (array(dim=0))
    }
    # vector
    parNames = grep(paste0('^', par, '[\\[\\.]'), sampleParNames, value = TRUE)
    assert_that(length(parNames) > 0, msg = sprintf('no samples found for %s', par))
    assert_that(!anyNA(parNames), msg='sampling param names not found')
    chainMats = lapply(stanmod@sim$samples, '[', parNames) %>%
      lapply(function(parList) {
        parVec = do.call(c, parList)
        assert_that(length(parVec) == nSamples * pardim)
        matrix(parVec, nrow=nSamples, ncol=pardim)})
    
    sampleArray = do.call(c, chainMats) %>%
      array(dim = c(nSamples, pardim, nChains),
        dimnames = list(Sample=NULL, Par=parNames, Chain=NULL)) %>%
      .[, , chains, drop = FALSE]
    
    if (permuteK && pardim[1] == nClusters(x)) { #TODO: proper param checking
      for (ic in seq_along(chains)) {
        icChain = chains[ic]
        sampleArray[,,ic] = sampleArray[, x@permutations[[icChain]], ic, drop=FALSE]
      }
    }
    
    y = sampleArray[samples,,, drop=FALSE]
    attr(y, 'parNames') = parNames
    return (y)
  }
  else if(length(pardim) == 2) {
    # matrix
    parNames = grep(paste0('^', par, '\\['), sampleParNames, value = TRUE)
    if (length(parNames) == 0) {
      parNames = grep(paste0('^', par, '\\.\\d+\\.\\d+'), sampleParNames, value = TRUE)
    }
    assert_that(!anyNA(parNames))

    chainArrays = lapply(stanmod@sim$samples, '[', parNames) %>%
      lapply(function(parList) {
        parVec = do.call(c, parList)
        assert_that(length(parVec) == nSamples * prod(pardim))
        array(parVec, dim = c(Sample = nSamples, Row = pardim[1], Col = pardim[2]))
      })
    
    sampleArray = do.call(c, chainArrays) %>% 
      array(dim = c(nSamples, pardim[1], pardim[2], nChains),
        dimnames = list(Sample=NULL, Row=NULL, Col=NULL, Chain=NULL)) %>%
      .[, , , chains, drop = FALSE]
    
    if (permuteK) {
      if (par %in% c('b_mu', 'beta_mu', 'beta_sigma') && 
          dim(sampleArray)[3] == nClusters(x)) {
        # column-wise permutation
        for (ic in seq_along(chains)) {
          icChain = chains[ic]
          sampleArray[,,,ic] = sampleArray[,,x@permutations[[icChain]],ic, drop=FALSE]
        }
      }
      else if (par %in% c('r0_mu', 'pp', 'sdz_mu') && 
          dim(sampleArray)[2] == nClusters(x)) {
        # row-wise permutation
        for (ic in seq_along(chains)) {
          icChain = chains[ic]
          sampleArray[,,,ic] = sampleArray[,x@permutations[[icChain]],,ic, drop=FALSE]
        }
      }
    }
    
    y = sampleArray[samples,,,, drop=FALSE]
    attr(y, 'parNames') = parNames
    return (y)
  } 
  else {
    stop('unsupported parameter dimensions')
  }
}

as.data.frame.lcModelStan = function(x, ...) {
  as.data.table(x, ...)
}

as.data.table.lcModelStan = function(x, pars = coefNames(x), ...) {
  makeParData = function(par) {
    sampleData = as.array(x, par, ...)
    parNames = attr(sampleData, 'parNames')
    assert_that(length(parNames) > 0)
    ndim = length(dim(sampleData))
    
    if(ndim == 2) {
      tt = as.table(sampleData)
      dtSample = as.data.table(tt)
      dtSample[, Sample := factor(Sample, levels = rownames(tt)) %>% as.integer()]
      dtSample[, Chain := factor(Chain, levels = colnames(tt)) %>% as.integer()]
      dtSample[, Par := par]
      setnames(dtSample, 'N', 'Value')
    }
    else if (ndim == 3) {
      dtSample = as.data.table(sampleData, value.name = 'Value')    
    }
    else if (ndim == 4) {
      dtSample = as.data.table(sampleData, value.name = 'Value')
      parIdx = matrix(seq_along(parNames), dim(sampleData)[c(-1, -ndim)])
      dtSample[, Par := parNames[parIdx[cbind(Row, Col)]]]
      dtSample[, Row := NULL]
      dtSample[, Col := NULL]
    } else {
      stop('unsupported dims')
    }
    setcolorder(dtSample, c('Par', 'Chain', 'Sample', 'Value'))
    setkey(dtSample, Par, Chain, Sample)
    
    dtSample[]
  }
 
  lapply(pars, makeParData) %>% rbindlist() 
}


#' @inheritDotParams as.array.lcModelStan
aggregate.lcModelStan = function(x, par, fun, ...) {
  assert_that(
    is.lcModel(x),
    is.character(par),
    par %in% names(x)
  )
  stanmod = x@model
  nChains = nChains(x)
  pardim = stanmod@par_dims[[par]]
  
  samples = as.array(x, par = par, ...)
  if(length(pardim) == 0) {
    # scalar
    outSc = fun(samples)
    names(outSc) = par
    return(outSc)
  }
  else if(length(pardim) == 1) {
    if(pardim == 0) {
      return(numeric())
    }
    # vector
    outVec = apply(samples, 2, fun)
    names(outVec) = paste0(par, '[', 1:pardim, ']')
    return(outVec)
  }
  else {
    # matrix
    outMat = apply(samples, c(2, 3), fun)
    attr(outMat, 'parNames') = attr(samples, 'parNames')
    dimnames(outMat) = NULL
    return(outMat)
  }
}

#' @inheritDotParams aggregate.lcModelStan
mean.lcModelStan = function(x, par, ...) {
  aggregate(x, par = par, fun = mean, ...)
}

#' @inheritDotParams aggregate.lcModelStan
sd.lcModelStan = function(x, par, ...) {
  aggregate(x, par = par, fun = sd, ...)
}

fitted.lcModelStan = function(object, what = 'mu', ..., clusters = trajectoryAssignments(object)) {
  what = match.arg(what, c('mu', 'sigma'))
  
  par = paste0(what, 's')
  
  if (hasName(object, par)) {
    pred = mean(object, par)
  } else {
    warning('Cannot compute fitted values: missing model parameter ', par)
    return (transformFitted(pred = NULL, model = object, clusters= clusters))
  }
  
  colnames(pred) = clusterNames(object)
  rownames(pred) = NULL
  
  transformFitted(pred = pred, model = object, clusters= clusters)
}

#. clusterProportions ####
setMethod('clusterProportions', signature('lcModelStan'), function(object, ...) {
  props = mean(object, par='theta', ...)
  names(props) = clusterNames(object)
  return(props)
})

#. postprob ####
setMethod('postprob', signature('lcModelStan'), function(object, ...) {
  if(hasName(object@model@.MISC, 'pp')) {
    pp = object@model@.MISC$pp
  } else if(not('pp' %in% names(object))) {
      warning('pp not available in model')
      return(callNextMethod())
  } else {
    pp = mean(object, 'pp', ...) %>% t()
  }
  
  return(pp)
})


validate.lcModelStan = function(model, ...) {
  tryCatch(is_valid_postprob(postprob(model), model), error = function(e) FALSE)
}

plotChains = function(x, ..., chains = 1:nChains(x)) {
  library(gridExtra)
  assert_that(is.lcModel(x))
  pList = lapply(chains, function(i) plot(x, ..., chains = i) + labs(title=NULL, subtitle = paste0('Chain ', i)))
  
  do.call(grid.arrange, c(pList, 
    top = list(textGrob('Cluster trajectories'))))
}

#. converged ####
setMethod('converged', 'lcModelStan', function(object) {
  if (length(object@model@sim) == 0) {
    return (FALSE)
  }
  
  all(parameterConvergence(object))
})

#' @description Determine for each relevant model parameter whether it has converged (as determined by Rhat)
parameterConvergence = function(
  object, 
  pars, 
  maxRhat = 1.1,
  minClusSize = 10
) {
  allPars = names(object)
  
  if (missing(pars)) {
    # determine which parameters should be considered in the model convergence
    # parameters of empty clusters are excluded
    pars = c('theta', 'alpha','b_mu', 'intercept_sigma', 'cv', 'sdz0_mu', 'sdz_mu')
    if (not('intercept_sigma' %in% allPars)) {
      pars = union(pars, 'sigma')
    }
    if (not('alpha' %in% allPars)) {
      pars = union(pars, c('intercept_mu', 'beta_mu'))
    }
  }

  pars = intersect(pars, allPars)
  
  # compute Rhat per parameter
  parConvs = lapply(pars, function(p) Rhat(object, p)) %>%
    unlist() %>%
    `<`(maxRhat) 
  flatPars = names(parConvs)
  
  # drop cluster-specific parameters based on cluster criteria
  invalidClusters = which(clusterSizes(object) < minClusSize)
  
  if (length(invalidClusters) > 0) {
    dropPars = sapply(invalidClusters, function(idx) sprintf(
      c('alpha\\[%d\\]', 'b_mu\\[\\d,%d\\]', 'beta_mu\\[\\d,%d\\]'), idx))
    matchedDropPars = grep(paste0(dropPars, collapse = '|'), flatPars, value = TRUE)
    message('NOTE: Dropping cluster parameters ', paste0(matchedDropPars, collapse = ', '), ' from parameter convergence assessment')
    parConvs = parConvs[setdiff(flatPars, matchedDropPars)]
  }
  
  return(parConvs)
}

names.lcModelStan = function(x) {
  allPars = x@model@sim$pars_oi
  
  # filter out pars to were removed from sim$samples
  flatPars = flatnames(x)
  validMask = vapply(
    allPars, 
    function(par) {
      any(grepl(paste0('^', par, '\\['), flatPars)) || par %in% flatPars
    },
    FUN.VALUE = TRUE
  )
  allPars[validMask]
}

flatnames = function(x) {
  UseMethod('flatnames', x)
}

flatnames.lcModelStan = function(x) {
  names(x@model@sim$samples[[1]])
}

hdi.lcModelStan = function(object, par, credMass = .95, ...) {
  if (!all(hasName(object, par))) {
    warning('model does not have parameter ', par)
    return(numeric())
  }
  
  if (length(par) > 1) {
    out = lapply(par, function(p) HDInterval::hdi(object, p, credMass = credMass, ...))
    return(do.call(cbind, out))
  }
  
  sampleArray = as.array(object, par = par)
  dims = which(not(names(dimnames(sampleArray)) %in% c('Sample', 'Chain')))
  
  if (length(dims) == 0) {
    # fix for models which did not define a param as a vector
    boundsArray = HDInterval::hdi(as.numeric(sampleArray), credMass = credMass, ...)
    parNames = paste0(par, '[1]')
  }
  else { 
    boundsArray = apply(sampleArray, dims, function(a) HDInterval::hdi(as.numeric(a), credMass = credMass, ...))
    parNames = attr(sampleArray, 'parNames')
  }
  
  boundsMat = matrix(boundsArray, nrow = 2)
  rownames(boundsMat) = c('Lower', 'Upper')
  colnames(boundsMat) = parNames
  
  return(boundsMat)
}


coef.lcModelStan = function(object, fun=mean, ...) {
  pars = coefNames(object, ...)
  
  coefs = lapply(pars, function(p) {
    out = fun(object, p, ...)
    if(is.matrix(out)) {
      idxStrings = expand.grid(1:nrow(out), 1:ncol(out)) %>%
        as.matrix() %>%
        apply(1, paste0, collapse=',')
      names(out) = paste0(p, '[', idxStrings, ']')
    }
    return(out)
  })
  
  unlist(coefs)
}

coefNames = function(object, maxLen = 100, ...) {
  assert_that(is.lcModel(object))
  
  parNames = names(object)
  
  suppressWarnings({ 
    parLens = sapply(object@model@par_dims, prod)[parNames] 
  })
  
  setdiff(parNames, union(parNames[parLens >= maxLen], c('pp', 'log_lik', 'log_lik_n', 'lp__', 'cfc', 'total_log_lik', 'total_log_lik_n')))
}

nChains = function(object) {
  assert_that(is.lcModel(object))
  object@model@sim$chains
}

#' @title Total number of samples, including warmup
nSamples = function(object) {
  assert_that(is.lcModel(object))
  allLengths = sapply(object@model@sim$samples, lengths)
  assert_that(all(allLengths == allLengths[1]), msg='Stan model parameters have mismatching number of samples')
  allLengths[1]
}

nWarmup = function(object) {
  assert_that(is.lcModel(object))
  # check if warmup samples are present
  nIter = object@model@sim$iter
  nWarm = object@model@sim$warmup
  nSamp = object@model@sim$samples[[1]][[1]] %>% as.matrix() %>% nrow()
  
  if (nSamp == nIter) {
    nWarm
  } else {
    assert_that(nIter - nWarm == nSamp, msg = 'invalid model state; some samples are missing')
    0L
  }
}

Rhat = function(object, par, ...) {
  if (length(par) > 1) {
    return(sapply(par, function(p) Rhat(object, p, ...)) %>% 
        unname() %>% 
        unlist())
  }
  
  if (!hasName(object, par)) {
    warning('model does not have parameter ', par)
    return(numeric())
  }

  sampleArray = as.array(object, par = par, ...)
  parNames = attr(sampleArray, 'parNames')
  
  if (par == 'theta' && nClusters(object) == 1) {
    rhats = 1
  } else if (is.null(parNames)) {
    # special case
    rhats = rstan::Rhat(sampleArray) %>% as.numeric()
    parNames = paste0(par, '[1]')
  } else {
    chainDim = length(dim(sampleArray))
    parDims = seq_along(dim(sampleArray))[c(-1, -chainDim)]
    
    # compute Rhat for each parameter
    rhats = apply(sampleArray, parDims, rstan::Rhat) %>% as.numeric()
  }
  
  names(rhats) = parNames
  if (anyNA(rhats)) {
    stop('Rhat is NA for par', par)
  }
  return(rhats)
}

permute = function(object, ...) {
  UseMethod('permute')
}

permute.lcModelStan = function(object, chain, permutations) {
  assert_that(
    is.lcModel(object),
    is.integer(permutations),
    is.count(chain),
    chain <= nChains(object),
    length(permutations) == nClusters(object)
  )
  
  if (length(object@permutations) == 0) {
    object@permutations = replicate(nChains(object), seq_len(nClusters(object)), simplify = FALSE)
  }
  
  object@permutations[[chain]] = object@permutations[[chain]][permutations]
  
  return(object)
}

autopermute = function(object, method = 'ECR-ITERATIVE-1') {
  library(label.switching)
  message('Auto-permuting Stan model using ', method, '...')
  if (nClusters(object) == 1) {
    message('Only 1 cluster. Nothing to permute')
    return(object)
  }
  
  data = model.data(object) %>% as.data.table()
  assert_that(!is.null(object@model@.MISC$classifications), msg='cannot autopermute stripped model: @model@.MISC$classifications is missing')
  
  if (hasName(data, 'Group')) {
    refClusters = data[, first(Group), keyby=Id]$V1 %>% as.integer()
    out1 = label.switching::label.switching(
      method = method, 
      z = object@model@.MISC$classifications[,,1], 
      K = nClusters(object),
      groundTruth = refClusters)
    perms = out1$permutations[[method]] %>% 
      matrixStats::colMedians() %>% 
      as.integer()
    object = permute(object, chain = 1, permutations = perms)
  } else {
    ref = label.switching::label.switching(
      method = method, 
      z = object@model@.MISC$classifications[,,1], 
      K = nClusters(object)
    )
    refClusters = as.integer(ref$clusters[method, ])  
  }

  if (nChains(object) > 1) {
    for (ic in seq.int(2, nChains(object), by=1)) {
      kOut = label.switching::label.switching(method = 'ECR-ITERATIVE-1', 
        z = object@model@.MISC$classifications[,,ic], 
        K = nClusters(object), 
        groundTruth = refClusters)
      perms = kOut$permutations[[method]] %>% 
        matrixStats::colMedians() %>% 
        as.integer()
      object = permute(object, chain = ic, permutations = perms)
    }
  }
  return(object)
}

#. predictForCluster ####
setMethod('predictForCluster', 'lcModelStan', function(object, newdata, cluster, what, transform = force, ...) {
  assert_that(
    !is.null(newdata),
    cluster %in% clusterNames(object),
    !has_name(newdata, 'Cluster')
  )
  
  clusIdx = match(cluster, clusterNames(object))
  
  f.par = formula(object, what=what) %>% 
    update(NULL ~ .)
  fMat = model.matrix(update(f.par, ~ . - Fit), data = newdata)
  
  alpha = numeric(nrow(fMat))
  betaNames = colnames(fMat)
  intPar = paste0('intercept_', what)
  betaPar = paste0('beta_', what)

  # special parameters
  if (what == 'sigma') {
    if (hasName(object, 'cv')) {
      assert_that(!has_name(newdata, 'Fit'))
      betaNames = setdiff(betaNames, 'Fit')
      mu = predictForCluster(object, newdata=newdata, cluster=cluster, what='mu', ...)
      assert_that(length(mu) == nrow(fMat))
      cv = mean(object, par='cv') %>%
        matrix(nrow=nClusters(object)) %>%
        set_colnames('cv')
      alpha = alpha + cv[clusIdx,] * mu
    } 
    else if (hasName(object, 'sigma')) {
      sigmas = mean(object, 'sigma', ...) %>% 
        log() %>%
        matrix(nrow = nClusters(object))
      alpha = alpha + sigmas[clusIdx, ]
    }
  }
  
  # intercept
  if (hasName(object, intPar)) {
    betaNames = setdiff(betaNames, '(Intercept)')
    intercept = mean(object, par=intPar, ...) %>% 
      matrix(nrow=nClusters(object)) %>%
      set_colnames('Intercept')
    alpha = alpha + intercept[clusIdx,]
  }
  
  # extract beta
  if (hasName(object, betaPar)) {
    betaMat = mean(object, par=betaPar, ...) %>%
      t() %>%
      set_rownames(clusterNames(object)) %>%
      set_colnames(betaNames)
    xMat = fMat[, betaNames]
    y = betaMat[clusIdx, , drop=FALSE] %*% t(xMat) + alpha
  }
  else {
    y = alpha
  }

  assert_that(length(y) == nrow(newdata))
  return(transform(y) %>% as.numeric())
})

waic.lcModelStan = function(object, ...) {
  mll = compute_marginal_loglik(object, ...)
  waic(mll)
}

loo.lcModelStan = function(object, par = 'log_lik') {
  if (par == 'log_lik' && hasName(object@model@.MISC, 'loo')) {
    object@model@.MISC$loo
  }
  else if (par == 'log_lik_n' && hasName(object@model@.MISC, 'loo_n')) {
    object@model@.MISC$loo_n
  }
  else {
    ll = loo::extract_log_lik(object@model, par, merge_chains = FALSE)
    loo::loo(ll, r_eff = loo::relative_eff(exp(ll)))
  }
}

defineInternalMetric('cLOO', function(m) {
  loo(m)$estimates['looic', 'Estimate']
})

psis.lcModelStan = function(object, par = 'log_lik') {
  if (par == 'log_lik' && hasName(object@model@.MISC, 'psis')) {
    object@model@.MISC$psis
  }
  else if (par == 'log_lik_n' && hasName(object@model@.MISC, 'psis_n')) {
    object@model@.MISC$psis_n
  }
  else {
    ll = loo::extract_log_lik(object@model, par, merge_chains = FALSE)
    loo::psis(-ll, r_eff = loo::relative_eff(exp(ll)))
  }
}

cWAIC = function(object, par, ...) {
  UseMethod('cWAIC')
}

cWAIC.lcModelStan = function(object, par = 'log_lik') {
  if (par == 'log_lik' && hasName(object@model@.MISC, 'waic')) {
    object@model@.MISC$waic
  }
  else {
    assert_that(has_name(object, par))
    
    ll = loo::extract_log_lik(object@model, par)
    
    if (nrow(ll) > nSamples(object) - nWarmup(object)) {
      warning('cWAIC() using warmup samples')
    }
    
    loo::waic(ll)
  }
}

defineInternalMetric('cWAIC', function(m) {
  cWAIC(m)$estimates['waic', 'Estimate']
})

waic2 = function(object, par = 'log_lik') {
  ll = as.array(object, par)
  training_error = -mean(log(colMeans(exp(ll))))
  fun_var_div = mean(colMeans(ll^2) - colMeans(ll)^2)
  c(waic = training_error + fun_var_div, p_waic = fun_var_div)
}

defineInternalMetric('cWAIC2', function(m) {
  waic2(m)['waic']
})

waic3 = function(object, par = 'log_lik') {
  ll = as.array(object, par)
  lppd = sum(log(mean(exp(ll))))
  p_waic = sum(var(ll))
  c(waic = -2 * lppd + 2 * p_waic, p_waic = p_waic)
}

defineInternalMetric('cWAIC3', function(m) {
  waic3(m)['waic']
})

mLogLik = function(object, ...) {
  UseMethod('mLogLik')
}

mLogLik.NULL = function(object, ...) {
  as.numeric(NA)
}

mLogLik.matrix = function(object, ..., df = NA) {
  N = ncol(object)
  S = nrow(object)
  mll = sum(matrixStats::colLogSumExps(object) - log(S))
  
  attr(mll, 'nobs') = N
  attr(mll, 'df') = df
  class(mll) = 'logLik'
  
  mll
}

#' @description Get the expected marginal likelihood 
mLogLik.lcModelStan = function(object, ...) {
  assert_that(inherits(object, 'lcModelStan'))
  df = length(coef(object)) - 1 #-1 for theta
  mll = compute_marginal_loglik(object, ...)
  mLogLik(mll, df = df)
}

defineInternalMetric('mLogLik', mLogLik)

logLik.lcModelStan = mLogLik.lcModelStan

#' @param M The number of MC samples for approximating the random effects integral
#' @param compute Whether to compute the marginal loglik, or only return cached results (if available)
#' @param recompute Whether to recompute the result when already cached. Overwrites the cache.
#' @param mixType The type of mixture model to assume. Default assumes individuals belong to each class.
#' 'fixed.alloc' allocates trajectories to their most likely class on average.
#' @param cache Whether to save the computed mll matrix inside the model
#' @return A `matrix` with the marginal loglik estimated per sample per individual
compute_marginal_loglik = function(
  object, 
  M = 5e3, 
  recompute = FALSE, 
  mixType = 'smax.alloc', 
  cache = TRUE,
  seed = 1
) {
  assert_that(
    is.count(M),
    is.flag(recompute),
    is.flag(cache)
  )
  
  mixType = match.arg(mixType[1], c('mix', 'max.alloc', 'smax.alloc', 'prop.alloc'))
  
  modelSuffix = switch(
    mixType,
    mix = '',
    max.alloc = '_z',
    smax.alloc = '_zs',
    prop.alloc = '_m'
  )
  
  mllKey = paste0('mll', modelSuffix, '_M', M)

  if (!recompute && hasName(object@model@.MISC, mllKey)) {
    message(
      sprintf(
        'Returning previous result for marginal likelihood (M=%d).', 
        object@model@.MISC$marginalLogLikIters
      )
    )
    return (object@model@.MISC[[mllKey]])
  }

  pars = coefNames(object)
  
  modelName = tools::file_path_sans_ext(object@model@model_name)
  
  loglikFile = file.path('loglik', paste0(modelName, '_mloglik', modelSuffix, '.stan'))
  message('Checking presence of ', loglikFile)
  assert_that(file.exists(loglikFile), msg = sprintf('mloglik Stan file not found: computation is not implemented for model %s', modelName))
  
  message('Compiling marginal likelihood computation model for mix type ', mixType, '...')
  cmod = stan_model(
    file = loglikFile, 
    allow_undefined = TRUE,
    auto_write = TRUE,
    includes = sprintf('\n#include "%s"\n', file.path(getwd(), 'loglik', 'get_iter.hpp'))
  )
  
  message('Generating standata...')
  S = nrow(as.array(object, pars[1]))
  standata = make.standata(getLcMethod(object), data = model.data(object))
  standata$S = S
  standata$M = M
  nSdz = 1L + standata$B_mu
  standata$z = rnorm(M * nSdz, mean = 0, sd = 1) %>% matrix(nrow = M)
  standata$g_i = as.integer(trajectoryAssignments(object))
  
  if (modelSuffix == '_zs') {
    message('Computing sample-specific class allocations...')
    standata$g_is = apply(as.array(object, 'pp'), c(1, 3), which.max) %>% t()
    assert_that(
      ncol(standata$g_is) == S,
      nrow(standata$g_is) == standata$I
    ) 
  }

  for (par in pars) {
    message(sprintf('\tAdding samples for parameter "%s"', par))
    a = as.array(object, par)
    chainDim = match('Chain', names(dimnames(a)))
    parSamples = abind::adrop(a, drop = chainDim)

    if (length(object@model@par_dims[[par]]) > 0) {
      standata[[par]] = parSamples
    } else {
      standata[[par]] = as.vector(parSamples)
    }
  }

  message(sprintf('Estimating marginal log-likelihood for %s samples using %s iterations...', scales::comma(standata$S), scales::comma(M)))
  stanmodel <<- sampling(
    cmod, 
    data = standata, 
    warmup = 0, 
    iter = standata$S - 1, 
    algorithm = 'Fixed_param', 
    chains = 1
  )
  
  assert_that(length(stanmodel@sim) > 0, msg = 'Stan model failed to complete')
  
  ll = as.array(stanmodel, 'log_lik')[, 1, ] # select chain
  colnames(ll) = ids(object)
  
  if (cache) {
    object@model@.MISC[[mllKey]] = ll
    message('Saved computed sample marginal likelihoods inside the model.')
  }
  
  ll
}

lppd = function(object, ...) {
  UseMethod('lppd')
}

lppd.NULL = function(object, ...) {
  as.numeric(NA)
}

lppd.matrix = function(object, ...) {
  S = nrow(object)
  sum(matrixStats::colLogSumExps(object) - log(S))
}

lppd.lcModelStan = function(object, ...) {
  mll = compute_marginal_loglik(object, ...)
  lppd(mll)
}

defineInternalMetric('lppd', lppd)


mWAIC = function(object, ...) {
  UseMethod('mWAIC')
}

mWAIC.NULL = function(object, ...) {
  as.numeric(NA)
}

mWAIC.matrix = function(object, ...) {
  elppd = lppd(object) - p_waic1(object)
  -2 * elppd
}

#' @inheritDotParams compute_marginal_loglik
mWAIC.lcModelStan = function(object, recompute = FALSE, ...) {
  if (!recompute && hasName(object@model@.MISC, 'mWAIC')) {
    object@model@.MISC$mWAIC
  } else {
    mll = compute_marginal_loglik(object, recompute = recompute, ...)
    mWAIC(mll, ...)  
  }
}

defineInternalMetric('mWAIC', mWAIC)

p_waic1 = function(object, ...) {
  UseMethod('p_waic1')
}

p_waic2 = function(object, ...) {
  UseMethod('p_waic2')
}

p_waic1.NULL = function(object, ...) {
  as.numeric(NA)
}

p_waic2.NULL = function(object, ...) {
  as.numeric(NA)
}

p_waic1.lcModelStan = function(object, ...) {
  mll = compute_marginal_loglik(object, ...)
  p_waic1(mll)
}

p_waic2.lcModelStan = function(object, ...) {
  mll = compute_marginal_loglik(object, ...)
  p_waic2(mll)
}

# see Gelman 2013
p_waic1.matrix = function(object, ...) {
  S = nrow(object)
  lppd_i = matrixStats::colLogSumExps(object) - log(S)
  mean_ll_i = colMeans(object)
  
  2 * sum(lppd_i - mean_ll_i)
}

# see Gelman 2013
p_waic2.matrix = function(object, ...) {
  sum(colVars(object))
}

defineInternalMetric('p_waic1', p_waic1)

defineInternalMetric('p_waic2', p_waic2)

#' @description Returns the highest correlation between cluster intercept samples
alphaCor = function(object, par, minClusSize = 10) {
  if (missing(par)) {
    par = intersect(c('alpha', 'intercept_mu', 'intercept'), names(object))[1]
  }
  
  a = as.array(object, par)[,,]
  
  if (nClusters(object) == 1) {
    return(0)
  }
  
  validClusterMask = clusterSizes(object) >= minClusSize
  a = a[, validClusterMask, drop = FALSE]
  
  if (ncol(a) == 1) {
    0
  } else if (ncol(a) == 2) {
    assert_that(is.matrix(a))
    cor(a[, 1], a[, 2])
  } else {
    assert_that(is.matrix(a))
    max(
      cor(a[, 1], a[, 2]),
      cor(a[, 1], a[, 3]),
      cor(a[, 2], a[, 3])
    )
  }
}

#. strip ####
setMethod('strip', signature('lcModelStan'), function(object) {
  newObject = callNextMethod()
  
  if (inherits(getLcMethod(object), 'lcFitConverged')) {
    newObject@method = strip(getLcMethod(object@method))
  }
  
  if(length(object@model@sim) == 0) {
    message('Model was not sampled. Nothing to strip')
    return(newObject)
  }
  
  stanmod = newObject@model
  nSamples = stanmod@sim$iter - stanmod@sim$warmup
  sampleIdx = seq(stanmod@sim$warmup + 1, stanmod@sim$iter)
  nChains = stanmod@sim$chains
  
  stanmod@stan_args[[1]]$warmup = 0
  stanmod@sim$permutation = NULL
  stanmod@sim$fnames_oi = NULL
  stanmod@.MISC$classifications = NULL
  stanmod@.MISC$stan_fit_instance = NULL
  slot(stanmod, 'stanmodel', check = FALSE) <- NULL
  
  excludePars = c(
    grep('^mus\\[', names(stanmod@sim$samples[[1]]), value=TRUE),
    grep('^sigmas\\[', names(stanmod@sim$samples[[1]]), value=TRUE),
    grep('^r0_mu\\[', names(stanmod@sim$samples[[1]]), value=TRUE),
    grep('^r0_sigma\\[', names(stanmod@sim$samples[[1]]), value=TRUE),
    grep('^r_mu\\[', names(stanmod@sim$samples[[1]]), value=TRUE),
    'lp__'
  ) %>% unique()
  
  # drop pp samples
  if('pp' %in% stanmod@sim$pars_oi && !hasName(stanmod@.MISC, 'pp')) {
    ppdim = stanmod@par_dims$pp
    assert_that(
      nClusters(newObject) == ppdim[1],
      nIds(newObject) == ppdim[2]
    )

    stanmod@.MISC$pp = mean(newObject, 'pp') %>% t()
    dimnames(stanmod@.MISC$pp) = NULL
    
    ppNames = grep('^pp\\[', names(stanmod@sim$samples[[1]]), value=TRUE)
    excludePars = union(excludePars, ppNames)
  }
  
  # drop log_lik samples
  if('log_lik' %in% stanmod@sim$pars_oi && !hasName(stanmod@.MISC, 'log_lik')) {
    stanmod@.MISC$log_lik = mean(newObject, 'log_lik') %>% sum()
    stanmod@.MISC$loo = loo(newObject, par = 'log_lik')
    stanmod@.MISC$waic = cWAIC(newObject)

    llNames = grep('^log_lik\\[', names(stanmod@sim$samples[[1]]), value=TRUE)
    excludePars = union(excludePars, llNames)
  }
  
  # drop log_lik_n samples
  if('log_lik_n' %in% stanmod@sim$pars_oi) {
    # not used anymore
    llnNames = grep('^log_lik_n\\[', names(stanmod@sim$samples[[1]]), value=TRUE)
    excludePars = union(excludePars, llnNames)
  }
  
  # drop log_lik_z samples
  if('log_lik_z' %in% stanmod@sim$pars_oi) {
    # not used anymore
    llnNames = grep('^log_lik_z\\[', names(stanmod@sim$samples[[1]]), value=TRUE)
    excludePars = union(excludePars, llnNames)
  }
  
  # drop parameters
  for (i in seq_len(nChains)) {
    stanmod@sim$samples[[i]][excludePars] = NULL
  }
  
  # drop warmup samples
  samplerParams = get_sampler_params(stanmod)
  
  for (i in seq_len(nChains)) {
    stanmod@sim$samples[[i]] = lapply(stanmod@sim$samples[[i]], '[', sampleIdx)
    # restore sampler param attributes..
    attr(stanmod@sim$samples[[i]], 'sampler_params') = list(samplerParams[[i]])
  }
  
  # drop inits
  stanmod@inits = list()
  
  newObject@model = stanmod
  newObject@tag = NA
  
  binModel = serialize(newObject, connection = NULL)
  
  # sanity check to ensure that stripped model is actually fully stripped
  assert_that(
    length(binModel) < 15 * 1024^2, 
    msg = 'stripped model remained large (> 15MB)'
  )
  
  return(newObject)
})

# to ensure compatbility with older models that dont have a @times slot
time.lcModelStan = function(x, ...) {
  if (.hasSlot(x, 'times') && length(x@times) > 0) {
    x@times
  } else {
    model.data(x)[[timeVariable(x)]] %>% unique()
  }
}

#. traceplot ####
setMethod('traceplot', signature('lcModelStan'), function(object, pars = coefNames(object), ...) {
  library(ggplot2)
  library(scales)
  pardata = as.data.table(object, pars = pars, ...) %>%
    .[, Chain := factor(Chain)]
  
  ggplot(pardata, aes(x = Sample, y = Value, color = Chain)) +
    scale_x_continuous(breaks = pretty_breaks(5), labels=comma) +
    scale_color_viridis_d() +
    geom_line(alpha = if(uniqueN(pardata$Chain) == 1) 1 else .5) +
    facet_wrap(~ Par, scales = 'free_y') +
    theme(axis.text.x = element_text(angle = 45, hjust=1)) +
    labs(title = 'Trace plot')
})

pairs.lcModelStan = function(x, pars, ...) {
  pardata = as.data.table(x, par = pars, ...)
  dtpairs = tidybayes::gather_pairs(pardata, key = 'Par', value = 'Value', triangle = 'lower only') %>%
    as.data.table() %>%
    .[, Chain := factor(Chain)]
    
  ggplot(dtpairs, aes(x = .x, y = .y)) +
    geom_hex() +
    scale_fill_viridis_c() +
    facet_grid(.row ~ .col, scales = 'free') +
    labs(x = 'Draw', y = 'Draw')
}
