library(latrend)
setClass('lcMethodStan', contains = 'lcMethod')

#' @param algorithm algorithm supported by sampling(), vb(), or optimizing()
lcMethodStan = function(
  model,
  formula = Value ~ 0 + Time + I(Time^2), # without intercept
  random = ~0,
  time = getOption('latrend.time'),
  id = getOption('latrend.id'),
  nClusters = 2,
  use_priors = TRUE,
  use_sigmas = FALSE,
  use_cor = TRUE,
  init = 0,
  chains = 1, 
  samples = 1e3,
  warmup = 1e3,
  thin = 1,
  algorithm = 'NUTS',
  adapt_delta = .6,
  max_treedepth = 10,
  exclude = character(0),
  seed = sample.int(.Machine$integer.max, 1),
  ...
) {
  mc = match.call.all()
  mc$Class = 'lcMethodStan'
  do.call(new, as.list(mc))
}

lcMethodStanGBTM = function(model = 'gbtm', use_sigmas = FALSE, ...) {
  mc = match.call.all()
  do.call(lcMethodStan, as.list(mc[-1]))
}

lcMethodStanGMM = function(
  random = ~1, 
  model = 'gmm_full_diag', 
  use_sigmas = FALSE, 
  use_sdzs = FALSE, 
  exclude = c('z_mu', 'r_mu', 'r0_mu', 'cfc'), 
  ...
) {
  mc = match.call.all()
  do.call(lcMethodStan, as.list(mc[-1]))
}

lcMethodStanGMMRV = function(
  random = ~ Time + I(Time^2), 
  model = 'gmm-rv_full_diag', 
  use_sdzs = FALSE,
  prior_sdzSigma_nu = 3,
  prior_sdzSigma_sigma = .1,
  exclude = c('r0_mu', 'r_mu', 'z_mu', 'cfc', 'r0_sigma'),
  ...
) {
  mc = match.call.all()
  do.call(lcMethodStan, as.list(mc[-1]))
}

lcMethodStanGMMMeanVar = function(
  random = ~ Time + I(Time^2), 
  model = 'gmm-mv_full_diag', 
  use_sdzs = FALSE,
  prior_cv_nu = 3,
  prior_cv_sigma = .1,
  exclude = c('r0_mu', 'r_mu', 'z_mu', 'cfc', 'mus', 'sigmas'),
  ...
) {
  mc = match.call.all()
  do.call(lcMethodStan, as.list(mc[-1]))
}


lcMethodStanGMMRMV = function(
  random = ~ Time + I(Time^2), 
  model = 'gmm-rmv_full_diag', 
  use_sdzs = FALSE,
  prior_sdzSigma_nu = 3,
  prior_sdzSigma_sigma = .1,
  prior_cv_nu = 3,
  prior_cv_sigma = .1,
  exclude = c('r0_mu', 'r_mu', 'r0_sigma', 'z_mu', 'cfc', 'mus', 'sigmas'),
  ...
) {
  mc = match.call.all()
  do.call(lcMethodStan, as.list(mc[-1]))
}

setValidity('lcMethod', function(object) { TRUE }) # disable to allow for character formulas
setValidity('lcMethodStan', function(object) { TRUE })

#. getName ####
setMethod('getName', signature('lcMethodStan'), function(object) 'Bayesian mixture model using Stan')

#. getShortName ####
setMethod('getShortName', signature('lcMethodStan'), function(object) 'stanmix')

#. getLabel ####
setMethod('getLabel', signature('lcMethodStan'), function(object) '')

# setMethod('responseVariable', signature('lcMethodStan'), function(object) latrend:::getResponse(object$formula))

#. preFit ####
setMethod('preFit', signature('lcMethodStan'), function(method, data, envir, verbose) {
  outEnvir = new.env()
  outEnvir$stanfile = file.path('models', paste0(method$model, '.stan'))
  outEnvir$standata = make.standata(method, data)
  outEnvir$control = list(
    adapt_delta = method$adapt_delta, 
    max_treedepth = method$max_treedepth
  )

  message('Compiling Stan model ', outEnvir$stanfile, '...')
  outEnvir$stanmodel = stan_model(
    file = outEnvir$stanfile, 
    model_name = basename(outEnvir$stanfile), 
    auto_write = TRUE
  )
  assign(basename(outEnvir$stanfile), value = outEnvir$stanmodel, envir = .GlobalEnv)
  
  
  
  if(length(outEnvir$standata$prior_betaSigma_nu) == 1 && G > 1) {
    outEnvir$standata$prior_betaSigma_nu = rep(outEnvir$standata$prior_betaSigma_nu, G)
  }
  
  if(length(outEnvir$standata$prior_betaSigma_sigma) == 1 && G > 1) {
    outEnvir$standata$prior_betaSigma_sigma = rep(outEnvir$standata$prior_betaSigma_sigma, G)
  }
  
  priordata = outEnvir$standata[grep('^prior_', names(outEnvir$standata), value = TRUE)]
  if(not(method$algorithm %in% c('LBFGS', 'BFGS', 'Newton'))) {
    message(paste0(capture.output(print(priordata)), collapse = '\n'))
  }
  
  times = data[[method$time]] %>% unique() %>% sort()
  message('Groups (I): ', comma(outEnvir$standata$I))
  message('Time points: ', length(times))
  message('Observations (N): ', comma(outEnvir$standata$N))
  
  outEnvir
})


#. fit ####
setMethod('fit', signature('lcMethodStan'), function(method, data, envir, verbose, ...) {
  if (method$algorithm %in% c('LBFGS', 'BFGS', 'Newton')) {
    message('Using EM optimization...')
    args = as.list(method, args = c('seed', 'init', 'algorithm', 'hessian', 'draws', 'constrained', 'importance_resampling'))
    args$as_vector = FALSE
    message('Note: priors are disabled')
    envir$standata$use_priors = FALSE
    
    fit = do.call(rstan::optimizing, c(envir$stanmodel, args, data=list(envir$standata)))
    new(
      'lcModelStanEM', 
      method = method,
      data = data,
      model = fit,
      clusterNames = make.clusterNames(method$nClusters),
      tag = list(model=envir$stanmodel)
    )
  } else {
    args = as.list(method, args = c(
      'chains', 'warmup', 'thin', 'seed', 'init', 'verbose', 'algorithm', 'control', 'cores', 'diagnostic_file'
    ))
    args$iter = method$warmup + method$samples
    args$control = envir$control
    args$par = method$exclude
    args$include = FALSE
    
    message('Excluding parameters ', paste0(method$exclude, collapse = ', '), '.')
    
    if (method$algorithm %in% c('meanfield', 'fullrank')) {
      message('Using variational Bayes...')
      vbArgs = args[setdiff(names(args), c('chains', 'warmup', 'thin', 'control'))]
      fit = do.call(rstan::vb, c(envir$stanmodel, vbArgs, data = list(envir$standata)))
    }
    else {
      message('Sampling...')
      fit = do.call(rstan::sampling, c(envir$stanmodel, args, data = list(envir$standata)))
    }
    
    new(
      'lcModelStan', 
      method = method,
      data = data,
      model = fit,
      permutations = list(),
      clusterNames = make.clusterNames(method$nClusters),
      tag = list(model = envir$stanmodel)
    )
  }
})

#. postFit ####
setMethod('postFit', signature('lcMethodStan'), function(method, data, model, envir, verbose, ...) {
  if(is(model, 'lcModelStanEM')) {
    return(model)
  }
  
  stanmod = model@model

  # compute latent class assignments per sample and chain
  model@model@.MISC$classifications = as.array(model, 'pp') %>% 
    apply(c(1, 3, 4), which.max)
  dimnames(model@model@.MISC$classifications) = list(Sample = NULL, Id = NULL, Chain = NULL)
  
  return(model)
})


make.standata = function(method, data) {
  library(stringr)
  library(R.utils)
  library(rstan)
  
  f = method$formula
  data = as.data.table(data)
  setkeyv(data, c(method$id, method$time))
  assert_that(
    attr(terms(f), 'response') == 1,
    !latrend:::hasIntercept(f),
    is.factor(data[[method$id]]) || is.integer(data[[method$id]])
  )
  response = all.vars(f)[1]
  mf = model.frame(f, data) 
  G = method$nClusters
  
  remZeros = function(x) {
    ifelse(x == 0, NA, x)
  }
  
  times = data[[method$time]] %>% unique() %>% sort()
  firstObsIdx = data[, which(get(method$time) == min(times))]
  standata = list(
    N = nrow(data),
    y = mf[[1]],
    G = G,
    I = uniqueN(data[[method$id]]),
    ii = as.integer(data[[method$id]]),
    offsets = cumsum(c(1, rle(as.integer(data[[method$id]]))$lengths)), 
    use_priors = method$use_priors,
    use_sigmas = method$use_sigmas,
    use_cor = method$use_cor,
    prior_theta = as.array(rep(1, G)),
    prior_cfc = 1,
    prior_nu_alpha = 2,
    prior_nu_beta = .1,
    prior_interceptMu_mean = quantile(mf[[1]][firstObsIdx], probs = seq(.1, .9, length.out=G)) %>% as.array(),
    prior_interceptMu_sigma = rep(sd(mf[[1]]), G) %>% as.array(),
    prior_interceptSigma_mean = data[, sd(diff(get(response)), na.rm=TRUE), by=c(method$id)]$V1 %>% remZeros() %>% log() %>% mean(na.rm=TRUE),
    prior_interceptSigma_sigma = 1,
    prior_sigma_nu = 3,
    prior_sigma_sigma = 1,
    prior_sdzMu_nu = 3,
    prior_sdzMu_sigma = 1
  )
  
  assert_that(
    length(standata$offsets) == 1 + standata$I,
    min(standata$ii) == 1,
    max(standata$ii) == standata$I,
    all(standata$ii %between% c(1, standata$I))
  )
  
  if(hasName(method, 'use_sdzs')) {
    standata$use_sdzs = method$use_sdzs
  }
  
  # fixed effects
  fNames = grep('^formula\\.', names(method), value = TRUE)
  fList = as.list(method, args = fNames)
  names(fList) = str_match(fNames, 'formula\\.(\\w+)')[, 2] %>% tolower()
  fList$mu = update(f, NULL ~ .)
  
  for(dpar in names(fList)) {
    f.par = fList[[dpar]] %>% update( ~ . - Fit)
    xmat = model.matrix(f.par, data=data)
    standata[[paste0('B_', dpar)]] = ncol(xmat)
    standata[[paste0('x_' , dpar)]] = as.array(xmat * 1)
    # priors
    standata[[paste0('prior_beta', capitalize(dpar), '_nu')]] = 3
    standata[[paste0('prior_beta', capitalize(dpar), '_mean')]] = 0
    standata[[paste0('prior_beta', capitalize(dpar), '_sigma')]] = 1
  }
  
  # random effects
  mfRan = model.matrix(method$random, data=data)
  standata$L_mu = ncol(mfRan)
  if(has_name(method, 'L_mu')) {
    message('Overriding L_mu for L_mu = ', method$L_mu)
    standata$L_mu = method$L_mu
  }
  standata$u_mu = mfRan
  assert_that(nrow(mfRan) == nrow(data))
  
  # priors
  priors = as.list(method, grep('^prior_', names(method), value=TRUE))
  
  modifyList(standata, priors)
}