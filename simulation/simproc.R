library(latrend)
library(data.table)
library(assertthat)
library(magrittr)

# Simulation study ####
#' @description Process simulation output of models into a results table
#' @param exp The experiment name
#' @param results A previously computed results table. Only cases not present yet will be computed, and added.
#' @param deleteMissingCases Whether to delete simulation cases which are missing the associated model file
#' @param dropInvalid Whether to discard invalid models from the results table
#' @param keepModels Whether to keep the model objects as part of the results table (requires a lot of memory!)
#' @param maxCases The maximum number of cases to process. By default, there is no limit.
#' @param ... Passed to [getExperimentResults]
processSim = function(
  exp, 
  ..., 
  results = NULL,
  deleteMissingCases = FALSE, 
  dropInvalid = TRUE,
  keepModels = FALSE,
  maxCases = Inf
) {
  assert_that(
    is.flag(keepModels),
    is.flag(dropInvalid),
    is.null(results) || is.data.table(results),
    is.count(maxCases)
  )
  library(foreach)
  message('Gathering cases...')
  expData = getExperimentResults(exp, ...)
  
  allCasesTable = lapply(expData, '[[', 'case') %>% rbindlist(fill = TRUE)
  
  message(sprintf('Found %d cases.', nrow(allCasesTable)))
  
  if (!is.null(results)) {
    assert_that(has_name(results, names(allCasesTable)))
    # setkeyv(allCasesTable, names(allCasesTable))
    # setkeyv(results, names(allCasesTable))
    newCasesTable = fsetdiff(allCasesTable, subset(results, select = names(allCasesTable)))
    
    if (nrow(newCasesTable) == nrow(allCasesTable)) {
      warning('No cases matched the results table', immediate. = TRUE)
    }
  } else {
    newCasesTable = allCasesTable
  }
  message(nrow(newCasesTable), ' cases to be processed')

  nCases = min(nrow(newCasesTable), maxCases)
  newResults = vector('list', nCases)
  for(i in seq_len(nCases)) {
    message(sprintf('Model %d/%d (%g%%)...', i, nCases, round(i / nCases * 100)))
    if (i %% 100 == 0) {
      gc()
    }
    case = newCasesTable[i, ]
    model = NULL
    tryCatch(
      {
        model = loadCaseModel(expData[[i]], deleteMissingCases = deleteMissingCases)
      }, 
      error = function(e) {
        warning(sprintf('loadCaseModel %s', e), immediate. = TRUE)
      }
    )
    
    out = NULL
    if (!is.null(model)) {
      tryCatch(
        {
          out = processSimModel(case, model)  
        },
        error = function(e) {
          warning(sprintf('processSimModel %s', e), immediate. = TRUE)
        }
      )
    }
    
    if (!is.null(out)) {
      if (keepModels) {
        setattr(out, 'model', model)
      }
      
      assert_that(all(lengths(out) == 1), msg = 'one or more metrics failed to compute')
      newResults[[i]] = as.data.frame(out)
    }
  }
  
  validMask = sapply(newResults, is.data.frame)
  if (!all(validMask)) {
    message(sprintf('Dropping %d results due to model processing errors', sum(!validMask)))
  }
  
  newResultsTable = rbindlist(newResults[validMask], fill = TRUE)
  if (keepModels) {
    models = lapply(newResults[validMask], attr, 'model')
    setattr(resultsTable, 'models', models)
  }
  
  gc()
  if (isFALSE(dropInvalid)) {
    newResultsTable = rbind(resultsTable, caseTable[!validMask], fill = TRUE)
  }
  
  if (is.null(results)) {
    newResultsTable
  } else {
    # append existing results
    if (length(intersect(names(results), names(newResultsTable))) != length(names(results))) {
      warning('results table columns do not match the new results')
    }
    
    rbind(results, newResultsTable, fill = TRUE)
  }
}

processSimModel = function(
  case, 
  model
) {
  assert_that(
    is.list(case),
    is(model, 'lcModelStan')
  )
  
  if (!isTRUE(validate.lcModelStan(model))) {
    message('Invalid model')
    return (NULL)
  }
  
  allParNames = coefNames(model)
  assert_that(length(allParNames) > 0L)
  
  argMap = c(
    theta = 'prop',
    sdz0_mu = 'random',
    sdz_mu = 'random',
    sdz0_sigma = 'randomSigma',
    sigma = 'sigma',
    cv = 'cv',
    intercept_mu = 'fixed',
    beta_mu = 'fixed'
  )
  
  # filter params which are not supported for comparison to the ref
  parNames = allParNames[which(!is.na(argMap[allParNames]))]
  # matching data args supported by the model
  dataArgs = na.exclude(argMap[allParNames]) %>% unique()
  
  # Compute results
  result = as.list(case)
  result$RunTime = metric(model, 'estimationTime')
  result$thetaMin = min(mean(model, 'theta'))
  result$hasEmptyCluster = any(clusterSizes(model) == 0)
  result$hasSolitaryCluster = any(clusterSizes(model) == 1)
  suppressMessages({
    # result$Entropy = metric(model, 'relativeEntropy')
    result$ConvProp = parameterConvergence(
      model, 
      pars = c('theta', 'intercept_mu', 'beta_mu', 'sdz0_mu', 'cv', 'sigma')
    ) %>% mean()
    result$AlphaCor = alphaCor(model)
  })
  
  # ARI
  refData = model.data(model) %>% as.data.table()
  assert_that(has_name(refData, 'Group'), msg = 'model.data has no reference Group column')
  refClus = refData[, first(Group), by=Id]$V1
  ref = lcModelPartition(refData, response = 'Value', trajectoryAssignments = refClus)
  
  result$AdjustedRand = externalMetric(model, ref, 'adjustedRand')
  #result$AdjustedRandB = (abs(result$AdjustedRand) * (.N-1) + .5) / .N # See Smithson and Verkuilen (2006)
  
  # Rhat
  rhats = Rhat(model, parNames) %>%
    set_names(., paste('Rhat', mapSimDataCoefs(names(.)), sep = '.'))
  result = c(result, rhats)
  
  # Parameter estimation and reference
  paramResults = computeParamRecovery(model, dataArgs)
  result = c(result, unlist(paramResults))
  
  # HDI
  hdiMat = hdi(model, parNames, credMass = .95)
  assert_that(is.matrix(hdiMat))
  hdiCoefs = mapSimDataCoefs(colnames(hdiMat))
  hdiVec = as.numeric(hdiMat) %>%
    set_names(
      paste(
        'HDI', 
        rep(hdiCoefs, each = 2), 
        tolower(rownames(hdiMat)), 
        sep = '.'
      )
    )
  result = c(result, hdiVec)
  
  result
}

# Class recovery study ####
#' @description Computes the relevant metrics for the given case model for 
#' class-recovery simulation scenarios
processClassRecoveryModel = function(
  case, 
  model, 
  simConds = c('data.numGroups', 'data.fixed', 'data.random', 'data.randomSigma', 'data.cv', 'data.seed', 'seed')
) {
  assert_that(
    is.list(case),
    is(model, 'lcModelStan'),
    is.character(simConds)
  )
  
  if (!isTRUE(validate.lcModelStan(model))) {
    message('Invalid model')
    return (NULL)
  }
  
  result = case
  result$RunTime = metric(model, 'estimationTime')
  result$thetaMin = min(mean(model, 'theta'))
  result$hasEmptyCluster = any(clusterSizes(model) == 0)
  result$hasSolitaryCluster = any(clusterSizes(model) == 1)
  suppressMessages({
    result$LogLik = logLik(model)
    result$BIC = BIC(model)
    result$mWAIC = metric(model, 'mWAIC')
    result$p_waic1 = p_waic1(model)
    result$p_waic2 = p_waic2(model)
    # result$Entropy = metric(model, 'relativeEntropy')
    result$ConvProp = mean(parameterConvergence(model))
    result$AlphaCor = alphaCor(model)
  })
  
  result
}

#' @param ... Arguments passed to [getExperimentResults()]
processClassRecovery = function(
  exp, 
  ..., 
  deleteMissingCases = FALSE, 
  dropInvalid = TRUE,
  simConds = c('data.numGroups', 'data.fixed', 'data.random', 'data.randomSigma', 'data.cv', 'data.seed', 'seed'),
  keepModels = FALSE
) {
  assert_that(
    is.flag(keepModels),
    is.flag(dropInvalid)
  )
  library(foreach)
  expData = getExperimentResults(exp, ...)
  
  allCasesTable = lapply(expData, '[[', 'case') %>% rbindlist(fill = TRUE)
  message(nrow(allCasesTable), ' cases to be processed')
  
  nCases = nrow(allCasesTable)
  results = foreach(i = seq_len(nCases), .errorhandling = 'pass') %do% {
    message(sprintf('Model %d/%d (%g%%)...', i, nCases, round(i / nCases * 100)))
    if (i %% 100 == 0) {
      gc()
    }
    case = allCasesTable[i, ]
    model = loadCaseModel(expData[[i]], deleteMissingCases = deleteMissingCases)
    out = processClassRecoveryModel(case, model)
    if (keepModels) {
      setattr(out, 'model', model)
    }
    model = NULL
    out
  }
  
  validMask = sapply(results, is.data.frame)
  if (!all(validMask)) {
    message(sprintf('Dropping %d results due to missing model', sum(!validMask)))
  }
  resultsTable = rbindlist(results[validMask], fill = TRUE)
  if (keepModels) {
    models = lapply(results[validMask], attr, 'model')
    setattr(resultsTable, 'models', models)
  }
  
  # check for completeness of simulation conditions
  table = copy(resultsTable)
  table[, i := 1:.N]
  dtSims = table[, c(simConds, 'nClusters'), with = FALSE] %>%
    .[, i := .GRP, by = simConds]
  assert_that(all(dtSims$nClusters == table$nClusters))
  if (!all(dtSims[, uniqueN(nClusters) >= 3, by = i]$V1 == TRUE)) {
    warning('some simulation conditions are not complete')
  }
  
  gc()
  if (dropInvalid) {
    return (resultsTable)
  } else {
    allResultsTable = rbindlist(list(resultsTable, caseTable[!validMask]), fill = TRUE)
    return (allResultsTable[])
  }
}

# Case study cases ####
processCaseStudy = function(
  exp, 
  ..., 
  results = NULL,
  deleteMissingCases = FALSE, 
  keepModels = FALSE,
  maxCases = Inf
) {
  assert_that(
    is.flag(keepModels),
    is.null(results) || is.data.table(results),
    is.count(maxCases)
  )
  message('Gathering cases...')
  expData = getExperimentResults(exp, ...)
  
  allCasesTable = lapply(expData, '[[', 'case') %>% rbindlist(fill = TRUE)
  
  message(sprintf('Found %d cases.', nrow(allCasesTable)))
  
  if (!is.null(results)) {
    assert_that(has_name(results, names(allCasesTable)))
    # setkeyv(allCasesTable, names(allCasesTable))
    # setkeyv(results, names(allCasesTable))
    newCasesTable = fsetdiff(allCasesTable, subset(results, select = names(allCasesTable)))
    
    if (nrow(newCasesTable) == nrow(allCasesTable)) {
      warning('No cases matched the results table', immediate. = TRUE)
    }
  } else {
    newCasesTable = allCasesTable
  }
  message(nrow(newCasesTable), ' cases to be processed')
  
  nCases = min(nrow(newCasesTable), maxCases)
  newResults = vector('list', nCases)

  for(i in seq_len(nCases)) {
    message(sprintf('Model %d/%d (%g%%)...', i, nCases, round(i / nCases * 100)))
    if (i %% 100 == 0) {
      gc()
    }
    case = newCasesTable[i, ]
    model = NULL
    tryCatch(
      {
        model = loadCaseModel(expData[[i]], deleteMissingCases = deleteMissingCases)
      }, 
      error = function(e) {
        warning(sprintf('loadCaseModel %s', e), immediate. = TRUE)
      }
    )
    out = NULL
    if (!is.null(model)) {
      tryCatch(
        {
          out = processCaseModel(case, model)  
        },
        error = function(e) {
          warning(sprintf('processCaseModel %s', e), immediate. = TRUE)
        }
      )
    }
    
    if (!is.null(out)) {
      if (keepModels) {
        setattr(out, 'model', model)
      }
      
      assert_that(all(lengths(out) == 1), msg = 'one or more metrics failed to compute')
      newResults[[i]] = as.data.frame(out)
    }
  }
  
  validMask = sapply(newResults, is.data.frame)
  if (!all(validMask)) {
    message(sprintf('Dropping %d results due to model processing errors', sum(!validMask)))
  }
  
  newResultsTable = rbindlist(newResults[validMask], fill = TRUE)
  
  if (keepModels) {
    models = lapply(newResults[validMask], attr, 'model')
    setattr(newResultsTable, 'models', models)
  }
  
  gc()
  
  if (is.null(results)) {
    newResultsTable
  } else {
    # append existing results
    if (length(intersect(names(results), names(newResultsTable))) != length(names(results))) {
      warning('results table columns do not match the new results')
    }
    
    rbind(results, newResultsTable, fill = TRUE)
  }
}

processCaseModel = function(case, model) {
  assert_that(
    is.list(case),
    is(model, 'lcModelStan')
  )

  if (!isTRUE(validate.lcModelStan(model))) {
    message('Invalid model')
    return (NULL)
  }
  
  result = case
  result$RunTime = metric(model, 'estimationTime')
  result$thetaMin = min(mean(model, 'theta'))
  result$hasEmptyCluster = any(clusterSizes(model) == 0)
  result$hasSolitaryCluster = any(clusterSizes(model) == 1)
  
  suppressMessages({
    result$mWAIC = metric(model, 'mWAIC')
    result$Entropy = metric(model, 'relativeEntropy')
    result$ConvProp = mean(parameterConvergence(model))
    result$AlphaCor = alphaCor(model)
  })
  
  result
}

# Helpers ####
#' @description dataCoef for all k and i
#' @return A named list of the dataCoef outputs
applyDataArgOptions = function(model, arg, fun = dataCoef) {
  assert_that(
    is.string(arg),
    is.function(fun)
  )
  dataArgs = model.data(model) %>% attr('args')
  
  g = dataArgs$numGroups
  assert_that(is.count(g))
  argOption = as.character(dataArgs[[arg]])
  
  if (arg == 'prop') {
    refValues = sqrt(1:g) / sum(sqrt(1:g))
  } else {
    assert_that(is.string(argOption))
    refValues = dataOptions[[g]][[arg]][[argOption]]
    assert_that(
      !is.null(refValues), 
      msg = sprintf('Could not find data option %s for argument %s', argOption, arg)
    )
  }
  
  if (is.list(refValues)) {
    # K x I matrix
    k = length(refValues)  
    out = lapply(seq_along(refValues),
      function(k) {
        lapply(seq_along(refValues[[k]]), function(i) fun(model, arg, i = i, k = k)) %>%
          set_names(paste0(arg, '[k=', k, ', i=', seq_along(refValues[[k]]), ']'))
      })
    unlist(out, recursive = FALSE)
  }
  else if (length(refValues) > 1) {
    if (arg == 'random') {
      refValues = diag(refValues)
    }
    
    if (arg %in% c('sigma', 'cv', 'prop')) {
      # vector of length K
      out = lapply(seq_along(refValues), function(k) fun(model, arg, k = k))
      names(out) = paste0(arg, '[k=', seq_along(refValues), ']')
    } else {
      # vector of length I
      out = lapply(seq_along(refValues), function(i) fun(model, arg, i = i))
      names(out) = paste0(arg, '[i=', seq_along(refValues), ']')
    }
    out
  }
  else {
    # scalar
    out = list(fun(model, arg))
    names(out) = arg
    out
  }
}

dataCoef = function(model, arg, i, k = 1, ...) {
  dataArgs = model.data(model) %>% attr('args')
  g = dataArgs$numGroups
  assert_that(
    !is.null(dataArgs),
    arg %in% c('prop', names(dataArgs)),
    missing(i) || is.count(i),
    missing(k) || is.count(k)
  )
  
  if (arg == 'prop') {
    if (nClusters(model) == 1) {
      y = 1
      yhat = 1
    } else {
      yhat = as.array(model, 'theta')[, k, ] %>% as.numeric()
      y = model.data(model) %>%
        as.data.table() %>%
        .[, first(Group), by=Id] %>% 
        `$`('V1') %>%
        table() %>%
        prop.table() %>% 
        as.numeric() %>% 
        .[k]
    }
  }
  else if (arg == 'fixed') {
    if (i == 1) {
      # intercept
      yhat = as.array(model, 'intercept_mu')[, k, ] %>% as.numeric()
      y = dataOptions[[g]]$fixed[[as.character(dataArgs$fixed)]][[k]][1]
    } else {
      yhat = as.array(model, 'beta_mu')[, i - 1, k, ] %>% as.numeric()
      y = dataOptions[[g]]$fixed[[as.character(dataArgs$fixed)]][[k]][i]
    }
  }
  else if (arg == 'random') {
    if (i == 1) {
      if ('sdz0_mu' %in% coefNames(model)) {
        yhat = as.array(model, 'sdz0_mu')[, k, ] %>% as.numeric()
      } else {
        warning('random-sdz0_mu not available for stan model')
        yhat = NA*0
      }
      y = dataOptions[[g]]$random[[as.character(dataArgs$random)]] %>% diag() %>% .[1] %>% sqrt()
    } else {
      if ('sdz_mu' %in% coefNames(model)) {
        yhat = as.array(model, 'sdz_mu')[, k, i - 1, ] %>% as.numeric()
      } else {
        warning('random-sdz_mu not available for stan model')
        yhat = NA*0
      }
      y = dataOptions[[g]]$random[[as.character(dataArgs$random)]] %>% diag() %>% .[i] %>% sqrt()
    }
    if (is.null(y) || all(is.na(y))) {
      y = 0
    }
  }
  else if (arg == 'sigma') {
    if('sigma' %in% coefNames(model)) {
      yhat = as.array(model, 'sigma') %>% as.numeric()
    } else {
      warning('sigma not available for stan model')
      yhat = NA*0
    }
    y = dataOptions[[g]]$sigma[[as.character(dataArgs$sigma)]]
  }
  else if (arg == 'randomSigma') {
    if('sdz0_sigma' %in% coefNames(model)) {
      yhat = as.array(model, 'sdz0_sigma') %>% as.numeric()
    } else {
      warning('sdz0_sigma not available for stan model')
      yhat = 0
    }
    y = dataOptions[[g]]$randomSigma[[as.character(dataArgs$randomSigma)]]
    if (is.null(y) || all(is.na(y))) {
      y = 0
    }
  }
  else if (arg == 'cv') {
    if('cv' %in% coefNames(model)) {
      yhat = as.array(model, 'cv')[, k, ] %>% as.numeric()
    } else 
    {
      warning('cv not available for stan model')
      yhat = 0
    }
    y = dataOptions[[g]]$cv[[as.character(dataArgs$cv)]][k]
    if (is.null(y) || all(is.na(y))) {
      y = 0
    }
  }
  else {
    stop('unsupported arg')
  }
  
  assert_that(is.numeric(y),
    length(y) == 1)
  return(list(yhat = yhat, y = y))
}

computeHDIs = function(model, pars, ...) {
  dtHDIs = lapply(models, hdi, pars, ...) %>% 
    lapply(as.data.table) %>%
    rbindlist(fill = TRUE) %>%
    .[, Bound := c('lower', 'upper') %>% rep(length(models))] %>%
    setcolorder('Bound')
  
  dtHDIs
}

computeParamEst = function(model, pars) {
  paramMeans = lapply(
    intersect(pars, names(model)), 
    function(par) {
      res = mean(model, par)
      if (is.vector(res)) {
        res
      } else {
        out = as.numeric(res)
        names(out) = attr(res, 'parNames')
        out
      }
    }
  )
  do.call(c, paramMeans)
}


computeParamRecovery = function(model, args) {
  paramRec = function(...) {
    out = dataCoef(...)
    list(ref = out$y, est = mean(out$yhat))
  }

  argResults = lapply(
    args, 
    function(arg) {
      applyDataArgOptions(model, arg, fun = paramRec) %>% 
        unlist() %>%
        as.list() %>% 
        as.data.table()    
    }
  )
  
  argResults
}


mapSimDataCoefs = function(params) {
  parMap = c(
    `theta[1]` = 'prop[k=1]',
    `theta[2]` = 'prop[k=2]',
    `theta[3]` = 'prop[k=3]',
    `intercept_mu[1]` = 'fixed[k=1, i=1]',
    `intercept_mu[2]` = 'fixed[k=2, i=1]',
    `intercept_mu[3]` = 'fixed[k=3, i=1]',
    `beta_mu[1,1]` = 'fixed[k=1, i=2]',
    `beta_mu[1,2]` = 'fixed[k=2, i=2]',
    `beta_mu[1,3]` = 'fixed[k=3, i=2]',
    `beta_mu[2,1]` = 'fixed[k=1, i=3]',
    `beta_mu[2,2]` = 'fixed[k=2, i=3]',
    `beta_mu[2,3]` = 'fixed[k=3, i=3]',
    `cv[1]` = 'cv[k=1]',
    `cv[2]` = 'cv[k=2]',
    `cv[3]` = 'cv[k=3]',
    `sdz0_mu[1]` = 'random[i=1]',
    `sdz_mu[1,1]` = 'random[i=2]',
    `sdz_mu[1,2]` = 'random[i=3]',
    `sdz0_sigma[1]` = 'randomSigma',
    `sigma[1]` = 'sigma')
  
  coefNames = parMap[params]
  assert_that(
    noNA(coefNames), 
    msg = sprintf(
      'unknown model parameters %s cannot be mapped to data coef', 
      paste(params[is.na(coefNames)], collapse = ', ')
    )
  )
  coefNames
}
