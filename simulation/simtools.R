generate_method_cases = function(methodClass, dataNames, model, ..., seed = 1, chains = 1) {
  assert_that(is.character(methodClass))
  mfun = get(methodClass)
  assert_that(is.function(mfun))
  stanfile = file.path('models', paste0(model, '.stan'))
  message('Compiling Stan model ', stanfile, '...')
  stan_model(file = stanfile, model_name = basename(stanfile), auto_write = TRUE)
  
  mllFiles = grep(model, list.files('loglik', pattern = 'stan$'), value = TRUE)
  for (file in mllFiles) {
    message('Compiling Stan mloglik model ', file, '...')
    stan_model(
      file = file.path('loglik', file), 
      model_name = basename(file), 
      allow_undefined = TRUE,
      auto_write = TRUE,
      includes = sprintf('\n#include "%s"\n', file.path(getwd(), 'loglik', 'get_iter.hpp'))
    )
  }
  
  message('Generating method cases...')
  mRef = mfun() %>% evaluate()
  methods = do.call(lcMethods, list(mRef, model = model, ..., seed = seed, chains = chains))
  
  methodCases = cbind(lcMethod = methodClass, as.data.frame(methods)) %>%
    as.data.table()
  methodCases[, .class := NULL]
  return(methodCases[])
}

generate_data_cases = function(
  numGroups = 2, 
  fixed = 'Isolated', 
  random = 'Low',
  sigma = 'None', 
  randomSigma = 'None',
  cv = 'None', 
  numTraj = 100,
  numObs = 10,
  tStart = 0,
  tEnd = 2,
  seed = 1
) {
  mc = stackoverflow::match.call.defaults()
  args = as.list(mc)[formalArgs(gen_gtsdata)]
  
  do.call(expand.grid, args) %>% 
    as.data.table() %>%
    setnames(paste0('data.', names(args))) %>%
    .[]
}

create_cases = function(methodCases, dataCases, methodFun = 'latrendFit', dataFun = 'caseData') {
  assert_that(
    is.data.frame(methodCases),
    is.data.frame(dataCases),
    nrow(methodCases) > 0,
    nrow(dataCases) > 0
  )
  
  tidyr::crossing(methodCases, dataCases) %>% 
    as.data.table() %>%
    .[, method := methodFun] %>%
    .[, data := dataFun] %>%
    setcolorder(c('method', 'lcMethod')) %>%
    .[]
}

print_cases = function(cases) {
  cols = names(cases)[sapply(cases, uniqueN) > 1]
  print(subset(cases, select=cols))
}

caseData = function(...) {
  args = list(...)
  
  dataArgs = args[sapply(names(args), startsWith, 'data.')]
  names(dataArgs) = sapply(names(dataArgs), substring, 6)
  
  data = do.call(gen_gtsdata, dataArgs)
  attr(data, 'args') = dataArgs
  return(data)
}

covidWeekData = function(data.state, data.minDate, data.maxDate, data.excludeZeroCounties, data.minNormNew) {
  args = match.call() %>% as.list() %>% `[`(-1)
  dataArgs = args[sapply(names(args[-1]), startsWith, 'data.')]
  
  alldata = load_csse_covid19_data() %>%
    .[Population > 0] %>%
    .[Lat != 0] %>%
    .[County != ''] %>%
    .[is.finite(FIPS)]
  
  if(data.state == 'all') {
    data.state = unique(alldata$State)
  }
  minDate = as.Date(data.minDate)
  maxDate = as.Date(data.maxDate)
  
  if (data.excludeZeroCounties) {
    fips = alldata[, max(CumConfirmed) > 0, keyby = FIPS][V1 == TRUE, FIPS]
    message('Dropping ', uniqueN(alldata$FIPS) - length(fips), ' counties without any cases.')
  } else {
    fips = unique(alldata$FIPS)
  }
  
  data = alldata %>%
    .[State %in% data.state] %>%
    .[FIPS %in% fips]

  allWeekData = data[, .(
      Confirmed = max(CumConfirmed), 
      Population = max(Population),
      FirstDate = min(Date),
      LastDate = max(Date)), 
      keyby=.(FIPS, State, County, Week)] %>%
    .[, NormConfirmed := Confirmed / Population * 1e5] %>%
    .[, NewConfirmed := c(0, diff(Confirmed)), by=.(FIPS)] %>%
    .[, NormNewConfirmed := NewConfirmed / Population * 1e5] %>%
    .[, NormNewConfirmed4 := NewConfirmed / Population * 1e4] %>%
    .[, NormNewConfirmed3 := NewConfirmed / Population * 1e3]
  
  weekData = allWeekData[FirstDate >= minDate & LastDate <= maxDate] %>%
    .[, Id := as.integer(factor(FIPS))] %>%
    .[, Time := (Week - min(Week)) / (max(Week) - min(Week))]

  message(comma(uniqueN(weekData$FIPS)), ' counties within date range with data')
  message(comma(nrow(weekData)), ' observations within date range')
  message(comma(weekData[, sum(NormNewConfirmed == 0)]), ' observations are zero')
  message(comma(weekData[NormNewConfirmed > 0, sum(NormNewConfirmed < data.minNormNew)]), ' observations are below the required minimum log-count')
  
  weekData[, LogNormNewConfirmed := ifelse(NormNewConfirmed == 0, log(data.minNormNew), log(NormNewConfirmed))]
  weekData[LogNormNewConfirmed < log(data.minNormNew), LogNormNewConfirmed := log(data.minNormNew)]
  
  nErrDays = weekData[, sum(diff(Confirmed) < 0), by=.(FIPS)][, sum(V1)]
  message('Days with negative new cases: ', comma(nErrDays), ' (', percent(nErrDays / nrow(data)), ')')

  attr(weekData, 'args') = dataArgs
  return(weekData)
}



latrendFit = function(data, lcMethod, ..., .experiment, .case) {
  library(splines)
  lcMethod = as.character(lcMethod)
  
  mfun = get(lcMethod)
  assert_that(is.function(mfun))
  args = list(...)
  args$formula = as.formula(args$formula)
  args$random = as.formula(args$random)
  
  if(isTRUE(is.na(args$formula.sigma))) {
    args$formula.sigma = NULL
  }
  
  dataArgNames = grep('^data\\.', names(args), value = TRUE)
  m = do.call(mfun, args[setdiff(names(args), dataArgNames)])
  
  assert_that(is.lcMethod(m))
  
  if (hasName(m, 'maxFits')) {
    n = m$maxFits
  } else {
    n = 10L
  }
  
  if (is.finite(n) && n > 1) {
    model = latrend(lcFitConverged(m, maxRep = n), data = data)
    model@method = m  
  } else {
    model = latrend(m, data = data)
  }
  
  model = autopermute(model)
  
  if (hasName(m, 'margLogLikIter') && m$margLogLikIter > 0) {
    message('Computing marginal log likelihood...')
    
    if (hasName(m, 'margLogLikCache')) {
      cacheMll = m$margLogLikCache
    } else {
      cacheMll = TRUE
    }
    
    mll = compute_marginal_loglik(
      model, 
      M = m$margLogLikIter,
      mixType = m$margLogLikMixType,
      cache = cacheMll
    )
    model@model@.MISC$mWAIC = mWAIC(mll)
    message('\tSuccess.')
  }
  
  message(sprintf('Model size: %s MiB', object.size(model) %>% format('Mb')))
  stripModel = strip(model)
  message(sprintf('Stripped model size: %s MiB', object.size(stripModel) %>% format('Mb')))

  dir.create(file.path(RESULTS_DIR, .experiment), showWarnings = FALSE)
  fname = paste0(UUIDgenerate(use.time=TRUE), '.rds')
  file = file.path(RESULTS_DIR, .experiment, fname)
  saveRDS(stripModel, file, compress = TRUE)
  
  return(file)
}

getCaseSettings = function(table) {
  table[, lapply(.SD, uniqueN)]
}

getCaseTable = function(exp, ..., deleteMissingCases = FALSE) {
  expData = getExperimentResults(exp, ...)
  
  allCasesTable = lapply(expData, '[[', 'case') %>% rbindlist(fill = TRUE)
  message(nrow(allCasesTable), ' results in this experiment')
  
  models = getCaseModels(
      exp = exp, 
      data = expData, 
      deleteMissingCases = deleteMissingCases
  )
  
  assert_that(
    length(models) == nrow(allCasesTable)
  )
  
  validMask = !sapply(models, is.null)
  if (!all(validMask)) {
    message(sprintf('Dropping %d results due to missing model', sum(!validMask)))
    models = models[validMask]
    casesTable = allCasesTable[validMask, ]
  } else {
    casesTable = allCasesTable
  }
  
  setattr(casesTable, 'models', models)
  
  casesTable[]
}



viewCaseTable = function(table, metric = 'AdjustedRand') {
  table[, .(Min = min(get(metric)), 
    Mean = mean(get(metric)), 
    Max = max(get(metric)), 
    N=.N, 
    i = paste0(i, collapse='; ')), 
    keyby=.(data.numGroups, data.fixed, data.random, data.sigma, data.randomSigma, data.cv, data.numObs)] %>% 
    View()
}

loadCaseModelFile = function(file, deleteMissingCases = FALSE) {
  assert_that(
    is.character(file),
    is.flag(deleteMissingCases)
  )
  
  # FIX
  if (RESULTS_DIR != 'results') {
    # correct path for old cases
    if (stringr::str_starts(file, 'results/')) {
      file = stringr::str_replace(
        string = caseFiles[fixMask], 
        pattern = 'results/', 
        replacement = paste0(RESULTS_DIR, '/')
      )
    }
  }
  
  if (file.exists(file)) {
    model = readRDS(file)
    model
  } else if (deleteMissingCases) {
    message('\tDeleting case.')
    missingCase = Hmisc::escapeRegex(caseKeys[i])
    deleteCases(exp = exp, pattern = missingCase, sim = FALSE)
    NULL
  } else {
    NULL
  }
}

loadCaseModel = function(data, ...) {
  assert_that(has_name(data, 'output'))
  caseFile = data[['output']]
  loadCaseModelFile(caseFile, ...)
}

getCaseModels = function(exp, data, deleteMissingCases = FALSE) {
  assert_that(is.character(exp))
  if (missing(data)) {
    data = getExperimentResults(exp)
  }
  caseKeys = names(data)
  assert_that(!is.null(caseKeys))
  caseFiles = sapply(data, '[[', 'output')
  message('To load ', length(caseFiles), ' model files...')
  
  models = vector('list', length = length(caseFiles))
  for (i in seq_along(caseFiles)) {
    if (i %% 10 == 0) {
      message('...', percent(i / length(caseFiles)))
    }
    mod = loadCaseModelFile(caseFiles[i], deleteMissingCases = deleteMissingCases)
    if (!is.null(mod)) {
      models[[i]] = mod
    }
  }
  message('Done.')
  return(models)
}

object.Mb = function(x, serialize = FALSE) {
  if (serialize) {
    x = base::serialize(x, connection = NULL)
  }
  as.numeric(object.size(x)) / 1024 / 1024
}
