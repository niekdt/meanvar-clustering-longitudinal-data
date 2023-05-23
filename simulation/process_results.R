sim_init()

exp = 'finalc'

# quick test
dtTest = processSim(
  exp, 
  lcMethod == 'lcMethodStanGMM' & data.fixed == 'Equal' & data.cv == 'HighDesc' & nClusters == 2 & data.tStart == 0
)

# model recovery (only for matching K=nClusters)
gmmFile = file.path(RESULTS_DIR, sprintf('sim_%s_gmm.rds', exp))
dtGmm = processSim(
  exp, 
  lcMethod == 'lcMethodStanGMM' & data.numGroups == nClusters &
    data.tStart == 0 &
    (data.fixed != 'Equal' | data.cv != 'None'),
  # results = readRDS(gmmFile),
  maxCases = Inf
)
saveRDS(dtGmm, gmmFile)

mvFile = file.path(RESULTS_DIR, sprintf('sim_%s_gmmMv.rds', exp))
dtGmmMv = processSim(
  exp, 
  lcMethod == 'lcMethodStanGMMMeanVar' & data.numGroups == nClusters &
    data.tStart == 0 &
    adapt_delta == .9 &
    (data.fixed != 'Equal' | data.cv != 'None'),
  # results = readRDS(mvFile),
  maxCases = Inf
)
saveRDS(dtGmmMv, mvFile)

rvFile = file.path(RESULTS_DIR, sprintf('sim_%s_gmmRv.rds', exp))
dtGmmRv = processSim(
  exp, 
  lcMethod == 'lcMethodStanGMMRV' & data.numGroups == nClusters &
    data.tStart == 0 &
    adapt_delta == .9 &
    (data.fixed != 'Equal' | data.cv != 'None'),
  # results = readRDS(rvFile),
  maxCases = Inf
)
saveRDS(dtGmmRv, rvFile)

rmvFile = file.path(RESULTS_DIR, sprintf('sim_%s_gmmRmv.rds', exp))
dtGmmRmv = processSim(
  exp, 
  lcMethod == 'lcMethodStanGMMRMV' & data.numGroups == nClusters &
    data.tStart == 0 &
    adapt_delta == .95 &
    (data.fixed != 'Equal' | data.cv != 'None'),
  # results = readRDS(rmvFile),
  maxCases = Inf
)
saveRDS(dtGmmRmv, rmvFile)

dtGmm = readRDS(gmmFile)
dtGmmMv = readRDS(mvFile)
dtGmmRv = readRDS(rvFile)
dtGmmRmv = readRDS(rmvFile)

dtFinal = rbindlist(
  list(
    GMM = dtGmm,
    `MV-GMM` = dtGmmMv,
    `RV-GMM` = dtGmmRv,
    `RMV-GMM` = dtGmmRmv
  ), 
  fill = TRUE, idcol = 'Model') %>%
  .[, Model := factor(Model, levels = c('GMM', 'MV-GMM', 'RV-GMM', 'RMV-GMM'))]

saveRDS(dtFinal, file.path(RESULTS_DIR, sprintf('sim_%s.rds', exp)))

# numclass recovery ####
# quick test
dtTest = processClassRecovery('final', lcMethod == 'lcMethodStanGMM' & nClusters == 2 & data.numGroups == 2 & data.seed == 1)

##
dtClassGmm = getCaseTable('final', lcMethod == 'lcMethodStanGMM') %>% processSimCaseTable(paramRecovery = FALSE)
setattr(dtClassGmm, 'models', NULL)
dtClassGmmMv = getCaseTable('final', lcMethod == 'lcMethodStanGMMMeanVar') %>% processSimCaseTable(paramRecovery = FALSE)
setattr(dtClassGmmMv, 'models', NULL)
dtClassGmmRmv = getCaseTable('final', lcMethod == 'lcMethodStanGMMRMV') %>% processSimCaseTable(paramRecovery = FALSE)
setattr(dtClassGmmRmv, 'models', NULL)

dtClass = rbindlist(list(
  GMM = dtClassGmm,
  `MV-GMM` = dtClassGmmMv,
  `RMV-GMM` = dtClassGmmRmv), 
  fill = TRUE, idcol = 'Model'
) %>%
  .[, Model := factor(Model, levels = c('GMM', 'MV-GMM', 'RMV-GMM'))]

saveRDS(dtClass, file.path(RESULTS_DIR, 'sim_class.rds'))


# class_mv
dtAll = getCaseTable('derp') %>% processSimCaseTable(paramRecovery = FALSE)
setattr(dtAll, 'models', NULL)
saveRDS(dtAll, file.path(RESULTS_DIR, 'sim_class_mv.rds'))


exp = 'class_gmmrmv_final2'
dtAll = processClassRecovery(exp, subset = data.numTraj == 250, deleteMissingCases = TRUE)
saveRDS(dtAll, file.path(RESULTS_DIR, paste0(exp, '.rds')))
# dtAll = readRDS(file.path(RESULTS_DIR, paste0(exp, '.rds')))

# testing
dtAll = processClassRecovery(
  'class_gmmv_final', 
  subset = data.randomSigma == 'High' & data.fixed == 'PartialLinear' & data.seed <= 2 & nClusters >= 2 & seed == 1, 
  keepModels = TRUE
)
models = attr(dtAll, 'models')
