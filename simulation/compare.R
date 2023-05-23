sim_init()
library(betareg)
options(try.outFile = stdout())

# Gather results ####
dtGbtm = getCaseTable('gbtm', data.fixed != 'Distinct') %>% processCaseTable()
gbtms = attr(dtGbtm, 'models')

dtGmm = getCaseTable('gmm', data.fixed != 'Distinct') %>% processCaseTable()
gmms = attr(dtGmm, 'models')

dtGmmMv = getCaseTable('gmmv', data.fixed != 'Distinct') %>% processCaseTable()
gmmMvs = attr(dtGmmMv, 'models')

gmmRvs = getCaseModels('gmmrv')
dtGmmRv = getCaseTable('gmmrv') %>% processCaseTable()

gmmRmvsResults = getExperimentResults('gmmrmv')
gmmRmvs = getCaseModels(gmmRmvsResults)
dtGmmRmv = getCaseTable(gmmRmvsResults) %>% processCaseTable()

dtAll = rbind(dtGmm, dtGmmMv, dtGmmRv, dtGmmRmv, fill = TRUE) %>%
  .[, lcMethodMV := lcMethod %in% c('lcMethodStanGMMMeanVar', 'lcMethodStanGMMMRMV')]

dataArgs = paste0('data.', c('fixed', 'random', 'sigma', 'cv', 'randomSigma'))
ss = quote(data.fixed %in% c('Partial', 'Equal') & data.numObs == 10 & data.numGroups == 2 & data.numTraj == 250 & adapt_delta == .7 & data.cv %in% c('None', 'LowDesc', 'HighDesc'))

# Model stability ####
# gbtm is incredibly stable across all datasets and conditions
dtGbtm[eval(ss), .(Mean = mean(AdjustedRand), SD = sd(AdjustedRand), .N), keyby=c(dataArgs)]
dtGmm[eval(ss), .(Mean = mean(AdjustedRand), SD = sd(AdjustedRand), .N), keyby=c(dataArgs)]
dtGmmMv[eval(ss), .(Mean = mean(AdjustedRand), SD = sd(AdjustedRand), .N), keyby=c(dataArgs)]
dtGmmRv[eval(ss), .(Mean = mean(AdjustedRand), SD = sd(AdjustedRand), .N), keyby=c(dataArgs)]
dtGmmRmv[eval(ss), .(Mean = mean(AdjustedRand), SD = sd(AdjustedRand), .N), keyby=c(dataArgs)]

# numTraj comparison
mod = lm(AdjustedRand ~ data.numTraj + data.cv + data.random + data.sigma + data.randomSigma, data = dtGmmRmv, subset = data.numGroups == 2); coef(mod)
mod = lm(Rhat.sigma ~ data.numTraj + data.cv + data.random + data.sigma + data.randomSigma, data = dtGmmRmv, subset = data.numGroups == 2); coef(mod)

# GBTM ####
dtGbtm[eval(ss), .(AdjustedRand = mean(AdjustedRand), WMAE = mean(WMAE)), keyby=.(data.fixed, data.sigma, data.cv)]
lm(AdjustedRand ~ data.fixed + data.random + data.sigma + data.cv, data = dtGbtm, subset = eval(ss)) %>% summary()

# GMM ####
dtGmm[eval(ss), .(AdjustedRand = mean(AdjustedRand), WMAE = mean(WMAE)), keyby=.(data.fixed, data.cv)]
lm(AdjustedRand ~ data.fixed + data.random + data.sigma * data.cv, data = dtGmm, subset = eval(ss)) %>% summary()

# Bias
dtGmm[eval(ss), .(int1 = mean(`Bias.fixed[k=1, i=1]`), int2 = mean(`Bias.fixed[k=2, i=1]`)), keyby=.(data.fixed, data.cv)]
lm(`Bias.fixed[k=1, i=1]` ~ data.fixed + data.random + data.cv + data.sigma, data = dtGmm, subset = eval(ss)) %>% summary()
lm(`Bias.fixed[k=2, i=1]` ~ data.fixed + data.random + data.cv + data.sigma, data = dtGmm, subset = eval(ss)) %>% summary()

# Rhat
dtGmm[eval(ss), .(int1 = median(Rhat.intercept)), keyby=.(data.fixed, data.cv)]
lm(`Bias.fixed[k=1, i=1]` ~ data.fixed + data.random + data.cv + data.sigma, data = dtGmm, subset = eval(ss)) %>% summary()
lm(`Bias.fixed[k=2, i=1]` ~ data.fixed + data.random + data.cv + data.sigma, data = dtGmm, subset = eval(ss)) %>% summary()


# GMM-MV ####
dtGmmMv[data.numGroups == 2, .(ARI = mean(AdjustedRand), Rhat.intercept = mean(Bias.cv)), keyby=.(data.fixed, data.sigma, data.cv)]
lm(AdjustedRand ~ data.fixed + data.random + data.fixed * data.cv, data = dtGmmMv, subset = eval(ss)) %>% summary()

# Bias intercept
dtGmmMv[eval(ss), .(int1 = mean(`Bias.fixed[k=1, i=1]`), int2 = mean(`Bias.fixed[k=2, i=1]`)), keyby=.(data.fixed, data.cv)]
lm(`Bias.fixed[k=1, i=1]` ~ data.fixed + data.random + data.cv + data.sigma, data = dtGmm, subset = eval(ss)) %>% summary()
lm(`Bias.fixed[k=2, i=1]` ~ data.fixed + data.random + data.cv + data.sigma, data = dtGmm, subset = eval(ss)) %>% summary()

# Bias CV
dtGmmMv[data.numGroups == 2, mean(`Bias.cv[k=1]`), keyby=.(data.fixed, data.cv)]
lm(`Bias.cv[k=1]` ~ data.fixed + data.random + data.cv + data.sigma, data = dtGmmMv, subset = eval(ss)) %>% summary()
lm(`Bias.cv[k=2]` ~ data.fixed + data.random + data.cv + data.sigma, data = dtGmmMv, subset = eval(ss)) %>% summary()



# GMM-RV ####
dtGmmRv[eval(ss), mean(AdjustedRand), keyby=.(data.fixed, data.sigma, data.cv)]
lm(AdjustedRand ~ data.fixed + data.random, data = dtGmmRv, subset = eval(ss)) %>% summary()

# Bias
dtGmmRv[eval(ss), .(int1 = mean(`Bias.fixed[k=1, i=1]`), int2 = mean(`Bias.fixed[k=2, i=1]`)), keyby=.(data.fixed, data.cv)]
dtGmmRv[eval(ss), mean(`Bias.sigma`), keyby=.(data.fixed, data.cv)] # sigma

# GMM-RMV ####
dtGmmRmv[, .(int1 = mean(`Bias.fixed[k=1, i=1]`), int2 = mean(`Bias.fixed[k=2, i=1]`)), keyby=.(data.fixed, data.cv)]
dtGmmRmv[eval(ss), mean(AdjustedRand), keyby=.(data.fixed, data.randomSigma, data.sigma, data.cv)]
lm(AdjustedRand ~ adapt_delta + data.numTraj + data.numObs + data.fixed + data.random, data = dtGmmRmv, subset = eval(ss)) %>% summary()

# Bias
# lower CV bias for higher numObs
dtGmmRmv[eval(ss), mean(`Bias.sigma`), keyby=.(data.fixed, data.cv)] # sigma
dtGmmRmv[eval(ss), mean(`Bias.cv[k=1]`), keyby=.(data.fixed, data.cv, data.numObs)] # cv

lm(`Bias.cv[k=1]` ~ data.random + data.sigma + data.cv, data = dtGmmRmv, subset = eval(ss)) %>% summary()
lm(`Bias.cv[k=2]` ~ data.random + data.sigma + data.cv, data = dtGmmRmv, subset = eval(ss)) %>% summary()

# Overall effects ####
lm(AdjustedRand ~ lcMethod + data.fixed + data.random + data.sigma + data.randomSigma + data.cv, data = dtAll, subset = data.numGroups == 2) %>% summary()


# Overall, dataset-specific ####
# Distinct
dtAll[data.fixed == 'Distinct' & data.numGroups == 2 & data.numObs == 10 & data.numTraj == 250] %>%
  .[, .(C=mean(AdjustedRand), E=mean(WMAE), .N), 
    keyby=.(data.random, data.sigma, data.randomSigma, data.cv, lcMethod)]

# No significant difference on cv, which makes sense because the group trajectories are identifiable on mu
lm(AdjustedRand ~ lcMethod * data.cv + data.random + data.sigma + data.randomSigma + data.cv, data = dtAll, 
  subset = data.fixed == 'Distinct' & data.numGroups == 2 & data.numObs == 10 & data.numTraj == 250) %>% summary()

# No practical differences in WMAE between models. Only notable difference in randomHigh
lm(WMAE ~ lcMethod + data.random + data.sigma + data.randomSigma + data.cv, data = dtAll, 
  subset = data.fixed == 'Distinct' & data.numGroups == 2 & data.numObs == 10 & data.numTraj == 250) %>% summary()

lm(`Bias.fixed[k=1, i=1]` ~ lcMethod + data.random + data.sigma + data.randomSigma + data.cv, data = dtAll, 
  subset = data.fixed == 'Distinct' & data.numGroups == 2 & data.numObs == 10 & data.numTraj == 250) %>% summary()

lm(`Bias.fixed[k=2, i=1]` ~ lcMethod + data.random + data.sigma + data.randomSigma + data.cv, data = dtAll, 
  subset = data.fixed == 'Distinct' & data.numGroups == 2 & data.numObs == 10 & data.numTraj == 250) %>% summary()

lm(`Bias.fixed[k=1, i=2]` ~ lcMethod + data.random + data.sigma + data.randomSigma + data.cv, data = dtAll, 
  subset = data.fixed == 'Distinct' & data.numGroups == 2 & data.numObs == 10 & data.numTraj == 250) %>% summary()

lm(`Bias.fixed[k=2, i=2]` ~ lcMethod * data.sigma + data.random + data.sigma + data.randomSigma + data.cv, data = dtAll, 
  subset = data.fixed == 'Distinct' & data.numGroups == 2 & data.numObs == 10 & data.numTraj == 250) %>% summary()


# Partial
# No significant difference on cv, which makes sense because the group trajectories are identifiable on mu
dtAll[data.fixed == 'Partial' & data.numGroups == 2 & data.numObs == 10 & data.numTraj == 250] %>%
  .[, .(C=mean(AdjustedRand), E=mean(WMAE)), 
    keyby=.(data.random, data.sigma, data.randomSigma, data.cv, lcMethod)]
lm(AdjustedRand ~ lcMethod * data.cv + data.random + data.sigma + data.randomSigma + data.cv, data = dtAll, 
  subset = data.fixed == 'Partial' & data.numGroups == 2 & data.numObs == 10 & data.numTraj == 250) %>% summary()

# no practical differences in WMAE between models. Only notable differences on randomHigh and sigmaHigh
lm(WMAE ~ lcMethod + data.random + data.sigma + data.randomSigma + data.cv, data = dtAll, 
  subset = data.fixed == 'Partial' & data.numGroups == 2 & data.numObs == 10 & data.numTraj == 250) %>% summary()

# Equal
# significantly better performance for CV models
lm(AdjustedRand ~ lcMethod * data.cv + data.random + data.sigma + data.randomSigma + data.cv, data = dtAll, 
  subset = data.fixed == 'Equal' & data.numGroups == 2 & data.numObs == 10 & data.numTraj == 250) %>% summary()

# same results for partial
lm(WMAE ~ lcMethod + data.random + data.sigma + data.randomSigma + data.cv, data = dtAll, 
  subset = data.fixed == 'Equal' & data.numGroups == 2 & data.numObs == 10 & data.numTraj == 250) %>% summary()

# Similar WMAE across models. Makes sense as the correct identification of the variance does not reduce WMAE
mod2e = lm(WMAE ~ lcMethod + data.fixed + data.random + data.sigma + data.cv, data = dtAll, subset = data.randomSigma == 'None' & data.numGroups == 2) #& lcMethod == 'lcMethodStanGMMRMV')
coef(mod2e)
summary(mod2e)

# Bias ####
mod2b = lm(`Bias.random[i=1]` ~ lcMethod + data.random + data.sigma, data = dtAll, subset = data.randomSigma == 'None' & data.numGroups == 2) #& lcMethod == 'lcMethodStanGMMRMV')
coef(mod2b)

mod2b = lm(`Bias.random[i=2]` ~ lcMethod + data.random + data.sigma, data = dtAll, subset = data.randomSigma == 'None' & data.numGroups == 2) #& lcMethod == 'lcMethodStanGMMRMV')
coef(mod2b)

mod2b = lm(`Bias.random[i=3]` ~ lcMethod + data.random + data.sigma, data = dtAll, subset = data.randomSigma == 'None' & data.numGroups == 2) #& lcMethod == 'lcMethodStanGMMRMV')
coef(mod2b)

mod2b = lm(`Bias.sigma` ~ lcMethod + data.random + data.sigma, data = dtAll, subset = data.randomSigma == 'None' & data.numGroups == 2) #& lcMethod == 'lcMethodStanGMMRMV')
coef(mod2b)

mod2b = lm(`Bias.fixed[k=1, i=1]` ~ lcMethod + data.random + data.sigma, data = dtAll, subset = data.randomSigma == 'None' & data.numGroups == 2) #& lcMethod == 'lcMethodStanGMMRMV')
coef(mod2b)

mod2b = lm(`Bias.fixed[k=1, i=2]` ~ lcMethod + data.random + data.sigma, data = dtAll, subset = data.randomSigma == 'None' & data.numGroups == 2) #& lcMethod == 'lcMethodStanGMMRMV')
coef(mod2b)

mod2b = lm(`Bias.fixed[k=2, i=1]` ~ lcMethod + data.random + data.sigma, data = dtAll, subset = data.randomSigma == 'None' & data.numGroups == 2) #& lcMethod == 'lcMethodStanGMMRMV')
coef(mod2b)

mod2b = lm(`Bias.cv[k=1]` ~ lcMethod + data.random + data.sigma, data = dtAll[is.finite(`Bias.cv[k=1]`) & data.randomSigma == 'None' & data.numGroups == 2] %>% droplevels, na.action = na.omit) #& lcMethod == 'lcMethodStanGMMRMV')
coef(mod2b)



gmmMod = lm(AdjustedRand ~ data.fixed * data.cv + data.random + data.sigma + data.cv, data = dtGmm)

gmmvMod = lm(AdjustedRand ~ data.fixed * data.cv + data.random + data.sigma + data.cv, data = dtGmmMv)
gmmrmvMod = lm(AdjustedRand ~ data.fixed * data.cv + data.random + data.sigma + data.cv, data = dtGmmRmv)

gmmMod = lm(WMAE ~ data.fixed * data.cv + data.random + data.sigma + data.cv, data = dtGmm, subset = data.numGroups == 2)
gmmvMod = lm(WMAE ~ data.fixed * data.cv + data.random + data.sigma + data.cv, data = dtGmmMv, subset = data.numGroups == 2)
gmmrmvMod = lm(WMAE ~ data.fixed * data.cv + data.random + data.sigma + data.cv, data = dtGmmRmv, subset = data.numGroups == 2)



# dtGmm[, AdjustedRandB := (ifelse(AdjustedRand < 0, 0, AdjustedRand) * (.N-1) + .5) / .N] # See Smithson and Verkuilen (2006)
gmmMod = betareg::betareg(AdjustedRandB ~ data.fixed * data.cv + data.random + data.sigma, data = dtGmm)
gmmvMod = betareg::betareg(AdjustedRandB ~ data.fixed + data.cv + data.random + data.sigma, data = dtGmmMv)
gmmrmvMod = betareg::betareg(AdjustedRandB ~ data.fixed + data.cv + data.random + data.randomSigma + data.sigma, data = dtGmmRmv)

# Distinct data performance ####
dtGmm[data.fixed == 'Distinct', mean(AdjustedRand), keyby=.(data.numGroups)]
