library(latrend)
# dataset mu trajectories are approximately between [-1, 2]
# what about high variance for low mu case?

plot_gtsdata = function(data) {
  pval = plotTrajectories(data, response = 'Value', cluster = 'Group', facet = FALSE) + ggtitle('Value trajectories')
  pmu = plotTrajectories(data, response = 'Mu', cluster = 'Group', facet = FALSE) + ggtitle('Mu trajectories')
  psigma = plotTrajectories(data, response = 'Sigma', cluster = 'Group', facet = FALSE) + ggtitle('Sigma trajectories')
  psigmadistr = ggplot(data, aes(x=Sigma, color=Group)) + geom_density() + ggtitle('Sigma distribution')
  gridExtra::grid.arrange(pval, pmu, psigma, psigmadistr, ncol=2)
}

dataOptions = list()
dataOptions$`1` = list(
  # intended to be similar to the multi-group datasets, so has larger variability to match multiple groups
  fixed = list(
    Constant=list(c(.5, 0, 0)), # no use
    Linear=list(c(0, 1, 0)), #no use
    Quad=list(c(0, 1.2, -.5))),
  random = list(
    None = diag(rep(1e-20, 3)), 
    Low = diag(.1 / c(1, 2, 10)), #no use
    LowInt = diag(.1 / c(1, 1e100, 1e100)), #no use
    High = diag(.2 / c(1, 2, 10))),
  sigma = list( # standard deviation
    None = .001, #no use
    Low = .05,
    Med = .1,
    High = .2),
  randomSigma = list( # variance
    None = 0,
    Low = .1,
    Med = .2,
    High = .3),
  cv = list(
    None = 0,
    Low = .3,
    Med = .45,
    High = .6)
)
    
dataOptions$`2` = list(
  fixed = list(
    Isolated = list(c(.5, 1, -.2), c(-.5, -1, .2)), # completely non-overlapping group trajectories
    Distinct = list(c(-.5, 1, -.2), c(1, -1, .2)), # clearly different coefficients, but overlapping group trajectories
    DistinctLinear = list(c(-.5, 1, 0), c(1, -1, 0)),
    Partial = list(c(-.25, 1, -.15), c(0, 0, 0)), # same intercept
    PartialLinear = list(c(-.25, 1, 0), c(0, 0, 0)),
    Equal = list(c(-.25, 1.1, -.2), c(0, .9, -.2)),
    EqualLinear = list(c(-.25, 1.1, 0), c(0., .9, 0))
    ), # very similar group trajectories
  random = list( # covariance matrix
    None = dataOptions$`1`$random[['None']],
    VeryLowInt = diag(c(.01, 1e-10, 1e-10)),
    VeryLowLinear = diag(c(.01, .01 / 5, 1e-10)),
    VeryLow = diag(c(.01 / c(1, 5, 25))),
    LowLinear = diag(c(.025, .025 / 5, 1e-10)),
    LowInt = diag(c(.025, 1e-10, 1e-10)),
    Low = diag(.025 / c(1, 5, 25)),
    MedInt = diag(c(.05, .05 / 1e-10, 1e-10)),
    MedLinear = diag(c(.05, .05 / 2.5, 1e-10)),
    Med = diag(.05 / c(1, 2.5, 10)),
    Med2 = diag(.05 / c(1, 2.5, 10)),
    HighInt = diag(c(.1, .1 / 1e-10, 1e-10)),
    HighLinear = diag(c(.1, .1 / 5, 1e-10)),
    High = diag(.1 / c(1, 5, 25))
  ),
  sigma = dataOptions$`1`$sigma,
  randomSigma = dataOptions$`1`$randomSigma,
  cv = list(
    None = c(0, 0),
    LowDesc = c(dataOptions$`1`$cv[['Low']], 0),
    MedDesc = c(dataOptions$`1`$cv[['Med']], 0),
    HighDesc = c(dataOptions$`1`$cv[['High']], 0)
  )
)
  
dataOptions$`3` = list(
  fixed = list(
    Isolated = list(c(.75, .5, -.1), c(0, 0, 0), c(-.75, -.5, .1)),
    Distinct = list(c(-1, 1, 0), c(0, -1.1, .2), c(1, 0, 0)),
    DistinctLinear = list(c(-1, 1, 0), c(0, -1.1, 0), c(1, 0, 0)),
    Partial = list(dataOptions$`2`$fixed$Partial[[1]], c(1.25, 0, -.35), dataOptions$`2`$fixed$Partial[[2]]),
    Equal = list(dataOptions$`2`$fixed$Equal[[1]], c(1.25, 0, -.35), dataOptions$`2`$fixed$Equal[[2]])), # first and last group are practically equal
  random = dataOptions$`2`$random,
  sigma = dataOptions$`2`$sigma,
  randomSigma = dataOptions$`2`$randomSigma,
  cv = list(
    None = c(0, 0, 0),
    LowDesc = c(dataOptions$`1`$cv[['Low']], dataOptions$`1`$cv[['Low']] / 2, 0),
    MedDesc = c(dataOptions$`1`$cv[['Med']], dataOptions$`1`$cv[['Med']] / 2, 0),
    HighDesc = c(dataOptions$`1`$cv[['High']], dataOptions$`1`$cv[['High']] / 2, 0))
)



if(sys.nframe() == 0L) {
  # Generate and visualize
  data = gen_gtsdata(numGroups = 2, 
                     numTraj = 100, 
                     numObs = 15, 
                     fixed = 'Equal', 
                     random = 'Med2', 
                     sigma = 'High', 
                     randomSigma='Med', 
                     cv='HighDesc',
                     seed = 1) %T>% 
    plot_gtsdata()
  # lmer(Value ~ Time + I(Time^2) + (Time | Id), data = data) %>% summary()
}
