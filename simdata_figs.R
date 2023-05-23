colors = c('grey70', 'black', 'grey50')

# Fixed effects ####
dataOpts = list(
  data.numTraj = 100, 
  data.numObs = 10, 
  data.random = 'Med2', 
  data.randomSigma = 'None', 
  data.sigma = 'Med',
  data.cv = 'None')

dataP2 = do.call(caseData, c(dataOpts, data.fixed = 'Partial', data.numGroups = 2))
dataP3 = do.call(caseData, c(dataOpts, data.fixed = 'Partial', data.numGroups = 3))
dataE2 = do.call(caseData, c(dataOpts, data.fixed = 'Equal', data.numGroups = 2))
dataE3 = do.call(caseData, c(dataOpts, data.fixed = 'Equal', data.numGroups = 3))

dataP = rbind(dataP2, dataP3, idcol = 'K')
dataE = rbind(dataE2, dataE3, idcol = 'K')
data = rbind(dataP, dataE, idcol = 'Dataset') %>%
  .[, Dataset := factor(Dataset, levels = 1:2, labels = c('Partial overlap', 'Full overlap'))] %>%
  .[, K := factor(K, levels = 1:2, labels = c('K = 2', 'K = 3'))]

trendsP = rbind(tsdata_trends(dataP2), tsdata_trends(dataP3), idcol = 'K')
trendsE = rbind(tsdata_trends(dataE2), tsdata_trends(dataE3), idcol = 'K')
trends = rbind(trendsP, trendsE, idcol = 'Dataset') %>%
  .[, Dataset := factor(Dataset, levels = 1:2, labels = c('Partial overlap', 'Full overlap'))] %>%
  .[, K := factor(K, levels = 1:2, labels = c('K = 2', 'K = 3'))]

ggplot(data, aes(x = Time, y = Value, group = Id, color = Group)) +
  scale_color_manual(name = 'Class', values = colors) +
  geom_line(size = .1) +
  geom_line(data = trends, mapping = aes(y = Mu, group = Group), color = 'white', size = 1.5) +
  geom_line(data = trends, mapping = aes(y = Mu, group = Group), color = 'black', size = 1) +
  guides(color = FALSE) +
  coord_cartesian(ylim = c(-.5, 1.5), expand = 0) +
  labs(y = expression(y['i,j'])) +
  #facet_grid(~ Dataset + K) +
  facet_grid(rows = vars(K), cols = vars(Dataset)) +
  theme(axis.title.y = element_text(vjust = .5, angle = 0),
    plot.margin=margin(),
    panel.spacing = unit(4, 'mm'),
    strip.text = element_text(margin=margin(5,5,5,5)))

ggsave(filename = file.path(FIG_DIR, 'simdata.pdf'), width = 7, height = 7, units = 'cm')

# Mean-variance ####
data = caseData(
  data.numGroups = 3,
  data.fixed = 'Equal',
  data.numTraj = 50, 
  data.numObs = 10, 
  data.random = 'Med2',
  data.randomSigma = 'None', 
  data.sigma = 'Med',
  data.cv = 'HighDesc')
trends = tsdata_trends(data)

ggplot(data, aes(x = Time, y = Value, group = Id, color = Group)) +
  scale_color_manual(values = colors) +
  geom_line(size = .1) +
  geom_line(data = trends, mapping = aes(y = Mu, group = Group), color = 'white', size = 1.5) +
  geom_line(data = trends, mapping = aes(y = Mu, group = Group), color = 'black', size = 1) +
  coord_cartesian(expand=0) +
  labs(y = expression(y['i,j'])) +
  guides(color = FALSE) +
  theme(axis.title.y = element_text(vjust = .5, angle = 0),
    plot.margin=margin(r=6))
ggsave(filename = file.path(FIG_DIR, 'simmv.pdf'), width = 3.7, height = 4, units = 'cm')

ggplot(data, aes(x = Time, y = Sigma, group = Id, color = Group)) +
  scale_color_manual(values = colors) +
  scale_y_continuous(breaks = pretty_breaks(4)) +
  geom_line(size = .1) +
  coord_cartesian(expand=0) +
  expand_limits(y = 0) +
  guides(color = FALSE) +
  labs(y = expression(sigma[epsilon][',i,j'])) +
  theme(axis.title.y = element_text(vjust = .5, angle = 0, margin = margin(r = -5)),
    plot.margin=margin(r=5))
ggsave(filename = file.path(FIG_DIR, 'simmv-sigma.pdf'), width = 4, height = 4, units = 'cm')
