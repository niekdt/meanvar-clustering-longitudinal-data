suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(assertthat)
  library(polynom)
  library(MASS)
  library(splines)
  library(magrittr)
  library(mvnfast)
  library(latrend)
  library(rstan)
  library(bayesplot)
  library(lme4)
  library(loo)
  library(matrixStats)
  library(scales)
  library(stringr)
  library(HDInterval)
  library(weights)
  
  eval(parse('util.R', encoding = 'UTF-8'))
  
  source(file.path('data', 'tsdata.R'))
  source(file.path('data', 'tsdata_builder.R'))
  source(file.path('data', 'gtsdata_builder.R'))
  source(file.path('data', 'util_datasets.R'))
  source(file.path('data', 'datasets_sim.R'))
  source(file.path('data', 'data_covid.R'))
  
  source(file.path('models', 'lcMethodStan.R'))
  source(file.path('models', 'lcModelStan.R'))
  source(file.path('models', 'lcModelStanEM.R'))

  source(file.path('simulation', 'simtools.R'))
  source(file.path('simulation', 'simproc.R'))
  source(file.path('redis', 'redis.R'))
})

rstan_options('auto_write' = TRUE)

theme_minimal(base_size = 9) %+replace%
  theme(
    plot.background = element_rect(colour = NA),
    plot.margin = unit(c(2,1,1,1), 'mm'),
    panel.background = element_rect(colour = NA),
    panel.spacing = unit(1, 'mm'),
    strip.background = element_rect(colour=NA, fill=NA),
    strip.text = element_text(face='plain', size=7, margin=margin()),
    legend.text = element_text(size=7),
    legend.title = element_text(size=9),
    legend.position = 'right',
    legend.spacing = unit(1, 'cm'),
    legend.margin = margin(0,0,0,0),
    legend.key.size = unit(9, 'pt'),
    legend.box.margin = margin(0,0,-5,0),
    axis.line.x = element_line(color='black', size = .1),
    axis.line.y = element_line(color='black', size = .1)
  ) %>%
  theme_set()
sigfig

message('Initialized scripts.')

message('Ensure that the Redis server is running (by starting redis/redis.bat). Then run sim_init() here to connect.')
