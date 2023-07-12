suppressPackageStartupMessages({
  library(R.utils)
  source('util.R')
  
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

sigfig

message('Initialized scripts.')

message('Ensure that the Redis server is running (by starting redis/redis.bat). Then run sim_init() here to connect.')
