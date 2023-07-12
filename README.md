# meanvar-clustering-longitudinal-data
Supplementary materials for the manuscript "Latent-class trajectory modeling with a heterogeneous mean-variance relation" by N. G. P. Den Teuling, F. Ungolo, S.C. Pauws, and E.R. van den Heuvel

This repository contains the source code for the conditional growth mixture Stan models, the estimation of the marginal loglikelihood thereof, and for running and analyzing the simulation study and case study.

## Setup
1. Install R.
2. Create an `.Rprofile` file with the following content, and fill in the placeholders:
```R
source("renv/activate.R")
FIG_DIR = 'figs'
RESULTS_DIR = 'results'
TABLES_DIR = 'tables'
COVID_DATA_DIR = '~/data/csse_covid_19_data' # set to correct folder
REDIS_HOST_FILE = file.path('redis', 'redis_host.txt') # used by worker.R to connect to the Redis server
options(
  redis.host = 'localhost',
  redis.port = 6379,
  redis.pwd = '', # set password if configured
  latrend.warnMetricOverride = FALSE,
  mc.cores = parallel::detectCores(logical = FALSE)
)

source('include.R')
```
3. Install the required packages via `renv::restore()` or manually.

## Setting up and using the batch job environment
In case you want to run the simulation or case study, proceed with the next steps.
1. Install Redis.
2. Configure Redis and the credentials in the `.Rprofile` file.
3. Run Redis
4. Test if you can connect from R by running `sim_init()`.  
5. Submit jobs, e.g., by running `sim_all.R`.
6. Run `worker.R` as one or more stand-alone processes, e.g., by executing `worker6.bat` on Windows.
7. Wait a long time for computations to finish.
8. Collect and process results (in case of the simulation study) by running `process_results.R`.

# Case study results
![image](https://github.com/niekdt/meanvar-clustering-longitudinal-data/assets/8193083/6531a62d-e94e-4e8b-aed3-29358014d081)

![image](https://github.com/niekdt/meanvar-clustering-longitudinal-data/assets/8193083/a22a21de-f961-4953-8d97-dd6d180cf52e)
