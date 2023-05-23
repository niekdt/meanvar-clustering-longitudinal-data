sim_init = function() {
  library(simworkr)
  library(assertthat)
  library(latrend)
  library(R.utils)
  library(stackoverflow)
  library(uuid)
  source(file.path('redis', 'redis.R'))
  
  assert_that(
    redis_isActive(), 
    msg = sprintf('redis is not running (host file "%s" not present', REDIS_HOST_FILE)
  )
	
  host_info = read.table(REDIS_HOST_FILE, stringsAsFactors=FALSE)
  
  options(redis.host = host_info$V1,
          redis.port = as.integer(host_info$V2),
          mc.cores = 1)
  
  message('Connecting to ', host_info$V1, ':', host_info$V2, '...')
  
  simworkr::connect(timeout = 10)
  
  if(!isConnected()) {
    stop('not connected to Redis server!')
  }
}

redis_isActive = function() {
  return(file.exists(REDIS_HOST_FILE))
}
