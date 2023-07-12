library(assertthat)

setClass('lcModelStanEM', contains='lcModelStan')

names.lcModelStanEM = function(x) {
  assert_that(is.list(x@model$par))
  names(x@model$par)
}

flatnames.lcModelStanEM = function(x) {
  colnames(x@model$theta_tilde)
}

mean.lcModelStanEM = function(x, par) {
  assert_that(is.list(x@model$par))
  yraw = x@model$par[[par]]
  
  if(length(dim(yraw)) %in% c(0, 1)) {
    y = as.numeric(yraw)
    if(length(y) == 0) {
      return(y)
    }
    else if(length(y) == 1) {
      # scalar
      names(y) = par
      return(y)
    }
    else {
      # vector
      names(y) = paste0(par, '[', seq_along(y), ']')
      return(y)
    }
  } else {
    return(yraw)
  }
}

coef.lcModelStanEM = function(object) {
  parNames = names(object)
  
  suppressWarnings({ 
    parLens = lapply(object@model$par, dim) %>% sapply(prod) 
  })
  
  pars = setdiff(parNames, union(parNames[parLens >= 100], c('pp', 'log_lik', 'lp__', 'cfc')))
  
  coefs = lapply(pars, function(p) {
    out = mean(object, p)
    if(is.matrix(out)) {
      idxStrings = expand.grid(1:nrow(out), 1:ncol(out)) %>%
        as.matrix() %>%
        apply(1, paste0, collapse=',')
      names(out) = paste0(p, '[', idxStrings, ']')
    }
    return(out)
  })
  
  unlist(coefs)
}

setMethod('converged', signature('lcModelStanEM'), function(object) {
  object@model$return_code == 0
})


setMethod('postprob', signature('lcModelStanEM'), function(object) {
  if(not('pp' %in% names(object))) {
    return(callNextMethod())
  }
  
  pp = mean(object, 'pp') %>% t()
  
  assert_that(nrow(pp) == nIds(object),
              ncol(pp) == nClusters(object))
  colnames(pp) = clusterNames(object)
  assert_that(is_valid_postprob(pp, object))
  return(pp)
})


logLik.lcModelStanEM = function(object) {
  ll = sum(object@model$par$log_lik)
  N = nIds(object)
  df = length(coef(object)) - 1 #-1 for theta
  attr(ll, 'nobs') = N
  attr(ll, 'df') = df
  class(ll) = 'logLik'
  return(ll)
}