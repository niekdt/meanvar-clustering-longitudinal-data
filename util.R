sigfig = function(x, digits = 2) {
  out = formatC(signif(x, digits = digits), digits = digits, big.mark = ',', format = 'fg', flag = '#')
  # remove leading zero
  out = sub('^(-)?0[.]', '\\1.', out)
  # removing ending zero
  out = sub('\\.$', '\\1', out)
  
  out[out == 'NaN'] = '-'
  
  out
}

pct = function(x, ...) {
  out = rep('-', length(x))
  out[is.finite(x)] = percent(x[is.finite(x)], ...)
  out
}

se = function(x) {
  sd(x, na.rm = TRUE) / sqrt(sum(is.finite(x)))
}

tsigfig = function(x, digits = 2, prec = -10) {
  if (log10(abs(x)) < prec) {
    sigfig(sign(x) * 10 ^ prec, 1)
  } else {
    sigfig(x, digits = digits)
  }
}

meanText = function(x, digits = 2, prec = -3) {
  m = mean(x, na.rm = TRUE)
  if (is.na(m)) {
    ''
  } else {
    tsigfig(m, digits = digits, prec = prec)
  }
}

meanSeText = function(x, mean.digits = 2, se.digits = 1, prec = -3) {
  mTxt = meanText(x, digits = mean.digits, prec = prec)
  s = se(x)
  
  if (is.na(s)) {
    mTxt
  } else {
    sTxt = tsigfig(s, digits = se.digits, prec = prec)
    paste0(mTxt, '±', sTxt)
  }
}

cohens_d = function(emmGrid, model) {
  assert_that(is(emmGrid, 'emmGrid'))
  
  if (missing(model)) {
    model = eval(attr(emmGrid, 'model.info')$call)
  }
  
  emmeans::eff_size(emmGrid, sigma = sqrt(mean(sigma(model)^2)), edf = df.residual(model))
}


emmTable = function(model, effects) {
  if (is.list(model) && !is(model, 'lm')) {
    lapply(model, emmTable, effects) %>%
      rbindlist(idcol = 'LM')
  } else {
    mapply(emmeans, list(model), effects, nesting = list(NULL)) %>% 
      lapply(as.data.table) %>%
      lapply(setnames, 1, 'Effect') %>%
      lapply(setnames, 2, 'Setting') %>%
      set_names(effects) %>%
      rbindlist(idcol = 'Effect')
  }
}

emm2Table = function(model, spec, effects, wide = TRUE) {
  if (is.list(model) && !is(model, 'lm')) {
    lapply(model, emm2Table, spec, effects, wide = wide) %>%
      rbindlist(idcol = 'LM', fill = TRUE)
  } else {
    dt = mapply(emmeans, list(model), spec, by = effects, nesting = list(NULL)) %>% 
      lapply(as.data.table) %>%
      lapply(setnames, 1, 'Model') %>%
      lapply(setnames, 2, 'Setting') %>%
      set_names(effects) %>%
      rbindlist(idcol = 'Effect')
    
    if (wide) {
      dt = dcast(dt, Effect + Setting ~ Model, value.var = 'emmean')
    }
    dt
  }
}


pairsTable = function(model, effects) {
  if (is.list(model) && !is(model, 'lm')) {
    lapply(model, pairsTable, effects) %>%
      rbindlist(idcol = '.model')
  } else {
    mapply(emmeans, list(model), effects, nesting = list(NULL)) %>% 
      lapply(pairs) %>%
      lapply(as.data.table) %>%
      set_names(effects) %>%
      rbindlist(idcol = 'Effect')
  }
}

pairsEffTable = function(model, effects) {
  if (is.list(model) && !is(model, 'lm')) {
    lapply(model, pairsEffTable, effects) %>%
      rbindlist(idcol = '.model')
  } else {
    mapply(emmeans, list(model), effects, nesting = list(NULL)) %>% 
      lapply(cohens_d) %>%
      lapply(as.data.table) %>%
      set_names(effects) %>%
      rbindlist(idcol = 'Effect')
  }
}
