library(data.table)
library(assertthat)
library(magrittr)

#' @description Get the group trend(s) table
tsdata_trends = function(data) {
  attr(data, 'trend')
}

#' @description Renames the standard tsdata columns
tsdata_standardize = function(data, time='Time', value='Value', id='Id', group='Group') {
  xdata = copy(data)
  
  newnames = c('Group', 'Id', 'Time', 'Value')
  oldnames = c(group, id, time, value)
  mask = oldnames %in% names(xdata) & newnames != oldnames
  
  # remove already existing columns
  delnames = intersect(newnames[mask], names(xdata))
  if(length(delnames) > 0) {
    xdata[, c(delnames) := NULL]
  }
  
  setnames(xdata, oldnames[mask], newnames[mask])
  
  # restore key because this tends to go wrong elsewhere... sometimes
  if('Id' %in% names(xdata)) {
    setkey(xdata, Id, Time)
  }
  
  return(xdata)
}

#' @description Combine a list of tsdata tables, including the trend tables
#' @param datasets list of tsdata tables
tsdata_merge = function(datasets) {
  stopifnot(!is.data.frame(datasets))
  stopifnot(is.list(datasets))
  
  # construct group trends table
  dt_trends = lapply(datasets, attr, 'trend') %>% 
    rbindlist() %>% 
    setkey(Group, Time)
  
  # create merged dataset
  dt_all = do.call(rbind, datasets) %>%
    setkey(Group, Id, Time) %>%
    setattr('trend', dt_trends)
  
  return(dt_all)
}

#' @description Get the vector of group assignment per trajectory
tsdata_groups = function(data) {
  if('Group' %in% names(data)) {
    data[, Group[1], by=Id]$V1
  } else {
    data[, NA*0, by=Id]$V1
  }
}

# sets the TIME column of trends to the column in tsdata specified by 'newtime'
tstrends_adopt_time = function(trends, data, newtime) {
  stopifnot(is.data.table(trends))
  stopifnot(is.data.table(data))
  stopifnot('Time' %in% names(trends))
  
  if(newtime == 'Time') {
    return(copy(trends))
  }
  
  stopifnot(c('Time', newtime) %in% names(data))
  xtrends = trends[unique(data[, c('Time', newtime), with=FALSE]), on='Time'] %>%
    .[, Time := NULL] %>%
    setnames(newtime, 'Time') %>%
    setkeyv(c(CLUSTER, 'Time'))
  return(xtrends)
}


tsdata_matrix = function(data, value='Value') {
  assert_that(is.data.frame(data), has_name(data, c('Id', 'Time', value)))
  dtWide = dcast(data, Id ~ Time, value.var=value)
  
  dataMat = as.matrix(dtWide[, -'Id'])
  assert_that(nrow(dataMat) == uniqueN(data$Id), ncol(dataMat) == uniqueN(data$Time))
  rownames(dataMat) = dtWide$Id
  colnames(dataMat) = names(dtWide)[-1]
  return(dataMat)
}


tsdata_plot = function(data, response='Value', trends=TRUE, facet=TRUE, alpha=1) {
  library(ggplot2)
  p = ggplot() +
    geom_line(data=data, aes_string(x='Time', y=response, color='Group', group='Id'), alpha=alpha) +
    labs(title=sprintf('Trajectories for %s', response))
  
  dt_trends = tsdata_trends(data)
  if(!isFALSE(trends) && has_name(dt_trends, response)) {
    p = p + geom_line(data=dt_trends, aes_string(x='Time', y=response, group='Group'), color='black', size=1)
  }
  
  if(facet) {
    p = p + facet_wrap(~Group)
  }
  
  return(p)
}