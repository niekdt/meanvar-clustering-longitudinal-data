library(data.table)
library(magrittr)

load_csse_covid19_data = function(dir = COVID_DATA_DIR, filter = TRUE) {
  dataDir = dir
  library(scales)
  library(zoo)
  library(usmap)
  
  usddata = parse_covid_data(file = file.path(dataDir, 'csse_covid_19_time_series', 'time_series_covid19_deaths_US.csv'), 'Deaths')
  uscdata = parse_covid_data(file = file.path(dataDir, 'csse_covid_19_time_series', 'time_series_covid19_confirmed_US.csv'), 'Confirmed')
  
  alldata = usddata[uscdata, Confirmed := i.Confirmed] %>%
    setnames('Province_State', 'State') %>%
    setnames('Admin2', 'County')
  
  alldata[, RevDay := as.numeric(max(Date) - Date)]
  alldata[, State := factor(State)]
  alldata[, County := factor(County)]
  
  if (filter) {
    message('Dropping US terroritories and virtual counties...')
    states = levels(alldata$State)[which(!is.na(usmap::fips(levels(alldata$State))))]
    alldata = alldata[State %in% states & Lat != 0]
    alldata[, State := factor(State)]
    alldata[, County := factor(County)]
  }
  
  # adjust Confirmed cases for future corrections
  message(sprintf('Updating cumulative confirmed cases in consideration of future corrections: %s out of %s...', 
    alldata[, sum(diff(Confirmed) < 0), by=Combined_Key]$V1 %>% sum %>% comma, 
    comma(nrow(alldata))))
  
  alldata[, MedConfirmed := rollmedian(c(-Inf, Confirmed, Inf), k = 3, align = 'center'), by=Combined_Key]
  
  alldata[, CumConfirmed := rev(zoo::rollapplyr(
    rev(MedConfirmed), 
    width = length(MedConfirmed), 
    FUN = min, 
    partial = TRUE)), by=Combined_Key]
  
  alldata[, NormCumConfirmed := CumConfirmed / Population * 1e5]
  alldata[, NormCumConfirmedAvg3 := movmeanc(NormCumConfirmed, 3), by=Combined_Key]
  alldata[, NormCumConfirmedAvg7 := movmeanc(NormCumConfirmed, 7), by=Combined_Key]
  
  # diff
  alldata[, NewConfirmed := c(0, diff(CumConfirmed)), by=Combined_Key]
  alldata[, NormNewConfirmed := c(0, diff(NormCumConfirmed)), by=Combined_Key]
  
  # helper columns
  alldata[, Week := lubridate::isoweek(Date)]
  alldata[, DayOfWeek := lubridate::wday(Date, label = TRUE, week_start = 1)]
  
  assert_that(min(alldata$NewConfirmed) >= 0)
  return(alldata)
}


parse_covid_data = function(file, value.name) {
  rawdata = fread(file)
  dateColumns = grep('\\d+/\\d+/\\d+', names(rawdata), value = TRUE)
  
  data = melt(rawdata, id.vars = setdiff(names(rawdata), dateColumns),
    variable.name = 'Date', 
    value.name = value.name)
  
  data[, Date := as.Date(Date, format = '%m/%d/%y')]
  setkey(data, Country_Region, Province_State, Admin2, Date)
  data[]
}

movsumr = function(x, n) {
  Reduce('+', shift(x, 0:(n-1), fill=0))
}

movmeanr = function(x, n) {
  Reduce('+', shift(x, 0:(n-1), fill=0)) / pmin(n, seq_along(x))
}

movmeanc = function(x, n) {
  nHalf = floor(n / 2)
  Reduce('+', shift(x, -nHalf:nHalf, fill=0)) / pmin(n, 
    c(seq(nHalf + 1, ceiling(length(x)/2) + nHalf), 
      seq(floor(length(x)/2) + nHalf, nHalf + 1)))
}
