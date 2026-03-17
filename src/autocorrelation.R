suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2) ; theme_set(theme_bw())
  library(lubridate)
  library(stringr)
  library(patchwork)
})

hosp = readRDS('../data/ari_hosp.rds')
ww = readRDS('../data/nml-database.RDS') |>
  mutate(date = ceiling_date(ymd(collectiondatetime),
                             unit = 'week'))


# --- Autocorrelation analysis

hosp_pacf <- function(p, h) {
  df = h |> 
    filter(prov == p, virus == 'COVID',
           source == 'DAD',
           date > ymd('2020-01-01'))
  
  y = df$count
  a = pacf(y, plot = FALSE)
  lags = 1:10
  b = Box.test(y, lag = max(lags), type = "Ljung-Box")
  
  res = data.frame(
    prov = p, 
    lag = lags,
    pacf = a$acf[lags],
    boxtest = unname(b$statistic)
  )
  return(res)
}

z = lapply(c('ON','BC','AB','SK','MB'), hosp_pacf, h = hosp) |> 
  bind_rows()

g.hosp = z |> ggplot(aes(x=lag, y=pacf)) + 
  facet_wrap(~prov) + 
  geom_col()+ 
  scale_x_continuous(breaks = seq(1,20,by=2))+
  labs(title = 'Hospital admissions', y = 'PACF')




ww_pacf <- function(city) {
  dfw = ww |> 
    filter(
      grepl(paste0('^',city,'[A-Z][A-Z]$'), siteid), 
      measureid == 'covN2', 
      fractionid == 'sol'
    ) |> 
    group_by(date) |> 
    summarise(v = mean(value, na.rm = TRUE), .groups = 'drop')
  
  # ggplot(dfw, aes(x=date, y=v)) + geom_line()
  
  y = dfw$v
  a = pacf(y, plot = FALSE)
  lags = 1:10
  b = Box.test(y, lag = max(lags), type = "Ljung-Box")
  
  res = data.frame(
    city = city, 
    lag = lags,
    pacf = a$acf[lags],
    boxtest = unname(b$statistic)
  )
  return(res)
}

zw = lapply(c('V', 'E', 'W', 'T'), ww_pacf) |> 
  bind_rows()

g.ww = zw |> ggplot(aes(x=lag, y=pacf)) + 
  facet_wrap(~city) + 
  geom_col()+ 
  scale_x_continuous(breaks = seq(1,20,by=2))+
  labs(title = 'Wastewater', y = 'PACF')


g = wrap_plots(g.hosp, g.ww)
pdf('../figs/plot-autocorrelation.pdf', width = 15)
plot(g)
dev.off()

# PACF drops sharply after lag 1 ==> AR(1) type. 
# ==> regression should include H(t) ~ H(t-1) + W(t-1)
# (unit of time is _week_)

