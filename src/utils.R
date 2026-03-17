## Utility functions ## 

# ---- Map city/province to dataframe----
#' @param dat dataframe
#' @param var mapping variable
map.prov <- function(dat, var = "city") {
  prov = read.csv("../prm/cityprov.csv")
  
  if(var == "city"){
    df = right_join(dat, prov, by = var) |>
      select(-province)
  }
  else if(var %in% c("prov", "province")){
    if("city" %in% names(dat)){
      dat = select(dat, -city)
    }
    df = right_join(dat, prov, by = var) |>
      select(-province)
  }
  return(df)
}

# ---- Apply date lag to dataframe ----
#' @param x Date
#' @param n lag (in weeks)
#' @return y lagged date

lag.date <- function(x, n){
  y = as_date(x) - n*7
  return(y)
}

# ---- Rename virus
#' @description Rename viruses from ww convention to clinical convention
#' @param x variable
  
rename.v <- function(x) {
  tmp = toupper(x)
  v = case_when(
    tmp == "COVN2" ~ "COVID",
    grepl("FLU*", tmp) ~ "FLU",
    tmp != "COVN2" | !grepl("FLU*", tmp) ~ tmp
  )
  return(v)
}

# ---- Aggregate WW data at city level ----
#' @param df Dataframe of wastewater data at site level
agg.ww <- function(df){
  pop = df %>%
    distinct(city, siteid, wwtp.pop) %>%
    group_by(city) %>%
    mutate(
      pct.pop = wwtp.pop/sum(wwtp.pop)) %>%
    ungroup() %>%
    mutate(across(pct.pop, ~na_if(.x, 0))) %>%
    filter(if_all(pct.pop, ~ !is.na(.x)))
  
  dfc = left_join(df, pop, by=c('siteid', 'city')) %>%
    group_by(city, date, virus, ww.lag) %>% 
    summarize(ww.conc = sum(value * pct.pop),  # weighted avg,
              cov = sum(pct.pop),         # total percentage pop coverage of
              # aggregated sites
              .groups='drop') %>%
    filter(cov >= 0.5) %>%
    select(-cov)
  
  return(dfc)
}

# ---- Prepare ww dataframe ----

#' @param ww Wastewater dataframe
#' @param lag Vector of time lags. Can be set to 0 if no lags required
#' 
prep.ww <- function(ww, lag){
  d = ww |>
    filter(measureid  == "covN2",
           fractionid == "sol") |>
    rename(virus = measureid,
           city = healthregion)
  
  dl = list()
  for(i in seq_along(lag)){
    tmp = d |> 
      mutate(virus = rename.v(virus),
             date = ceiling_date(collectiondatetime, "week")) |>
      # Calculate the weekly average
      group_by(siteid, date, virus, city, wwtp.pop) |> 
      summarize(value = mean(value, na.rm = TRUE), # take weekly conc avg.
                ww.lag = lag[i],
                .groups = 'drop')
    
    # Flow normalization option
    norm = as.logical(fetch_model_params()[["norm.ww"]])
    d = tmp
    
    if(norm){
      d = norm.ww(tmp)
    }
    
    dl[[i]] = d|>
      # Aggregate concentration 
      # at city level
      agg.ww() |>
      # Add province by mapping city 
      # to associated province.
      map.prov()
  }
  
  dj = bind_rows(dl)
  
  saveRDS(dj, '../out/ww.rds')
  
  return(dj)
}

# ---- Calculate city and total catchment population (for ms) ----
calc.pop <- function(){
  cities = read.csv(file = '../prm/cityprov.csv')[["city"]]
  ww = readRDS(file = '../data/nml-database.RDS')
  
  wwtp.pop = ww |>
    rename(city = "healthregion") |>
    filter(city %in% cities) |>
    distinct(siteid, city, wwtp.pop)
  
  city.pop = wwtp.pop |>
    group_by(city) |>
    summarise(
      pop = sum(wwtp.pop),
      .groups = 'drop'
    )
  
  total = summarise(city.pop,
                    total = sum(pop)) |>
    as.numeric()
  
  pct = total/36991981 * 100 # 2021 Canadian census estimate
  
  return(list(
    wwtp.pop = wwtp.pop,
    city.pop = city.pop,
    total = total,
    pct = pct
  ))
}

# ---- Prepare hosp dataframe ----
#' @param hosp Dataframe containing hospital admissions data
#' @param lag Vector of time lags. Can be set to 0.
#' @param src Source of data (e.g. DAD, CNISP, etc.)
prep.hosp <- function(hosp, lag, src = 'DAD'){
  dat = hosp |>
    filter(source %in% src,
           virus == "COVID",
           !is.na(city)) |>
    transmute(
      date = ceiling_date(date, "week"),
      source,
      virus,
      count,
      prov,
      city = case_when(
        city == "Montréal" ~ "Montreal",
        city == "StJohns" ~ "St. John's",
        !city %in% c("Montréal", "StJohns") ~ city
      ),
      geo
    ) |>
    # Fix DAD hosp data stitching issue around Marc/Arpil 2022
    group_by(date, source, virus, prov, city, geo) |>
    summarise(
      hosp.count = mean(count, na.rm = TRUE),
      .groups = 'drop')
    
  
  dl = list()
  for(i in seq_along(lag)){
    dl[[i]] = dat |>
      mutate(
        cl.lag = lag[i],
        date = lag.date(date, n = lag[i])
      )
  }
  dj = bind_rows(dl)
  
  # Option to normalize hosp by population
  norm = as.logical(fetch_model_params()[["norm.hosp"]])
  
  d = dj
  
  if(norm){
    d = norm.hosp(dj)
  }
  
  saveRDS(d, '../out/hosp.rds')
  return(d)
}

# ---- Retrieve city pop ----
retrieve.city.pop <- function(){
  pop = readRDS('../data/nml-database.RDS') |>
    rename(city = "healthregion") |>
    select(siteid, city, wwtp.pop) |>
    distinct() |>
    group_by(city) |>
    summarise(
      pop = sum(wwtp.pop),
      .groups = 'drop'
    )
  return(pop)
}

# ---- Normalize hosp to population ----
#' @param df Dataframe of hospital admissions
norm.hosp <- function(df){
  # retrieve pop of each city
  pop = retrieve.city.pop()
  
  # normalize hosp.count by population
  d = df |>
    left_join(pop, by = "city") |>
    mutate(hosp.count = hosp.count/pop * 100000)
  
  return(d)
}

# ---- Normalize ww concentration to flow ----
#' @param df Dataframe of hospital admissions
norm.ww <- function(df){
  # Retrieve flow data
  
  flow = readRDS('../data/ww-flow.rds') |>
    rename(siteid = "location") |>
    mutate(date_sampled = as.Date(date_sampled, format = "%m/%d/%Y"),
           date = ceiling_date(date_sampled, "week")) |>
    group_by(siteid, date) |>
    summarize(flow = mean(dinflvol, na.rm = TRUE),
              .groups = 'drop')
  
  # Merge with ww.conc
  d = df |>
    left_join(flow, by = c("siteid", "date")) |>
    mutate(value = case_when(
      city %in% c("Winnipeg", "Regina", "St. John's") ~ value, # currently no flow data
      !city %in% c("Winnipeg", "Regina", "St. John's") ~ value/flow 
    )) |>
    drop_na(value)
  
  return(d)
  
}

# ---- Prepare variant dataframe ----
#' @param var Dataframe of variant counts data
#' @param merge.vars Logical flag to merge variant lineages
prep.var <- function(var, merge.vars = TRUE){
  counts = get.var.count()
  
  dat = var |>
    bind_rows() |>
    transmute(
      date,
      province,
      virus = "COVID",
      variant = name,
      var.prop = value,
      dom.variant = dominant
    ) |>
    drop_na(var.prop, dom.variant) |>
    filter(var.prop > 0) |>
    left_join(counts, by = c("date", "variant", "province")) |>
    map.prov(var = "province")
  
  if(merge.vars){
    dat = merge.var(dat)
  }
  
  saveRDS(dat, '../out/var.rds')
  return(dat)
}

# ---- Get count data ----
get.var.count <- function(){
  counts = readRDS('../data/counts-by-prov.rds')
  
  df = lapply(counts, pivot_longer, cols = -c(province, epiweek),
              names_to = "variant", values_to = "count") |>
    bind_rows() |>
    rename(date = "epiweek") |>
    mutate(date = as_date(date)) |>
    filter(province != "Canada")
  return(df)
}

# ---- Summarise count data (for ms) ----
#' @param var.data Dataframe of variant counts
summarise.counts <- function(var.data){
  counts = var.data |>
    group_by(prov, date) |>
    summarise(counts = sum(count),
              .groups = 'drop') |>
    group_by(prov) |>
    summarise(min = round(min(counts, na.rm = TRUE)),
              mean = round(mean(counts, na.rm = TRUE)),
              max = round(max(counts, na.rm = TRUE)),
              .groups = 'drop')
  return(counts)
}

# ---- Retrieve pango lineage information ----
retrieve.pango <- function(){
  message("Retrieving Pango lineages..")
  
  network = CCT::checkCache(network = TRUE)
  vars = CCT::getVariants() |>
    CCT::expandChildren(network) |>
    rename(variant = "child")
  
  # Save pango lineage information (for ms)
  saveRDS(vars, file = '../out/pango-lineages.rds')
  
  message("Pango lineages retrieved.")
  return(vars)
}

# ---- Combine variants and their counts/proportions by pango lineage ----
#' @param df Variant dataframe
merge.var <- function(df){
  merge = retrieve.pango()
  d = df |>
    mutate_at("variant",
              sapply, replace.var, m = merge) |>
    group_by(city, prov, date, virus, variant) |>
    summarise(var.prop = sum(var.prop, na.rm = TRUE),
              count = sum(count),
              .groups = 'drop') |>
    group_by(city, prov, date, virus) |>
    mutate(dom.variant = variant[which.max(var.prop)]) |>
    ungroup()
  
  return(d)
}

# ---- Function to replace variant name ----
#' @param x variant name
#' @param m Dataframe
replace.var <- function(x, m){
  if(x %in% unique(m$variant)){
    y = m[m$variant == x, "lineage"]
  }
  else{
    y = x
  }
  return(y)
}

# ---- Filter data to start on emergence of specified variant ----
#' @param d Dataframe
#' @param var.start String. Name of starting variant
filter.var <- function(d,var.start){
  df = d |>
    arrange(prov, date) |>
    group_by(prov) |>
    mutate(start = first(date[dom.variant == var.start])) |>
    filter(date >= start) |>
    ungroup() |>
    select(-start)
  return(df)
}

# ---- Fetch variant parameters ----
fetch_var_params <- function(){
  params = fetch_model_params() |>
    select(contains("var")) |>
    mutate_all(as.numeric)
  return(params)
}

# ---- Filter hosp/ww data based on variant start and end dates ----
#' @param df Dataframe of hospital/ww data
#' @param dates Dataframe of dates for each variant
#' @param v Variant name
#' @param sensitivity Logical flag for use in sensitivity analysis
filter.data <- function(df, dates, v, sensitivity = FALSE){
  var = filter(dates, variant == v)
  message(paste0("Variant: ", v, "\n", var[["start"]], " to ", var[["end"]]))
  d = df |>
    filter(date >= var[["start"]] & date <= var[["end"]]) |>
    mutate(variant = v)
  if(sensitivity){
    d = d |>
      mutate(threshold = var$threshold,
             window = var$min_weeks,
             pct.above = var$pct.above)
  }
  return(d)
}

# ---- Create model object ---
#' @param df Dataframe of hospital/ww data
#' @param dates Dataframe of dates for each variant
#' @param sensitivity Logical flag for use in sensitivity analysis
create.mdl.obj <- function(df, dates, sensitivity = FALSE){
  message("Creating model object..")
  cities = unique(df$city)
  d = list()
  for(c in cities){
    message(paste0("Filtering data for: ", c))
    tmp.df = filter(df, city == c)
    tmp.dates = filter(dates, city == c)
    vars = unique(tmp.dates$variant)
    y = lapply(vars, filter.data,
               dates = tmp.dates, df = tmp.df,
               sensitivity)
    z = purrr::keep(y, ~nrow(.) >= 6)
    d[[c]] = bind_rows(z)
  }
  dat = bind_rows(d)
  if(!sensitivity){
    saveRDS(dat, file = "../out/obj-model.rds")
  }
  message("Complete.")
  return(dat)
}

# ---- Prepare data before MCMC inference
#' @param data.source Hospital admissions data source
#' @param logscale Logical flag whether or not to log-transform data
#' @param variant.ignore String of variant to exclude from data
#' @param sensitivity Logical flag if function being used during sensitivity
#'  analysis
#' @param obj Optional. Model object if function being used in sensitivity 
#'  analysis  
prepare_data <- function(data.source, logscale, variant.ignore,
                         sensitivity = FALSE,
                         obj = NULL) {
  if(!sensitivity){
    obs = readRDS("../out/obj-model.rds")
  }
  
  else if(sensitivity){
    obs = obj
  }
  
  dat  = obs |> 
    filter(
      source == data.source, 
      !variant %in% variant.ignore,
      ww.conc > 0,
      hosp.count > 0) |> 
    mutate(
      w = logscale * log(ww.conc) + (1-logscale) * ww.conc,
      h = logscale * log(hosp.count) + (1-logscale) * hosp.count,
      # Associate an index to variants and locations
      variant.idx = as.integer(as.factor(variant)),
      geo.idx     = as.integer(as.factor(geo))
    ) |> 
    # insert time-lagged values
    group_by(geo, variant) |> 
    arrange(date, .by_group = TRUE) |> 
    mutate(w_1 = lag(w, n=1),
           h_1 = lag(h, n=1)) |> 
    filter(! is.na(w_1), ! is.na(w_1))
  
  if(!sensitivity){
    # Save index and name link 
    vidx = dat |> select(variant, variant.idx) |>
      distinct() |> arrange(variant.idx)
    write.csv(vidx, file = '../out/idx_variant.csv', quote = F, row.names = F)
    
    vgeo = dat |> select(geo, geo.idx) |>
      distinct() |> arrange(geo.idx)
    write.csv(vgeo, file = '../out/idx_geo.csv', quote = F, row.names = F)
  }
  
  return(dat)
}

# ---- Create simulated ts series data (for model optimization) ----
#' @param obs dataframe of real data
#'
create_sim_data <- function(obs){
  
  if(1) set.seed(123)
  
  message("Creating simulated data from observed time-series...\n\n")
  # Retrieve sim parameters
  params = read.csv("../prm/sim-params.csv", strip.white = TRUE)
  
  # universal level
  mu   = params[params$param == "mu", "value"] 
  mu1  = params[params$param == "mu1", "value"] 
  beta = params[params$param == "beta", "value"] 
  
  # variant level
  n.var = params[params$param == "n.var", "value"]
  if(n.var == 0){
    n.var = length(unique(obs$variant)) # number of variants (based on obs data)
  }
  sigma_M  = params[params$param == "sigma_M", "value"]  
  sigma_M1 = params[params$param == "sigma_M1", "value"]  
  sigma_B  = params[params$param == "sigma_B", "value"]  # sd of beta
  
  M  = rnorm(n = n.var, mean = mu,   sd = sigma_M)
  M1 = rnorm(n = n.var, mean = mu1,  sd = sigma_M1)
  B  = rnorm(n = n.var, mean = beta, sd = sigma_B)
  
  # location level
  n.loc = params[params$param == "n.loc", "value"]
  if(n.loc == 0){
    n.loc = length(unique(obs$geo)) # number of locations/prov (based on obs data)
  }
  sigma_m  = params[params$param == "sigma_m", "value"] 
  sigma_m1 = params[params$param == "sigma_m1", "value"] 
  sigma_b  = params[params$param == "sigma_b", "value"] 
  
  m  = sapply(M, rnorm,  n = n.loc, sd = sigma_m) |> t()
  m1 = sapply(M1, rnorm, n = n.loc, sd = sigma_m1) |> t()
  b  = sapply(B, rnorm,  n = n.loc, sd = sigma_b) |> t()
  
  message('Simulated data:\n number of variants: ', n.var,
          '\n number of locations: ', n.loc)
  
  # construct time series
  nmin = params[params$param == "nmin", "value"] 
  nmax = params[params$param == "nmax", "value"] 
 
  # number of data points for each location/variant 
  nt = matrix(sample(nmin:nmax, size = n.var*n.loc, replace = TRUE),
              nrow = n.var)
  
  sigma = params[params$param == "sigma", "value"] 
  
  d = list()
  x = 0
  
  for(v in 1:n.var){
    for(l in 1:n.loc){
      x = x + 1
      ndata = nt[v,l]
      
      # simulated ww concentration
      int_ww   = runif(n=1, min = 0.1, max = 5)
      rand_ww  = rnorm(n = ndata, mean = 0, sd = 10)
      slope_ww = runif(n=1, min = 1, max = 3)
      w = slope_ww * c(1:ndata) + int_ww + rand_ww
      w[w<=0] = 1
      log.w = log(w)
        if(0) plot(log.w, typ = 'o')
      
      # Simulated hosp
      log.h = numeric(length = ndata)
      # stochastic component of hospital counts
      u = rnorm(n = ndata, mean = 0, sd = 1) 
      
      log.h[1] =  m[v,l] * log.w[1] + b[v,l] + sigma * u[1]
      
      for(t in 2:ndata){
        log.h[t] = m[v,l] * log.w[t] + m1[v,l] * log.w[t-1] + 
          b[v,l] + sigma * u[t]
      }
      # log.h = m[v,l] * log(w) + b[v,l] + sigma * u
      h = exp(log.h)
      
      if(0) {
        plot(h,typ='o') ; grid()
      }
      
      d[[x]] = data.frame(
        date       = ymd('2020-01-01') + 7 * c(0:(ndata-1)),
        ww.conc    = w,
        hosp.count = round(h),
        variant    = paste0("variant_", v),
        geo        = paste0("loc_", l)
      )
    }
  }
  
  # factor variables (for model construction)
  dat = bind_rows(d) |>
    mutate(variant = factor(variant),
           geo = factor(geo))
  
  if(0){
    g.dat = dat |> ggplot(aes(x=ww.conc, y = hosp.count)) + 
      geom_abline(intercept = beta, slope = mu) + 
      geom_smooth(method = 'lm', formula = 'y~x')+
      geom_point() + facet_grid(variant ~ geo, scales = 'fixed')
    g.dat+ scale_x_log10() + scale_y_log10()
  }

  ## ---- Construct dataframe of true values (for posterior analysis) ----
  df.scl = data.frame(
    mu   = mu,
    mu1  = mu1,
    beta = beta,
    sigma_M  = sigma_M,
    sigma_M1 = sigma_M1,
    sigma_B  = sigma_B,
    sigma_m  = sigma_m,
    sigma_m1 = sigma_m1,
    sigma_b  = sigma_b
  ) |>
    pivot_longer(cols = everything()) |>
    as.data.frame()
  
  df.M  = data.frame(name = paste0('M[',1:n.var,']'),  value = M)
  df.M1 = data.frame(name = paste0('M1[',1:n.var,']'), value = M1)
  df.B  = data.frame(name = paste0('B[',1:n.var,']'),  value = B)
  
  idx = expand.grid(1:n.var, 1:n.loc)
  vec.m  = matrix(m, ncol = 1, byrow = TRUE)
  vec.m1 = matrix(m1, ncol = 1, byrow = TRUE)
  vec.b  = matrix(b, ncol = 1, byrow = TRUE)
  
  df.m = data.frame(
    name = paste0('m[', idx[,1], ', ', idx[,2], ']'),
    value = vec.m
  )
  df.m1 = data.frame(
    name = paste0('m1[', idx[,1], ', ', idx[,2], ']'),
    value = vec.m1
  )
  
  df.b = data.frame(
    name = paste0('b[', idx[,1], ', ', idx[,2], ']'),
    value = vec.b
  )
  
  df.sigma = data.frame(
    name = "sigma",
    value = sigma
  )
  
  truth = rbind(df.scl, df.M, df.M1,
                df.B, 
                df.m, df.m1,
                df.b, df.sigma) |>
    rename(truth.val = "value")
  
  # Store drawn slopes and intercepts
  draws = list(
    mu = mu, # mu and # beta are not drawn, but stored for comparison
    beta = beta,
    M = M,
    M1 = M1,
    B = B,
    m = m,
    m1 = m1, 
    b = b
  )
 message('Simulated data created.\n') 
  return(
    list(
      draws = draws,
      dat   = dat,
      truth = truth
    )
  )
}

# ---- Fetch model parameters ----
fetch_model_params <- function(){
  params = read.csv('../prm/model-params.csv', strip.white = TRUE) |>
    select(-description) |>
    pivot_wider(names_from = param)
  return(params)
}

# ---- Fetch mcmc parameters ----
fetch_mcmc_params <- function(){
  params = fetch_model_params() |>
    select(iter, burnin, chains) |>
    mutate_all(as.numeric)
  return(params)
}

# ---- Process posterior output (from MCMC algorithm of simulated data) ----
#' @param x MCMC object
#' @param gelman.psrf Gelman-Rubin potential scale factor value
process_post <- function(x, gelman.psrf) {
  
  post    = x$post
  postsum = x$postsum
  
  names.univ = c("mu", "mu1",
                 "beta", 
                 "sigma_beta", "sigma_mu",  
                 "sigma_mu1",
                 "mu_mean")
  names.variant = c('^M', '^B')
  names.variant.sigma =  c("sigma_M", "sigma_M1", "sigma_B")
  names.location = c("sigma_m", "sigma_m1", "sigma_b", "sigma")
  
  post.plot = post |>
    mutate(
      level = case_when(
        name %in% names.univ ~ "universal",
        grepl(names.variant[1], name) | 
          grepl(names.variant[2], name) | 
          grepl(names.variant[3],name) ~ "variant",
        name %in% names.variant.sigma ~ "variant",
        grepl("^m", name) | 
          grepl("^b", name) ~ "location",
        name %in% names.location ~ "location"
      )
    )
  
  df.post = reduce(list(post.plot, postsum,
                        gelman.psrf, truth), 
                   full_join, by = "name")
  return(df.post)
}

# ---- Filter for coefficients where truth is outside 95% CI or Gelman PSRF ----
# is above threshold
#' @param mcmcobj MCMC object
filter_error_truth <- function(mcmcobj, gelman.psrf, truth) {
  
  df.postsum = reduce(list(mcmcobj$postsum,
                           gelman.psrf,
                           truth), 
                      full_join, by = "name")
  
  df.error = df.postsum |>
    mutate(
      truth.err = case_when(
        truth.val < lo | truth.val > hi ~ TRUE,
        truth.val >= lo | truth.val <= hi ~ FALSE
      ),
      hi.psrf = case_when(
        psrf > 1.01 ~ TRUE,
        psrf <= 1.01 ~ FALSE
      )
    ) |>
    filter(
      truth.err == TRUE | hi.psrf == TRUE
    )
  return(list(postsum = df.postsum, error = df.error)) 
}

# ---- Print diagnostic statistics of MCMC execution with simulated data ----
#' @param df.postsum Dataframe of posterior distributions
#' @param df.error Dataframe of coefficients with high gelman PSRF or outside CI
#' @param gelman.psrf Gelman-Rubin potential scale factor value
diagnostic_simdata <- function(df.postsum, df.error, gelman.mpsrf) {
  n.coeff = nrow(df.postsum)
  n.trutherr = nrow(df.error[df.error$truth.err == TRUE, ])
  n.hipsrf = nrow(df.error[df.error$hi.psrf == TRUE, ])
  
  message(paste0(
    "Number of coeffecients: ", n.coeff,
    "\nNumber of true values outside 95% CI: ", n.trutherr,
    "\nNumber of coefficients with high Gelman shrink factor: ", n.hipsrf,
    "\nGelman multivariate shrink factor: ", round(gelman.mpsrf,4)))
}

# ---- Create model inits for MCMC execution ----
#' @param dat Dataframe to be used in MCMC run
create_mdl_inits <- function(dat){
  
  message("Creating model inits...", appendLF = FALSE)
  
  nchains = fetch_mcmc_params()[["chains"]]
  
  inits = list()
  V = length(unique(dat$variant))
  L = length(unique(dat$geo))
  
  for(i in 1:nchains){
    inits[[i]] = list(
        muu       = runif(1, -1, 3),
        mu1u      = runif(1, -1, 3),
        mu_mean   = runif(1, -1, 3),
        mu1_mean  = runif(1, -1, 3),
        betau     = runif(1, -5,5),
        beta_mean = runif(1, -5,5),
        Mu       = rep(runif(1, -1,3), V),
        M1u      = rep(runif(1, -1,3), V),
        Bu       = rep(runif(1, -5,5), V),
        sigma_M  = rexp(1, 1),
        sigma_M1 = rexp(1, 1),
        sigma_B  = rexp(1, 1),
        m_u     = matrix(runif(1,-1,3),nrow = V, ncol = L),
        m1_u    = matrix(runif(1,-1,3),nrow = V, ncol = L),
        b_u     = matrix(runif(1,-5,5),nrow = V, ncol = L),
        sigma_m    = rexp(1, 1),
        sigma_m1   = rexp(1, 1),
        sigma_b    = rexp(1, 1),
        sigma_mu   = rexp(1, 1),
        sigma_mu1  = rexp(1, 1),
        sigma_beta = rexp(1, 1),
        sigma      = rexp(1, 1/20)
    )
  }

  message("complete.\n")
  return(inits)
}

# Get final start and end dates (for ms)
#' @param df Dataframe to be used in model
get_mdl_dates <- function(df){
  dates = df |>
    group_by(city, prov, variant) |>
    summarise(
      start = first(date),
      end = last(date),
      .groups = 'drop'
    )
  saveRDS(dates, file = '../out/mdl-dates.rds')
  return(dates)
}

# Create summary table (for ms)
#' @param df Dataframe of model output
create_sum_table <- function(df){
  # retrieve variant idx table
  var.idx = read.csv('../out/idx_variant.csv')
  
  # merge tables and mutate values
  tmp = df |>
    filter(grepl("^B|^M", name) | name %in% c("mu", "beta")) |>
    mutate(variant.idx = as.numeric(str_extract_all(name, "\\d+"))) |>
    left_join(var.idx, by = 'variant.idx') |>
    mutate(variant = case_when(
      name %in% c("mu", "beta") ~ "Universal",
      !name %in% c("mu", "beta") ~ variant
    ),
    value = paste0(round(m, 3), " [", round(lo, 3), "; ", round(hi,3), "]"))
  
  # rearrange and polish table
  tbl = tmp |>
    select(name, variant, value) |>
    mutate(name = case_when(
      grepl("^M", name) | name == "mu" ~ "m",
      grepl("^B", name) | name == "beta" ~ "b"
    )) |>
    pivot_wider(names_from = name, values_from = value) |>
    rename(Variant = "variant") |>
    relocate(Variant, m, b)
  
  colnames(tbl)[1:3] <- c(
    '\\textbf{Variant}',
    '\\textbf{slope $m$}', 
    '\\textbf{intercept $b$}')
  
  return(tbl)
}

# ---- Sensitivity functions ---- 
#' @description Create of model object for sensitivity analysis
#' @param thresholds Vector of variant thresholds
#' @param windows Vector of time windows
#' @param pct.above Vector of minimum percentage of days above threshold
#' 
create_sens_mdl_obj <- function(thresholds,
                                windows,
                                pct.above){
  message("Creating model objects for sensitivity analysis...\n\n")
  # Retrieve data
  wwhosp = readRDS('../out/df.merged.rds')
  domvars = readRDS('../out/domvars.rds')
  
  mdl.obj = list()
  x = 1
  for(t in thresholds){
    for(w in windows){
      for(p in pct.above){
        message(paste0("Creating model object for:\nthreshold: ", t,
                       "\nwindow: ", w,
                       "\npct.above: ", p))
        dates = get_var_dates(domvars,
                              sensitivity = TRUE,
                              var.prop = t,
                              var.window = w,
                              var.pctabove = p)
        mdl.obj[[x]] = create.mdl.obj(df = wwhosp,
                                      dates = dates,
                                      sensitivity = TRUE)
        
        x = x + 1
      }
    }
  }
  message("Model objects complete.")
  return(mdl.obj)
}

#' @description Create summary of sensitivity analysis
#' @param mcmcobj Output of MCMC execution
create_sens_sum_table <- function(mcmcobj){
  params = mcmcobj[["params"]]
  df = mcmcobj[["postsum"]] |>
    filter(name %in% c("mu", "beta")) |>
    mutate(
      value = paste0(round(m, 3), " [", round(lo, 3), "; ", round(hi,3), "]"),
      threshold = (params$threshold - 0.6) * 100,
      window = params$window - 8,
      pct.above = (params$pct.above - 0.75) * 100
    ) |>
    mutate_at(
      c("threshold", "window", "pct.above"),
      sign_mutate
    ) |>
    select(name, value, threshold, window, pct.above) |>
    pivot_wider(names_from = name, values_from = value) |>
    mutate_at(
      c("threshold", "window", "pct.above"),
      reference_mutate
    ) |>
    relocate(threshold, window, pct.above, mu, beta)
    
    colnames(df)[1:5] <- c(
      '\\textbf{Threshold (\\%)}',
      '\\textbf{Window (weeks)}',
      '\\textbf{\\% above}',
      '\\textbf{slope $m$}', 
      '\\textbf{intercept $b$}')
  return(df)
}

#' @description Adding + or - to values for sensitivity table
#' @param x Numeric column
sign_mutate <- function(x){
  y = case_when(
    sign(x) == 1 ~ paste0("+", x),
    sign(x) == -1 | sign(x) == 0 ~ as.character(x)
  )
  return(y)
}

#' @description
#' Renaming values to Reference if they are the values used in ms
#' @param x Numeric column
reference_mutate <- function(x){
  y = case_when(
    x == "0" ~ "Reference",
    x != "0" ~ x
  )
  return(y)
}

# Translate variant index to variant name
#' @param d Dataframe of variant index with their corresponding name.
helper_variantname <- function(d) {
  
  v = read.csv('../out/idx_variant.csv')
  
  res = d |>
    mutate(variant.idx = str_extract(name, '\\[\\d+\\]') |> 
             str_remove_all(pattern = '\\W') |> 
             as.integer()) |>
    left_join(v, by = 'variant.idx') |> 
    mutate(
      name2 = str_extract(name, '^\\w+'),
      name3 = paste(name2,variant, sep = '__') |> 
        str_remove('__NA$')) 
  return(res)
}

# Retrieve sublineage merging information (for ms appendix)
retrieve.sublineage <- function(){
  # Retrieve pango lineage information
  pango = readRDS('../out/pango-lineages.rds')
  
  # Retrieve variants included in model
  vars = unique(readRDS('../out/var-analysis.rds')[['df.dates']][['variant']])
  
  df = pango |>
    filter(lineage %in% vars) |>
    rename(variant = 'lineage',
           sublineage = 'variant',
           name = 'who') |>
    relocate(variant, .after = sublineage)
  
  return(df)
  
}
