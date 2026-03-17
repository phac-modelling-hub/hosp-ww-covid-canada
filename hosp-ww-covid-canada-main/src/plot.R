source('lm.R')

# ---- Plot functions ----
# Plotting the dominant periods for each location
#' @param dates.model Dataframe containing start and end dates of dominant
#'  variants in each province.
plot_domvariant_periods <- function(dates.model) {
  g = dates.model |> 
    ggplot() + 
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) + 
    geom_segment(aes(x = start, xend = end, y=1, yend=1, 
                     colour = variant), 
                 linewidth = 3) + 
    guides(colour = 'none', fill = 'none')+
    labs(title = 'Dominant variants periods', y='', x='')
  
  g1 = g + facet_grid(prov~ variant, scales = 'free_x')
  g2 = g + 
    facet_wrap(~prov, ncol = 1) +
    geom_label(aes(x = start, y=2, label = variant, fill = variant),
               hjust = 0, vjust = 1, size=2, fontface = 'bold')
  return(list(g1,g2))
}

# Plotting a variable posterior distribution
#' @param dat Dataframe of posteriors.
#' @param thename String. Variable name.
plot_single_post <- function(dat, thename){
  
  message('plotting posterior: ', thename)
  
  if(0){  #DEBUG
    thename = 'M[1]' 
  }
  
  datplot = dat |>
    filter(name == thename, 
           iter%%5 == 0) # thin dataframe to speed up `geom_density`
  
  truth.val = datplot$truth.val[1]
  lo = datplot$lo[1]
  hi = datplot$hi[1]
  m  = datplot$m[1]
  psrf = unique(datplot$psrf)
  
  g = ggplot(datplot, aes(x = value)) +
    geom_density(fill = 'azure2') +
    geom_vline(xintercept = truth.val,
               linewidth = 1.3, color = 'red2', alpha = 0.8) +
    annotate(geom='segment', x = lo, xend = hi, y = 0, yend = 0, linewidth = 2) +
    annotate(geom = 'point', x = m, y = 0, size=3, shape = 21, 
             fill = 'white', stroke = 2)+
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) +
    labs(title = thename,x='',y='',
         subtitle = paste("Gelman shrink factor:", round(psrf, 4)))
  g
}

# Plotting posterior distributions of variables at a given hierarchical level
#' @param df.post Dataframe of posteriors
#' @param level String. Hierarchical level.
plot_post_level <- function(df.post, level){
  
  dat = df.post |>
    filter(level == !!level)
  
  names = unique(dat$name)
  
  if(level == 'universal') names = c('beta', 'mu', 'mu1')
  
  g = lapply(names, plot_single_post, dat = dat)
  
  gg = wrap_plots(g, guides = 'collect') + 
    plot_annotation(title = paste("Posterior distribution on simulated data:", 
                                  level, "level"))
  
  return(gg)
}

# Plotting 2D correlations between different posterior distributions
#' @param post Dataframe of posteriors of different variables
#' @param pattern String. Regex pattern of variables to include in plotting.
plot_post_2d <- function(post, pattern) {
  
  d = filter(post, grepl(pattern, name)) |> 
    pivot_wider(names_from = name, values_from = value) |> 
    filter(iter %% 5 == 0) |>
    select(-iter, -chain) 
  nc = ncol(d)
  tmp = tmpcor = list() ; k=1
  for(i in 1:(nc-1)){
    message('plot post correl ',i,'/',nc-1)
    for(j in (i+1):nc){
      tmpdf = d[,c(i,j)]
      tmpcor[[k]] = data.frame(
        pair  = paste(names(tmpdf), collapse = '_'),
        correl = cor(tmpdf)[1,2]
      )
      names(tmpdf) = c('x','y')
      tmp[[k]] =  ggplot(tmpdf, aes(x=x,y=y))+
        geom_density_2d_filled() + guides(fill = 'none')+
        labs(title= paste(names(d)[c(i,j)], collapse = ' ; '),
             x='',y='')+
        theme(panel.grid = element_blank())
      tmp[[k]]
      k = k+1
    }
  }
  correls = bind_rows(tmpcor)
  g.cor = correls |> ggplot(aes(x=pair, y=correl))+
    geom_col() + 
    theme(axis.text.x = element_text(angle=30, hjust=1, size = rel(0.7)))+
    labs(x='', y = '', title =  paste('Posterior correlations', pattern))
  # g.cor
  g = wrap_plots(tmp)
  return(list(density=g, correl = g.cor))
}

# Plotting density plots of posterior universal variables against their 
#  associated sigma
#' @param post Dataframe of posteriors.
plot_universal_sigma_post <- function(post) {
  
  q.b = post |> filter(grepl('beta',name)) |>
    pivot_wider(names_from = name, values_from = value)
  
  q.m = post |> filter(grepl('mu',name)) |>
    pivot_wider(names_from = name, values_from = value)
  
  gb = gm = list() 
  
  gb[[1]] = q.b |> 
    ggplot(aes(x=beta, y=sigma_beta)) + geom_point(alpha = 0.1)
  gb[[2]] = ggplot(q.b, aes(x=sigma_beta)) + geom_density()
  gb[[3]] = ggplot(q.b, aes(x=beta)) + geom_density()
  
  gm[[1]] = q.m |> 
    ggplot(aes(x=mu, y=sigma_mu)) + geom_point(alpha = 0.1)
  gm[[2]] = ggplot(q.m, aes(x=sigma_mu)) + geom_density()
  gm[[3]] = ggplot(q.m, aes(x=mu)) + geom_density()
  
  
  ggb = wrap_plots(gb)
  ggm = wrap_plots(gm)
  
  return(ggb/ggm)
}

# Plot simulated data
#' @param dat Dataframe of simulated data.
plot_simulated_data <- function(dat) {
  g = dat |> 
    ggplot(aes(x=ww.conc, y = hosp.count)) + 
    geom_point(color = 'steelblue4', shape = 1) + 
    facet_grid(variant ~ geo, scales = 'free') + 
    theme(panel.grid.minor = element_blank())+
    scale_x_log10() + scale_y_log10()+
    labs(title = 'Simulated data',
         x = 'concentration in wastewater',
         y = 'hospital admissions') 
  g
  
  g.ts.w = dat |> 
    ggplot(aes(x = date, y = ww.conc))+
    facet_grid(variant ~ geo)+
    geom_line(aes(y = w_1), color = 'grey')  + 
    geom_line() + geom_point() +
    labs(title = 'WW time series')
  g.ts.w
  
  g.ts.h = dat |> 
    ggplot(aes(x = date, y = hosp.count))+
    facet_grid(variant ~ geo)+
    geom_line(aes(y = h_1), color = 'grey')  + 
    geom_line() + geom_point() +
    labs(title = 'Hospitalization time series')
  g.ts.h
  
  return(list(
    xy = g, 
    ts.w = g.ts.w,
    ts.h = g.ts.h
  ))
}

# Plotting observed data prior to transformation
#' @param obs Dataframe of observed data
#' @param data.source String. Source of hospital admissions data.
#' @param logscale Logical. Flag to observe data on log-scale.
plot_obs <- function(obs, data.source, logscale) {
  
  dplot = obs |> 
    filter(source == data.source) 
  
  g.wh = dplot |>
    ggplot(aes(x = ww.conc, y = hosp.count)) + 
    geom_smooth(formula = 'y~x', method = 'lm')+
    geom_point(shape = 21 , alpha = 0.6, color = 'red3')+
    facet_grid(variant ~ geo, scales = 'free')+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = 'grey97')) +
    labs(title = 'Real observations', 
         subtitle = paste('source hosp data:', data.source))
  
  # Time series
  
  g.ts.h = dplot |> 
    ggplot(aes(x=date, y = hosp.count))+
    facet_grid(variant ~ geo, scales = 'free') + 
    geom_line(linewidth = 0.8, colour = 'seagreen') + 
    labs(title= 'Time series of hospital admissions', x='')
  g.ts.h
  
  g.ts.w = dplot |> 
    ggplot(aes(x=date, y = ww.conc))+
    facet_grid(variant ~ geo, scales = 'free') + 
    geom_line(linewidth = 0.8, colour = 'chocolate3') + 
    labs(title= 'Time series of wastewater concentrations', x='')
  g.ts.w
  
  if(logscale) {
    g.wh = g.wh + scale_x_log10()+ scale_y_log10()
    g.ts.h = g.ts.h + scale_y_log10()
    g.ts.w = g.ts.w + scale_y_log10()
  }
  
  return(list(
    wh = g.wh,
    ts.h = g.ts.h, 
    ts.w = g.ts.w
  ))
}

# Plotting slopes of naive independent linear regressions
#' @param dflm Dataframe of naive linear regression outputs
#' @param xx String. Column to be used as x-axis of plot.
#' @param ff String. Column to be used for the facetting of plot.
hlp_plotslope <- function(dflm, xx, ff) {
  g = dflm|> 
    ggplot(aes(x = .data[[xx]])) + 
    geom_segment(aes(xend =  .data[[xx]], y = slope.lo, yend = slope.hi),
                 linewidth = 2) + 
    geom_point(aes( y = slope.m), shape = 21, fill='white') + 
    facet_grid(~ .data[[ff]]) + 
    labs(title = 'Slope indep. linear regressions',
         x='', y = 'slope') 
  g
}

# Plotting intercepts of naive independent linear regressions
#' @param dflm Dataframe of naive linear regression outputs
#' @param xx String. Column to be used as x-axis of plot.
#' @param ff String. Column to be used for the facetting of plot.
hlp_plotint <- function(dflm, xx, ff) {
  g = dflm|> 
    ggplot(aes(x = .data[[xx]],)) + 
    geom_segment(aes(xend = .data[[xx]], 
                     y = intercept.lo, yend = intercept.hi),
                 linewidth = 2) + 
    geom_point(aes( y = intercept.m), shape = 21, fill='white') + 
    facet_grid( ~ .data[[ff]], scales = 'free_y' )+
    labs(title = 'intercept indep. linear regressions',
         x='', y = 'intercept')
  g
}



#' Plot linear regressions' slope and intercept estimates
#' @param dat Dataframe to be used in linear regressions.
plot_lm_data <- function(dat) {
  
  dflm = lm_prov_var(dat)
  
  g.lm.slope  = hlp_plotslope(dflm, xx = 'prov', ff='variant')+
    theme(axis.text.x = element_text( size = rel(0.8)))
  
  g.lm.slope2 = hlp_plotslope(dflm, xx = 'variant', ff='prov')+
    theme(axis.text.x = element_text(angle=45, hjust = 1, size = rel(0.8)))
  
  g.lm.int  = hlp_plotint(dflm, xx='prov', ff='variant')+
    theme(axis.text.x = element_text( size = rel(0.8)))
  
  g.lm.int2 = hlp_plotint(dflm, xx='variant', ff='prov')+
    theme(axis.text.x = element_text(angle=45, hjust = 1, size = rel(0.8)))
  
  return(list(
    lm.slope  = g.lm.slope,
    lm.slope2 = g.lm.slope2,
    lm.int    = g.lm.int,
    lm.int2   = g.lm.int2
  ))
}

# Plotting observed data after transformations
#' @param dat Dataframe of observed data to be used in model.
#' @param ms Logical. Flag if plot object is to be used in manuscript.
#' @param size Numeric. Size of text in plot.
#' @param show.lm Logical. Flag to show linear regression in each panel of plot.
plot_data <- function(dat, ms = FALSE, size = 8,
                      show.lm = FALSE) {
  pt.size = ifelse(ms, 3, 2)
  if(size > 20) pt.size = 4
  pt.alpha = ifelse(ms, 0.8, 0.4)
  
  g = dat |> ggplot(aes(x = w, y=h))+
    facet_grid(variant ~ city, scales = 'free') + 
    geom_point(size=pt.size, alpha = pt.alpha, shape = 1)+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = 'grey97'))+
    labs(title = 'Data used in inference')
  
  if(!ms){
    if(show.lm) g = g + geom_smooth(method = 'lm', formula = 'y~x')
    
    g
  }
  if(ms){
    g = g +
      theme(text = element_text(size = size),
            strip.background = element_rect(fill ='steelblue4' ), 
            strip.text = element_text(color= 'white', face = 'bold')) +
      labs(title = "",
           x = "log viral wastewater concentration (cp/mL)",
           y = "log hosp. admissions per 100,000")
  if(size > 20)
    g = g + theme(axis.text = element_text(size = rel(0.5)),
                  axis.title = element_text(size = rel(1.5)))
  # g  
  }
  return(g)
}


plot_data_tally <- function(dat) {
  g = dat |>
    group_by(variant, geo) |> 
    tally() |>
    ggplot(aes(x=geo, y = n))+
    geom_col()+
    facet_wrap(~variant) +
    theme(axis.text.x = element_text(angle = 45, hjust=1, size = rel(1)),
          strip.background = element_rect(fill = 'steelblue4'),
          strip.text = element_text(color = 'white', face = 'bold'),
          panel.grid.minor.y = element_blank(),
          panel.grid.major = element_line(color = 'grey97'))+
    labs(title = 'Number of data points', y = 'count', x = 'location')
  g
}


# Plotting the traces of the MCMC execution
#' @param mcmcobj MCMC object
#' @param iter.displayed Numeric. Number of iterations to display.
#' @param subtitle Logical. Flag to display total # of iterations and # displayed
plot_trace <- function(mcmcobj, iter.displayed = 1e3, subtitle = TRUE) {
  
  d = mcmcobj$post
  # Reduce the size of the plot 
  # (else very slow to display)
  n = max(d$iter)
  m = iter.displayed # The smaller, the less data displayed
  if(n > m){
    q = as.integer(n / m)
    d = filter(d, iter %% q == 0)
  }
  
  subt = NULL
  if(subtitle) subt = paste('Total MCMC iterations:', n, 
                                ' (total displayed: ',m,')')
  
  g = helper_variantname(d) |>
    ggplot(aes (x= iter, y = value, color = chain))+
    geom_line(alpha = 0.5) + 
    facet_wrap(~name3, scales = 'free_y')+
    theme(panel.grid = element_blank(), 
          axis.text.x = element_text(size = rel(0.6), angle = 45, hjust=1)) + 
    scale_color_brewer(palette = 'Dark2') + 
    labs(title = 'Trace plots', 
         subtitle = subt, 
         x = 'iteration', y = 'value')
  g
}

# Plotting the Gelman Rubin statistic for each posterior variable
#' @param gelman.psrf Object containing Gelman-Rubin statistic.
plot_gelman <- function(gelman.psrf) {
  g = gelman.psrf |>
    helper_variantname() |>
    ggplot(aes(x = name3))+
    geom_segment(aes(y = psrf, yend = psrf.upper, xend=name3),
                 arrow = arrow(angle=90, length = unit(0.01, "npc"))) + 
    geom_point(aes(y=psrf)) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = rel(0.9)),
          panel.grid.minor = element_blank()) + 
    labs(title = 'Gelman-Rubin statistic', x='parameter',
         y = 'potential scale reduction factor (PSRF)')
  g
}

# Plotting all posterior distributions from the MCMC execution
#' @param mcmcobj MCMC object output.
plot_post_all <- function(mcmcobj) {
  g.post = mcmcobj$post |> 
    ggplot(aes(x = value)) +
    facet_wrap(~name, scales = 'free') + 
    geom_density(fill = 'steelblue3', alpha = 0.3) + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(title = 'Posterior distributions')
  g.post
}

# Plotting posterior distribution of variable from MCMC execution
#' @param mcmcobj MCMC object
#' @param pattern.prm String. Regex pattern of variable to be displayed.
#' @param subtitle Subtitle to include in plot.
plot_post <- function(mcmcobj, pattern.prm, subtitle = 'pattern') {
  
  if(subtitle == 'pattern') subtitle = pattern.prm
  if(subtitle == 'none') subtitle = NULL
  
  v = read.csv('../out/idx_variant.csv')
  
  d = mcmcobj$post |>
    slice(seq(1, n(), by = 50)) |>
    helper_variantname()    
  
  g = d |> 
    filter(grepl(pattern.prm, name)) |> 
    ggplot(aes(x = value)) +
    facet_wrap(name3 ~ ., scales = 'free') + 
    geom_density(fill = 'steelblue3') + 
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(title = 'Posterior distributions', 
         subtitle = subtitle,
         x='', y = '')
  g
}

# Comparison with naive linear regression
#' @param dflm Dataframe of independent linear regressions.
#' @param mcmcobj MCMC object
#' @param dat Dataframe of observed data used in model.
#' @param ms Logical. Flag if plot to be used in manuscript.
plot_compare_lm_mcmc <- function(dflm, mcmcobj, dat, ms = FALSE) {
  
  # Link to variant names
  vv = select(dat, starts_with('variant')) |> 
    distinct() |> 
    mutate(name = paste0('M[',variant.idx,']'),
           name2 = paste0('B[',variant.idx,']'))
  
  # Posteriors
  ps.slope = mcmcobj$postsum |> filter(grepl('M|^mu$',name)) |> 
    left_join(vv, by = 'name') |>
    filter(!name %in% c("mu", "sigma_M")) |>
    mutate(name.plot = paste(name,variant))
  
  ps.int = mcmcobj$postsum |>
    rename(name2 = name) |>
    filter(grepl('^B',name2)) |> 
    left_join(vv, by = 'name2') |>
    mutate(name.plot = paste(name2,variant))
  
  
  df = left_join(dflm, vv, by='variant') |>
    mutate(name = paste0('M[',variant.idx,']'),
           name.plot = paste(name, variant),
           cityshort = substr(city,1,3)) 
  
  mu.m = as.numeric(mcmcobj$postsum[mcmcobj$postsum$name == "mu", "m"])
  beta.m = as.numeric(mcmcobj$postsum[mcmcobj$postsum$name == "beta", "m"])
  
  col.mcmc = 'indianred'
  x.nudge = 0.19
  
  if(!ms){
    gg = ggplot(data = df, aes(x=name.plot))
  }
  
  if(ms){
    gg = ggplot(data = df)
  }
  
  g.slope = gg +
    geom_hline(yintercept = 0, color = 'grey80') + 
    #
    # -- simple linear regression
    geom_segment(aes(x=variant, y= slope.lo, yend=slope.hi), alpha = 0.6)+
    geom_point(aes(x=variant, y=slope.m, size = n),
               shape = 21, fill='white', alpha = 0.6) +
    # -- hierarchical bayesian model
    geom_hline(yintercept = mu.m, color = col.mcmc, linetype = 'dashed')+
    geom_segment(data = ps.slope, 
                 aes(x=variant, xend=variant, y=lo, yend=hi),
                 position = position_nudge(x=x.nudge),
                 color = col.mcmc, alpha = 0.5,
                 linewidth = 2) + 
    geom_point(data = ps.slope, aes(x=variant, y=m), 
               position = position_nudge(x=x.nudge),
               color = col.mcmc, size = 2, stroke = 1,
               shape = 15) +
    coord_cartesian(ylim = c(-0.6,1.5)) +  
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = 'grey97'),
          axis.text.x = element_text(angle = 0, hjust = 0.5)) + 
    labs(title = 'Regression slope (parameter m)',
         x = '', y = '')
  
  g.int = gg +
    geom_hline(yintercept = 0, color = 'grey80') + 
    #
    # -- simple linear regression
    geom_segment(aes(x=variant, y= intercept.lo, yend=intercept.hi), 
                 alpha = 0.6)+
    geom_point(aes(x=variant, y=intercept.m, size = n),
               shape = 21, fill='white', alpha = 0.6) +
    #
    # -- hierarchical bayesian model
    geom_hline(yintercept = beta.m, color = col.mcmc, linetype = 'dashed')+
    geom_segment(data = ps.int, 
                 aes(x=variant, xend=variant, y=lo, yend=hi),
                 position = position_nudge(x=x.nudge),
                 color = col.mcmc, alpha = 0.5,
                 linewidth = 2) + 
    geom_point(data = ps.int, aes(x=variant, y=m), 
               position = position_nudge(x=x.nudge),
               color = col.mcmc, size = 2, stroke = 1,
               shape = 15) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = 'grey97'),
          axis.text.x = element_text(angle = 0, hjust = 0.5)) + 
    labs(title = 'Regression intercept (parameter b)',
         x = '', y = '')
  
  g = wrap_plots(g.slope, g.int, ncol = 1, guides = 'collect') + 
    plot_annotation(title = "Comparison inference Linear Model v.s. Hamiltonian MCMC")

  return(g)
}

plot_leaveoneout_gelman <- function(gelman.loo) {
  g = gelman.loo |> 
    ggplot(aes(x = vl.left, y=value))+
    facet_wrap(~name, ncol = 2, scales = 'free_y')+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.7)),
          panel.grid.minor = element_blank())+
    geom_hline(yintercept = 1.02, color = 'indianred', linetype = 'dashed')+
    geom_point()+
    labs(title = 'Gelman-Rubin PSRF', 
         x = 'variant / location left out',
         y = 'PSRF')
  return(g)
}

plot_leaveoneout_post <- function(post.loo){
  col = 'steelblue'
  
  g = post.loo |> ggplot(aes(x=vl.left)) +
    geom_segment(aes(xend = vl.left, y = lo, yend = hi), 
                 linewidth = 2, alpha = 0.2, color = col) + 
    geom_point(aes(y = m), color = col) + 
    theme(axis.text.x = element_text(angle = 30, hjust = 1, size = rel(0.7)),
          panel.grid.minor = element_blank())+
    facet_wrap(~name, scales = 'free_y')+
    labs(title = 'Posterior values', 
         x = 'variant / location left out',
         y = 'PSRF')
  return(g)
}
