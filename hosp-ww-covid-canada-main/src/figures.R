## Functions to generate figures for ms ##
library(tidyverse)
library(patchwork)
library(kableExtra)

source('utils.R')
source('model-hmc.R')
source('plot.R')

#' Create synthetic data to overlay to the `main_result` figure.
#' @param mu Numeric. Value of slope parameter
#' @param beta Numeric. Value of intercept parameter
#' @param wmin Numeric. Minimum wastewater concentration
#' @param wmax Numeric. Maximum wastewater concentration
#' @param n Numeric. Number of datapoints.
#' @param sd Numeric. Standard deviation when estimating hosp admissions.
synthetic_data <- function(mu, beta, wmin, wmax, n, sd ) {
  w = runif(n = n, min = wmin, max=wmax)
  h = ff(w = w, mu = mu, beta = beta)*rnorm(n=n, mean=1, sd=sd)
  res = data.frame(w = w, h = h)
}

# Function that translates the linear equation
# of the log variables into "natural" scale, ie,
# log(h) = m * log(w) + b
# -->  h = exp(b) * W^m
#' @param w Numeric. Value of wastewater concentration.
#' @param mu Numeric. Value of slope parameter.
#' @param beta Numeric. Value of intercept parameter.
ff = function(w, mu, beta) return(exp(beta)*w^mu)

#' Main result figure.
#' Show categorization of future variant waves
#'  based on inference on historical data.
#' @param file.mcmc String. Filepath of MCMC object
#' @param add.synth.data Logical. Flag to include hypothetical observations.   
fig_main_result <- function(file.mcmc, add.synth.data,
                            xlim = c(0,1000),
                            ylim = c(0,50),
                            legend.version = FALSE) {
  
  a = readRDS(file.mcmc)  
  
  # Retrieve the posterior samples for the 
  # universal level slope and intercept
  df = a$post |> 
    filter(name %in% c('mu', 'beta')) |>
    mutate(idx= paste(chain, iter)) |>
    select(name, value, idx) |> 
    pivot_wider(names_from = 'name', values_from = 'value')
  
  # We only need a small sub-sample
  subsample.size = 100
  dfs = df[sample(1:nrow(df), size = subsample.size, replace = F),]
  
  # wastewater concentration values for the plot
  w = c(0, 2, 5, seq(10, 90, by = 10), seq(100,1000, by=100))
  w
  
  # Loop over all posterior (sub)samples for the
  # slope and intercept (mu, beta) and calculate
  # corresponding hosp admission (h) for a given 
  # wastewater concentration (w):
  tmp = list()
  for(i in 1:nrow(dfs)){
    mu   = dfs$mu[i]
    beta = dfs$beta[i] 
    tmp[[i]] = data.frame(
      h = ff(w, mu = mu, beta = beta),
      w = w,
      mu = mu, beta = beta, 
      prmset = factor(i)
    )
  }
  dplot = bind_rows(tmp) 
  
  # Credible intervals that defines
  # the "severity classes"
  ci = c(0.5, 0.8, 0.95)
  
  cibrks = c(0.5, c(1 - ci, 1 + ci) / 2 )
  cibrks = sort(cibrks) * 100
  
  # Calculate quantiles for severity classes
  dsum = dplot |> 
    group_by(w) |>
    summarise(m = mean(h),
              lo = quantile(h, probs = 0.5 - ci[1]/2), 
              hi = quantile(h, probs = 0.5 + ci[1]/2), 
              vlo = quantile(h, probs = 0.5 - ci[2]/2), 
              vhi = quantile(h, probs = 0.5 + ci[2]/2), 
              vvlo = quantile(h, probs = 0.5 - ci[3]/2), 
              vvhi = quantile(h, probs = 0.5 + ci[3]/2), 
    )
  
  # Plot
  
  alpha = 0.7
  
  if(legend.version){  # asked by reviewer
    
    dsum2 = dplot |> 
      group_by(w) |>
      summarise(
        a_lo = quantile(h, probs = 0.5 - ci[3]/2), 
        a_hi = quantile(h, probs = 0.5 - ci[2]/2), 
        b_lo = quantile(h, probs = 0.5 - ci[2]/2), 
        b_hi = quantile(h, probs = 0.5 - ci[1]/2), 
        c_lo = quantile(h, probs = 0.5 - ci[1]/2), 
        c_hi = median(h),
        d_lo = median(h),
        d_hi = quantile(h, probs = 0.5 + ci[1]/2), 
        e_lo = quantile(h, probs = 0.5 + ci[1]/2), 
        e_hi = quantile(h, probs = 0.5 + ci[2]/2), 
        f_lo = quantile(h, probs = 0.5 + ci[2]/2), 
        f_hi = quantile(h, probs = 0.5 + ci[3]/2), 
      ) |> 
      pivot_longer(-w) |> 
      separate(col = 'name', into = c('cat','lvl'), sep = '_' ) |>
      pivot_wider(names_from = 'lvl', values_from = 'value') |> 
      mutate(quantiles = case_when(
        cat == 'a' ~ paste0(cibrks[1],'-',cibrks[2],'%'),
        cat == 'b' ~ paste0(cibrks[2],'-',cibrks[3],'%'),
        cat == 'c' ~ paste0(cibrks[3],'-',cibrks[4],'%'),
        cat == 'd' ~ paste0(cibrks[4],'-',cibrks[5],'%'),
        cat == 'e' ~ paste0(cibrks[5],'-',cibrks[6],'%'),
        cat == 'f' ~ paste0(cibrks[6],'-',cibrks[7],'%'),
        
      ))
    
    unique(dsum2$quantiles)
    ql = c("2.5-10%" , 
           "10-25%" ,  
           "25-50%"  , 
           "50-75%" ,  
           "75-90%" ,  
           "90-97.5%")
    dsum2$quantiles = factor(dsum2$quantiles, 
                             levels = ql )
    
    g = dsum2 |> 
      ggplot(aes(x=w)) +
      geom_ribbon(aes(fill = quantiles, ymin=lo, ymax=hi), alpha = 0.7) + 
      scale_fill_brewer(palette= 'BuPu')+
      theme_bw() + 
      coord_cartesian(xlim = xlim, ylim = ylim) + 
      theme(panel.grid.minor = element_blank()) + 
      labs(x = 'viral concentration in wastewater\n(gene copies per mL)', 
           y = 'hosp. adm. per 100,000')
  }
  
  
  if(!legend.version){
    
    # color palette (from https://colorbrewer2.org/#type=sequential&scheme=BuPu&n=6)
    cols = c('#edf8fb','#bfd3e6','#9ebcda',
             '#8c96c6','#8856a7','#810f7c')
    
    g = ggplot(dsum, aes(x=w)) + 
      geom_line(aes(y = m), alpha=alpha)+
      geom_line(aes(y = vvlo), colour = 'grey80')+
      geom_ribbon(aes(ymin=vhi, ymax=vvhi), alpha = alpha, fill = cols[6])+
      geom_ribbon(aes(ymin=hi, ymax=vhi), alpha = alpha, fill =cols[5]) + 
      geom_ribbon(aes(ymin=m, ymax=hi), alpha = alpha, fill = cols[4]) + 
      geom_ribbon(aes(ymin=lo, ymax=m), alpha = alpha, fill = cols[3]) + 
      geom_ribbon(aes(ymin=vlo, ymax=lo), alpha = alpha, fill = cols[2]) + 
      geom_ribbon(aes(ymin=vvlo, ymax=vlo), alpha = alpha, fill = cols[1]) + 
      theme_bw() + 
      coord_cartesian(xlim = xlim, ylim = ylim) + 
      theme(panel.grid.minor = element_blank()) + 
      labs(x = 'viral concentration in wastewater\n(gene copies per mL)', 
           y = 'hosp. adm. per 100,000')
  }
  
  res  = g
  
  # Add hypothetical observations
  if(add.synth.data){
    df.A = synthetic_data(mu = 0.45, beta = 0.40, wmin = 10, wmax = 250, n = 12, sd = 0.05)
    df.B = synthetic_data(mu = 0.30, beta = 0.35, wmin = 10, wmax = 550, n = 12, sd = 0.09)
    res = g + 
      geom_point(data = df.A, aes(x=w, y=h), shape = 1, color = 'grey20')+ 
      geom_point(data = df.B, aes(x=w, y=h), shape = 4, color = 'grey20')
  }
  return(res)
}

# Generate appendix figures
figure_appendix_mcmcoutputs <- function() {
  
  fname =  '../out/mcmcobj.rds'
  mcmcobj = readRDS(fname)
  
  # Traces
  
  g.trace = plot_trace(mcmcobj, iter.displayed = 1e3, subtitle = FALSE)
  pdf('../ms/figs/appendix-mcmcout-trace.pdf', width = 9, height = 6)
  plot(g.trace)
  dev.off()
  
  # Gelman-Rubin statistics
  
  gelman_diag  = run_gelman(mcmcobj$mcmc)
  gelman.mpsrf = gelman_diag[["gelman.mpsrf"]]
  gelman.psrf  = gelman_diag[["gelman.psrf"]]
  g.gelman = plot_gelman(gelman.psrf)
  pdf('../ms/figs/appendix-mcmcout-gelman.pdf', 
      width = 7, height = 4)
  plot(g.gelman)
  dev.off() 
  ggsave('../ms/figs/appendix-mcmcout-gelman.png',
         g.gelman, width = 14, height = 8)
  
  # Posterior distributions
  post.u = plot_post(mcmcobj, pattern.prm = 'beta|mu', subtitle = 'universal level')
  post.v = plot_post(mcmcobj, pattern.prm = 'B|M', subtitle = 'variant level')
  
  pdf('../ms/figs/appendix-mcmcout-postu.pdf', width=5, height=3)
  plot(post.u)
  dev.off()
  
  pdf('../ms/figs/appendix-mcmcout-postv.pdf', width=9, height=6)
  plot(post.v)
  dev.off()
  
  # LM vs MCMC
  data.source    = 'DAD'
  logscale       = TRUE
  variant.ignore = c('B.1.438.1') 
  
  dat = prepare_data(data.source, logscale, variant.ignore) 
  dflm = lm_prov_var(dat)
  g.compare.lm.mcmc.ms = plot_compare_lm_mcmc(dflm, mcmcobj, dat, ms = TRUE)
  path.gcomp = '../ms/figs/lm-vs-mcmc'
  pdf(paste0(path.gcomp, '.pdf'), height = 7, width = 9)
  plot(g.compare.lm.mcmc.ms)
  dev.off()
  ggsave(paste0(path.gcomp, '.png'), g.compare.lm.mcmc.ms,
         width = 9, height = 7)
  
}
