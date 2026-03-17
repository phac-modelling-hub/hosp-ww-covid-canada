source('utils.R')

digest_data <- function(data.source, 
                        logscale, 
                        variant.ignore, 
                        do.plot = TRUE) {
  
  # All (raw) data
  obs = read_rds('../out/obj-model.rds')
  
  if(do.plot){
    message('plotting observations... ', appendLF = F)
    g.obs = plot_obs(obs, data.source, logscale = 1)
    pdf('../figs/plot-obs.pdf', width = 12, height = 10)
    plot(g.obs$wh)
    plot(g.obs$ts.h)
    plot(g.obs$ts.w)
    dev.off()
    message('done.')
  }
  
  # Prepare data: retrieve all observations, 
  # filter, transform and format.
  message('preparing data... ', appendLF = F)
  dat = prepare_data(data.source, logscale, variant.ignore)
  saveRDS(dat, '../out/data-used.rds')
  message('done.')
  if(do.plot){
    message('plotting data used... ', appendLF = F)
    pdf('../figs/plot-data.pdf', width = 12, height = 10)
    plot(plot_data(dat))
    dev.off()
    path.datams = '../ms/figs/plot-data'
    pdf(paste0(path.datams, '.pdf'), width = 6, height = 5)
    plot(plot_data(dat, ms = TRUE))
    dev.off()
    ggsave(paste0(path.datams, '.png'), 
           plot_data(dat, ms = TRUE,
                     size = 14.5),
           width = 12, height = 10)
    
    pdf('../ms/figs/appendix-data-tally.pdf', width = 12)
    plot(plot_data_tally(dat))
    dev.off()
    message('done.')
  }
  
  # Perform naive independent linear regressions
  message('naive linear regression... ', appendLF = F)
  dflm = lm_prov_var(dat)
  message('done.')
  
  if(do.plot){
    message('plotting naive linear regressions... ', appendLF = F)
    g.lm = plot_lm_data(dat)
    pdf('../figs/plot-lm.pdf', width = 15)
    for(i in seq_along(g.lm)) plot(g.lm[[i]])
    dev.off()
    message('done.')
  }
  
  out = list(
    data = dat, 
    lm = dflm
  )
}
