suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2) ; theme_set(theme_bw())
  library(lubridate)
  library(stringr)
  library(patchwork)
})

source('plot.R')
source('utils.R')
source('func-data.R')
source('figures.R')

# Retrieve upstream objects
a = read_rds('../out/datvlpairs.rds')
dat = a$dat
vl.pairs = a$vl.pairs

# -- Loop on diagnostics

# Process only the successful iterations
fout = list.files(path = '../out', pattern = 'mcmcobj-leave-', full.names = TRUE)
foutnum = grep('\\d+', fout)

gel.mu = gel.beta = numeric(length(foutnum))
tmp.post = tmp.gel = list() 

message('Running diagnostics...')

for(i in seq_along(fout)){
  mcmcobj = readRDS(fout[i])
  loo.i   =  vl.pairs[foutnum[i]]
  message(i, ' - Extract posterior and Gelman diagnostic: ', loo.i)
  # Gelman-Rubin statistics
  try({
    gelman_diag  = run_gelman(mcmcobj$mcmc)
    gelman.psrf  = gelman_diag[["gelman.psrf"]]
    tmp.gel[[i]] = data.frame(
      mu   = gelman.psrf$psrf[gelman.psrf$name == 'mu'],
      beta = gelman.psrf$psrf[gelman.psrf$name == 'beta'],
      sigma_mu   = gelman.psrf$psrf[gelman.psrf$name == 'sigma_mu'],
      sigma_beta = gelman.psrf$psrf[gelman.psrf$name == 'sigma_beta'],
      vl.left = loo.i
    )
  })
  tmp.post[[i]] = mcmcobj$postsum |> 
    filter(name %in% c('mu', 'beta', 'sigma_mu', 'sigma_beta')) |>
    mutate(vl.left = loo.i)
}  

post.loo   = bind_rows(tmp.post)
gelman.loo = bind_rows(tmp.gel) |> pivot_longer(-vl.left)

g.gelman = plot_leaveoneout_gelman(gelman.loo)
g.post   = plot_leaveoneout_post(post.loo)

pdf('../ms/figs/appendix-loo-gelman.pdf', width = 14)
plot(g.gelman)
dev.off()

pdf('../ms/figs/appendix-loo-post.pdf', width = 14)
plot(g.post)
dev.off()

# Generate main text figure

fig.tmp = list()
mhvec = numeric(length(foutnum))

for(j in seq_along(foutnum)){
  loo = vl.pairs[foutnum[j]]
  message('Generating leave-one-out figures: ',loo)
  dat.only.j = filter(dat, vl == loo)
  loo = dat.only.j$vl[1] 
  mhvec[j] = mean(dat.only.j$hosp.count)
  fig.tmp[[j]] = fig_main_result(fout[j], 
                                 add.synth.data = F,
                                 xlim = c(0,400), 
                                 ylim = c(0,45)) +
    geom_point(data = dat.only.j, aes(x=ww.conc, y = hosp.count), 
               shape = 16, size = 2) + 
    labs(title = paste('Leave-one-out:', loo))
}


fig.tmp2 = fig.tmp[order(mhvec)]

# Final figure comparing the "new" data (the one left out)
# to the inferred range of severities.

fig.main = wrap_plots(fig.tmp2)

file.fig = paste0('../ms/figs/appendix-loo-postHW.pdf')
pdf(file = file.fig, width=24, height = 18)
plot(fig.main)
dev.off()
