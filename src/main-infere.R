###
###  MAIN SCRIPT 
###

message('\n', rep('=',50), 
        '\n\tMAIN INFERENCE WITH OBSERVATIONS\n',
        rep('=',50),' \n\n')

suppressMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2) ; theme_set(theme_bw())
  library(lubridate)
  library(stringr)
  library(patchwork)
  library(parallel)
})

source('lm.R')
source('plot.R')
source('utils.R')
source('model-hmc.R')
source('func-data.R')

set.seed(12345)

# Main global parameters
data.source    = 'DAD'
logscale       = TRUE
variant.ignore = c('B.1.438.1') 

# ---- Data ----

dd = digest_data(data.source = data.source, 
                 logscale = logscale,
                 variant.ignore = variant.ignore,
                 do.plot = TRUE)

dat  = dd$data
dflm = dd$lm


# ---- H-MCMC ----

message('\nH-MCMC launched.')
mcmc.params = fetch_mcmc_params()
nburnin     = mcmc.params$burnin
niter       = mcmc.params$iter
nchains     = mcmc.params$chains

message('\n  iter: ', niter,
        '\n  burn: ', nburnin,
        '\n  chains: ', nchains)

stopifnot(niter > nburnin + 10)

# HMC-NUTS tweaking to improve convergence
#
#  *** WARNING: FOR EXPERT USERS ***
#  
#  KEEP DEFAULT VALUES UNLESS YOU KNOW WHAT YOU ARE DOING
#
#  - https://mc-stan.org/docs/2_18/reference-manual/hmc-algorithm-parameters.html
#  - help('NUTS') (library(nimbleHMC)) 
#    `delta`: target acceptance probability used 
# during the initial period of step-size adaptation. (default = 0.8)
#    `maxTreeDepth`: The maximum allowable depth of the binary leapfrog
# search tree for generating candidate transitions. (default = 10)
#
ctrl.nuts = list(
  delta = 0.90, # default = 0.80 
  maxTreeDepth = 13  # default = 10
  )  

model.type = 'simple'  # simple, lagW, ARlagW, cauchy

message('Model type: ', model.type)

tmp = paste(names(ctrl.nuts), ctrl.nuts, sep = ' = ')
message('H-MCMC-NUTS controls:\n',
        paste(tmp, collapse = '\n'))

t1 = Sys.time()
cluster.mcmc = parallel::makeCluster(nchains)
mcmcp = parLapply(cl = cluster.mcmc,
                  X = 1:nchains,
                  fun = run_one_parallel,
                  dat = dat, 
                  mcmc.params = mcmc.params,
                  unique.init.chain = TRUE,
                  ctrl.nuts = ctrl.nuts,
                  model.type = model.type)
stopCluster(cluster.mcmc)
saveRDS(mcmcp, file = '../out/mcmcp.rds')

dt = difftime(Sys.time(), t1,units = 'mins') |> round(2)
message('>> Parallel H-MCMC done in ',dt,' mins\n')

mcmcobj = digest_inference(mcmcp, nchains)
saveRDS(mcmcobj, file = "../out/mcmcobj.rds")
saveRDS(mcmcobj, file = paste0("../out/mcmcobj-",model.type,".rds"))

# ---- Diagnostics ----

message('Running diagnostics...')

# Gelman-Rubin statistics
try({
  gelman_diag  = run_gelman(mcmcobj$mcmc)
  gelman.mpsrf = gelman_diag[["gelman.mpsrf"]]
  gelman.psrf  = gelman_diag[["gelman.psrf"]]
  g.gelman = plot_gelman(gelman.psrf)
  pdf('../figs/diagn-gelman.pdf', width = 12)
  plot(g.gelman)
  dev.off()
})
message('  Gelman-Rubin done.')

# Traces
g.trace = plot_trace(mcmcobj)
pdf('../figs/diagn-trace.pdf', width = 18, height = 12)
plot(g.trace)
dev.off()

message('  Traces done.')

# Correlations
g.2d   = plot_post_2d(post = mcmcobj$post, pattern = 'M|mu')
g.2d.B = plot_post_2d(post = mcmcobj$post, pattern = 'B|beta')

pdf('../figs/diagn-correl2D.pdf', 
    width = 30, height = 20)
plot(g.2d$density)
plot(g.2d$correl)
plot(g.2d.B$density)
plot(g.2d.B$correl)
dev.off()

message('  Correlations done.')

# ---- inference ----

message('Plotting inferences')
g.postall = plot_post_all(mcmcobj)
pdf(paste0('../figs/post-all-',model.type,'.pdf'), 
    width = 30, height = 20)
plot(g.postall)
dev.off()

g.post.slope = plot_post(mcmcobj, pattern.prm = 'M|^mu')
g.post.autor = plot_post(mcmcobj, pattern.prm = 'P|^rho')
g.post.inter = plot_post(mcmcobj, pattern.prm = 'B|^be')
g.post.univ  = plot_post(mcmcobj, pattern.prm = 'mu|beta|^sigma|^rho')

pdf(paste0('../figs/post-prms-',model.type,'.pdf'), 
    width = 20, height = 12)
plot(g.post.slope)
plot(g.post.inter)
try({plot(g.post.autor)}, silent = TRUE)
plot(g.post.univ)
dev.off()

# Comparison with naive linear regression
g.compare.lm.mcmc = plot_compare_lm_mcmc(dflm, mcmcobj, dat)
pdf('../figs/lm-vs-mcmc.pdf', width = 14)
plot(g.compare.lm.mcmc)
dev.off()

message('\n=== H-MCMC inference completed ===\n')
