#
# This script is very similar to the "main" one. 
# One variant/location pair is left out of the 
# fitting dataset to subsequently check the inference.
#

message('\n', rep('=',50), 
        '\n\t LEAVE-ONE-OUT INFERENCES \n',
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
source('figures.R')


set.seed(12345)

# Main global parameters
data.source    = 'DAD'
logscale       = TRUE
variant.ignore = c('B.1.438.1') 

# Data 

dd = digest_data(data.source = data.source, 
                 logscale = logscale,
                 variant.ignore = variant.ignore,
                 do.plot = FALSE)

dat  = dd$data |> mutate(vl = paste(variant, city))
dflm = dd$lm


g.dat = dat |> 
  group_by(vl) |> 
  summarise(m = mean(hosp.count)) |> 
  ggplot(aes(x=reorder(vl,m), y = m)) + 
  geom_col()+
  coord_flip()+
  theme(panel.grid = element_blank())+
  labs(title = 'Sorted mean hospital admissions',
       y = 'mean hospital admissions per 100,000 population',
       x = 'Variant / City')
pdf('../ms/figs/appendix-data-hosp-sorted.pdf', height = 12)
plot(g.dat)
dev.off()

# HMC parameters

mcmc.params = fetch_mcmc_params()
nburnin     = mcmc.params$burnin
niter       = mcmc.params$iter
nchains     = mcmc.params$chains

message('\n  iter: ', niter,
        '\n  burn: ', nburnin,
        '\n  chains: ', nchains)

stopifnot(niter > nburnin + 10)
ctrl.nuts = list(
  delta = 0.90, # default = 0.80 
  maxTreeDepth = 13  # default = 10
)  

model.type = 'simple'  # simple, lagW, ARlagW, cauchy

message('Model type: ', model.type)


vl.pairs = unique(dat$vl)

# DEBUG
#i = which(vl.pairs == 'BA.5 Regina')
#ii = c(1,12, 33, 41)
ii = seq_along(vl.pairs)
print(paste(seq_along(ii),vl.pairs[ii], sep = '--'))

# For downstream analysis (other script)
saveRDS(list(dat=dat, vl.pairs=vl.pairs), 
        file = '../out/datvlpairs.rds')

# -- Loop on inferences


for(i in ii){
  try({
    message(" Leaving out < ",vl.pairs[i]," > from inference.\n")
    dat.i = filter(dat, vl != vl.pairs[i])
    
    message('\nH-MCMC launched.')
    
    tmp = paste(names(ctrl.nuts), ctrl.nuts, sep = ' = ')
    message('H-MCMC-NUTS controls:\n',
            paste(tmp, collapse = '\n'))
    
    t1 = Sys.time()
    cluster.mcmc = parallel::makeCluster(nchains)
    mcmcp = parLapply(cl  = cluster.mcmc,
                      X   = 1:nchains,
                      fun = run_one_parallel,
                      dat = dat.i, 
                      mcmc.params = mcmc.params,
                      unique.init.chain = TRUE,
                      ctrl.nuts = ctrl.nuts,
                      model.type = model.type)
    stopCluster(cluster.mcmc)
    saveRDS(mcmcp, file = paste0('../out/mcmcp-leaveout-',i,'.rds'))
    
    dt = difftime(Sys.time(), t1,units = 'mins') |> round(2)
    message('\nHMCMC ', vl.pairs[i],' done in ',dt,' mins\n')
    
    mcmcobj = digest_inference(mcmcp, nchains)
    file.mcmc.i =  paste0("../out/mcmcobj-leave-",i,".rds")
    saveRDS(mcmcobj, file = file.mcmc.i)
  }) 
}

message('\n=== H-MCMC LEAVE ONE OUT inference completed ===\n')
