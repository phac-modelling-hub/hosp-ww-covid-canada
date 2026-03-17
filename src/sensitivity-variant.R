# Sensitivity analysis at different variant threshold/pct above levels

suppressMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2) ; theme_set(theme_bw())
  library(lubridate)
  library(stringr)
  library(patchwork)
  library(parallel)
  
})

source('plot.R')
source('utils.R')
source('model-hmc.R')
source('var.R')

set.seed(12345)

# Sensitivity parameters
thresholds = c(0.6, 0.75, 0.45)
windows = c(8, 10, 6)
pct.above = c(0.75, 0.95, 0.55)

# Generate various model objects
mdlobj.sens = create_sens_mdl_obj(thresholds, windows, pct.above)
dat.sens = lapply(mdlobj.sens, prepare_data,
                  data.source = 'DAD', logscale = TRUE,
                  variant.ignore = c('B.1,438,1'), sensitivity = TRUE)

# Run HMC on model objects
iter = 1e2
nburnin = 5e1
nchains = 4

mcmc.sens = run_sens_mcmc(dat.sens,
                          iter, nburnin, nchains) 

saveRDS(mcmc.sens, file = '../out/sensitivity-mcmc.rds')

message("\nGenerating summary table...")
# Generate summary table
table = lapply(mcmc.sens, create_sens_sum_table) |>
  bind_rows()

# Convert to latex format and save to ms folder
tbl.latex = knitr::kable(table, 
                         format = 'latex',
                         align = c('l','l', 'l', 'c','c'),
                         linesep = '', 
                         vline= '',
                         escape = FALSE)
writeLines(tbl.latex, '../ms/tables/sens-tbl.tex')

message("\nSensitivity analysis complete.")
