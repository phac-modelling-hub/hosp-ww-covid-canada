## Generate figures for LaTeX document ##
library(tidyverse)
library(patchwork)
library(kableExtra)

source('figures.R')

set.seed(12345)

# Retrieve files
file.mcmc = '../out/mcmcobj.rds'

# Generate main text figure
fig.main = fig_main_result(file.mcmc, 
                           add.synth.data = TRUE, 
                           legend.version = TRUE) 

fig.main.path = '../ms/figs/main'

pdf(paste0(fig.main.path, '.pdf'), height = 3, width = 5)
plot(fig.main)
dev.off()

ggsave(paste0(fig.main.path, '.png'), fig.main)

# Appendix figures
figure_appendix_mcmcoutputs()
