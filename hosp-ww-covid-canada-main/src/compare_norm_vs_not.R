###
### This script is a late addition following 
### a request by a reviewer of the manuscript. 
###


suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2) ; theme_set(theme_bw())
  library(lubridate)
  library(stringr)
  library(patchwork)
})


# * * * WARNING * * *

# Run twice the main script, with and without normalization,
# save the RDS files somewhere, and update their path below (manual!)
# The switch to turn on/off normalization is in:
# howareg/prm/model-param.csv --> `norm.ww = TRUE/FALSE`

obj.raw = readRDS('../out/mcmcobj.rds')
obj.nrm = readRDS('../out/mcmcobj-norm.rds')


post.raw = obj.raw$postsum |> mutate(type = 'raw')
post.nrm = obj.nrm$postsum |> mutate(type = 'normalized')

post = bind_rows(post.nrm, post.raw) |> 
  mutate(prmtype = str_extract(name, '^\\w+'))

post.plot = filter(post, prmtype %in% c('M','B', 'mu', 'beta'))

dodgewidth = 0.7

g = post.plot |> ggplot(aes(x = name, color = type, fill= type))+
  facet_wrap(~prmtype, scales = 'free', ncol = 2) +
  geom_linerange(aes(ymin = lo, ymax = hi),
                 linewidth = 1.2,
                 alpha =0.7,
               position = position_dodge(width=dodgewidth))+
  geom_point(aes(y=m),
             shape = 23, fill = 'white', size = 2,
             position = position_dodge(width=dodgewidth))+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank()) + 
  labs(title = 'Posteriors comparison: raw vs. normalized wastewater',
       x = 'parameter', y = 'value', 
       color = 'Wastewater\nconcentration\ndata type')

g

pdf('../ms/figs/appendix-comp-ww-norm.pdf',
    height = 4, width = 9 )
plot(g)
dev.off()
