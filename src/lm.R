###
###   SIMPLE LINEAR REGRESSION ON DATA
###


library(tidyr)
library(dplyr)
library(ggplot2) ; theme_set(theme_bw())
library(lubridate)
library(stringr)
library(patchwork)

#' Perform linear regressions independently
#' for each province and variant
#' @param dat Dataframe containing hospital admissions and wastewater data
lm_prov_var <- function(dat) {
  tmp = list() ; k=1
  
  for(p in unique(dat$prov)){
    for(v in unique(dat$variant)){ # p = 'ON'; v = 'BA.1'
      
      d = filter(dat, prov == p, variant == v)
      
      if(nrow(d)>2){
        m  = lm(data = d, formula = 'h ~ w')
        ci = confint(object = m)
        s = summary(m)
       
        if(0){  # DIAGNOSTICS for debug
          r = m$residuals
          
          acf(r)
          plot(m)
          lmtest::dwtest(m)
        }
         
        tmp[[k]] = data.frame(
          prov     = p,
          city     = unique(d$city),
          variant  = v,
          slope.m  = m$coefficients[2],
          slope.lo = ci[2,1],
          slope.hi = ci[2,2],
          intercept.m  = m$coefficients[1],
          intercept.lo = ci[1,1],
          intercept.hi = ci[1,2],
          n            = nrow(d),
          Rsqradj      = s$adj.r.squared)
        k = k+1
      }
    }
  }
  
  dflm = bind_rows(tmp)
  return(dflm)
}




