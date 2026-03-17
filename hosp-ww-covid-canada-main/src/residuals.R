###
###   CALCULATE (ERROR) RESIDUALS FROM POSTERIOR SAMPLES
###
###   Residual := Hhat(variant, geo) - Hobs(variant, geo)
###   This can only be calculated for observations, hence
###   at the city/variant level (not for latent levels, eg prov and universal)


suppressPackageStartupMessages({
  library(tidyr)
  library(dplyr)
  library(ggplot2) ; theme_set(theme_bw())
  library(lubridate)
  library(stringr)
  library(patchwork)
  library(snowfall)
  library(xtable)
})




calc_residuals_t <- function(t, dat.vl, m.vec, b.vec) {
  # t = 1
  
  w.t = dat.vl$w[t]
  h.t = dat.vl$h[t]
  
  h.t.hat = m.vec * w.t + b.vec
  
  residual.t = h.t.hat - h.t
  
  if(0){
    plot(density(h.t.hat))
    abline(v = h.t)
    plot(density(residual.t))
  }
  
  res = data.frame(
    replicate = seq_along(residual.t),
    time = t,
    residual = residual.t,
    variant.idx = dat.vl$variant.idx[1],
    geo.idx = dat.vl$geo.idx[1],
    geo = dat.vl$geo[1],
    variant = dat.vl$variant[1]
  )
  return(res)
}


calc_acf_r <- function(r, dfres) {
  # r = 33
  dfres.r = dfres[dfres$replicate == r,]
  acf.r   = acf(dfres.r$residual, lag.max = 10, plot = FALSE)
  acf.r.1 = acf.r$acf[-1]
  res = data.frame(repl = r, 
                   correl = acf.r.1, 
                   lag = seq_along(acf.r.1))
  res
}

#' Calculate error residuals and
#' autocorrelations of errors
calc_residual_vargeo <- function(vari, 
                                 geoi, 
                                 post, 
                                 do.plot = FALSE) {
  # vari =3 ; geoi = 1
  
  dat.vl = filter(dat, 
                  geo.idx == geoi, 
                  variant.idx == vari)
  if(nrow(dat.vl) == 0){
    return(list(res = NULL, acf = NULL))
  }
  
  post.vl = filter(post, 
                   geo.idx == geoi,
                   variant.idx == vari)
  
  b.vec = post.vl$value[post.vl$prm == 'b']
  m.vec = post.vl$value[post.vl$prm == 'm']
  
  var.geo = paste(dat.vl$variant[1], dat.vl$geo[1], sep = ' - ')
  
  g.line = NULL
  if(do.plot){
    g.line = ggplot(dat.vl, aes(x=w,y=h)) + 
      geom_smooth(method = 'lm', formula = 'y~x', se = T,
                  alpha = 0.15, 
                  color = 'steelblue', 
                  fill = 'steelblue')+
      geom_abline(slope = mean(m.vec), 
                  intercept = mean(b.vec), 
                  color = 'indianred', linewidth = 2)+
      geom_point(size=2) +
      labs(title = 'Linear regression', subtitle =var.geo,
           caption = 'blue: naive lm() ; red: Bayesian model')
    g.line
  }
  
  # Calculate residuals at all times
  dfres = lapply(1:nrow(dat.vl), 
                 calc_residuals_t, 
                 m.vec = m.vec, 
                 b.vec = b.vec,
                 dat.vl = dat.vl) |>
    bind_rows()
  
  # autocorrelation (in time) of residuals
  sfInit(parallel = TRUE, cpus = 4)
  sfExport('dfres', 'calc_acf_r')
  df.acf = sfLapply(
    unique(dfres$replicate), 
    calc_acf_r,
    dfres = dfres) |>
    bind_rows() |> 
    mutate(variant = dat.vl$variant[1], 
           geo = dat.vl$geo[1])
  sfStop()
  
  df.acf$lag = as.factor(df.acf$lag)
  
  if(0){
    g.acf = df.acf |> 
      ggplot(aes(x=lag, y=correl)) +
      geom_hline(yintercept = 0)+
      geom_boxplot()+
      scale_y_continuous(limits = c(-1,1))+
      labs(title = )
    g.acf
  }
  
  out = list(
    res = dfres,
    acf = df.acf,
    plot = g.line
  )
  return(out)
}

# ---- RUN ----

do.plot.fits = TRUE

# Retrieve data and fit objects
dat      = readRDS('../out/data-used.rds')
mcmcobj  = readRDS('../out/mcmcobj.rds')
post.all = mcmcobj$post

# Extract the coefficients used at the observations level
# and extract the variant and location indices

post = filter(post.all, grepl('^[mb]\\[', name)) |> 
  separate(col = name, into = c('tmp1','tmp2'), sep = ',',
           remove = FALSE) |>
  mutate(
    variant.idx = str_extract(tmp1, '\\d+'),
    geo.idx = str_extract(tmp2, '\\d+'),
    prm = str_extract(name, '^\\w')
  ) |>
  select(-starts_with('tmp'))

# Loop through variants and locations
df.loop = dat |> 
  select(starts_with('variant'), 
         starts_with('geo')) |>
  distinct()

tmp = list()
for(i in 1:nrow(df.loop)){
  message(i, '/', nrow(df.loop))
  tmp[[i]] = calc_residual_vargeo(
    vari = df.loop$variant.idx[i],
    geoi = df.loop$geo.idx[i], 
    post = post,
    do.plot = do.plot.fits)  
}

residuals = lapply(tmp, '[[', 1) |> bind_rows()
resid.acf = lapply(tmp, '[[', 2) |> bind_rows()

saveRDS(residuals, '../out/residuals.rds')
saveRDS(residuals, '../out/residuals-acf.rds')

# For table of residuals
tab.res = residuals |>
  group_by(geo, variant) |> 
  summarise(m = median(residual, na.rm = TRUE), 
            qlo = quantile(residual, probs = 0.25, na.rm = TRUE),
            qhi = quantile(residual, probs = 0.75, na.rm = TRUE),
            rmse = sqrt(mean(residual^2)),
            mae = mean(abs(residual)),
            .groups = 'drop') |>
  mutate(tmp = paste0(round(m,2), ' (',
                      round(qlo,2), ' ; ',
                      round(qhi,2), ')')) |> 
  select(geo, variant, tmp, rmse, mae)|> 
  rename(`median residual (IQR)` = tmp,
         location = geo)

summary(tab.res)

# Export to a LaTeX table
xtable::xtable(tab.res, label = 'tab:residuals',
               caption = 'Residual errors of the fitted hierarchical model.') |>
  print(file = '../ms/tables/residuals-tbl.tex')

# ---- PLOTS

q = 0.95   # quantile interval for plots

if(do.plot.fits){
  g.fits = lapply(tmp, '[[', 3) |> wrap_plots(ncol = 8)
  pdf('../figs/diagn-fits.pdf', width = 25, height = 15)
  plot(g.fits)
  dev.off()
}


# Residuals with time dimension

g.res.t = residuals |> 
  group_by(time, variant, geo) |>
  summarise(
    m = mean(residual),
    lo = quantile(residual, probs = 0.5 - q/2),
    hi = quantile(residual, probs = 0.5 + q/2),
    .groups = 'drop' ) |>
  ggplot(aes(x=time, y = m)) + 
  geom_hline(yintercept = 0, color = 'darkgreen') + 
  geom_segment(aes(xend=time,y=lo, yend =hi))+
  geom_point()+
  theme(panel.grid.minor = element_blank())+
  facet_grid(variant ~ geo)+
  labs(title = 'posterior residuals')
g.res.t

# Summarized residuals

g.res = residuals |> 
  ggplot(aes(x=variant, y = residual)) + 
  geom_hline(yintercept = 0, color = 'grey40', 
             linewidth = 0.5, linetype = 'dashed') + 
  geom_boxplot(color = 'black', outliers = FALSE)+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1, size = rel(0.8)),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_rect(fill = 'steelblue4'),
        strip.text = element_text(color = 'white', face = 'bold'))+
  facet_wrap( ~ geo)+
  guides(fill = 'none', color = 'none')+
  labs(title = 'posterior residuals', subtitle = 'averaged over time',
       y = 'residual\nlog(H_estimated / H_observed)')
g.res  

# For manuscript
pdf(file = '../ms/figs/appendix-residuals.pdf', 
    width = 10)
plot(g.res)
dev.off()

message('residuals:',
        '\n   median = ', median(residuals$residual),
        '\n   mean   = ', mean(residuals$residual))

# Autocorrelation of residuals

g.acf = resid.acf |> 
  ggplot(aes(x = factor(lag), y = correl)) + 
  geom_hline(yintercept = 0, color = 'darkgreen') + 
  geom_boxplot(outliers = FALSE) + 
  facet_grid(variant ~ geo)+
  theme(panel.grid.minor = element_line(color = 'grey97'),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(size = rel(0.7))) + 
  scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,by=0.5))+
  labs(title = 'Autocorrelation of residuals', 
       x = 'lag', y = 'correlation')
g.acf

now = format(Sys.time(),'%Y-%m-%dT%H%M%S')
fname = paste0('../figs/plot-residuals-post-',now,'.pdf')
pdf(file = fname, width = 20, height = 10)
plot(g.res)
plot(g.res.t)
plot(g.acf)
dev.off()




