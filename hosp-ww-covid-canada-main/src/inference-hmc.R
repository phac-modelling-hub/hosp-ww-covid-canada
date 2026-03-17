
message('\n\n==== HMC inference on Simulated Data ====\n\n')

suppressMessages({
  library(tidyverse)
  library(ggplot2) ; theme_set(theme_bw())
  library(lubridate)
  library(stringr)
  library(patchwork)
  library(coda)
  library(nimble)
  library(nimbleHMC)
  library(parallel)

  source("plot.R")
  source("model-hmc.R")
  source("utils.R")
})

set.seed(12345)

t1 = Sys.time()
plot.correl = 0


## ----- observation data -----

message("Reading time-series observation data...\n")
obs = readRDS("../out/obj-model.rds")

g.obs = plot_obs(obs, data.source = 'DAD', logscale = 1)

tally.obs = obs|> group_by(geo,variant) |> tally() |> arrange(-n)
pdf('../figs/plot-obs-sim.pdf', width = 16, height = 12)
plot(g.obs$wh)
plot(g.obs$ts.h)
plot(g.obs$ts.w)
dev.off()

## ---- simulated data ----

df    = create_sim_data(obs)
draws = df$draws
dat   = df$dat |> 
  mutate(h = hosp.count, w = ww.conc,
         variant.idx = as.integer(as.factor(variant)),
         geo.idx = as.integer(as.factor(geo))) |> 
  group_by(geo, variant) |> 
  arrange(date, .by_group = TRUE) |> 
  mutate(w_1 = lag(w, n=1),
         h_1 = lag(h, n=1)) |> 
  filter(! is.na(w_1), ! is.na(w_1))

truth = df$truth
g.simdata  = plot_simulated_data(dat)

pdf(file = "../ms/figs/appendix-sim-data.pdf", 
    width = 7, height = 7)
plot(g.simdata$xy)
dev.off()
pdf(file = "../figs/simulated-data.pdf", 
    width = 7, height = 7)
plot(g.simdata$xy)
plot(g.simdata$ts.w)
plot(g.simdata$ts.h)
dev.off()


# ----- Inference ----

mcmc.params = fetch_mcmc_params()
nburnin     = mcmc.params$burnin
niter       = mcmc.params$iter
nchains     = mcmc.params$chains

stopifnot(niter > nburnin + 10)

ctrl.nuts = list(
  delta = 0.90, # default = 0.80 
  maxTreeDepth = 13  # default = 10
)

message('HMC algorithm started with:\n',
        ' n.iter = ', niter,
        '\n n.burn = ', nburnin, 
        '\n n.chains = ', nchains)


t1 = Sys.time()
cluster.mcmc = parallel::makeCluster(nchains)
mcmcp = parLapply(cl = cluster.mcmc,
                 X = 1:nchains,
                 fun = run_one_parallel,
                 dat = dat, 
                 mcmc.params = mcmc.params,
                 unique.init.chain = TRUE,
                 ctrl.nuts = ctrl.nuts)
stopCluster(cluster.mcmc)
dt = difftime(Sys.time(), t1,units = 'mins') |> round(2)
message('>> Parallel HMC done in ',dt,' mins\n')

mcmcobj = digest_inference(mcmcp, nchains = nchains)
saveRDS(mcmcobj, file = "../out/mcmcobj-sim.rds")

## ---- Run Gelman diagnostic ---

gelman_diag  = run_gelman(mcmcobj$mcmc)
gelman.mpsrf = gelman_diag[["gelman.mpsrf"]]
gelman.psrf  = gelman_diag[["gelman.psrf"]]

# save summarized posterior dataframe w/ Gelman and truth data 
err = filter_error_truth(mcmcobj, gelman.psrf, truth)
df.postsum = err$postsum
df.error   = err$error
save(draws,df.postsum, df.error, file = "../out/sim-model-post-sim.RData")

diagnostic_simdata(df.postsum, df.error, gelman.mpsrf)


## ---- plot model output ----

# plot posterior distributions (at universal level,
# variant and location level have rendering issues due to number of coeff.)

df.post  = process_post(mcmcobj, gelman.psrf)
g.post.u = plot_post_level(df.post, level = "universal")
g.post.v = plot_post_level(df.post, level = "variant")

chains.prms = c('mu_mean', 'beta_mean', 'rho1_mean')

g.chains = df.post |> 
  filter(name %in% chains.prms, iter < 1e4) |> 
  ggplot(aes(x=iter, y=value, color = chain)) + 
  facet_grid(name ~ chain, scales = 'free_y')+
  # geom_line() + 
  geom_point(alpha = 0.3)+
  labs(title = 'Trace plot')
g.chains

if(plot.correl){
  g.post.c1 = plot_post_2d(mcmcobj$post, pattern = '^(M|m.1|mu)')
  g.post.c2 = plot_post_2d(mcmcobj$post, pattern = '^(B|be|b.1)')
}

g.univ.sigma = plot_universal_sigma_post(mcmcobj$post)

# --- save plots


# For manuscript appendix

pdf('../ms/figs/appendix-sim-inference-univ.pdf', 
    width = 7, height = 3)
plot(g.post.u)
dev.off()

pdf('../ms/figs/appendix-sim-inference-variant.pdf', 
    width = 13, height = 9)
plot(g.post.v)
dev.off()



fname = paste0("../figs/sim-model-outputs-burn",
               nburnin,"-post",niter-nburnin,".pdf")
pdf(file = fname, width = 12, height = 12)
plot(g.post.u)
plot(g.post.v)
plot(g.univ.sigma)
plot(g.chains)
if(plot.correl){
  plot(g.post.c1$density)
  plot(g.post.c1$correl)
  plot(g.post.c2$density)
  plot(g.post.c2$correl)
}
dev.off()

t2 = Sys.time()
diff = difftime(t2, t1, units = "mins")
message("\n--- Inference on simulated data completed in ", 
              round(diff,2), " mins.\n")

