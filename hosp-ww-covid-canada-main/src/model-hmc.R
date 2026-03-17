## ---- Functions for Hamiltonian MCMC execution ---- ## 

library(tidyverse)
library(ggplot2) ; theme_set(theme_bw())
library(lubridate)
library(stringr)
library(patchwork)
library(coda)
library(nimble)
library(nimbleHMC)
library(parallel)

# Retrieve model constants
#' @param dat Dataframe to be used in model
get_constants <- function(dat) {
  return(list(
    N = nrow(dat),
    V = max(dat$variant.idx),
    L = max(dat$geo.idx)
  ))
}
# Generate model definition
#' @param dat Dataframe to be used in model
define_model <- function(dat) {
  
  csts = get_constants(dat)
  N = csts[['N']]
  L = csts[['L']]
  V = csts[['V']]
  
  res = nimbleCode(
    {
      # --- priors
      
      # universal slope
      mu <- mu_mean + sigma_mu * muu
      muu ~ dnorm(0,1)
      mu_mean ~ dnorm(1,0.5)
      
      # universal intercept
      beta <- beta_mean + sigma_beta * betau
      betau ~ dnorm(0,1)
      beta_mean ~ dnorm(0,1)
      
      sigma_mu  ~ dexp(1)
      sigma_beta ~ dexp(1)
      sigma_M  ~ dexp(1)
      sigma_B  ~ dexp(1)
      sigma_m  ~ dexp(1)
      sigma_b  ~ dexp(1)
      sigma    ~ dexp(4)
      
      # ---
      # Define hierarchy for mean slope & intercept 
      # at variant and location levels
      
      for(v in 1:V){
        Mu[v]  ~ dnorm(0,1)
        Bu[v]  ~ dnorm(0,1)
        
        M[v]  <- mu   + sigma_M  * Mu[v]
        B[v]  <- beta + sigma_B  * Bu[v]
        
        for(l in 1:L){
          m_u[v,l]  ~ dnorm(0,1)
          b_u[v,l]  ~ dnorm(0,1)
          
          m[v,l]  <- M[v]  + sigma_m  * m_u[v,l] 
          b[v,l]  <- B[v]  + sigma_b  * b_u[v,l]
        }
      }  
      
      # hypothetical and likelihood estimation
      for (i in 1:N) {
        # Hypothetical / latent levels:
        # universal
        h_univ[i] <- mu  * w[i] + beta
        
        # variant level
        h_var[i] <- M[var[i]] * w[i] + B[var[i]] 
        
        # Likelihood at location level
        h[i] <- m[var[i],loc[i]] * w[i] + b[var[i],loc[i]] 
        
        y[i] ~ dnorm(h[i], sd = sigma)
      }
    }
  )
  return(res) 
}



#' Generate model definition using Cauchy distribution
#' instead of normal for likelihood
#' @param dat Dataframe to be used in model
#' 
define_model_cauchy <- function(dat) {
  
  csts = get_constants(dat)
  N = csts[['N']]
  L = csts[['L']]
  V = csts[['V']]
  
  res = nimbleCode(
    {
      # --- priors
      
      # universal slope
      mu <- mu_mean + sigma_mu * muu
      muu ~ dnorm(0,1)
      mu_mean ~ dnorm(1,0.5)
      
      # universal intercept
      beta <- beta_mean + sigma_beta * betau
      betau ~ dnorm(0,1)
      beta_mean ~ dnorm(0,1)
      
      sigma_mu  ~ dexp(1)
      sigma_beta ~ dexp(1)
      sigma_M  ~ dexp(1)
      sigma_B  ~ dexp(1)
      sigma_m  ~ dexp(1)
      sigma_b  ~ dexp(1)
      sigma    ~ dexp(4)
      
      # ---
      # Define hierarchy for mean slope & intercept 
      # at variant and location levels
      
      for(v in 1:V){
        Mu[v]  ~ dnorm(0,1)
        Bu[v]  ~ dnorm(0,1)
        
        M[v]  <- mu   + sigma_M  * Mu[v]
        B[v]  <- beta + sigma_B  * Bu[v]
        
        for(l in 1:L){
          m_u[v,l]  ~ dnorm(0,1)
          b_u[v,l]  ~ dnorm(0,1)
          
          m[v,l]  <- M[v]  + sigma_m  * m_u[v,l] 
          b[v,l]  <- B[v]  + sigma_b  * b_u[v,l]
        }
      }  
      
      # hypothetical and likelihood estimation
      for (i in 1:N) {
        # Hypothetical / latent levels:
        # universal
        h_univ[i] <- mu  * w[i] + beta
        
        # variant level
        h_var[i] <- M[var[i]] * w[i] + B[var[i]] 
        
        # Likelihood at location level
        h[i] <- m[var[i],loc[i]] * w[i] + b[var[i],loc[i]] 
        
        # Cauchy is Students distribution ("dt")
        # with 1 degree of freedom ("df=1")
        y[i] ~ dt(h[i], sigma = sigma, df = 1)
      }
    }
  )
  return(res) 
}




#' Generate model definition. 
#' Use lagged W, ie: 
#' H(t) ~ W(t) + W(t-1)
#' @param dat Dataframe to be used in model
define_model_lagW <- function(dat) {
  
  csts = get_constants(dat)
  N = csts[['N']]
  L = csts[['L']]
  V = csts[['V']]
  
  res = nimbleCode(
    {
      # --- priors
      
      # universal slope
      mu <- mu_mean + sigma_mu * muu
      muu ~ dnorm(0,1)
      mu_mean ~ dnorm(1,0.5)
      
      # For lagged ww:
      mu1 <- mu1_mean + sigma_mu1 * mu1u
      mu1u ~ dnorm(0,1)
      mu1_mean ~ dnorm(1,0.5)
      
      # universal intercept
      beta <- beta_mean + sigma_beta * betau
      betau ~ dnorm(0,1)
      beta_mean ~ dnorm(0,1)
      
      sigma_mu  ~ dexp(1)
      sigma_mu1 ~ dexp(1)
      sigma_beta ~ dexp(1)
      sigma_M  ~ dexp(1)
      sigma_M1 ~ dexp(1)
      sigma_B  ~ dexp(1)
      sigma_m  ~ dexp(1)
      sigma_m1 ~ dexp(1)
      sigma_b  ~ dexp(1)
      sigma    ~ dexp(4)
      
      # ---
      # Define hierarchy for mean slope & intercept 
      # at variant and location levels
      
      for(v in 1:V){
        Mu[v]  ~ dnorm(0,1)
        M1u[v] ~ dnorm(0,1)
        Bu[v]  ~ dnorm(0,1)
        
        M[v]  <- mu   + sigma_M  * Mu[v]
        M1[v] <- mu1  + sigma_M1 * M1u[v]
        B[v]  <- beta + sigma_B  * Bu[v]
        
        for(l in 1:L){
          m_u[v,l]  ~ dnorm(0,1)
          m1_u[v,l] ~ dnorm(0,1)
          b_u[v,l]  ~ dnorm(0,1)
          
          m[v,l]  <- M[v]  + sigma_m  * m_u[v,l] 
          m1[v,l] <- M1[v] + sigma_m1 * m1_u[v,l] 
          b[v,l]  <- B[v]  + sigma_b  * b_u[v,l]
        }
      }  
      
      # hypothetical and likelihood estimation
      for (i in 1:N) {
        # Hypothetical / latent levels:
        # universal
        h_univ[i] <- mu  * w[i] + mu1 * w_1[i] + beta
        
        # variant level
        h_var[i] <- M[var[i]]  * w[i] + M1[var[i]] * w_1[i] + 
          B[var[i]] 
        
        # Likelihood at location level
        h[i] <- m[var[i], loc[i]]  * w[i] + 
          m1[var[i], loc[i]] * w_1[i] + 
          b[var[i],loc[i]] 
        
        y[i] ~ dnorm(h[i], sd = sigma)
      }
    }
  )
  return(res) 
}


#' Generate model definition
#' autoregression on H and W
#' 
#' @param dat Dataframe to be used in model
define_model_AR_H_W <- function(dat) {
  
  csts = get_constants(dat)
  N = csts[['N']]
  L = csts[['L']]
  V = csts[['V']]
  
  res = nimbleCode(
    {
      # --- priors
      
      # universal slope
      mu <- mu_mean + sigma_mu * muu
      muu ~ dnorm(0,1)
      mu_mean ~ dnorm(1,0.5)
      
      # For lagged ww:
      mu1 <- mu1_mean + sigma_mu1 * mu1u
      mu1u ~ dnorm(0,1)
      mu1_mean ~ dnorm(1,0.5)
      
      # universal intercept
      beta <- beta_mean + sigma_beta * betau
      betau ~ dnorm(0,1)
      beta_mean ~ dnorm(0,2)
     
      # universal lagged hosp
      rho1 <- rho1_mean + sigma_rho1 * rho1u
      rho1u ~ dnorm(0,1)
      rho1_mean ~ dnorm(0,2)      
       
      sigma_mu  ~ dexp(3)
      sigma_mu1 ~ dexp(3)
      sigma_beta ~ dexp(3)
      sigma_M  ~ dexp(2)
      sigma_M1 ~ dexp(2)
      sigma_P1 ~ dexp(2)
      sigma_B  ~ dexp(2)
      sigma_m  ~ dexp(2)
      sigma_m1 ~ dexp(2)
      sigma_p1 ~ dexp(2)
      sigma_b  ~ dexp(2)
      sigma_rho1 ~ dexp(2)
      sigma    ~ dexp(10)
      
      # ---
      # Define hierarchy for mean slope & intercept 
      # at variant and location levels
      
      for(v in 1:V){
        Mu[v]  ~ dnorm(0,1)
        M1u[v] ~ dnorm(0,1)
        Bu[v]  ~ dnorm(0,1)
        P1u[v] ~ dnorm(0,1)
        
        M[v]  <- mu   + sigma_M  * Mu[v]
        M1[v] <- mu1  + sigma_M1 * M1u[v]
        P1[v] <- rho1 + sigma_P1 * P1u[v]
        B[v]  <- beta + sigma_B  * Bu[v]
        
        for(l in 1:L){
          m_u[v,l]  ~ dnorm(0,1)
          m1_u[v,l] ~ dnorm(0,1)
          p1_u[v,l] ~ dnorm(0,1)
          b_u[v,l]  ~ dnorm(0,1)
          
          m[v,l]  <- M[v]  + sigma_m  * m_u[v,l] 
          m1[v,l] <- M1[v] + sigma_m1 * m1_u[v,l] 
          p1[v,l] <- P1[v] + sigma_p1 * p1_u[v,l] 
          b[v,l]  <- B[v]  + sigma_b  * b_u[v,l]
        }
      }  
      
      # hypothetical and likelihood estimation
      for (i in 1:N) {
        # Hypothetical / latent levels:
          # universal
        h_univ[i] <- rho1 * h_1[i] + 
          mu  * w[i] + 
          mu1 * w_1[i] + 
          beta
          
          # variant level
        h_var[i] <- P1[var[i]] * h_1[i] +  
          M[var[i]]  * w[i] + 
          M1[var[i]] * w_1[i] + 
          B[var[i]] 
        
        # Likelihood at location level
        h[i] <-  p1[var[i], loc[i]] * h_1[i] + 
          m[var[i], loc[i]]  * w[i] + 
          m1[var[i], loc[i]] * w_1[i] + 
          b[var[i],loc[i]] 
        
        y[i] ~ dnorm(h[i], sd = sigma)
      }
    }
  )
  return(res) 
}

# Construction of model definition
#' @param model.definition Nimble model code
#' @param dat Dataframe used in model
#' @param seed Numeric. Value of chain number 
#' (if unique.init.chain set TRUE)
#' @param unique.init.chain Logical. Flag if each chain
#'  uses a unique set of model inits
build_model <- function(model.definition, 
                        dat,
                        seed = NULL,
                        unique.init.chain = TRUE) {
  
  if(0) model.definition = define_model(dat)
  
  # Initialize all stochastic variables ("nodes")
    
  csts = get_constants(dat)
  N = csts[['N']]
  L = csts[['L']]
  V = csts[['V']]
  
  # Chains initialization
  
  if(unique.init.chain){
    inits.model = create_mdl_inits(dat)[[seed]]
  }
  
  if(!unique.init.chain){
    inits.model = list(
      muu = 1, 
      mu1u = 1, 
      mu_mean = 1,
      mu1_mean = 1,
      # rho_mean = 1,
      # rho1_mean = 1,
      sigma_mu = 1,
      sigma_mu1 = 1,
      sigma_rho1 = 1,
      sigma_beta = 1,
      betau = 1, 
      beta_mean = 0,
      Mu = rep(1, V),
      M1u = rep(1, V),
      # P1u = rep(1, V),
      Bu = rep(0, V),
      sigma_M = 2,  
      sigma_M1 = 2,  
      # sigma_P1 = 2,  
      sigma_B = 1,
      m_u  = matrix(1, nrow = V, ncol = L),
      m1_u = matrix(1, nrow = V, ncol = L),
      # p1_u = matrix(1, nrow = V, ncol = L),
      b_u = matrix(0, nrow = V, ncol = L),
      sigma_m = 1, 
      sigma_m1 = 1, 
      # sigma_p1 = 1, 
      sigma_b = 1,
      sigma = 20
    ) 
  }
  
  # Define data and constants
  
  data.model = list(
    y   = dat$h,
    w   = dat$w,
    w_1 = dat$w_1
  )
  
  constants.model = list(
    N = N,
    V = V,
    L = L,
    var = dat$variant.idx,
    loc = dat$geo.idx
  )
  
  model = nimbleModel(
    code      = model.definition, 
    name      = paste('howareg', seed, sep='_'), 
    constants = constants.model,
    data      = data.model, 
    inits     = inits.model, 
    calculate = FALSE, 
    buildDerivs = TRUE)
  
  return(model)
}

# Compilation of model
#' @param model Nimble model object
#' @param ctrl.nuts HMC-NUTS parameters to modify 
#' control argument to HMC sampler
compile_model_HMC <- function(model, ctrl.nuts, model.type) {
  
  Cmodel = compileNimble(model)
  
  mon = c(
    "mu", "mu_mean", 
    # "rho1", "rho1_mean",
    "beta", "beta_mean",
    "sigma",
    "sigma_mu", 
    # "sigma_rho1",
    "sigma_beta",
    "m",
    "b",
    "sigma_M", 
    # "sigma_P1", 
    "sigma_B",
    "M",
    # "P1",
    "B")
  
  if(model.type == 'lagW') mon = c(mon, 
                                   "mu1", 
                                   "mu1_mean", 
                                   "sigma_mu1",
                                   "sigma_M1", 
                                   "M1")

  HMC = buildHMC(
    model    = model, 
    control  = ctrl.nuts,
    monitors = mon
    )
  
  res = compileNimble(HMC, 
                      project = model,
                      showCompilerOutput = 0)
  return(res)
}

# Execution of one MCMC chain
#' @param seed Numeric (if unique.init.chain = TRUE). Chain number. 
#' Can be NULL if unique.init.chain = FALSE.
#' @param dat Dataframe for model.
#' @param mcmc.params List. Parameters for MCMC execution. 
#' Required parameters
#'  include niter (number of iterations) and burnin
#'  (number of iterations to be burned in)
#' @param unique.init.chain Logical. Flag if each chain uses 
#' a unique set of model inits
#' @param ctrl.nuts HMC-NUTS parameters to modify 
#' control argument to HMC sampler
run_one_parallel <- function(seed, dat, mcmc.params, 
                             unique.init.chain,
                             ctrl.nuts,
                             model.type = 'simple') {
  
  # Needs these library and source
  # calls here so that they are 
  # loaded on each core
  library(nimble)
  library(nimbleHMC)
  source("model-hmc.R")
  source("utils.R")
  
  nburnin     = mcmc.params$burnin
  niter       = mcmc.params$iter
 
  mf = NULL 
  if(model.type == 'simple') mf = define_model(dat)
  if(model.type == 'lagW')   mf = define_model_lagW(dat)
  if(model.type == 'ARlagW') mf = define_model_AR_H_W(dat)
  if(model.type == 'cauchy') mf = define_model_cauchy(dat)

  if(is.null(mf)) stop('Unknown `model.type`: ', model.type)
  
  model = mf |> 
    build_model(dat = dat, seed = seed, 
                unique.init.chain = unique.init.chain) 
  
  chmc = compile_model_HMC(model, ctrl.nuts, model.type)
  
  res = run_inference(chmc = chmc, 
                      nburnin = nburnin, 
                      niter = niter, 
                      nchains = 1) # <-- must be 1!!
  return(res)
}

# Execution of MCMC inference
#' @param chmc Compiled HMC model
#' @param nburnin Numeric. Number of iterations to be burned in
#' @param niter Numeric. Number of iterations
#' @param nchains Numeric. Number of chains
run_inference <- function(chmc, nburnin, niter, nchains){
    
  message(paste("Running MCMC run...",
                "\nburn in    :", nburnin,
                "\niterations :", niter,
                "\nchains     :", nchains, 
                "\n"))
  
  res = nimble::runMCMC(chmc, 
                        nburnin = nburnin, 
                        niter   = niter,
                        nchains = nchains )
  
  return(res)
}

# Processing of output of Hamiltonian MCMC execution
#' @param x MCMC model output
#' @param nchains Numeric. Number of chains
#' @param ci Numeric. Size of confidence interval (0.95 corresponds to 95% CI)
digest_inference <- function(x, nchains, ci = 0.95) {
  # x = mcmcp
  # Reformat output
  if(nchains == 1){
    post.tmp =list(as.data.frame(x))
  }
  
  if(nchains > 1){
    post.tmp = lapply(x, as.data.frame)
  }
  
  post = lapply(post.tmp, mutate, iter = row_number()) |>
    bind_rows(.id = 'chain') |>
    pivot_longer(cols = -c(chain, iter))
  
  # Summary statistics
  postsum = post |>
    group_by(name) |>
    summarise(m = mean(value),
              lo = quantile(value, probs = 0.5 - ci/2),
              hi = quantile(value, probs = 0.5 + ci/2),
              min = min(value),
              max = max(value))
  return(
    list(
      mcmc = x,
      post = post,
      postsum = postsum
    )
  )
}

# Retrieve Gelman-Rubin statistic to test for convergence of chains
#' @param x MCMC model output
run_gelman <- function(x){
  
  # x = mcmcobj$mcmc
  
  message("\nRunning Gelman diagnostic on MCMC object ...",
          appendLF = F)
  mcmc_list   = coda::mcmc.list(lapply(x, as.mcmc))
  gelman_diag = coda::gelman.diag(mcmc_list)
  
  gelman.psrf = gelman_diag[[1]] |>
    as.data.frame() |>
    rownames_to_column(var = "name") |>
    rename(
      psrf = "Point est.",
      psrf.upper = "Upper C.I."
    )
  
  gelman.mpsrf = gelman_diag[[2]]
  message(' done.')
  return(
    list(
      gelman.psrf = gelman.psrf,
      gelman.mpsrf = gelman.mpsrf
    )
  )
}

# Run MCMC on sensitivity data
#' @param sens Dataframe to be used in sensitivity analysis
#' @param iter Numeric. Number of iterations in MCMC run
#' @param burnin Numeric. Number of iterations to be burned in.
#' @param nchains Numeric. Number of MCMC chains
run_sens_mcmc <- function(sens, iter, burnin, nchains){
  t1 = Sys.time()
  message("\n\nRunning sensitivity MCMC runs across model objects...\n\n")
  ctrl.nuts = list(
    delta = 0.90, # default = 0.80 
    maxTreeDepth = 13  # default = 10
  )
  
  mcmc.params = tibble(
    iter = iter,
    burnin = burnin,
    nchains = nchains
  )
  
  mcmcobj = list()
  for(i in seq_along(sens)){
    t = unique(sens[[i]][["threshold"]])
    w = unique(sens[[i]][["window"]])
    p = unique(sens[[i]][["pct.above"]])
    message(paste0("MCMC run for:\nthreshold: ", t,
                   "\nwindow: ", w,
                   "\npct.above: ", p))
    cluster.mcmc = parallel::makeCluster(nchains)
    mcmc = parLapply(cl = cluster.mcmc,
                     X = 1:nchains,
                     fun = run_one_parallel,
                     dat = sens[[i]], 
                     mcmc.params = mcmc.params,
                     unique.init.chain = TRUE,
                     ctrl.nuts = ctrl.nuts)
    stopCluster(cluster.mcmc)
    mcmcobj[[i]] = digest_inference(mcmc, nchains)
    mcmcobj[[i]][["params"]] = list(
      threshold = t,
      window = w,
      pct.above = p
    )
    message("complete.\n\n")
  }
  t2 = Sys.time()
  dtime = difftime(t2, t1, units = "mins")
  message(paste0("\n\nSensitivity MCMC runs completed in ", round(dtime,2), " mins."))
  return(mcmcobj)
}
