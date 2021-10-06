# packages ---------------------------------------------------------------------
library(tidyverse)
library(mgcv)
library(furrr)
library(patchwork)
library(zoo)
library(doParallel)

# loading in functions from other files ----------------------------------------

# MAKE SURE TO SET WD
# MAKE SURE TO SET BOTH FILE PATHS IN run.holford.sim() 

source('data_generation.R')
source('holford_models.R')
source('analysis.R')

## make sure in the run.holford.sim() function the file source() file paths are correct!!!

# simulation study example  ----------------------------------------------------

# clear enviroment of all but functions
rm(list = setdiff(ls(), lsf.str()))

# parameters for knots and different M
modParams <- list(M = 1, knots = list(age = 15, period = 5, cohort = 20)) # equal aggregtaion results
modParams <- list(M = 5, knots = list(age = 5, period = 5, cohort = 8)) # unequal aggretation results

# fixed parameters for all distributions
fixedParams <- list(A = 60, P = 20, N = 150, nSim = 100, mod = "apc", scenario = c("apc", "pc", 'ac', 'ap'))

# specifc distribution parameters
distParams <- list(distribution = "binomial", FUNage = age.fun.binomial, FUNperiod = per.fun.binomial, FUNcohort = coh.fun.binomial) # binomial results
# distParams <- list(distribution = "normal", FUNage = age.fun.normal, FUNperiod = per.fun.normal, FUNcohort = coh.fun.normal) # normal results
# distParams <- list(distribution = "poisson", FUNage = age.fun.poisson, FUNperiod = per.fun.poisson, FUNcohort = coh.fun.poisson) # poisson results

# run the simulation with a timer
ptm <- proc.time() # Start the clock!
simRun <- run.holford.sim(fixedParams = fixedParams, modParams = modParams, distParams = distParams)
proc.time()-ptm # Stop the clock

# extract results for each scenario 
res <- get.final.plot(simRun, 'apc')

# show the plots
res$p
res$simPlots[[1]]
res$simPlots[[2]]
res$simPlots[[3]]
res$individualEffectPlots[[1]]
res$individualEffectPlots[[2]]
res$individualEffectPlots[[3]]

