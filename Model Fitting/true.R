# packages ----

library(tidyverse)

# Simulation Study  ----

## Folder to all simulation study ----

# simulationFolder <- 'C:/Users/cg863/OneDrive - University of Bath/Bath PhD/Year 3/Simulation Study/'
simulationFolder <- '/beegfs/scratch/user/r/cg863/code/Simulation Study/'

## import functions ----

source(paste0(simulationFolder, 'Code/functions.R'))

# command line ---

## actual arguments ----

## First read in the arguments listed at the command line
args <- commandArgs(trailingOnly = TRUE)

## args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args) == 0){
  print("WARNING: No arguments supplied.")
  family <- c('normal', 'poisson', 'binomial')[2]
  dataMod <- c('apc', 'ap', 'ac', 'a')[1]
  M <- c(1, 3, 5)[1]
} else {
  family <- c('normal', 'poisson', 'binomial')[as.numeric(args[1])]
  dataMod <- c('apc', 'ap', 'ac', 'a')[as.numeric(args[2])]
  M <- c(1, 3, 5)[as.numeric(args[3])]
}

## dependent on command line arguments ----

### what data generating functions to use

if(family == 'normal'){
  family = 'normal'; FUNage = age.fun.normal;  FUNperiod = per.fun.normal; FUNcohort = coh.fun.normal
} else if(family == 'poisson'){
  FUNage = age.fun.poisson;  FUNperiod = per.fun.poisson; FUNcohort = coh.fun.poisson
} else if(family == 'binomial'){
  FUNage = age.fun.binomial;  FUNperiod = per.fun.binomial; FUNcohort = coh.fun.binomial
}

# model fitting ----

## fixed arguments ----

A <- 60; P <- 20

## model fit ----

trueFit <- true.values(A = A, P = P, M = M, mod = dataMod, FUNage = FUNage, FUNperiod = FUNperiod, FUNcohort = FUNcohort)

## saving ----

### create directory (if needed) ----

# results folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results'))
}

# family folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/', family))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/', family))
}

# data model folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod))
}

# M folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M))
}

# factor model folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/True'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/True'))
}

### model saves ----

saveRDS(trueFit , file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/True/fit.rds'))
