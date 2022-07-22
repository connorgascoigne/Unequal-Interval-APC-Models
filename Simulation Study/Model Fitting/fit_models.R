# packages ----

library(tidyverse)
library(mgcv)

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
  i <- 1
} else {
  family <- c('normal', 'poisson', 'binomial')[as.numeric(args[1])]
  dataMod <- c('apc', 'ap', 'ac', 'a')[as.numeric(args[2])]
  M <- c(1, 3, 5)[as.numeric(args[3])]
  i <- as.numeric(args[4])
}

## dependent on command line arguments ----

### data load ----

dataList <- readRDS(file = paste0(simulationFolder, 'Data/', family, '/', dataMod, '/', M, '/allData.rds'))

### picking only one data set ----

data <- dataList[[i]]

### spline parameters for different M ----

if(M == 1){
  knots <- list(age = 15, period = 5, cohort = 20)
} else if(M  ==  3){
  knots <- list(age = 5, period = 5, cohort = 10)
} else if(M == 5){
  knots <- list(age = 5, period = 5, cohort = 8)
}

# on less than the distinct number of data points for each 
addKnots <- 
  list(age = (data$a %>% unique %>% length) - 1,
       period = (data$p %>% unique %>% length) - 1,
       cohort = (data$c %>% unique %>% length) - 1)

# model fitting ----

## model fit ----

### normal simulation study ----

facFit <-
  hol.factor.extract(data = data, family = family, # dependent on command line
                     # fixed arguments
                     slopeDrop = 'c', mod = 'apc')

rssFit <-
  hol.spline.extract(data = data, family = family, knots = knots, # dependent on command line
                     # fixed arguments
                     bs = 'cr', include.periodic = FALSE, 
                     slopeDrop = 'c', mod = 'apc',
                     fixed = list(age = TRUE, period = TRUE, cohort = TRUE))

pssFit <-
  hol.spline.extract(data = data, family = family, knots = knots, # dependent on command line
                     # fixed arguments
                     bs = 'cr', include.periodic = FALSE, 
                     slopeDrop = 'c', mod = 'apc',
                     fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

### unnatural basis ----

#### additional knots ----

rssAddKnotFit <-
  hol.spline.extract(data = data, family = family, knots = addKnots, # dependent on command line
                     # fixed arguments
                     bs = 'cr', include.periodic = FALSE, 
                     slopeDrop = 'c', mod = 'apc',
                     fixed = list(age = TRUE, period = TRUE, cohort = TRUE))

pssAddKnotFit <-
  hol.spline.extract(data = data, family = family, knots = addKnots, # dependent on command line
                     # fixed arguments
                     bs = 'cr', include.periodic = FALSE, 
                     slopeDrop = 'c', mod = 'apc',
                     fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

#### periodic ----

rssPeriodicFit <- 
  hol.spline.extract(data = data, family = family, knots = knots, # dependent on command line
                     # fixed arguments
                     bs = 'cr', include.periodic = TRUE, 
                     slopeDrop = 'c', mod = 'apc',
                     fixed = list(age = TRUE, period = TRUE, cohort = TRUE))


pssPeriodicFit <- 
  hol.spline.extract(data = data, family = family, knots = knots, # dependent on command line
                     # fixed arguments
                     bs = 'cr', include.periodic = TRUE, 
                     slopeDrop = 'c', mod = 'apc',
                     fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

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
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/FA'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/FA'))
}

# RSS model folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/RSS'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/RSS'))
}

# PSS model folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/PSS'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/PSS'))
}

# RSS additional knots model folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/RSS Additonal Knots'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/RSS Additonal Knots'))
}

# PSS additional knots model folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/PSS Additonal Knots'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/PSS Additonal Knots'))
}

# RSS periodic model folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/RSS Periodic'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/RSS Periodic'))
}

# PSS periodic model folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/PSS Periodic'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/PSS Periodic'))
}

### model saves ----

saveRDS(facFit , file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/FA/fit_', i, '.rds'))
saveRDS(rssFit , file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/RSS/fit_', i, '.rds'))
saveRDS(pssFit , file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/PSS/fit_', i, '.rds'))
saveRDS(rssAddKnotFit , file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/RSS Additonal Knots/fit_', i, '.rds'))
saveRDS(pssAddKnotFit , file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/PSS Additonal Knots/fit_', i, '.rds'))
saveRDS(rssPeriodicFit , file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/RSS Periodic/fit_', i, '.rds'))
saveRDS(pssPeriodicFit , file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/PSS Periodic/fit_', i, '.rds'))
