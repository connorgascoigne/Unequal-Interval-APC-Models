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

# collecting results -----

## model fits ----

nSims <- 100

### defining results dataframe ----

# can use one of the model fits to define results data frame
## all models fit to the same number (and explcit groups) of age, period and cohort groups
## all models fit to the same number of data sets
temp <- readRDS(file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/FA/fit_1.rds'))
resultsColNames <- c('Effect', 'Group', 'Model', c(sprintf("s%05d",1:nSims)))

# dimensions
## estimate
facEstRes <- rssEstRes <- pssEstRes <- 
  rssAddKnotEstRes <- pssAddKnotEstRes <- rssPeriodicEstRes <- pssPeriodicEstRes <- 
  ## curvature
  facCurvRes <- rssCurvRes <- pssCurvRes <- 
  rssAddKnotCurvRes <- pssAddKnotCurvRes <- rssPeriodicCurvRes <- pssPeriodicCurvRes <- 
  ## temp
  array(NA, dim = c(nrow(temp$estData), ncol(temp$estData)+nSims)) %>% as.data.frame()

# Effect column
## estimate
facEstRes[,1] <- rssEstRes[,1] <- pssEstRes[,1] <- 
  rssAddKnotEstRes[,1] <- pssAddKnotEstRes[,1] <- rssPeriodicEstRes[,1] <- pssPeriodicEstRes[,1] <- 
  ## curvature
  facCurvRes[,1] <- rssCurvRes[,1] <- pssCurvRes[,1] <- 
  rssAddKnotCurvRes[,1] <- pssAddKnotCurvRes[,1] <- rssPeriodicCurvRes[,1] <- pssPeriodicCurvRes[,1] <- 
  ## temp
  temp$estData$Effect


# Group column
## estimate
facEstRes[,2] <- rssEstRes[,2] <- pssEstRes[,2] <- 
  rssAddKnotEstRes[,2] <- pssAddKnotEstRes[,2] <- rssPeriodicEstRes[,2] <- pssPeriodicEstRes[,2] <- 
  ## curvature
  facCurvRes[,2] <- rssCurvRes[,2] <- pssCurvRes[,2] <- 
  rssAddKnotCurvRes[,2] <- pssAddKnotCurvRes[,2] <- rssPeriodicCurvRes[,2] <- pssPeriodicCurvRes[,2] <- 
  ## temp
  temp$estData$Group


## column names
## esimate
colnames(facEstRes) <- colnames(rssEstRes) <- colnames(pssEstRes) <- 
  colnames(rssAddKnotEstRes) <- colnames(pssAddKnotEstRes) <- colnames(rssPeriodicEstRes) <- colnames(pssPeriodicEstRes) <-
  ## curvature
  colnames(facCurvRes) <- colnames(rssCurvRes) <- colnames(pssCurvRes) <- 
  colnames(rssAddKnotCurvRes) <- colnames(pssAddKnotCurvRes) <- colnames(rssPeriodicCurvRes) <- colnames(pssPeriodicCurvRes) <-
  # temp
  resultsColNames

## model type
facEstRes[,3] <- facCurvRes[,3] <- 'FA'
rssEstRes[,3] <- rssCurvRes[,3] <- 'RSS'
pssEstRes[,3] <- pssCurvRes[,3] <- 'PSS'
rssAddKnotEstRes[,3] <- rssAddKnotCurvRes[,3] <- 'RSS'
pssAddKnotEstRes[,3] <- pssAddKnotCurvRes[,3] <- 'PSS'
rssPeriodicEstRes[,3] <- rssPeriodicCurvRes[,3] <- 'RSS'
pssPeriodicEstRes[,3] <- pssPeriodicCurvRes[,3] <- 'PSS'

### collecting all fits into one data frame ----

for(i in 1:nSims){
  # i <- 1
  
  # read in result fits
  facFit <- readRDS(file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/FA/fit_', i, '.rds'))
  rssFit <- readRDS( file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/RSS/fit_', i, '.rds'))
  pssFit <- readRDS(file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/PSS/fit_', i, '.rds'))
  rssAddKnotFit <- readRDS(file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/RSS Additonal Knots/fit_', i, '.rds'))
  pssAddKnotFit <- readRDS(file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/PSS Additonal Knots/fit_', i, '.rds'))
  rssPeriodicFit <- readRDS(file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/RSS Periodic/fit_', i, '.rds'))
  pssPeriodicFit <- readRDS(file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/PSS Periodic/fit_', i, '.rds'))
  
  # update the effect estimates
  facEstRes[,3+i] <- facFit$estData$Estimate
  rssEstRes[,3+i] <- rssFit$estData$Estimate
  pssEstRes[,3+i] <- pssFit$estData$Estimate
  rssAddKnotEstRes[,3+i] <- rssAddKnotFit$estData$Estimate
  pssAddKnotEstRes[,3+i] <- pssAddKnotFit$estData$Estimate
  rssPeriodicEstRes[,3+i] <- rssPeriodicFit$estData$Estimate
  pssPeriodicEstRes[,3+i] <- pssPeriodicFit$estData$Estimate
  
  # update the curvature estimates
  facCurvRes[,3+i] <- facFit$curvData$Estimate
  rssCurvRes[,3+i] <- rssFit$curvData$Estimate
  pssCurvRes[,3+i] <- pssFit$curvData$Estimate
  rssAddKnotCurvRes[,3+i] <- rssAddKnotFit$curvData$Estimate
  pssAddKnotCurvRes[,3+i] <- pssAddKnotFit$curvData$Estimate
  rssPeriodicCurvRes[,3+i] <- rssPeriodicFit$curvData$Estimate
  pssPeriodicCurvRes[,3+i] <- pssPeriodicFit$curvData$Estimate
  
}

### true values ----

trueVals <- readRDS(file = paste0(simulationFolder, 'Code/Results/', family, '/', dataMod, '/', M, '/True/fit.rds'))

### all results together -----

normalBasisResults <- 
  list(simEstData = rbind(facEstRes, rssEstRes, pssEstRes),
       simCurvData = rbind(facCurvRes, rssCurvRes, pssCurvRes),
       trueEffData = trueVals$trueEffData,
       trueCurvData = trueVals$trueCurvData)

additionalKnotsBasisResults <- 
  list(simEstData = rbind(rssAddKnotEstRes, pssAddKnotEstRes),
       simCurvData = rbind(rssAddKnotCurvRes, pssAddKnotCurvRes),
       trueEffData = trueVals$trueEffData,
       trueCurvData = trueVals$trueCurvData)

periodicBasisResults <- 
  list(simEstData = rbind(rssPeriodicEstRes, pssPeriodicEstRes),
       simCurvData = rbind(rssPeriodicCurvRes, pssPeriodicCurvRes),
       trueEffData = trueVals$trueEffData,
       trueCurvData = trueVals$trueCurvData)

### saving ----

#### creating directories -----

# results folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results'))
}

# Final results folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Final'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Final'))
}

# family folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Final/', family))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Final/', family))
}

# data model folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Final/', family, '/', dataMod))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Final/', family, '/', dataMod))
}

# M folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Final/', family, '/', dataMod, '/', M))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Final/', family, '/', dataMod, '/', M))
}


#### save final ----

saveRDS(normalBasisResults , file = paste0(simulationFolder, 'Code/Results/Final/', family, '/', dataMod, '/', M, '/normalBasisResults', str_to_title(family), toupper(dataMod), M, '.rds'))
saveRDS(additionalKnotsBasisResults , file = paste0(simulationFolder, 'Code/Results/Final/', family, '/', dataMod, '/', M, '/additionalKnotsBasisResults', str_to_title(family), toupper(dataMod), M, '.rds'))
saveRDS(periodicBasisResults , file = paste0(simulationFolder, 'Code/Results/Final/', family, '/', dataMod, '/', M, '/periodicBasisResults', str_to_title(family), toupper(dataMod), M, '.rds'))
