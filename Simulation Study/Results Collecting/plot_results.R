# packages ----

library(tidyverse)

# Simulation Study  ----

## Folder to all simulation study ----

simulationFolder <- 'C:/Users/cg863/OneDrive - University of Bath/Bath PhD/Year 3/Simulation Study/'

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
  family <- c('normal', 'poisson', 'binomial')[3]
  dataMod <- c('apc', 'ap', 'ac', 'a')[1]
  M <- c(1, 3, 5)[3]
} else {
  family <- c('normal', 'poisson', 'binomial')[as.numeric(args[1])]
  dataMod <- c('apc', 'ap', 'ac', 'a')[as.numeric(args[2])]
  M <- c(1, 3, 5)[as.numeric(args[3])]
}

# load final results ----

normalBasisResults <- readRDS(file = paste0(simulationFolder, 'Code/Results/Final/', family, '/', dataMod, '/', M, '/normalBasisResults', str_to_title(family), toupper(dataMod), M, '.rds'))
additionalKnotsBasisResults <- readRDS(file = paste0(simulationFolder, 'Code/Results/Final/', family, '/', dataMod, '/', M, '/additionalKnotsBasisResults', str_to_title(family), toupper(dataMod), M, '.rds'))
periodicBasisResults <- readRDS(file = paste0(simulationFolder, 'Code/Results/Final/', family, '/', dataMod, '/', M, '/periodicBasisResults', str_to_title(family), toupper(dataMod), M, '.rds'))

# plots ----

## creating directories (if needed) ----

# results folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results'))
}

# plots folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Plots'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Plots'))
}

# family folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Plots/', family))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Plots/', family))
}

# data model folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod))
}

# M folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M))
}

# final results folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Final Results'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Final Results'))
}

# indvidual simulations folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Simulations'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Simulations'))
}

# indvidual effect folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Effect'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Effect'))
}

## making the plots ----

normalBasisPlots <- get.final.plot(normalBasisResults, text = element_text(size = 20))
additionalKnotsBasisPlots <- get.final.plot(additionalKnotsBasisResults, text = element_text(size = 20), boxplots = FALSE)
periodicBasisPlots <- get.final.plot(periodicBasisResults, text = element_text(size = 20), boxplots = FALSE)

width = height = 10

## saving ----

# final plot
## normal basis
ggsave(normalBasisPlots$p, 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Final Results/normalBasis', str_to_title(family), toupper(dataMod), M, 'Final.png'),
       width = width, height = height)
## additional knots basis
ggsave(additionalKnotsBasisPlots$p, 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Final Results/additionalKnotsBasis', str_to_title(family), toupper(dataMod), M, 'Final.png'),
       width = width, height = height)
## periodic basis
ggsave(periodicBasisPlots$p, 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Final Results/periodicBasis', str_to_title(family), toupper(dataMod), M, 'Final.png'),
       width = width, height = height)

# individual simualtions plots
## normal basis
ggsave(normalBasisPlots$simPlots[[1]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Simulations/normalBasis', str_to_title(family), toupper(dataMod), M, 'FA.png'),
       width = width, height = height)
ggsave(normalBasisPlots$simPlots[[2]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Simulations/normalBasis', str_to_title(family), toupper(dataMod), M, 'RSS.png'),
       width = width, height = height)
ggsave(normalBasisPlots$simPlots[[3]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Simulations/normalBasis', str_to_title(family), toupper(dataMod), M, 'PSS.png'),
       width = width, height = height)
## additional knots basis
ggsave(additionalKnotsBasisPlots$simPlots[[1]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Simulations/additionalKnotsBasis', str_to_title(family), toupper(dataMod), M, 'RSS.png'),
       width = width, height = height)
ggsave(additionalKnotsBasisPlots$simPlots[[2]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Simulations/additionalKnotsBasis', str_to_title(family), toupper(dataMod), M, 'PSS.png'),
       width = width, height = height)
## periodic basis
ggsave(periodicBasisPlots$simPlots[[1]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Simulations/periodicBasis', str_to_title(family), toupper(dataMod), M, 'RSS.png'),
       width = width, height = height)
ggsave(periodicBasisPlots$simPlots[[2]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Simulations/periodicBasis', str_to_title(family), toupper(dataMod), M, 'PSS.png'),
       width = width, height = height)

# individual effects plots
## normal basis
ggsave(normalBasisPlots$individualEffectPlots[[1]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Effect/normalBasis', str_to_title(family), toupper(dataMod), M, 'Age.png'),
       width = width, height = height)
ggsave(normalBasisPlots$individualEffectPlots[[2]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Effect/normalBasis', str_to_title(family), toupper(dataMod), M, 'Period.png'),
       width = width, height = height)
ggsave(normalBasisPlots$individualEffectPlots[[3]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Effect/normalBasis', str_to_title(family), toupper(dataMod), M, 'Cohort.png'),
       width = width, height = height)
## additional knots basis
ggsave(additionalKnotsBasisPlots$individualEffectPlots[[1]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Effect/additionalKnotsBasis', str_to_title(family), toupper(dataMod), M, 'Age.png'),
       width = width, height = height)
ggsave(additionalKnotsBasisPlots$individualEffectPlots[[2]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Effect/additionalKnotsBasis', str_to_title(family), toupper(dataMod), M, 'Period.png'),
       width = width, height = height)
ggsave(additionalKnotsBasisPlots$individualEffectPlots[[3]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Effect/additionalKnotsBasis', str_to_title(family), toupper(dataMod), M, 'Cohort.png'),
       width = width, height = height)
## periodic basis
ggsave(periodicBasisPlots$individualEffectPlots[[1]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Effect/periodicBasis', str_to_title(family), toupper(dataMod), M, 'Age.png'),
       width = width, height = height)
ggsave(periodicBasisPlots$individualEffectPlots[[2]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Effect/periodicBasis', str_to_title(family), toupper(dataMod), M, 'Period.png'),
       width = width, height = height)
ggsave(periodicBasisPlots$individualEffectPlots[[3]], 
       file = paste0(simulationFolder, 'Code/Results/Plots/', family, '/', dataMod, '/', M, '/Individual Effect/periodicBasis', str_to_title(family), toupper(dataMod), M, 'Cohort.png'),
       width = width, height = height)
