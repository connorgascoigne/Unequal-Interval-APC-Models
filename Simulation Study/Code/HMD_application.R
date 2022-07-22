# packages ----

library(demography) # for data
library(tidyverse) # for data wrangling
library(mgcv) # for fitting models
library(patchwork) # for ggplot+ggplot
library(xtable) # r tabels in latex

# path ----

simulationFolder <- 'C:/Users/cg863/OneDrive - University of Bath/Bath PhD/Year 3/Simulation Study/'
# simulationFolder <- '/beegfs/scratch/user/r/cg863/code/Simulation Study/'

## functions ----

source(paste0(simulationFolder, 'Code/functions.R'))

# data import ----

uk <- hmd.mx('GBR_NP', 'xxx ENTER USER NAME xxx', 'xxx ENTER PASSWORD xxx', 'UK')

# rate
ukRateAll <- uk$rate$total
ukRate <- ukRateAll[!rownames(ukRateAll) %in% c(as.character(100:109), '110+'),
                  !colnames(ukRateAll) %in% as.character(c(1922:1925, 2016:2018))]

# population
ukPopAll <- uk$pop$total
ukPop <- ukPopAll[!rownames(ukPopAll) %in% c(as.character(100:109), '110+'),
                !colnames(ukPopAll) %in% as.character(c(1922:1925, 2016:2018))]

# directory for plots ----

# results folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results'))
}

# plots folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Plots'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Plots'))
}

# HMD application folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Plots/HMD Application'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Plots/HMD Application'))
}


# fitting models for different M ----

## data ----

# round data to ensure we have counts
dataA <- aggregate.data(ukRate, ukPop, ageAgg = 1, perAgg = 1) %>% mutate(N = round(N), y = round(y))
dataB <- aggregate.data(ukRate, ukPop, ageAgg = 5, perAgg = 1) %>% mutate(N = round(N), y = round(y))
dataC <- aggregate.data(ukRate, ukPop, ageAgg = 5, perAgg = 5) %>% mutate(N = round(N), y = round(y))

# # tables for paper
# print(xtable(dataA[c(1:3, (nrow(dataA)-2):nrow(dataA)),4:7], type = 'latex'), include.rownames = FALSE)
# print(xtable(dataB[c(1:3, (nrow(dataB)-2):nrow(dataB)),4:7], type = 'latex'), include.rownames = FALSE)
# print(xtable(dataC[c(1:3, (nrow(dataC)-2):nrow(dataC)),4:7], type = 'latex'), include.rownames = FALSE)

## gradient plot of data ----

dataAug <-
  dataA %>%
  mutate(rate = logit(y/N))

trueLexis <- 
  ggplot(dataAug, aes(x = p, y = a, fill = rate)) +
  geom_tile() +
  scale_fill_gradientn(colors = colorRamps::matlab.like2(10000))  +
  labs(x = 'Year', y = 'Age', fill = 'Logit \nMortality \nRate') +
  my.theme(text = element_text(size = 20))

ggsave(trueLexis, 
       file = paste0(simulationFolder, 'Code/Results/Plots/HMD Application/trueLexis1x1.png'),
       width = 10, height = 10)

## model fit ----

# smallest number of unique points to base the knots off of
A <- length(unique(dataC$a))
P <- length(unique(dataC$p))
C <- length(unique(dataC$c))
A;P;C

# age 50% and period and cohort are 25% number of unique data points for smallest 
# partition of data.
# Keep knots constant for comparison
knots <- list(age = 10,period = 10,cohort = 20)
# want to penalise to fixed = F
fixed <- list(age = FALSE, period = FALSE, cohort = FALSE)

fitA <- hol.spline.fun(data = dataA, mod = 'apc', knots = knots, fixed = fixed, slopeDrop = 'c', family = 'binomial')
fitB <- hol.spline.fun(data = dataB, mod = 'apc', knots = knots, fixed = fixed, slopeDrop = 'c', family = 'binomial')
fitC <- hol.spline.fun(data = dataC, mod = 'apc', knots = knots, fixed = fixed, slopeDrop = 'c', family = 'binomial')

modA <- fitA$mod
modB <- fitB$mod
modC <- fitC$mod

# getting plots for smooth terms for all models ----

predA <- plot.terms(gamObj = modA, data = dataA)
predB <- plot.terms(gamObj = modB, data = dataB)
predC <- plot.terms(gamObj = modC, data = dataC)

# pred.terms.1$plot + pred.terms.2$plot + pred.terms.3$plot

all.mod.term <- 
  rbind(predA$termsPred %>% mutate(M = '1X1'),
        predB$termsPred %>% mutate( M = '5X1'),
        predC$termsPred %>% mutate( M = '5X5')) %>% 
  mutate(Effect = factor(Effect,
                       levels = c('s(a)', 's(p)', 's(c)'),
                       labels = c(latex2exp::TeX(r'($\hat{f}_{A_C}(a)$)'), 
                                latex2exp::TeX(r'($\hat{f}_{P_C}(p)$)'), 
                                latex2exp::TeX(r'($\hat{f}_{C_C}(c)$)'))))

smoothCurvaturePlot <- 
  ggplot(all.mod.term, aes(x = Index, y = Estimate, group = M, col = M))+
  geom_line()+
  geom_point(aes(shape = M))+
  facet_wrap(~Effect, scales = 'free', labeller = label_parsed)+
  labs(x = 'Years',
       y = 'Estimate')+
  my.theme(text = element_text(size = 20),
           legend.title=element_blank())

smoothCurvaturePlot

ggsave(smoothCurvaturePlot, 
       file = paste0(simulationFolder, 'Code/Results/Plots/HMD Application/smoothCurvaturePlot.png'),
       width = 10, height = 10)

# expected response plots ----

yPredA <- my.predict(gamObj = fitA$mod, type = 'response')
yPredB <- my.predict(gamObj = fitB$mod, type = 'response')
yPredC <- my.predict(gamObj = fitC$mod, type = 'response')

dataPredA <-
  dataA %>%
  mutate(rateHat = yPredA$est)
dataPredB <-
  dataB %>%
  mutate(rateHat = yPredB$est)
dataPredC <-
  dataC %>%
  mutate(rateHat = yPredC$est)

predLexisA <- 
  ggplot2::ggplot(dataPredA, aes(x = p, y = a, fill = rateHat)) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradientn(colors = colorRamps::matlab.like2(10000))  +
  ggplot2::labs(x = 'Year', y = 'Age', fill = 'Predicted \nLogit \nMortality \nRate') +
  my.theme(text = element_text(size = 20))

predLexisB <- 
  ggplot2::ggplot(dataPredB, aes(x = p, y = a, fill = rateHat)) +
  ggplot2::geom_tile(height = 5) +
  ggplot2::scale_fill_gradientn(colors = colorRamps::matlab.like2(10000))  +
  ggplot2::labs(x = 'Year', y = 'Age', fill = 'Predicted \nLogit \nMortality \nRate') +
  my.theme(text = element_text(size = 20))

predLexisC <- 
  ggplot2::ggplot(dataPredC, aes(x = p, y = a, fill = rateHat)) +
  ggplot2::geom_tile(width = 5, height = 5) +
  # ggplot2::geom_tile() +
  ggplot2::scale_fill_gradientn(colors = colorRamps::matlab.like2(10000))  +
  ggplot2::labs(x = 'Year', y = 'Age', fill = 'Predicted \nLogit \nMortality \nRate') +
  my.theme(text = element_text(size = 20))

(trueLexis + predLexisA) / (predLexisB + predLexisC)


ggsave(predLexisA, 
       file = paste0(simulationFolder, 'Code/Results/Plots/HMD Application/predLexis1x1.png'),
       width = 10, height = 10)
ggsave(predLexisB, 
       file = paste0(simulationFolder, 'Code/Results/Plots/HMD Application/predLexis5x1.png'),
       width = 10, height = 10)
ggsave(predLexisC, 
       file = paste0(simulationFolder, 'Code/Results/Plots/HMD Application/predLexis5x5.png'),
       width = 10, height = 10)
