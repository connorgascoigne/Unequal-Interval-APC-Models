# packages ---------------------------------------------------------------------

library(demography) # for data
library(tidyverse) # for data wrangling
library(mgcv) # for fitting models
library(patchwork) # for ggplot+ggplot
library(xtable) # r tabels in latex

# source functions -------------------------------------------------------------

source('holford_models.R') # for model fit
source('analysis.R') # for my.theme()

# functions --------------------------------------------------------------------

# aggregates data by finding an age and period bin using a and p group only
# then something like agg(~ ageBin+periodBin)
# then finds c = p-a where p and 1 are midpoints of bins
aggregate.data <- function(rate, population, ageAgg, perAgg){
  
  # collecting rate and population into one df and finding counts
  p <- as.data.frame(as.table(population)); colnames(p) <- c('age', 'year', 'N')
  r <- as.data.frame(as.table(rate)); colnames(r) <- c('age', 'year', 'R')
  data <- left_join(p, r, by = c('age', 'year'), keep = FALSE)
  data$age <- as.numeric(levels(data$age))[data$age]
  data$year <- as.numeric(levels(data$year))[data$year]
  data$y <- data$N*data$R
  
  # define a sequence where last element is included in from:to !%% by
  my.seq <- function (from, to, by){
    vec <- do.call(what = seq, args = list(from, to, by))
    if ( tail(vec, 1) !=  to ) {
      return(c(vec, to))
    } else {
      return(vec)
    }
  }
  
  # bin observations from my.seq
  my.cut <- function(x, from, to, by){
    brks <- my.seq(from = from, to = to, by = by)
    x = cut(x = x, breaks = brks, right = F, include.lowest = T, dig.lab = 4)
    x
  }
  
  # vector of the mid points from my.seq
  mid.seq <- function(x, from, to, by){
    brks <- my.seq(from = from, to = to, by = by)
    x = cut(x = x, breaks = brks, include.lowest = T, dig.lab = 4)
    vec = my.seq(from = from+by/2, to = to-by/2, by = by)
    names(vec) = levels(x)
    vec
  }
  
  # define the bins and add column to data
  data$ageBin <- my.cut(x = data$age, from = min(data$age), to = max(data$age), by = ageAgg)
  data$yearBin <- my.cut(x = data$year, from = min(data$year), to = max(data$year), by = perAgg)
  
  # define midpoints of the bins
  a <- mid.seq(data$age, from = min(data$age), to = max(data$age), by = ageAgg)
  p <- mid.seq(data$year, from = min(data$year), to = max(data$year), by = perAgg)
  
  # aggregate data over the bin and add vector of midpoints
  agg <- aggregate(cbind(N,y) ~ ageBin + yearBin, data = data, sum, na.rm = T) %>% arrange(ageBin, yearBin)
  agg$a <- a[as.numeric(agg$ageBin)]
  agg$p <- p[as.numeric(agg$yearBin)]
  agg$c <- agg$p - agg$a
  
  dataAgg <- agg %>% relocate(c('a', 'p', 'c'), .before = ageBin) %>% arrange(a, p, c)
  
  dataAgg
}

# make a gradient plot
gradient.plot <- function(x, is.log = F){
  
  x1 <- as.data.frame(as.table(x))
  colnames(x1) <- c('age', 'year', 'value')
  x1$age <- as.numeric(levels(x1$age))[x1$age]
  x1$year <- as.numeric(levels(x1$year))[x1$year]
  
  if(is.log == T){
    x1$value <- log(x1$value)
  }
  
  ggplot(x1, aes(x = year, y = age, fill = value))+ 
    geom_raster(hjust = 0.5, vjust = 0.5, interpolate = FALSE)+
    scale_fill_gradientn(colors = colorRamps::matlab.like2(10000), name = 'Log \nMortality \nRate')+
    labs(x = 'Year', y = 'Age')+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          # legend.title = element_blank(),
          legend.text.align = 0,
          legend.key = element_rect(fill = NA))
}

# takes and model and plots the smooth functions (identifiable in APC model)
plot.terms <- function(gamObj, data){
  
  pred <- my.predict(gamObj = gamObj,type = 'terms')
  pred.est <- pred$est
  pred.se <- pred$se
  
  n.pterms <- length(attr(gamObj$pterms, 'term.labels'))
  n.smooth <- length(gamObj$smooth)
  names <- colnames(pred.est)
  
  allTermsPred <- c()
  effect.levels <- c()
  for(i in 1:n.smooth){
    # labelling
    label <- gamObj$smooth[[i]]$label
    term <- gamObj$smooth[[i]]$term
    
    # dataframe for estimate and 95% CI
    all <- data.frame(Effect = label, Index = data[,term], Estimate = pred.est[,label])
    all$lower <- all$Estimate + stats::qnorm(0.025)*pred.se[,label]
    all$upper <- all$Estimate + stats::qnorm(0.975)*pred.se[,label]
    
    # unique values
    final <- aggregate(cbind(Estimate, lower, upper) ~ Effect + Index, FUN = mean, data = all)
    
    allTermsPred <- rbind(allTermsPred, final)
    effect.levels <- c(effect.levels, label)
  }
  
  termsPred <- 
    allTermsPred %>% 
    mutate(Effect = factor(Effect, levels = unique(effect.levels, fromLast = TRUE)))
  
  plot <- 
    ggplot(termsPred, aes(x = Index, y = Estimate, group = Effect, col = Effect)) + 
    geom_line() + 
    geom_point() + 
    facet_grid(~Effect,scales = 'free') + 
    geom_ribbon(aes(ymin = lower,ymax = upper,fill = Effect), alpha = .2, linetype = 0, show.legend = F) + 
    my.theme()
  
  return(list(plot = plot, termsPred = termsPred)) 
}


# data import ------------------------------------------------------------------

rm(list = setdiff(ls(), lsf.str()))

uk <- hmd.mx('GBR_NP', 'xxx ENTER USER NAME xxx', 'xxx ENTER PASSWORD xxx', 'UK')

# rate
ukRateAll <- uk$rate$total
ukRate <- ukRateAll[!rownames(ukRateAll) %in% c(as.character(100:109), '110+'),
                  !colnames(ukRateAll) %in% as.character(c(1922:1925, 2016:2018))]

# population
ukPopAll <- uk$pop$total
ukPop <- ukPopAll[!rownames(ukPopAll) %in% c(as.character(100:109), '110+'),
                !colnames(ukPopAll) %in% as.character(c(1922:1925, 2016:2018))]


# exploratory plots ------------------------------------------------------------

# gradient plot
p1 <- gradient.plot(x = ukRate,is.log = T);p1

# fitting models for different M -----------------------------------------------

# round data to ensure we have counts
dataA <- aggregate.data(ukRate, ukPop, ageAgg = 1, perAgg = 1) %>% mutate(N = round(N), y = round(y))
dataB <- aggregate.data(ukRate, ukPop, ageAgg = 5, perAgg = 1) %>% mutate(N = round(N), y = round(y))
dataC <- aggregate.data(ukRate, ukPop, ageAgg = 5, perAgg = 5) %>% mutate(N = round(N), y = round(y))

# # tables for paper
# print(xtable(dataA[c(1:3, (nrow(dataA)-2):nrow(dataA)),4:7], type = 'latex'), include.rownames = FALSE)
# print(xtable(dataB[c(1:3, (nrow(dataB)-2):nrow(dataB)),4:7], type = 'latex'), include.rownames = FALSE)
# print(xtable(dataC[c(1:3, (nrow(dataC)-2):nrow(dataC)),4:7], type = 'latex'), include.rownames = FALSE)

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

fitA <- hol.spline.fun(data = dataA, mod = 'apc', knots = knots, fixed = fixed, colDrop = 'c', distribution = 'binomial')
fitB <- hol.spline.fun(data = dataB, mod = 'apc', knots = knots, fixed = fixed, colDrop = 'c', distribution = 'binomial')
fitC <- hol.spline.fun(data = dataC, mod = 'apc', knots = knots, fixed = fixed, colDrop = 'c', distribution = 'binomial')

modA <- fitA$mod
modB <- fitB$mod
modC <- fitC$mod

# getting plots for smooth terms for all models --------------------------------

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

p1 <- 
  ggplot(all.mod.term, aes(x = Index, y = Estimate, group = M, col = M))+
  geom_line()+
  geom_point(aes(shape = M))+
  facet_wrap(~Effect,scales = 'free',labeller = label_parsed)+
  labs(x = 'Years',
       y = 'Estimate')+
  my.theme();p1
