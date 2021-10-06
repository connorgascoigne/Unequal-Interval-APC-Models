# data generation functions ----------------------------------------------------

## true functions --------------------------------------------------------------

### normal ----------------------------------------------------------------------
age.fun.normal <- function(a){
  0.3*a - 0.01*(a^2)
}

per.fun.normal <- function(p){
  -0.04*p + 0.02*(p^2)
}

coh.fun.normal <- function(c){
  0.35*c - 0.0015*(c^2)
}

### binomial --------------------------------------------------------------------
age.fun.binomial <- function(a){
  0.4 + (1/50)*(0.3*a - 0.01*(a^2))
}

per.fun.binomial <- function(p){
  (1/50)*(-0.04*p + 0.02*(p^2))
}

coh.fun.binomial <- function(c){
  (1/50)*(0.35*c - 0.0015*(c^2))
}

### poisson ---------------------------------------------------------------------
age.fun.poisson <- function(a){
  -1.5 + (1/50)*(0.3*a - 0.01*(a^2))
}

per.fun.poisson <- function(p){
  (1/50)*(-0.04*p + 0.02*(p^2))
}

coh.fun.poisson <- function(c){
  (1/50)*(0.35*c - 0.0015*(c^2))
}

## data frame simulation -------------------------------------------------------

data.sim <- function(A, P, M, N, distribution = c('normal', 'poisson', 'binomial'),
                   mod = c('apc','ac','ap','pc','a','p','c'),
                   FUNage, FUNperiod, FUNcohort){
  
  # a and p continuous and then find the midpoints between relative parts
  # i.e. if a= 1, 2, 3,.. aMid = 0.5, 1.5, 2.5,...
  a <- seq(from = 0,to = A,by = 1); aMid <- zoo::rollapply(a, width = 1 + 1, by = 1, mean, align = 'left') # want to average over 1 group in a
  p <- seq(from = 0,to = P,by = 1); pMid <- zoo::rollapply(p, width = 1 + 1, by = 1, mean, align = 'left') # want to average over 1 group in p
  data <- data.frame(a = rep(aMid, each = P),
                   p = rep(pMid, times = A))
  data$c <- data$p-data$a # cohort found on single years with c = p-a
  
  # generate the expected value of y
  # notice that input for cohort function is A + c to get c starting at 1
  if(mod == 'apc'){
    data$expec <- FUNage(data$a) + FUNperiod(data$p) + FUNcohort(A + data$c)
  }else if(mod == 'ac'){
    data$expec <- FUNage(data$a) + FUNcohort(A + data$c)
  }else if(mod == 'pc'){
    data$expec <- FUNperiod(data$p) + FUNcohort(A + data$c)
  }else if(mod == 'ap'){
    data$expec <- FUNage(data$a) + FUNperiod(data$p)
  }else if(mod == 'c'){
    data$expec <- FUNcohort(A + data$c)
  }else if(mod == 'p'){
    data$expec <- FUNperiod(data$p)
  }else if(mod == 'a'){
    data$expec <- FUNage(data$a)
  }
  
  # define a new vector for a with M = M aggregation
  data$aFactor <- as.factor(data$a) # a as a factor to make relabeling easier
  Aprime <- ceiling(A/M) # the new number of age groups
  aNew <- zoo::rollapply(a, width = M + 1, by = M, mean, align = 'left') # want to average over M groups of a
  aTmat <- kronecker(diag(Aprime), matrix(1, ncol = M , nrow = 1))[,1:A] # matrix to get aNew the correct length of A
  levels(data$aFactor) <- aNew %*% aTmat # renaming levels of a to new sequence
  data$aNew <- as.numeric(levels(data$aFactor))[data$aFactor]
  
  # normal distribution
  ## generate y from M = 1 a, p and c
  ## remove M = 1 a, use M = M a to define c
  if(distribution == 'normal'){
    
    data <- purrr::map_dfr(seq_len(N), ~data) # repeat each data set N times
    data$N <- rep(1:N, each = A*P) # N is a label for the replication in this case
    for(i in 1:nrow(data)){
      data$y[i] <- rnorm(n = 1, mean = data$expec[i], sd = 1)
    }
    
    dataAgg <- 
      data %>% 
      select(aNew, p, N, y) %>%  # remove unwanted columns (incl M = 1 a)
      rename(a = aNew) %>% 
      mutate(c = p-a) %>%  # define c based off of M = M a
      relocate(a, p, c, .before = N)
    
  }else if(distribution %in% c('poisson', 'binomial')){
    
    # binomial and poisson
    ## generate y from M = 1 a,p and c
    ## aggregate (sum) N and y over p and  M = M a then define new c 
    
    if(distribution == 'poisson'){
      data$N <- rep(N, times = nrow(data)) # N is population at risk per a-p combination
      data$lambda <- (data$N * exp(data$expec))
      # data$lambda <- (exp(data$expec))
      for(i in 1:nrow(data)){
        data$y[i] <- rpois(n = 1, lambda = data$lambda[i])
      }
    } else if (distribution == 'binomial'){
      expit <- function(x){
        exp(x)/(1 + exp(x))
      }
      data$N <- rep(N, times = nrow(data)) # N is the number of observations per a-p combination
      data$prob <- expit(data$expec)
      for(i in 1:nrow(data)){
        data$y[i] <- rbinom(n = 1, size = data$N[i], prob = data$prob[i])
      }
    }
    
    # aggregation over aNew (where M = M) and p
    agg <- aggregate(cbind(N,y) ~ aNew + p, data = data, sum, na.rm = T) # summing N and y over p adn M = M a
    agg$c <- agg$p-agg$aNew # defining c using p and M = M a
    
    dataAgg <- 
      agg %>% 
      relocate(aNew,p,c, .before = N) %>% 
      rename(a = aNew) %>% 
      arrange(N,a,p,c)
    
  }
  
  dataAgg # display the aggregated data
  
}

## function to find the estimates of effects and curvatures given y ------------

find.all.estimates <- function(data, y){
  
  data$y <- y
  
  # max value of each temproal effect
  A <- length(unique(data$a));P <- length(unique(data$p));C <- length(unique(data$c))
  
  # aggregated by finding the mean over each age, period and cohort index
  f_a <- aggregate(y ~ a,FUN = mean,data = data)[,2]-mean(data$y)
  f_p <- aggregate(y ~ p,FUN = mean,data = data)[,2]-mean(data$y)
  f_c <- aggregate(y ~ c,FUN = mean,data = data)[,2]-mean(data$y)
  
  # hat matrices for a, p and c
  H_a <- cbind(rep(1, times = A), 1:A)
  H_p <- cbind(rep(1, times = P), 1:P)
  H_c <- cbind(rep(1, times = C), 1:C)
  
  # detrended temporal terms
  f_a_curv <- c(f_a-(H_a%*%qr.solve(H_a, f_a)))
  f_p_curv <- c(f_p-(H_p%*%qr.solve(H_p, f_p)))
  f_c_curv <- c(f_c-(H_c%*%qr.solve(H_c, f_c)))
  
  full <- list(age = f_a, period = f_p, cohort = f_c)
  curvature <- list(age = f_a_curv, period = f_p_curv, cohort = f_c_curv)
  
  return(list(full = full,curvature = curvature))
}

## true value generation for continuous data -----------------------------------

true.values <- function(A, P, M, mod, FUNage, FUNperiod, FUNcohort){
  
  # a and p continuous
  ## for true points, we only want for aggregated A, not for single A then aggregate after
  a <- seq(from = 0,to = A,by = 1); aMid <- zoo::rollapply(a, width = 1 + 1, by = 1, mean, align = 'left') # want to average over 1 group in a
  p <- seq(from = 0,to = P,by = 1); pMid <- zoo::rollapply(p, width = 1 + 1, by = 1, mean, align = 'left') # want to average over 1 group in p
  dataTemp <- data.frame(a = rep(aMid, each = length(pMid)),
                       p = rep(pMid, times = length(aMid)))
  dataTemp$c <- dataTemp$p-dataTemp$a
  
  # adding new columns to data temp for the aggregated age and cohort 
  ## age
  # new way to define a new vector for a with M = M aggregation
  dataTemp$aFactor <- as.factor(dataTemp$a) # a as a factor to make relabeling easier
  Aprime <- ceiling(A/M) # the new number of age groups
  aNew <- zoo::rollapply(aMid, width = M, by = M, mean, align = 'left', partial = TRUE) # want to average over M groups of aMid
  aTmat <- kronecker(diag(Aprime), matrix(1, ncol = M , nrow = 1))[,1:A] # matrix to get aNew the correct length of A
  levels(dataTemp$aFactor) <- aNew %*% aTmat # renaming levels of a to new sequence
  dataTemp$aNew <- as.numeric(levels(dataTemp$aFactor))[dataTemp$aFactor]
  ## cohort
  dataTemp$cNew <- dataTemp$p-dataTemp$aNew
  
  # the balanced data frame
  data <- expand.grid(a = aMid, p = pMid, c = sort(unique(dataTemp$c)))
  # data <- dataTemp; data$c <- sort(dataTemp$c) # to check values are same as in data.sim
  
  # getting the expected value of y given model
  if(mod == 'apc'){
    data$ytrue <- FUNage(data$a) + FUNperiod(data$p) + FUNcohort(A + data$c)
  }else if(mod == 'ac'){
    data$ytrue <- FUNage(data$a) + FUNcohort(A + data$c)
  }else if(mod == 'pc'){
    data$ytrue <- FUNperiod(data$p) + FUNcohort(A + data$c)
  }else if(mod == 'ap'){
    data$ytrue <- FUNage(data$a) + FUNperiod(data$p)
  }else if(mod == 'c'){
    data$ytrue <- FUNcohort(A + data$c)
  }else if(mod == 'p'){
    data$ytrue <- FUNperiod(data$p)
  }else if(mod == 'a'){
    data$ytrue <- FUNage(data$a)
  }
  
  estimates <- find.all.estimates(data = data, y = data$ytrue)
  
  # the full (M = 1) versions of age and cohort estimates
  ## age
  f_a_full <- data.frame(Group = unique(dataTemp[,c('a','aNew')])$aNew,
                       Estimate = estimates$full$age)
  f_a_curv_full <- data.frame(Group = unique(dataTemp[,c('a','aNew')])$aNew,
                            Estimate = estimates$curvature$age)
  ## cohort
  f_c_full <- data.frame(Group = sort(unique(dataTemp[,'c'])), Estimate = estimates$full$cohort)
  f_c_curv_full <- data.frame(Group = sort(unique(dataTemp$c)), Estimate = estimates$curvature$cohort)
  
  # taking the reduced (M = M) versions of age and cohort estimates
  ## age: average over every M groups with M groups gap between
  f_a <- zoo::rollapply(f_a_full, width = M, by = M, FUN = mean)[,2]
  f_a_curv <- zoo::rollapply(f_a_curv_full, width = M, by = M, FUN = mean)[,2]
  ## cohort: average over every M groups with 1 group gap between
  f_c <- zoo::rollapply(f_c_full, width = M, by = 1, FUN = mean)[,2]
  f_c_curv <- zoo::rollapply(f_c_curv_full, width = M, by = 1, FUN = mean)[,2]
  
  # period needs no extra attention
  f_p <- estimates$full$period
  f_p_curv <- estimates$curvature$period
  
  # data frame
  trueEffData <- 
    data.frame(Effect = factor(c(rep('Age', times = length(f_a)),
                               rep('Period', times = length(f_p)),
                               rep('Cohort', times = length(f_c))),
                             levels = c('Age', 'Period', 'Cohort')))
  # data frame of the true model curvatures
  trueCurvData <- 
    data.frame(Effect = factor(c(rep('Age', times = length(f_a_curv)),
                               rep('Period', times = length(f_p_curv)),
                               rep('Cohort', times = length(f_c_curv))),
                             levels = c('Age', 'Period', 'Cohort')))
  
  trueEffData$Group <- c(sort(unique(dataTemp$aNew)),sort(unique(dataTemp$p)),sort(unique(dataTemp$cNew)))
  trueEffData$Estimate <- c(f_a, f_p, f_c)
  trueEffData$Model <- 'True'
  
  trueCurvData$Group <- c(sort(unique(dataTemp$aNew)),sort(unique(dataTemp$p)),sort(unique(dataTemp$cNew)))
  trueCurvData$Estimate <- c(f_a_curv, f_p_curv, f_c_curv)
  trueCurvData$Model <- 'True'
  
  # plot of the true model estimates
  trueEffPlot <- 
    ggplot(trueEffData, aes(x = Group, y = Estimate))  + 
    geom_point(aes(col = Effect, shape = Effect)) + 
    geom_line(aes(col = Effect)) + 
    facet_wrap(. ~Effect, scales = 'free')  + 
    guides(col = 'none')  + 
    my.theme(legend.position = c(0.95, 0.95),
             legend.justification = c(1, 1))
  # plot of the true model curvatures
  trueCurvPlot <- 
    ggplot(trueCurvData, aes(x = Group, y = Estimate))  + 
    geom_point(aes(col = Effect, shape = Effect)) + 
    geom_line(aes(col = Effect)) + 
    facet_wrap(. ~Effect, scales = 'free')  + 
    guides(col = 'none')  + 
    my.theme(legend.position = c(0.95, 0.95),
             legend.justification = c(1, 1))
  
  return(list(trueEffData = trueEffData, trueCurvData = trueCurvData, 
              trueEffPlot = trueEffPlot, trueCurvPlot = trueCurvPlot))
}

