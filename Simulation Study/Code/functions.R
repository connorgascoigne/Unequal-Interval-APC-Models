# data generation functions ----

## normal ----
age.fun.normal <- function(a){
  0.3*a - 0.01*(a^2)
}

per.fun.normal <- function(p){
  -0.04*p + 0.02*(p^2)
}

coh.fun.normal <- function(c){
  0.35*c - 0.0015*(c^2)
}

## binomial ----
age.fun.binomial <- function(a){
  0.4 + (1/50)*(0.3*a - 0.01*(a^2))
}

per.fun.binomial <- function(p){
  (1/50)*(-0.04*p + 0.02*(p^2))
}

coh.fun.binomial <- function(c){
  (1/50)*(0.35*c - 0.0015*(c^2))
}

## poisson ----
age.fun.poisson <- function(a){
  -1.5 + (1/50)*(0.3*a - 0.01*(a^2))
}

per.fun.poisson <- function(p){
  (1/50)*(-0.04*p + 0.02*(p^2))
}

coh.fun.poisson <- function(c){
  (1/50)*(0.35*c - 0.0015*(c^2))
}

## data frame simulation ----

# this is for continuous data (previous version is for factor)
# this is for continuous data (previous version is for factor)
# works for all M not just A//M
data.sim <- function(A, P, M, N, family = c('normal', 'poisson', 'binomial'),
                     mod = c('apc','ac','ap','pc','a','p','c'),
                     FUNage, FUNperiod, FUNcohort){
  
  # # inputs to run all parts of the function
  # # clear and reset params
  # rm(list = setdiff(ls(), lsf.str()))
  # A = 20; P = 20; N = 25; mod = 'apc';
  # M = 1;
  # # family = 'normal'; FUNage = age.fun.normal;  FUNperiod = per.fun.normal; FUNcohort = coh.fun.normal
  # # family = 'poisson'; FUNage = age.fun.poisson;  FUNperiod = per.fun.poisson; FUNcohort = coh.fun.poisson
  # family = 'binomial'; FUNage = age.fun.binomial;  FUNperiod = per.fun.binomial; FUNcohort = coh.fun.binomial
  
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
  
  # new way to define a new vector for a with M = M aggregation
  data$aFactor <- as.factor(data$a) # a as a factor to make relabeling easier
  Aprime <- ceiling(A/M) # the new number of age groups
  aNew <- zoo::rollapply(a, width = M + 1, by = M, mean, align = 'left') # want to average over M groups of a
  aTmat <- kronecker(diag(Aprime), matrix(1, ncol = M , nrow = 1))[,1:A] # matrix to get aNew the correct length of A
  levels(data$aFactor) <- aNew %*% aTmat # renaming levels of a to new sequence
  data$aNew <- as.numeric(levels(data$aFactor))[data$aFactor]
  
  # normal family
  ## generate y from M = 1 a, p and c
  ## remove M = 1 a, use M = M a to define c
  if(family == 'normal'){
    
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
    
  }else if(family %in% c('poisson', 'binomial')){
    
    # binomial and poisson
    ## generate y from M = 1 a,p and c
    ## aggregate (sum) N and y over p and  M = M a then define new c 
    
    if(family == 'poisson'){
      data$N <- rep(N, times = nrow(data)) # N is population at risk per a-p combination
      data$lambda <- (data$N * exp(data$expec))
      # data$lambda <- (exp(data$expec))
      for(i in 1:nrow(data)){
        data$y[i] <- rpois(n = 1, lambda = data$lambda[i])
      }
    } else if (family == 'binomial'){
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
  # data
  # all.equal(dataAgg,data[,c('a','p','c','N','y')]) # check if true or nah
  
}

# function to find the estimates of effects and curvatures given y
## now for continuous
find.all.estimates <- function(data, y){
  
  # data = data; y = data$ytrue
  
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

# true value generation for continuous data
true.values <- function(A, P, M, mod, FUNage, FUNperiod, FUNcohort){
  
  # # clear and reset params
  # rm(list = setdiff(ls(), lsf.str()))
  # A = 60; P = 20; N = 25; mod = 'apc';
  # M = 3
  
  # a and p continuous
  ## for true points, we only want for aggregated A, not for single A then aggregate after
  # a <- seq(from = 0,to = A,by = M); aMid <- a[-length(a)]  +  diff(a)/2 # using M = M
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

# model fitting ----

## for spline model ----

# general linear constraint <- - for holford M = linear term
gen.lin.constr <- function(M, X){
  # goal t(M) %*% (X %*% Z)  =  0
  C <- t(M) %*% X
  qrc <- qr(t(C))
  Z <- qr.Q(qrc, complete = TRUE)[,(nrow(C) + 1):ncol(C)] # define Z the last n-1 columns of find the Q
  Z
}

# updated to have re-parameterisation with cbind(1,linear)
my.PredictMat <- function(smObj, data, n = nrow(data), include.periodic = FALSE){
  # smObj = gamObj$smooth[[k]]; data = data
  
  # get the correct term and the basis dimensions
  term <- smObj$term
  bs.dim <- smObj$bs.dim
  data$x <- data[,term] # needed to relabel to fit gam
  
  # if periodic, alter the period and cohort terms ONLY
  if(include.periodic == TRUE && term %in% c('p', 'c')){
    
    # retrieve number of knots for periodic and non periodic parts
    nk <- bs.dim + 1-5
    
    # getting a cylic vector for x
    xFact <- factor(data$x)
    data$xVec <- as.numeric(xFact)
    data$xCC <- (data$xVec %% 5) + 1
    
    # smooth bases matrix with intercept constraint
    smNonPer <- mgcv::smoothCon(s(x, bs = 'cr', k = nk + 1),data = data, absorb.cons = F)[[1]] # non periodic: increase dim by 1 to orthogomalise correctly
    smPer <- mgcv::smoothCon(s(xCC, bs = 'cc', k = 5), knots = list(x = c(1,5)), data = data, absorb.cons = F)[[1]] # periodic part additionally add one dimension for end constraints
    
    X <- cbind(lab = data[,term], smNonPer$X, smPer$X)
    Xuni <- unique(X)[,-1]
    smS <- as.matrix(Matrix::bdiag(smNonPer$S[[1]], smPer$S[[1]]))
    
    
  } else {
    
    # smooth basis matrix with intercept constraint
    # need to increase the basis dimension by 1
    # to then orthogonalise correctly
    sm <- mgcv::smoothCon(s(x, bs = 'cr',k = bs.dim + 1), data = data, absorb.cons = F)[[1]]
    X <- cbind(lab = data[,term], sm$X)
    Xuni <- unique(X)[,-1]
    smS <- sm$S[[1]]
    
  }
  
  # data SHOULD be centered (orthogonal to intercept) linear column 
  lin <- cbind(lab = data[,term], lin = data[,term])
  linUni <- as.matrix(unique(lin)[,-1], ncol = 1)
  
  # perform re-parameterisation to unique columns of X
  Z <- gen.lin.constr(M = cbind(1, linUni), X = Xuni)
  curUni <- Xuni %*% Z
  S <- t(Z) %*% smS %*% Z
  
  # t(linUni) %*% curUni # should be zeros!
  
  # expand the curvature terms to be full number of rows in correct order
  lab <- data.frame(lab = data[,term]);
  curUni2 <- data.frame(lab = unique(data[,term]), curUni)
  cur <- left_join(lab, curUni2, by = 'lab') %>% select(-lab) %>% as.matrix
  
  # t(lin[,-1]) %*% cur # might change
  
  list(cur = cur, S = S)
}

# my prediction function fior gamObjs
my.predict <- function(gamObj, include.periodic = FALSE, newdata = NULL, type = 'response', exclude = NULL, inverse.link = FALSE){
  
  # gamObj <- fit$m; newdata = data; exclude = NULL; response.scale = F
  
  # use the new data if not use the data present in gam fit
  if(!is.null(newdata)){
    data <- newdata
  } else {
    data <- gamObj$model
  }
  
  # define matrix to store information in
  X <- matrix(0, nrow = nrow(data), ncol = ncol(model.matrix(gamObj)))
  
  # parametric parts
  Terms <- list(delete.response(gamObj$pterms))
  term.labels <- attr(gamObj$pterms, 'term.labels')
  n.pterms <- length(term.labels)
  
  # find columns for intercept + parameteric terms
  Xp <- model.matrix(Terms[[1]], data)
  # populating intercept and parametric columns of X
  X[,1:(1 + n.pterms)] <- Xp # nsdf: number of parametric, non-smooth, model terms including the intercept.
  
  # colnames for labelling
  colNames <- colnames(Xp); termNames <- term.labels
  
  # populating the smooth columns of X
  n.smooth <- length(gamObj$smooth)
  for (k in 1:n.smooth) {
    klab <- gamObj$smooth[[k]]$label
    if ((is.null(exclude) || !(klab %in% exclude))) { # do this for only null
      term <- gamObj$smooth[[k]]$term
      if(term %in% c('a', 'p', 'c')){ # do this for only temporal terms
        predMats <- my.PredictMat(gamObj$smooth[[k]], data, include.periodic = include.periodic) # my.PredictMat function defines smooth then reparameterises it
        Xfrag <- predMats$cur
      } else {
        Xfrag <- mgcv::PredictMat(gamObj$smooth[[k]], data) # for non temproal smooth terms
      }
      first.last <- gamObj$smooth[[k]]$first.para:gamObj$smooth[[k]]$last.para
      kname <- paste0(klab, '.', 1:length(first.last))
      X[, first.last] <- Xfrag
      colNames[first.last] <- kname
      termNames[n.pterms + k] <- klab
    }
  }
  
  # for labeling
  colnames(X) <- colNames
  
  # model family for inverse and delta function
  fam <- gamObj$family
  linkinv <- fam$linkinv 
  dmu.deta <- fam$mu.eta # for delta function
  
  # parameters and variance-covariance matrix after model fitting and 
  beta <- gamObj$coefficients
  Vp <- gamObj$Vp; colnames(Vp) <- rownames(Vp) <- colNames
  
  if(type == 'lpmatrix'){
    # return the design matrix
    H <- X
  }else if(type == 'response'){
    
    est2 <- X %*% beta # linear predictions with X defined above
    se2 <- sqrt(pmax(0, rowSums((X  %*%  Vp) * X))) # standard errors 
    se <- se2*abs(dmu.deta(est2)) # standard errors found via delta method
    
    # inverse the link function
    if(inverse.link == T){
      H <- data.frame(est = linkinv(est2), se = se)
    }else {
      H <- data.frame(est = est2, se = se)
    }
  }else if(type == 'terms'){
    
    est <- matrix(0, nrow = nrow(data), ncol = length(termNames))
    colnames(est) <- termNames
    se <- est
    
    # update terms for parameteric components
    for(j in 1:n.pterms){
      jlab <- term.labels[j]
      est2 <- X[,jlab,drop = F] %*% beta[jlab]
      
      se3 <- sqrt(pmax(0, rowSums((X[,jlab,drop = F]  %*%  Vp[jlab,jlab]) * X[,jlab,drop = F])))
      se2 <- se3*abs(dmu.deta(est2))
      
      est[,jlab] <- est2
      se[,jlab] <- se2
    }
    
    # update terns for smooth components
    for(k in 1:n.smooth){
      klab <- gamObj$smooth[[k]]$label
      if ((is.null(exclude) || !(klab %in% exclude))) { # do this for only null
        first.last <- gamObj$smooth[[k]]$first.para:gamObj$smooth[[k]]$last.para
        est2 <- X[,first.last,drop = F] %*% beta[first.last]
        
        se3 <- sqrt(pmax(0, rowSums((X[,first.last,drop = F]  %*%  Vp[first.last,first.last]) * X[,first.last,drop = F])))
        se2 <- se3*abs(dmu.deta(est2))
        
        est[,n.pterms + k] <- est2
        se[,n.pterms + k] <- se2
      }
    }
    H <- list(est = est,se = se)
  }
  
  H # return H
}

# fitting the spline model
## has the option to choose whether the spline is fixed (regression) or not (penalised)
## has the option to choose the basis function for the spline (basis from mgcv)
hol.spline.fun <- function(data, mod = c('apc','ac','ap','pc','a','p','c'), bs = 'cr', include.periodic = FALSE,
                           knots = NULL, fixed = NULL, slopeDrop = NULL, family = c('normal','poisson','binomial')){
  
  # # rm(list = setdiff(ls(), lsf.str()))
  # A <- 60; P <- 20; N <- 150; M = 5;
  # slopeDrop = 'c'; mod = 'apc';
  # knots <- list(age = 5, period = 5, cohort = 8)
  # fixed <- list(age = F,period = F,cohort = F)
  # bs <- 'cr'
  # include.periodic = FALSE
  # # family = 'normal'; FUNage = age.fun.normal;  FUNperiod = per.fun.normal; FUNcohort = coh.fun.normal
  # family = 'poisson'; FUNage = age.fun.poisson;  FUNperiod = per.fun.poisson; FUNcohort = coh.fun.poisson
  # # family = 'binomial'; FUNage = age.fun.binomial;  FUNperiod = per.fun.binomial; FUNcohort = coh.fun.binomial
  # data <- data.sim(A = A, P = P, M = M, N = N, mod = mod, family = family,
  #                FUNage = FUNage,FUNperiod = FUNperiod,FUNcohort = FUNcohort)
  
  if(is.null(knots)|!is.null(knots) & !is.list(knots)) stop('Warning: Need a list for number of knots for each effect labelled age, period, cohort.')
  if(is.null(fixed)|!is.null(fixed) & !is.list(fixed)) stop('Warning: Need a list for whether an effect is penalised or not.')
  if(missing(family)) stop('Warning: Need to specify the family (normal, poisson or binomial)')
  if(mod == 'apc'&&is.null(slopeDrop)) stop('Warning: When fitting an APC model need to drop linear column of one of age (a), period (p) or cohort (c).')
  if(mod == 'apc'&&!(slopeDrop %in% c('a', 'p', 'c'))) stop('slopeDrop needs to be one of a, p or c')
  
  # if data needs additional changes
  ## center data
  dataAug <- 
    data %>%
    mutate(a = a %>% scale(., scale = FALSE),
           p = p %>% scale(., scale = FALSE),
           c = c %>% scale(., scale = FALSE))
  
  ak <- knots$age; pk <- knots$period; ck <- knots$cohort
  afx <- fixed$age; pfx <- fixed$period; cfx <- fixed$cohort # default is to have fixed = FALSE for penalisation
  
  if(include.periodic == TRUE){
    formula <- 
      as.formula(y ~ 
                   a + p + c + 
                   s(a,bs = bs,k = ak-1,fx = afx) + 
                   s(p,bs = bs,k = (pk-1) + (5),fx = pfx) + 
                   s(c,bs = bs,k = (ck-1) + (5),fx = cfx))
    # update pre and post forumla for the temporal models
    # # do not remove any linear terms from pre fit as they are all needed later in 
    # # re-parameterisation
    if (mod == 'apc'){
      if(slopeDrop == 'a'){
        formula <- update(formula, ~.-a)
      } else if (slopeDrop == 'p'){
        formula <- update(formula, ~.-p)
      } else if (slopeDrop == 'c'){
        formula <- update(formula, ~.-c)
      }
    } else if (mod == 'ac'){
      formula <- update(formula, ~.-p-s(p, bs = bs, k = (pk-1) + (5) ,fx = pfx))
    } else if (mod == 'pc'){
      formula <- update(formula, ~.-a-s(a, bs = bs, k = ak-1, fx = afx))
    } else if (mod == 'ap'){
      formula <- update(formula, ~.-c-s(c, bs = bs, k = (ck-1) + (5), fx = cfx))
    } else if (mod == 'c'){
      formula <- update(formula, ~.-a-s(a, bs = bs, k = ak-1, fx = afx)-p-s(p, bs = bs, k = (pk-1) + (5), fx = pfx))
    } else if (mod == 'a'){
      formula <- update(formula, ~.-p-s(p, bs = bs, k = (pk-1) + (5), fx = pfx)-c-s(c, bs = bs, k = (ck-1) + (5), fx = cfx))
    } else if (mod == 'p'){
      formula <- update(formula, ~.-a-s(a, bs = bs, k = ak-1, fx = afx)-c-s(c, bs = bs, k = (ck-1) + (5), fx = cfx))
    }
    
    # use all cores avaliable
    ctrl <- list(nthreads = parallel::detectCores())
    
    # the extra information depending on the model
    if(family == 'normal'){
      fit <- mgcv::gam(formula, family = 'gaussian', data = dataAug, method = 'REML', control = ctrl, fit = F)
    } else if(family == 'poisson'){
      fit <- mgcv::gam(formula, offset = log(N), family = 'poisson', data = dataAug, method = 'REML', control = ctrl, fit = F)
    } else if(family == 'binomial'){
      formula <- update(formula, cbind(y, N-y)~.)
      fit <- mgcv::gam(formula, family = 'binomial', data = dataAug, method = 'REML', control = ctrl, fit = F)
    }
    
    # only need to redefine smooth terms for temporal terms
    # do not want to re-parameterise the cluster term
    nsmooth <- length(fit$smooth)
    for(i in 1:nsmooth){
      # the correct smooth term and label
      smObj <- fit$smooth[[i]]
      lab <- smObj$term
      if(lab %in% c('a', 'p', 'c')){
        mats <- my.PredictMat(smObj = smObj, data = dataAug, include.periodic = TRUE)
        # updating the fit with re-parameterised TEMPROAL terms
        Xpos <- smObj$first.para:smObj$last.para
        fit$X[,Xpos] <- mats$cur
        if(smObj$fixed == FALSE){
          fit$S[[i]] <- mats$S
        }
      }
    }
    
  } else {
    
    formula <- 
      as.formula(y ~ 
                   a + p + c + 
                   s(a, bs = bs, k = ak-1, fx = afx) + 
                   s(p, bs = bs, k = pk-1, fx = pfx) + 
                   s(c, bs = bs, k = ck-1, fx = cfx))
    # update pre and post forumla for the temporal models
    # # do not remove any linear terms from pre fit as they are all needed later in 
    # # re-parameterisation
    if (mod == 'apc'){
      if(slopeDrop == 'a'){
        formula <- update(formula, ~.-a)
      } else if (slopeDrop == 'p'){
        formula <- update(formula, ~.-p)
      } else if (slopeDrop == 'c'){
        formula <- update(formula, ~.-c)
      }
    } else if (mod == 'ac'){
      formula <- update(formula, ~.-p-s(p, bs = bs, k = pk-1, fx = pfx))
    } else if (mod == 'pc'){
      formula <- update(formula, ~.-a-s(a, bs = bs, k = ak-1, fx = afx))
    } else if (mod == 'ap'){
      formula <- update(formula, ~.-c-s(c, bs = bs, k = ck-1, fx = cfx))
    } else if (mod == 'c'){
      formula <- update(formula, ~.-a-s(a, bs = bs, k = ak-1, fx = afx)-p-s(p, bs = bs, k = pk-1, fx = pfx))
    } else if (mod == 'a'){
      formula <- update(formula, ~.-p-s(p, bs = bs, k = pk-1, fx = pfx)-c-s(c, bs = bs, k = ck-1, fx = cfx))
    } else if (mod == 'p'){
      formula <- update(formula, ~.-a-s(a, bs = bs, k = ak-1, fx = afx)-c-s(c, bs = bs, k = ck-1, fx = cfx))
    }
    
    # use all cores avaliable
    ctrl <- list(nthreads = parallel::detectCores())
    
    # the extra information depending on the model
    if(family == 'normal'){
      fit <- mgcv::gam(formula, family = 'gaussian', data = dataAug, method = 'REML', control = ctrl, fit = F)
    } else if(family == 'poisson'){
      fit <- mgcv::gam(formula, offset = log(N), family = 'poisson', data = dataAug, method = 'REML', control = ctrl, fit = F)
    } else if(family == 'binomial'){
      formula <- update(formula, cbind(y, N-y)~.)
      fit <- mgcv::gam(formula, family = 'binomial', data = dataAug, method = 'REML', control = ctrl, fit = F)
    }
    
    # only need to redefine smooth terms for temporal terms
    # do not want to re-parameterise the cluster term
    nsmooth <- length(fit$smooth)
    for(i in 1:nsmooth){
      # the correct smooth term and label
      smObj <- fit$smooth[[i]]
      lab <- smObj$term
      if(lab %in% c('a', 'p', 'c')){
        mats <- my.PredictMat(smObj = smObj,data = dataAug, include.periodic = FALSE)
        # updating the fit with re-parameterised TEMPROAL terms
        Xpos <- smObj$first.para:smObj$last.para
        fit$X[,Xpos] <- mats$cur
        if(smObj$fixed == FALSE){
          fit$S[[i]] <- mats$S
        }
      }
    }
    
  }
  
  mod = mgcv::gam(G = fit)
  
  return(list(mod = mod, X = fit$X, S = fit$S))
}

# extracts the spline function
## has the fixed argument now
## has the type of basis argument now
hol.spline.extract <- function(data, mod = c('apc','ac','ap','pc','a','p','c'), bs = 'cr', include.periodic = FALSE,
                               knots = NULL, fixed = NULL, slopeDrop = NULL, family = c('normal','poisson','binomial')){
  
  # # rm(list = setdiff(ls(), lsf.str()))
  # A <- 60; P <- 20; N <- 150; M = 5;
  # slopeDrop = 'c'; mod = 'apc';
  # knots <- list(age = 5, period = 5, cohort = 8)
  # fixed <- list(age = F,period = F,cohort = F)
  # bs <- 'cr'
  # include.periodic = FALSE
  # # family = 'normal'; FUNage = age.fun.normal;  FUNperiod = per.fun.normal; FUNcohort = coh.fun.normal
  # family = 'poisson'; FUNage = age.fun.poisson;  FUNperiod = per.fun.poisson; FUNcohort = coh.fun.poisson
  # # family = 'binomial'; FUNage = age.fun.binomial;  FUNperiod = per.fun.binomial; FUNcohort = coh.fun.binomial
  # data <- data.sim(A = A, P = P, M = M, N = N, mod = mod, family = family,
  #                FUNage = FUNage,FUNperiod = FUNperiod,FUNcohort = FUNcohort)
  
  # the reparameterised fit
  fit <- hol.spline.fun(data = data, mod = mod, bs = bs, include.periodic = include.periodic, knots = knots, fixed = fixed, slopeDrop = slopeDrop, family = family)
  
  # model fit
  m <- fit$mod
  
  # defining balanced data frame to define a new X
  ## the unique values of a, p and c used in the model
  a <- sort(unique(data$a)); p <- sort(unique(data$p)); c <- sort(unique(data$c))
  ## a balanced data frame of the a, p and c values used in the model
  data2 <- expand.grid(a = a, p = p, c = c)
  ## new X defined at the same values as the model but with balanced data
  X <- my.predict(gamObj = m, newdata = data2, include.periodic = include.periodic, type = 'lpmatrix')
  ## coefficients//parameters
  beta <- coef(m)
  ## yhat
  yhat <- X %*% beta
  
  # find estimate for age, period and cohort full and curvature
  estimates <- find.all.estimates(data = data2, y = yhat)
  
  hat_f_a <- estimates$full$age
  hat_f_p <- estimates$full$period
  hat_f_c <- estimates$full$cohort
  
  hat_f_a_curv <- estimates$curvature$age
  hat_f_p_curv <- estimates$curvature$period
  hat_f_c_curv <- estimates$curvature$cohort
  
  estData <- 
    data.frame(Effect = 
                 factor(c(rep('Age', times = length(hat_f_a)),
                          rep('Period', times = length(hat_f_p)),
                          rep('Cohort', times = length(hat_f_c))),
                        levels = c('Age', 'Period', 'Cohort')))
  curvData <- 
    data.frame(Effect = 
                 factor(c(rep('Age', times = length(hat_f_a_curv)),
                          rep('Period', times = length(hat_f_p_curv)),
                          rep('Cohort', times = length(hat_f_c_curv))),
                        levels = c('Age', 'Period', 'Cohort')))
  
  estData$Group <- c(a, p, c)
  estData$Estimate <- c(hat_f_a, hat_f_p, hat_f_c)
  rownames(estData) <- 1:sum(length(hat_f_a), length(hat_f_p), length(hat_f_c))
  
  curvData$Group <- c(a, p, c)
  curvData$Estimate <- c(hat_f_a_curv, hat_f_p_curv, hat_f_c_curv)
  rownames(curvData) <- 1:sum(length(hat_f_a_curv), length(hat_f_p_curv), length(hat_f_c_curv))
  
  # plotting the results
  # the smooth functions of the estimates plot
  estPlot <- 
    ggplot(estData, aes(x = Group, y = Estimate)) + 
    geom_point(aes(col = Effect), alpha = 1.2) + 
    geom_line(aes(col = Effect), alpha = 1.2) + 
    facet_wrap(.~Effect, scales = 'free') + 
    guides(col = 'none') + 
    
    my.theme()
  # the un-centered smooth functions of curvatures plot
  curvPlot <- 
    ggplot(curvData, aes(x = Group, y = Estimate)) + 
    geom_point(aes(col = Effect), alpha = 1.2) + 
    geom_line(aes(col = Effect), alpha = 1.2) + 
    facet_wrap(.~Effect, scales = 'free') + 
    guides(col = 'none') + 
    my.theme()
  
  
  return(list(estData = estData, curvData = curvData, 
              estPlot = estPlot, curvPlot = curvPlot, 
              yhat = yhat))
  
}

## for factor model ----

# function to make the factor matrix based on model
hol.factor.X <- function(data, mod, slopeDrop = NULL){
  
  # data = dataAug; mod = mod; slopeDrop = slopeDrop
  
  aX <- model.matrix(~factor(a), data = data)
  pX <- model.matrix(~factor(p), data = data)
  cX <- model.matrix(~factor(c), data = data)
  
  # linear terms SHOULD BE CENTERED
  aLin <- data$a
  pLin <- data$p
  cLin <- data$c
  
  # Z matrices
  aZ <- gen.lin.constr(M = cbind(1,aLin), X = aX)
  pZ <- gen.lin.constr(M = cbind(1,pLin), X = pX)
  cZ <- gen.lin.constr(M = cbind(1,cLin), X = cX)
  
  # curvature terms
  aCur <- aX %*% aZ
  pCur <- pX %*% pZ
  cCur <- cX %*% cZ
  
  # names
  colnames(aCur) <- paste0('aCur', 1:ncol(aCur))
  colnames(pCur) <- paste0('pCur', 1:ncol(pCur))
  colnames(cCur) <- paste0('cCur', 1:ncol(cCur))
  
  # # check orthogonal
  # t(aLin) %*% aCur
  # t(pLin) %*% pCur
  # t(cLin) %*% cCur
  
  # set apx X as the first and update accordingly
  X <- cbind(1, aLin, pLin, cLin, aCur, pCur, cCur)
  if(mod == 'apc'){
    if(slopeDrop == 'a'){
      X <- cbind(1, pLin, cLin, aCur, pCur, cCur)
    } else if(slopeDrop == 'p'){
      X <- cbind(1, aLin, cLin, aCur, pCur, cCur)
    } else if(slopeDrop == 'c'){
      X <- cbind(1, aLin, pLin, aCur, pCur, cCur)
    }
  } else if(mod == 'ac'){
    X <- cbind(1, aLin, cLin, aCur, cCur)
  } else if(mod == 'pc'){
    X <- cbind(1, pLin, cLin, pCur, cCur)
  } else if(mod == 'ap'){
    X <- cbind(1, aLin, pLin, aCur, pCur)
  } else if(mod == 'c'){
    X <- cbind(1, cLin, cCur)
  } else if(mod == 'a'){
    X <- cbind(1, aLin, aCur)
  } else if(mod == 'p'){
    X <- cbind(1, pLin, pCur)
  }
  
  X
  
}

# fit a factor based model
hol.factor.fun <- function(data, mod = c('apc','ac','ap','pc','a','p','c'),
                           slopeDrop = NULL, family = c('normal','poisson','binomial')){
  
  # # rm(list = setdiff(ls(), lsf.str()))
  # A <- 60; P <- 20; N <- 150; M = 5;
  # slopeDrop = 'c'; mod = 'apc';
  # # family = 'normal'; FUNage = age.fun.normal;  FUNperiod = per.fun.normal; FUNcohort = coh.fun.normal
  # family = 'poisson'; FUNage = age.fun.poisson;  FUNperiod = per.fun.poisson; FUNcohort = coh.fun.poisson
  # # family = 'binomial'; FUNage = age.fun.binomial;  FUNperiod = per.fun.binomial; FUNcohort = coh.fun.binomial
  # data <- data.sim(A = A, P = P, M = M, N = N, mod = mod, family = family,
  #                FUNage = FUNage,FUNperiod = FUNperiod,FUNcohort = FUNcohort)
  
  
  if(missing(family)) stop('Warning: Need to specify the family (normal, poisson or binomial)')
  if(mod == 'apc'&&is.null(slopeDrop)) stop('Warning: When fitting an APC model need to drop linear column of one of age (a), period (p) or cohort (c).')
  if(mod == 'apc'&&!(slopeDrop %in% c('a', 'p', 'c'))) stop('slopeDrop needs to be one of a, p or c')
  
  # if data needs additional changes
  ## center data
  dataAug <- 
    data %>%
    mutate(a = a %>% scale(., scale = FALSE),
           p = p %>% scale(., scale = FALSE),
           c = c %>% scale(., scale = FALSE))
  
  # define the factor model matrix 
  X <- hol.factor.X(data = dataAug, mod = mod, slopeDrop = slopeDrop)
  
  # use all cores avaliable
  ctrl <- list(nthreads = parallel::detectCores())
  
  # the extra information depending on the model
  if(family == 'normal'){
    mod <- mgcv::gam(y ~ -1 + X, family = 'gaussian', data = dataAug, method = 'REML', control = ctrl)
  } else if(family == 'poisson'){
    mod <- mgcv::gam(y ~ -1 + X, offset = log(N), family = 'poisson', data = dataAug, method = 'REML', control = ctrl)
  } else if(family == 'binomial'){
    mod <- mgcv::gam(cbind(y, N-y) ~ -1 + X, family = 'binomial', data = dataAug, method = 'REML', control = ctrl)
  }
  
  return(list(mod = mod, X = X))
}

# extract the results from a factor based model
hol.factor.extract <- function(data, slopeDrop = NULL, family = c('normal','poisson','binomial'),
                               mod = c('apc', 'ac', 'ap', 'pc', 'a', 'p', 'c')){
  
  # # rm(list = setdiff(ls(), lsf.str()))
  # A <- 60; P <- 20; N <- 150; M = 5;
  # slopeDrop = 'c'; mod = 'apc';
  # # family = 'normal'; FUNage = age.fun.normal;  FUNperiod = per.fun.normal; FUNcohort = coh.fun.normal
  # family = 'poisson'; FUNage = age.fun.poisson;  FUNperiod = per.fun.poisson; FUNcohort = coh.fun.poisson
  # # family = 'binomial'; FUNage = age.fun.binomial;  FUNperiod = per.fun.binomial; FUNcohort = coh.fun.binomial
  # data <- data.sim(A = A, P = P, M = M, N = N, mod = mod, family = family,
  #                FUNage = FUNage,FUNperiod = FUNperiod,FUNcohort = FUNcohort)
  
  # the reparameterised fit
  fit <- hol.factor.fun(data = data, mod = mod, slopeDrop = slopeDrop, family = family)
  
  # model fit
  m <- fit$mod
  
  # defining balanced data frame to define a new X
  ## the unique values of a, p and c used in the model
  a <- sort(unique(data$a)); p <- sort(unique(data$p)); c <- sort(unique(data$c))
  ## a balanced data frame of the a, p and c values used in the model
  data2 <- expand.grid(a = a, p = p, c = c)
  ## new X defined at the same values as the model but with balanced data
  X <- hol.factor.X(data = data2, mod = mod, slopeDrop = slopeDrop)
  ## coefficients//parameters
  beta <- coef(m)
  ## yhat
  yhat <- X %*% beta
  
  # find estimate for age, period and cohort full and curvature
  estimates <- find.all.estimates(data = data2, y = yhat)
  
  hat_f_a <- estimates$full$age
  hat_f_p <- estimates$full$period
  hat_f_c <- estimates$full$cohort
  
  hat_f_a_curv <- estimates$curvature$age
  hat_f_p_curv <- estimates$curvature$period
  hat_f_c_curv <- estimates$curvature$cohort
  
  estData <- 
    data.frame(Effect = 
                 factor(c(rep('Age', times = length(hat_f_a)),
                          rep('Period', times = length(hat_f_p)),
                          rep('Cohort', times = length(hat_f_c))),
                        levels = c('Age', 'Period', 'Cohort')))
  curvData <- 
    data.frame(Effect = 
                 factor(c(rep('Age', times = length(hat_f_a_curv)),
                          rep('Period', times = length(hat_f_p_curv)),
                          rep('Cohort', times = length(hat_f_c_curv))),
                        levels = c('Age', 'Period', 'Cohort')))
  
  estData$Group <- c(sort(unique(data$a)), sort(unique(data$p)), sort(unique(data$c)))
  estData$Estimate <- c(hat_f_a,hat_f_p, hat_f_c)
  rownames(estData) <- 1:sum(length(hat_f_a), length(hat_f_p), length(hat_f_c))
  
  curvData$Group <- c(sort(unique(data$a)), sort(unique(data$p)), sort(unique(data$c)))
  curvData$Estimate <- c(hat_f_a_curv, hat_f_p_curv, hat_f_c_curv)
  rownames(curvData) <- 1:sum(length(hat_f_a_curv), length(hat_f_p_curv), length(hat_f_c_curv))
  
  # plotting the results
  # the smooth functions of the estimates plot
  estPlot <- 
    ggplot(estData, aes(x = Group, y = Estimate)) + 
    geom_point(aes(col = Effect), alpha = 1.2) + 
    geom_line(aes(col = Effect), alpha = 1.2) + 
    facet_wrap(.~Effect, scales = 'free') + 
    guides(col = 'none') + 
    my.theme()
  
  # the un-centered smooth functions of curvatures plot
  curvPlot <- 
    ggplot(curvData, aes(x = Group, y = Estimate)) + 
    geom_point(aes(col = Effect), alpha = 1.2) + 
    geom_line(aes(col = Effect), alpha = 1.2) + 
    facet_wrap(.~Effect, scales = 'free') + 
    guides(col = 'none') + 
    my.theme()
  
  return(list(estData = estData, curvData = curvData, 
              estPlot = estPlot, curvPlot = curvPlot))
  
}

# analysis ----

# theme for plots
my.theme<-function(...){
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black'),
        # legend.title=element_blank(),
        legend.text.align=0,
        legend.key=element_rect(fill=NA),
        ...)
}

# get the plots from simulation runs
get.final.plot <- function(simRun, boxplots = TRUE, ...){
  
  # to use plot/plot
  library(patchwork)
  
  # simRun = normalBasisResults
  
  # extracting the data frames from the simulation output
  simEst <- simRun$simEstData
  simCurv <- simRun$simCurvData
  truEff <- simRun$trueEffData
  trucurv <- simRun$trueCurvData
  
  # binding the estimate and curvature data frames together
  ## simulation
  simRes <- rbind(simEst %>% mutate(Type = 'est'),
                  simCurv %>% mutate(Type = 'curv')) %>%
    relocate(Type, .before = 's00001')
  
  ## true values
  true <- rbind(truEff %>% mutate(Type = 'est'),
                trucurv %>% mutate(Type = 'curv'))
  
  trueTemp <-
    true %>%
    mutate(truth = Estimate) %>%
    select(Effect, Group, Type, truth)
  
  summary <- 
    simRes %>%
    left_join(., trueTemp, by = c('Effect', 'Group', 'Type')) %>%
    mutate(Estimate = select(., starts_with('s0')) %>% rowMeans,
           Bias = (select(., starts_with('s0')) - truth) %>% rowMeans,
           MSE = ((select(., starts_with('s0')) - truth)^2) %>% rowMeans,
           Model = Model %>% haven::as_factor()) %>%
    select(-starts_with('s0'), -truth)
  
  # combining the simulation results and the true results together for specific scenario
  ## simulation results in long format
  simLong <- 
    summary %>% 
    pivot_longer(cols = c('Estimate', 'Bias', 'MSE'),
                 names_to = 'Statistic',
                 values_to = 'Value') %>% 
    filter(!(Type == 'est' & Statistic %in% c('Bias', 'MSE')))
  ## true results in long format
  trueLong <- 
    true %>% 
    pivot_longer(cols = 'Estimate',
                 names_to = 'Statistic',
                 values_to = 'Value')
  # ## combinging long simulation and true data frames
  finalDat <- 
    rbind(simLong, trueLong) %>% 
    mutate(int = factor(interaction(Type, Statistic),
                        levels = c('est.Estimate', 'curv.Estimate', 'curv.Bias', 'curv.MSE'),
                        labels = c('Effect', 'Curvature', 'Curvature Bias', 'Curvature MSE')))
  
  # creating the plots
  ## plots of continuous variables: Effects and curvatures
  p1 <- 
    ggplot()+
    geom_point(data = finalDat %>% filter(!(Statistic %in% c('Bias', 'MSE'))),
               aes(x = Group, y = Value, col = Model, shape = Model))+
    geom_line(data = finalDat %>% filter(!(Statistic %in% c('Bias', 'MSE'))),
              aes(x = Group, y = Value, col = Model, linetype = Model))+
    facet_grid(int~Effect, scales = 'free')+
    labs(x = 'Relative Years', y = '')+
    my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
             panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
             legend.title=element_blank(), # no legend title
             ...)
  ## plots of discrete variables: curvature bias and MSE box plots
  p2 <- 
    ggplot()+
    geom_boxplot(data = finalDat %>% filter(Statistic %in% c('Bias')),
                 aes(x = Model, y = Value, col = Model, fill = Model),
                 alpha = 0.2)+
    geom_jitter(data = finalDat %>% filter(Statistic %in% c('Bias')),
                aes(x = Model, y = Value, col = Model, shape = Model),
                alpha = 0.5)+
    geom_boxplot(data = finalDat %>% filter(Statistic %in% c('MSE')),
                 aes(x = Model, y = Value, col = Model, fill = Model),
                 alpha = 0.2)+
    geom_jitter(data = finalDat %>% filter(Statistic %in% c('MSE')),
                aes(x = Model, y = Value, col = Model, shape = Model),
                alpha = 0.5)+
    facet_grid(int~Effect, scales = 'free')+
    labs(y = '')+
    my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
             panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
             strip.text.x = element_blank(), # removing column (x-axis) facet labels as same as continuous
             legend.title=element_blank(), # no legend title
             ...)
  if(boxplots == TRUE){
    ## combining plots vertically
    p <- p1/p2
  } else {
    p <- p1
  }
  ## return combined plot
  # p
  
  # have plots for each model with their simulations over
  ## get the simulation data frame in order
  simFitsLong <- 
    simRes %>% 
    pivot_longer(cols = starts_with('s0'),
                 names_to = 'Simulation',
                 values_to = 'Estimate') %>% 
    mutate(int = factor(Type,
                        levels = c('est', 'curv'),
                        labels = c('Effect', 'Curvature'))) %>% 
    select(-Type)
  # levels(simFitsLong$Type) <- c('Effect', 'Curvature')
  
  # set up loop to creat plot for each model and store as a list
  simPlots <- list()
  models <- unique(simFitsLong$Model)
  for(i in 1:length(models)){
    dataSim <- simFitsLong %>% filter(Model == models[i])
    dataFinal <- finalDat %>% filter(!(Statistic %in% c('Bias', 'MSE')), Model %in% c('True', models[i]))
    p3 <- 
      ggplot()+
      geom_line(data = dataSim, alpha = 0.05,
                aes(x = Group, y = Estimate, col = Model, linetype = Model, group = interaction(Model,Simulation)))+
      geom_point(data = dataSim, alpha = 0.05, 
                 aes(x = Group, y = Estimate, col = Model, shape = Model, group = interaction(Model,Simulation)))+
      geom_point(data = dataFinal,
                 aes(x = Group, y = Value, col = Model, shape = Model))+
      geom_line(data = dataFinal,
                aes(x = Group, y = Value, col = Model, linetype = Model))+
      facet_grid(int~Effect, scales = 'free')+
      labs(x = 'Relative Years', y = '')+
      my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
               panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
               legend.title=element_blank(), # no legend title
               ...)
    simPlots[[i]] <- p3
  }
  
  individualEffectPlots <- list()
  effects <- c('Age', 'Period', 'Cohort')
  for(i in 1:length(effects)){
    individualEffectDat <- 
      finalDat %>% 
      filter(Effect == effects[i])
    
    # creating the plots
    ## plots of continuous variables: Effects and curvatures
    p1.a <- 
      ggplot()+
      geom_point(data = individualEffectDat %>% filter(!(Statistic %in% c('Bias', 'MSE'))),
                 aes(x = Group, y = Value, col = Model, shape = Model))+
      geom_line(data = individualEffectDat %>% filter(!(Statistic %in% c('Bias', 'MSE'))),
                aes(x = Group, y = Value, col = Model, linetype = Model))+
      facet_grid(int~Effect, scales = 'free')+
      labs(x = 'Relative Years', y = '')+
      my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
               panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
               legend.title=element_blank(), # no legend title
               ...)
    ## plots of discrete variables: curvature bias and MSE box plots
    p2.a <- 
      ggplot()+
      geom_boxplot(data = individualEffectDat %>% filter(Statistic %in% c('Bias')),
                   aes(x = Model, y = Value, col = Model, fill = Model),
                   alpha = 0.2)+
      geom_jitter(data = individualEffectDat %>% filter(Statistic %in% c('Bias')),
                  aes(x = Model, y = Value, col = Model, shape = Model),
                  alpha = 0.5)+
      geom_boxplot(data = individualEffectDat %>% filter(Statistic %in% c('MSE')),
                   aes(x = Model, y = Value, col = Model, fill = Model),
                   alpha = 0.2)+
      geom_jitter(data = individualEffectDat %>% filter(Statistic %in% c('MSE')),
                  aes(x = Model, y = Value, col = Model, shape = Model),
                  alpha = 0.5)+
      facet_grid(int~Effect, scales = 'free')+
      labs(y = '')+
      my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
               panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
               strip.text.x = element_blank(), # removing column (x-axis) facet labels as same as continuous
               legend.title=element_blank(), # no legend title
               ...)
    
    if(boxplots == TRUE){
      ## combining plots vertically
      individualEffectPlots[[i]] <- p1.a/p2.a
    } else {
      individualEffectPlots[[i]] <- p1.a
    }
    
    
  }
  
  return(list(p = p, simPlots = simPlots, individualEffectPlots = individualEffectPlots))
  
}

# Human mortality database application ----

# logit function
logit <- function(x){
  log(x/(1 - x))
}

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
gradient.plot <- function(data, is.log = F){
  
  # data = dataA; is.log = TRUE
  
  dataAug <-
    data %>%
    mutate(rate = y/N,
           logRate = log(y/N))
  
  if(isTRUE(is.log)){
    p <- 
      ggplot(dataAug, aes(x = p, y = a, fill = logRate)) +
      geom_tile() +
      scale_fill_gradientn(colors = colorRamps::matlab.like2(10000))  +
      labs(x = 'Year', y = 'Age')
  } else {
    p <- 
      ggplot(dataAug, aes(x = p, y = a, fill = rate)) +
      geom_tile() +
      scale_fill_gradientn(colors = colorRamps::matlab.like2(10000))  +
      labs(x = 'Year', y = 'Age')
  }
  
  p
  
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
    geom_ribbon(aes(ymin = lower,ymax = upper,fill = Effect), alpha = .2, linetype = 0, show.legend = F)
  
  return(list(plot = plot, termsPred = termsPred)) 
}