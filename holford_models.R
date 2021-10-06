# for spline model -------------------------------------------------------------

## general linear constraint <- - for holford M = linear term ------------------
gen.lin.constr <- function(M, X){
  
  # goal t(M) %*% (X %*% Z)  =  0
  C <- t(M) %*% X
  qrc <- qr(t(C))
  Z <- qr.Q(qrc, complete = TRUE)[,(nrow(C) + 1):ncol(C)] # define Z the last n-1 columns of find the Q
  Z
}

## updated to have re-parameterisation with cbind(1,linear) --------------------
my.PredictMat <- function(smObj, data, n = nrow(data)){
  
  # get the correct term and the basis dimensions
  term <- smObj$term
  bs.dim <- smObj$bs.dim
  data$x <- data[,term] # needed to relabel to fit gam
  
  # smooth basis matrix with intercept constraint
  # need to increase the basis dimension by 1
  # to then orthogonalise correctly
  sm <- mgcv::smoothCon(s(x, bs = 'cr', k = bs.dim + 1), data = data, absorb.cons = F)[[1]]
  X <- cbind(lab = data[,term], sm$X)
  Xuni <- unique(X)[,-1]
  smS <- sm$S[[1]]
  
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
  
  list(cur = cur, S = S)
}

## my prediction function for gamObjs ------------------------------------------
my.predict <- function(gamObj, newdata = NULL, type = 'response', exclude = NULL, inverse.link = FALSE){
  
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
        predMats <- my.PredictMat(gamObj$smooth[[k]], data) # my.PredictMat function defines smooth then reparameterises it
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

## fitting the spline model -----------------------------------------------------

hol.spline.fun <- function(data, mod = c('apc','ac','ap','pc','a','p','c'),
                         knots = NULL, fixed = NULL, colDrop = NULL, 
                         distribution = c('normal','poisson','binomial')){
  
  if(is.null(knots)|!is.null(knots) & !is.list(knots)) stop('Warning: Need a list for number of knots for each effect labelled age, period, cohort.')
  if(is.null(fixed)|!is.null(fixed) & !is.list(fixed)) stop('Warning: Need a list for whether an effect is penalised or not.')
  if(missing(distribution)) stop('Warning: Need to specify the distributions (normal or poisson)')
  if(mod == 'apc'&&is.null(colDrop)) stop('Warning: When fitting an APC model need to drop linear column of one of age (a), period (p) or cohort (c).')
  if(mod == 'apc'&&!(colDrop %in% c('a', 'p', 'c'))) stop('colDrop needs to be one of a, p or c')
  
  # data augmented 
  dataAug <- data
  # knots
  ak <- knots$age; pk <- knots$period; ck <- knots$cohort
  # fixed for penalty
  afx <- fixed$age; pfx <- fixed$period; cfx <- fixed$cohort # default is to have fixed = FALSE for penalisation
  
  formula <- 
    as.formula(y ~ 
                 a + p + c + 
                 s(a, bs = 'cr', k = ak-1, fx = afx) + 
                 s(p, bs = 'cr', k = pk-1, fx = pfx) + 
                 s(c, bs = 'cr', k = ck-1, fx = cfx))
  # update pre and post forumla for the temporal models
  # # do not remove any linear terms from pre fit as they are all needed later in 
  # # re-parameterisation
  if (mod == 'apc'){
    if(colDrop == 'a'){
      formula <- update(formula, ~.-a)
    } else if (colDrop == 'p'){
      formula <- update(formula, ~.-p)
    } else if (colDrop == 'c'){
      formula <- update(formula, ~.-c)
    }
  } else if (mod == 'ac'){
    formula <- update(formula, ~.-p-s(p, bs = 'cr', k = pk-1, fx = pfx))
  } else if (mod == 'pc'){
    formula <- update(formula, ~.-a-s(a, bs = 'cr', k = ak-1, fx = afx))
  } else if (mod == 'ap'){
    formula <- update(formula, ~.-c-s(c, bs = 'cr', k = ck-1, fx = cfx))
  } else if (mod == 'c'){
    formula <- update(formula, ~.-a-s(a, bs = 'cr', k = ak-1, fx = afx)-p-s(p, bs = 'cr', k = pk-1, fx = pfx))
  } else if (mod == 'a'){
    formula <- update(formula, ~.-p-s(p, bs = 'cr', k = pk-1, fx = pfx)-c-s(c, bs = 'cr', k = ck-1, fx = cfx))
  } else if (mod == 'p'){
    formula <- update(formula, ~.-a-s(a, bs = 'cr', k = ak-1, fx = afx)-c-s(c, bs = 'cr', k = ck-1, fx = cfx))
  }
  
  # use all cores avaliable
  ctrl <- list(nthreads = parallel::detectCores())
  
  # the extra information depending on the model
  if(distribution == 'normal'){
    fit <- mgcv::gam(formula, family = 'gaussian', data = dataAug, method = 'REML', control = ctrl, fit = F)
  } else if(distribution == 'poisson'){
    fit <- mgcv::gam(formula, offset = log(N), family = 'poisson', data = dataAug, method = 'REML', control = ctrl, fit = F)
  } else if(distribution == 'binomial'){
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
      mats <- my.PredictMat(smObj = smObj,data = dataAug)
      # updating the fit with re-parameterised TEMPROAL terms
      Xpos <- smObj$first.para:smObj$last.para
      fit$X[,Xpos] <- mats$cur
      if(smObj$fixed == FALSE){
        fit$S[[i]] <- mats$S
      }
    }
  }
  
  mod = mgcv::gam(G = fit)
  
  return(list(mod = mod, X = fit$X, S = fit$S))
}

## extract results from the spline function ------------------------------------

hol.spline.extract <- function(data, mod = c('apc','ac','ap','pc','a','p','c'),
                               knots = NULL, fixed = NULL, colDrop = NULL, 
                               distribution = c('normal','poisson','binomial')){
  
  # the reparameterised fit
  fit <- hol.spline.fun(data = data, mod = mod, knots = knots, fixed = fixed, colDrop = colDrop, distribution = distribution)
  
  # model fit
  m <- fit$mod
  
  # defining balanced data frame to define a new X
  ## the unique values of a, p and c used in the model
  a <- sort(unique(data$a)); p <- sort(unique(data$p)); c <- sort(unique(data$c))
  ## a balanced data frame of the a, p and c values used in the model
  data2 <- expand.grid(a = a, p = p, c = c)
  ## new X defined at the same values as the model but with balanced data
  X <- my.predict(gamObj = m,newdata = data2,type = 'lpmatrix')
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

# for factor model -------------------------------------------------------------

## function to make the factor matrix based on model ---------------------------

hol.factor.X <- function(data, mod, colDrop = NULL){
  
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
    if(colDrop == 'a'){
      X <- cbind(1, pLin, cLin, aCur, pCur, cCur)
    } else if(colDrop == 'p'){
      X <- cbind(1, aLin, cLin, aCur, pCur, cCur)
    } else if(colDrop == 'c'){
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

## fit a factor based model ----------------------------------------------------

hol.factor.fun <- function(data, mod = c('apc','ac','ap','pc','a','p','c'),
                         colDrop = NULL, distribution = c('normal','poisson','binomial')){
  
  if(missing(distribution)) stop('Warning: Need to specify the distributions (normal or poisson)')
  if(mod == 'apc'&&is.null(colDrop)) stop('Warning: When fitting an APC model need to drop linear column of one of age (a), period (p) or cohort (c).')
  if(mod == 'apc'&&!(colDrop %in% c('a', 'p', 'c'))) stop('colDrop needs to be one of a, p or c')
  
  # data augmentation
  dataAug <- data
  
  # define the factor model matrix 
  X <- hol.factor.X(data = dataAug, mod = mod, colDrop = colDrop)
  
  # use all cores avaliable
  ctrl <- list(nthreads = parallel::detectCores())
  
  # the extra information depending on the model
  if(distribution == 'normal'){
    mod <- mgcv::gam(y ~ -1 + X, family = 'gaussian', data = dataAug, method = 'REML', control = ctrl)
  } else if(distribution == 'poisson'){
    mod <- mgcv::gam(y ~ -1 + X, offset = log(N), family = 'poisson', data = dataAug, method = 'REML', control = ctrl)
  } else if(distribution == 'binomial'){
    mod <- mgcv::gam(cbind(y, N-y) ~ -1 + X, family = 'binomial', data = dataAug, method = 'REML', control = ctrl)
  }
  
  return(list(mod = mod, X = X))
}

## extract the results from a factor based model -------------------------------

hol.factor.extract <- function(data, colDrop = NULL, 
                               distribution = c('normal','poisson','binomial'),
                               mod = c('apc', 'ac', 'ap', 'pc', 'a', 'p', 'c')){
  
  # the reparameterised fit
  fit <- hol.factor.fun(data = data, mod = mod, colDrop = colDrop, distribution = distribution)
  
  # model fit
  m <- fit$mod
  
  # defining balanced data frame to define a new X
  ## the unique values of a, p and c used in the model
  a <- sort(unique(data$a)); p <- sort(unique(data$p)); c <- sort(unique(data$c))
  ## a balanced data frame of the a, p and c values used in the model
  data2 <- expand.grid(a = a, p = p, c = c)
  ## new X defined at the same values as the model but with balanced data
  X <- hol.factor.X(data = data2, mod = mod, colDrop = colDrop)
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

# run all the models -----------------------------------------------------------

# run a simulation for all models related to holford
## all variants of the holford model (FA, RS and PRS)
run.holford.sim <- function(fixedParams, modParams, distParams){
  
  # Fixed Effects for all simulations
  list2env(fixedParams,globalenv())
  
  # Fixed Effects for all simulations
  list2env(modParams,globalenv())
  
  # functions for data gen
  list2env(distParams,globalenv())
  
  # setting up environment for parallel computing
  ncores <- parallel::detectCores()
  cl <- parallel::makeCluster(ncores-1)
  doParallel::registerDoParallel(cl)
  
  comb_c <- function(...) {
    mapply('cbind', ..., SIMPLIFY = FALSE)
  }
  
  comb_r <- function(...) {
    mapply('rbind', ..., SIMPLIFY = FALSE)
  }
  
  simresults <- foreach::foreach(i = scenario,
                               .packages = c('tidyverse', 'mgcv', 'zoo'),
                               .combine = 'comb_r', 
                               .multicombine = TRUE) %:% 
    foreach::foreach(j = 1:nSim,
                     .packages = c('tidyverse', 'mgcv', 'zoo'),
                     .combine = 'comb_c', 
                     .multicombine = TRUE) %dopar% {
                       
                       # # working directiory where all the files are stored
                       # wd <- getwd()
                       # 
                       # # CHANGE THESE TO RELEVENT FILE PATH
                       # source(paste0(wd,'/data_generation_v2.R'))
                       # source(paste0(wd,'/holford_all_models.R'))
                       # source(paste0(wd,'/analysis_v2.R'))
                       
                       # CHANGE THESE TO RELEVENT FILE PATH
                       source('data_generation.R')
                       source('holford_models.R')
                       source('analysis.R')
                       
                       # source('C:\Users\cg863\OneDrive - University of Bath\Bath PhD\Year 3\Simulation Study\Final Code\Simulation Study\GitHub\data_generation_v2.R')
                       # source('C:\Users\cg863\OneDrive - University of Bath\Bath PhD\Year 3\Simulation Study\Final Code\Simulation Study\GitHub\holford_all_models.R')
                       # source('C:\Users\cg863\OneDrive - University of Bath\Bath PhD\Year 3\Simulation Study\Final Code\Simulation Study\GitHub\analysis_v2.R')
                       
                       # Fixed Effects for all simulations
                       list2env(fixedParams,globalenv())
                       
                       # Fixed Effects for all simulations
                       list2env(modParams,globalenv())
                       
                       # functions for data gen
                       list2env(distParams,globalenv())
                       
                       fixedT = list(age = TRUE, period = TRUE, cohort = TRUE)
                       fixedF = list(age = FALSE, period = FALSE, cohort = FALSE)
                       
                       # setting colDrop
                       # as specified in paper
                       if(i == 'apc'|i == 'ap'){
                         colDrop = 'c'
                       } else if (i == 'pc'){
                         colDrop = 'a'
                       } else if (i == 'ac'){
                         colDrop = 'p'
                       }
                       
                       # geneate the data
                       data <- data.sim(A = A, P = P, M = M, N = N, distribution = distribution, mod = i, 
                                      FUNage = FUNage, FUNperiod = FUNperiod, FUNcohort = FUNcohort)
                       
                       # fit the models
                       facMod <- hol.factor.extract(data = data, colDrop = colDrop, distribution = distribution, mod = mod)
                       rssMod <- hol.spline.extract(data = data, mod = mod, knots = knots, fixed = fixedT, colDrop = colDrop, distribution = distribution)
                       pssMod <- hol.spline.extract(data = data, mod = mod, knots = knots, fixed = fixedF, colDrop = colDrop, distribution = distribution)
                       
                       # estimate results
                       ## preparing
                       facEstRes <- cbind(facMod$estData[,1:2],'FA', i, facMod$estData[,3])
                       rssEstRes <- cbind(rssMod$estData[,1:2], 'RSS', i, rssMod$estData[,3])
                       pssEstRes <- cbind(pssMod$estData[,1:2], 'PSS', i, pssMod$estData[,3])
                       ## naming
                       colnames(facEstRes) <- c('Effect', 'Group', 'Model', 'Scenario', c(sprintf('s%05d',j)))
                       colnames(rssEstRes) <- c('Effect', 'Group', 'Model', 'Scenario', c(sprintf('s%05d',j)))
                       colnames(pssEstRes) <- c('Effect', 'Group', 'Model', 'Scenario', c(sprintf('s%05d',j)))
                       
                       # curvature results
                       ## preparing
                       facCurvRes <- cbind(facMod$curvData[,1:2],'FA', i, facMod$curvData[,3])
                       rssCurvRes <- cbind(rssMod$curvData[,1:2], 'RSS', i, rssMod$curvData[,3])
                       pssCurvRes <- cbind(pssMod$curvData[,1:2], 'PSS', i, pssMod$curvData[,3])
                       ## naming
                       colnames(facCurvRes) <- c('Effect', 'Group', 'Model', 'Scenario', c(sprintf('s%05d',j)))
                       colnames(rssCurvRes) <- c('Effect', 'Group', 'Model', 'Scenario', c(sprintf('s%05d',j)))
                       colnames(pssCurvRes) <- c('Effect', 'Group', 'Model', 'Scenario', c(sprintf('s%05d',j)))
                       
                       # combine results
                       list(
                         est = rbind(facEstRes, rssEstRes, pssEstRes),
                         curv = rbind(facCurvRes, rssCurvRes, pssCurvRes)
                       )
                       
                     }
  
  
  
  trueresults <- foreach::foreach(i = scenario,
                                .packages = c('tidyverse', 'mgcv', 'zoo'),
                                .combine = 'comb_r', 
                                .multicombine = TRUE) %dopar% {
                                  
                                  # CHANGE THESE TO RELEVENT FILE PATH
                                  source('data_generation.R')
                                  source('holford_models.R')
                                  source('analysis.R')
                                  
                                  # Fixed Effects for all simulations
                                  list2env(fixedParams,globalenv())
                                  
                                  # Fixed Effects for all simulations
                                  list2env(modParams,globalenv())
                                  
                                  # functions for data gen
                                  list2env(distParams,globalenv())
                                  
                                  # true values for data model
                                  trueVals <- true.values(A = A, P = P, M = M, mod = i, FUNage = FUNage, FUNperiod = FUNperiod, FUNcohort = FUNcohort)
                                  
                                  list(
                                    eff = trueVals$trueEffData %>% mutate(Scenario = as.factor(i)),
                                    curv = trueVals$trueCurvData %>% mutate(Scenario = as.factor(i))
                                  )
                                  
                                }
  
  # Stop the clusters
  stopCluster(cl)
  
  # collecting results for outputs
  simEstData <- cbind(simresults$est[,1:4], simresults$est %>% select(starts_with('s0'))) %>% mutate(Scenario = as.factor(Scenario))
  simCurvData <- cbind(simresults$curv[,1:4], simresults$curv %>% select(starts_with('s0'))) %>% mutate(Scenario = as.factor(Scenario))
  
  trueEffData <- trueresults$eff
  trueCurvData <- trueresults$curv
  
  allResults <- list(simEstData = simEstData,
                   simCurvData = simCurvData,
                   trueEffData = trueEffData,
                   trueCurvData = trueCurvData)
  
  allResults
  
  
  
}


