# packages ----

library(demography) # for data
library(tidyverse) # for data wrangling
library(mgcv) # for fitting models
library(patchwork) # for ggplot+ggplot
library(xtable) # r tabels in latex

# path ----

simulationFolder <- 'C:/Users/cg863/OneDrive - University of Bath/Bath PhD/Year 3/Simulation Study/'

## functions ----

source(paste0(simulationFolder, 'Code/functions.R'))

## directories for plots ----

# results folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results'))
}

# plots folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Plots'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Plots'))
}

# Additional Plots folder
if(!dir.exists(paths = paste0(simulationFolder, 'Code/Results/Plots/Additional Plots'))) {
  dir.create(path = paste0(simulationFolder, 'Code/Results/Plots/Additional Plots'))
}

# cyclic pattern plot ----

## arguments for what data ----

family <- c('normal', 'poisson', 'binomial')[1]
dataMod <- c('apc', 'ap', 'ac', 'a')[1]
M <- c(1, 3, 5)[3]
i <- 1

## data load ----

dataList <- readRDS(file = paste0(simulationFolder, 'Data/', family, '/', dataMod, '/', M, '/allData.rds'))
data <- dataList[[i]]

## fit model ----

fit <- hol.factor.extract(data = data, mod = dataMod, slopeDrop = 'c', family = family)

## extract period results only ----

cohData <-
  fit$curvData %>%
  filter(Effect == 'Cohort') %>%
  mutate(Group = Group + abs(min(fit$curvData$Group)))

## plot ----

cyclicPatternPlot <-
  ggplot2::ggplot(cohData, aes(x = Group, y = Estimate)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::geom_vline(xintercept = seq(from = min(cohData$Group), to = max(cohData$Group), by = 5), linetype = 'dotted') +
  ggplot2::labs(x = 'Cohort Group', y = 'Curvature') +
  my.theme(text = element_text(size = 20))

cyclicPatternPlot
  
## saving plot ----

ggsave(cyclicPatternPlot, 
       file = paste0(simulationFolder, 'Code/Results/Plots/Additional Plots/cyclicPatternPlot.png'),
       width = 10, height = 10)


# basis constraint plot ----

## generate data ----

A<-20; P<-1; N<-1; M=1; mod='apc'
family="binomial"; FUNage=age.fun.binomial;  FUNperiod=per.fun.binomial; FUNcohort=coh.fun.binomial
data<-data.sim(A=A, P=P, M=M, N=N, mod=mod, family=family,
               FUNage=FUNage,FUNperiod=FUNperiod,FUNcohort=FUNcohort)

## initial basis ----

# with intercept
sm<-mgcv::smoothCon(s(a, k=5), data=data, absorb.cons=F)[[1]]
Xcheck<-sm$X
int<-rep(1, length=nrow(Xcheck))
Zcheck<-gen.lin.constr(M=int,X=Xcheck)
X<-Xcheck%*%Zcheck
lin<-data$a
Z<-gen.lin.constr(M=cbind(1,lin),X=Xcheck)
Xtilde<-Xcheck%*%Z

Xcheck_df<-data.frame(b0=Xcheck[,4], b1=Xcheck[,5], b2=Xcheck[,1], b3=Xcheck[,2], b4=Xcheck[,3])
X_df<-data.frame(b1=X[,4], b2=X[,3], b3=X[,1], b4=X[,2])
Xtilde_df<-data.frame(b2=Xtilde[,2], b3=Xtilde[,3], b4=Xtilde[,1])

## data frame of results ----

b1<-as.data.frame(Xcheck_df)%>%mutate(id = 1:nrow(data))%>%reshape2::melt(.,id.var="id")%>%mutate(type="Smooth")
b2<-as.data.frame(X_df)%>%mutate(id = 1:nrow(data))%>%reshape2::melt(.,id.var="id")%>%mutate(type="Orthogonal to intercept")
b3<-as.data.frame(Xtilde_df)%>%mutate(id = 1:nrow(data))%>%reshape2::melt(.,id.var="id")%>%mutate(type="Orthogonal to intercept and slope")
b<-rbind(b1,b2,b3)
b$type<-factor(b$type, levels=c('Smooth','Orthogonal to intercept','Orthogonal to intercept and slope'))

## plot ----

baisConstraintPlot<-
  ggplot2::ggplot(b, aes(x=id,y=value,group=interaction(type,variable),colour=variable))+
  ggplot2::geom_line(size=1)+
  ggplot2::facet_grid(variable~type)+
  ggplot2::labs(x='Age', y='Bases')+
  my.theme(legend.position = 'none',
           text = element_text(size = 15))

baisConstraintPlot

## saving plot ----

ggsave(baisConstraintPlot, 
       file = paste0(simulationFolder, 'Code/Results/Plots/Additional Plots/baisConstraintPlot.png'),
       width = 10, height = 10)


# Different orthogonality Experiment ----

## generate data ----

A<-20; P<-1; N<-150; M=1; mod='apc'
family="binomial"; FUNage=age.fun.binomial;  FUNperiod=per.fun.binomial; FUNcohort=coh.fun.binomial
data<-data.sim(A=A, P=P, M=M, N=N, mod=mod, family=family,
               FUNage=FUNage,FUNperiod=FUNperiod,FUNcohort=FUNcohort)

proj.ip <- function(M, X , orth = FALSE, weight = rep(1, nrow(M))){
  Pp <- solve( crossprod(M * sqrt(weight)), t(M*weight) ) %*% X
  PX <- M %*% Pp
  if(orth){
    PX <- X - PX
  } else {
    PX
  }
}

Thin.col <- function(M, tol = 1e-06){
  QR <- qr(M, tol = tol, LAPACK= FALSE)
  M[, QR$pivot[seq(length = QR$rank)], drop = FALSE]
}

## initial basis ----

# with intercept
sm<-mgcv::smoothCon(s(a, k=5), data=data, absorb.cons=F)[[1]]
X<-sm$X

Z<-gen.lin.constr(M = cbind(1,data$a),  X = X)
X1<-X%*%Z

# new projections
X2 <- Thin.col(proj.ip(M = cbind(1, data$a), X = X, orth = TRUE, weight = rep(1, nrow(data))))
X3 <- Thin.col(proj.ip(M = cbind(1, data$a), X = X, orth = TRUE, weight = c(data$y)))

# check
# cbind(1,data$a) %>% t %*% X1
# cbind(1,data$a) %>% t %*% X2
# cbind(1,data$a) %>% t %*% X3

X_df <- data.frame(b0=X[,4], b1=X[,5], b2=X[,1], b3=X[,2], b4=X[,3])
X1_df <- data.frame(b2=X1[,2], b3=X1[,3], b4=X1[,1])
X2_df <- data.frame(b2=X2[,1], b3=X2[,2], b4=X2[,3])
X3_df <- data.frame(b2=X3[,1], b3=X3[,2], b4=X3[,3])

## data frame of results ----

b1<-as.data.frame(X_df)%>%mutate(id = 1:nrow(data))%>%reshape2::melt(.,id.var="id")%>%mutate(type="Smooth")
b2<-as.data.frame(X1_df)%>%mutate(id = 1:nrow(data))%>%reshape2::melt(.,id.var="id")%>%mutate(type="Paper: W = diag(0)")
b3<-as.data.frame(X2_df)%>%mutate(id = 1:nrow(data))%>%reshape2::melt(.,id.var="id")%>%mutate(type="New: W = diag(0)")
b4<-as.data.frame(X3_df)%>%mutate(id = 1:nrow(data))%>%reshape2::melt(.,id.var="id")%>%mutate(type="New: W = diag(y)")
b<-rbind(b1,b2,b3, b4)
b$type<-factor(b$type, levels=c('Smooth', 'Paper: W = diag(0)', 'New: W = diag(0)', 'New: W = diag(y)'))

## plot ----

baisConstraintPlot<-
  ggplot2::ggplot(b, aes(x=id,y=value,group=interaction(type,variable),colour=variable))+
  ggplot2::geom_line(size=1)+
  ggplot2::facet_grid(variable~type)+
  ggplot2::labs(x='Age', y='Bases')+
  my.theme(legend.position = 'none',
           text = element_text(size = 15))

baisConstraintPlot

## saving plot ----

ggsave(baisConstraintPlot, 
       file = paste0(simulationFolder, 'Code/Results/Plots/Additional Plots/differentWeightOrthogonalisation.png'),
       width = 10, height = 10)



