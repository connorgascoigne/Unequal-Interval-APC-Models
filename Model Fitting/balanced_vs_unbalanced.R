pssFit <-
  hol.spline.extract(data = data, family = family, knots = knots, # dependent on command line
                     # fixed arguments
                     bs = 'cr', include.periodic = FALSE, 
                     slopeDrop = 'c', mod = 'apc',
                     fixed = list(age = FALSE, period = FALSE, cohort = FALSE))

data = data; family = family; knots = knots; # dependent on command line
# fixed arguments
bs = 'cr'; include.periodic = FALSE; 
slopeDrop = 'c'; mod = 'apc';
fixed = list(age = FALSE, period = FALSE, cohort = FALSE)

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





X1 <- my.predict(gamObj = m, include.periodic = include.periodic, type = 'lpmatrix')
yhat1 <- X1 %*% beta
estimates1 <- find.all.estimates(data = data, y = yhat1)
hat_f_a1 <- estimates1$full$age
hat_f_p1 <- estimates1$full$period
hat_f_c1 <- estimates1$full$cohort
hat_f_a_curv1 <- estimates1$curvature$age
hat_f_p_curv1 <- estimates1$curvature$period
hat_f_c_curv1 <- estimates1$curvature$cohort

estData1 <- 
  data.frame(Effect = 
               factor(c(rep('Age', times = length(hat_f_a1)),
                        rep('Period', times = length(hat_f_p1)),
                        rep('Cohort', times = length(hat_f_c1))),
                      levels = c('Age', 'Period', 'Cohort')))
curvData1 <- 
  data.frame(Effect = 
               factor(c(rep('Age', times = length(hat_f_a_curv1)),
                        rep('Period', times = length(hat_f_p_curv1)),
                        rep('Cohort', times = length(hat_f_c_curv1))),
                      levels = c('Age', 'Period', 'Cohort')))

estData1$Group <- c(a, p, c)
estData1$Estimate <- c(hat_f_a1, hat_f_p1, hat_f_c1)
rownames(estData1) <- 1:sum(length(hat_f_a1), length(hat_f_p1), length(hat_f_c1))

curvData1$Group <- c(a, p, c)
curvData1$Estimate <- c(hat_f_a_curv1, hat_f_p_curv1, hat_f_c_curv1)
rownames(curvData1) <- 1:sum(length(hat_f_a_curv1), length(hat_f_p_curv1), length(hat_f_c_curv1))



estDataFinal <- rbind(estData %>% mutate(data = 'balanced'),
                      estData1 %>% mutate(data = 'unbalanced'))

curvDataFinal <- rbind(curvData %>% mutate(data = 'balanced'),
                       curvData1 %>% mutate(data = 'unbalanced'))

# plotting the results
# the smooth functions of the estimates plot
estPlot <- 
  ggplot(estDataFinal, aes(x = Group, y = Estimate, group = data, color = data)) + 
  geom_point(alpha = 1.2) + 
  geom_line(alpha = 1.2) + 
  facet_wrap(.~Effect, scales = 'free') + 
  # guides(col = 'none') + 
  my.theme()
# the un-centered smooth functions of curvatures plot
curvPlot <- 
  ggplot(curvDataFinal, aes(x = Group, y = Estimate, group = data, color = data)) + 
  geom_point(alpha = 1.2) + 
  geom_line(alpha = 1.2) + 
  facet_wrap(.~Effect, scales = 'free') + 
  # guides(col = 'none') + 
  my.theme()

library(patchwork)
estPlot / curvPlot
