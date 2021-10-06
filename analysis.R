# analysis functions -----------------------------------------------------------

## theme for plots -------------------------------------------------------------
my.theme <- function(...){
  theme(...,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.text.align = 0,
        legend.key = element_rect(fill = NA))
}

## get the plots from simulation runs ------------------------------------------
get.final.plot <- function(simRun, scenarioChoice){
  
  # extracting the data frames from the simulation output
  simEst <- simRun$simEstData
  simCurv <- simRun$simCurvData
  truEff <- simRun$trueEffData
  trucurv <- simRun$trueCurvData
  
  # binding the estimate and curvature data frames together
  ## simulation
  simRes <- rbind(simEst %>% mutate(Type = 'est'),
                simCurv %>% mutate(Type = 'curv'))
  ## true values
  true <- rbind(truEff %>% mutate(Type = 'est'),
              trucurv %>% mutate(Type = 'curv'))
  
  
  # summary data frame
  ## extract only simulation results for summaries
  eachSim <- simRes %>% select(starts_with("s0"))
  ## data frame
  summary <- 
    simRes %>% 
    select(-starts_with('s0'))
  summary$Estimate <- rowMeans(eachSim)
  summary$Bias <- rowMeans(eachSim-true$Estimate)
  summary$MSE <- rowMeans((eachSim-true$Estimate)^2)
  
  # combining the simulation results and the true results together for specific scenario
  ## simulation results in long format
  simLong <- 
    summary %>% 
    filter(Scenario == scenarioChoice) %>%  # filter for scenario of choice
    pivot_longer(cols = c('Estimate', 'Bias', 'MSE'),
                 names_to = 'Statistic',
                 values_to = 'Value') %>% 
    filter(!(Type == 'est' & Statistic %in% c('Bias', 'MSE')))
  ## true results in long format
  trueLong <- 
    true %>% 
    filter(Scenario == scenarioChoice) %>%  # filter for scenario of choice
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
    labs(x = 'Relative Years')+
    my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
             panel.spacing.y = unit(0.125,'lines') # gap between facets in rows (y-axis)
    )
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
    my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
             panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
             strip.text.x = element_blank() # removing column (x-axis) facet labels as same as continuous
    )
  ## combining plots vertically
  p <- p1/p2
  ## return combined plot
  # p
  
  # have plots for each model with their simulations over
  ## get the simulation data frame in order
  simFitsLong <- 
    simRes %>% 
    filter(Scenario == scenarioChoice) %>%  # filter for scenario of choice
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
      labs(x = 'Relative Years')+
      my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
               panel.spacing.y = unit(0.125,'lines') # gap between facets in rows (y-axis)
      )
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
      labs(x = 'Relative Years')+
      my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
               panel.spacing.y = unit(0.125,'lines') # gap between facets in rows (y-axis)
      )
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
      my.theme(panel.spacing.x = unit(1,'lines'), # gap between facets in columns (x-axis)
               panel.spacing.y = unit(0.125,'lines'), # gap between facets in rows (y-axis)
               strip.text.x = element_blank() # removing column (x-axis) facet labels as same as continuous
      )
    individualEffectPlots[[i]] <- p1.a/p2.a
  }
  
  return(list(p = p, simPlots = simPlots, individualEffectPlots = individualEffectPlots))
  
}
