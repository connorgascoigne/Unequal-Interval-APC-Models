# Unequal Interval APC Models

Here are the files to recreate the data and plots from the Manuscript. 

## Folder Structure

### Simulation Study

#### Simulation Study/Data

- This folder is created in save_data.R
- This is where the data is stored once generated 

#### Simulation Study/Code

- functions.R
  - These are all the functions needed to create the data, fit the models and create the plots
  - `age.fun.binomial`, `per.fun.binomial` and `coh.fun.binomial` - the true age, period and cohort functions used to generate binomial data. There are analogous functions for the Gaussian and Possion data
  - `data.sim()` - simulates equal age-period data and (if chosen) aggregates appropriately
  - `find.all.estimates()` - takes a data frame with the age, period and cohort combinations and a vector of respones and returns each temporal effect and curavture
  - `true.values()` - returns a data frame and plots of true values from the temporal functions
  - `gen.lin.const()` - prduces the Z constraint matrix based off of the matrix to constrain and the columsn to orthogonalise against
  - `my.predictMat()` - produces the model matrix for a spline model version of Holfords reparameterisation
  - `my.predict()` - produces predictions for a give data frame (or uses the data supplied to fit the model) using a spline model and Holfords reparameterisation
  - `hol.spline.fun()` - fits the spline model based on Holford reparameterisation
    - `data` - data frame of data including age, period, cohort, y (response) and N (total) columns
    - `mod`  - what type of model to fit, can be APC or any sub model
    - `bs` - the choice of basis to use (fixed on cubic regression spline)
    - `include.periodic` - include additional cylic regression bases in the period and cohort terms
    - `knots` - list of number of knots for temporal functions
    - `fixed` - are the splines penalised or not
    - `colDrop` - if fitting an APC model, what slope is to be dropped. The sub models are reparameterised in the same manor but include all linear slopes as they do not have the structural link identification problem
    - `family` - what is the distribution used to generate the data
    - `hol.spline.extract()` - produces a data frame and plots of the effect and curvature estimates from the spline (penalised or not) model fit
    - `hol.factor.X()` - produces the model matrix for a factor model version of Holfords reparameterisation
    - `hol.factor.fun()` - fits the factor model based on Holfords reparameterisation.
    - `hol.factor.extract()` - produces a data frame and plots of the effect and curvature estimates from the factor model fit
  - `my.theme()` - produces the same style for all plots
  - `get.final.plot()` - collect the results of siumulation study and the make the plots seen in the Manuscript and Supplementary Material
  - `logit` - logit link function
  - `aggregate.data()` - aggregated HMD data
  - `gradient.plot()` - replicates the gradient plot from the supplementary material
  - `plot.terms()` - plots the smooth functions of curvature for each temporal term from the model fit
 
- additional_plots.R
  - This contains code to plot additional plots not related to the simulation study or HMD application in the Manuscript
  - Figure 1 - Cylic pattern in cohort curvature  
  - Figure 2 - how the basis changes after different orthogonalization
  - Experimenting how the basis changes with different types of orthogonalization procedures
 
- HMD_application.R
  - This loads in the HMD data using `demography::hmd.mx` but needs a username and password which are giving after registration
  - Runs all the code relating to the Application in both the Manuscript and Supplementary Material
  
- save_data.R
  - This is an R script to generate all the data needed for the simualtion study for specifc `family`, `datMod` and `M`
  - This requries a large amount of computation so it is run over a cluster using two additional scripts
  - Explicity uses the `data.sim` function
  - Implictly uses the `age.fun.binomial`, `per.fun.binomial` and `coh.fun.binomial` (and the Gaussian and Poisson evivalent) functions
  
- run_save_data.slm
  - code to run the save_data.R file for specifc `family`, `datMod` and `M`
  
- batch_run_save_data.slm
  - Loops over run_save_data.slm inputs for `family`, `datMod` and `M`

##### Simulation Study/Code/Model Fitting

- true.R
  - Generates the true data for specifc `family`, `datMod` and `M`
- run_true.slm
  - code to run the true.R file for specifc `family`, `datMod` and `M`
- batch_run_true.slm
  - Loops over run_true.slm inputs for `family`, `datMod` and `M`
- fit_models.R
  - Fits all models for specifc `family`, `datMod`, `M` and `dataSet`
- run_fit_models.slm
  - code to run the fit_models.R file for specifc `family`, `datMod`, `M` and `dataSet`
- batch_run_fit_models.slm
  - Loops over run_fit_models.slm inputs for `family`, `datMod`, `M` and `dataSet`

##### Simulation Study/Code/Results Collecting

- collect_results.R
  - Collects the results from each of the 100 `dataSet` simulations and combines them together by `family`, `datMod` and `M`
  - Stores in final data sets in Simulation Study/Code/Results/Final/...
- run_collect_results.slm
  - Code to run the collect_results.R file for specifc `family`, `datMod` and `M`
- batch_run_collect_results.slm
  - Loops over run_collect_results.slm inputs for `family`, `datMod` and `M`
- plot_results
  - Plots the results for each `family`, `datMod` and `M` 
  - Stores the plots in Simulation Study/Code/Results/Plots/...
- run_plot_results.sh
  - Code to run the plot_results.R file for specifc `family`, `datMod` and `M` 

##### Simulation Study/Code/Results

- This folder and subfolders are created in various functions within Model Fitting and Results Collecting

###### Simulation Study/Code/Results/Final

- The final data

###### Simulation Study/Code/Results/Plots

- The final plots

## Use of Cluster

- The simulation study is implemented over a cluster
- The format is as follows
    - XXX.R - perform an action (e.g., fit a model or generate data). This requires a number of arguments (e.g., type of APC model, family for distribution, aggregation type)
    - run_XXX.slm - passes one set of arguments to XXX.R to run one instance (e.g., fit one model for the given arguments)
    - batch_run_XXX.slm - loops over all combinations of the arguments and supplies this to run_XXX.slm which, in turn, supplies this to XXX.R   

## Practical Guide

- Each .R file follows a similar pattern with imported library's in the first lines followed by manually setting the working directory to a variable
- Setting the working directory to this variable is paramount to ensuring the rest of the code can be run in a ‘press and play’ fashion
- The working directory will change depending upon whether running off of local machine or cluster
- Files can run on a local machine for a small number of simualtions or combination of different arguments, but we suggest using a cluster to fully run the simulation
- Any subfolder needed for saving (data, images, outputs, etc, ...) will be created within the code


