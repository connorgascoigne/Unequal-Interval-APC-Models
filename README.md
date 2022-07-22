# Unequal Interval APC Models

Here are the files to recreate the data and plots from the Manuscript. 

## Folder Structure

```
Simulation Study
│   README.md
│   file001.txt    
│
└───Code
│   │   functions.R
│   │   additional_plots.R
│   │   save_data.R
│   │   run_save_data.slm
│   │   batch_run_save_data.slm
│   │
│   └───Model Fitting
│   |   │   true.R
│   |   │   run_true.slm
│   |   │   batch_run_true.slm
│   |   │   fit_models.R
│   |   │   run_fit_models.slm
│   |   │   batch_run_fit_models.slm
│   └───Model Fitting
│       │   true.R
│       │   run_true.slm
│       │   batch_run_true.slm
│       │   fit_models.R
│       │   run_fit_models.slm
│       │   batch_run_fit_models.slm
│   
└───Data
    │   file021.txt
    │   file022.txt
```

## Use of Cluster

## Practical Guide

- Each .R file follows a similar pattern with imported library's in the first lines followed by manually setting the working directory to a variable
- Setting the working directory to this variable is paramount to ensuring the rest of the code can be run in a ‘press and play’ fashion
- Any subfolder needed for savin (data, images, outputs, etc, ...) will be created within the code

## R Files

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
  - This is an R script to generate all the data needed for the simualtion study
  - This requries a large amount of computation so it is run over a cluster using two additional scripts
  - Explicity uses the `data.sim` function
  - Implictly uses the `age.fun.binomial`, `per.fun.binomial` and `coh.fun.binomial` (and the Gaussian and Poisson evivalent) functions
  
 - run_save_data.slm
  - code to run the save_data.R file for specifc arguments
   - `family` - family for the distribution
   - `datMod` - the type of model used to generate the data 
   - `M` - aggregation of age
  - This is run using slurm
  
 - batch_run_save_data.slm
  - Loops over inputs for `family`, `datMod` and `M`
  


## R Files

- data_generation.R
  - Funtions to create the curves for each temporal function for each case of Gaussian, binomial and Poisson data.
  - `data.sim()` - simulates equal age-period data and (if chosen) aggregates appropriately.
  -  `find.all.estimates()` - takes a data frame with the age, period and cohort combinations and a vector of respones and returns each temporal effect and curavture curve.
  -  `true.values()` - returns a data frame and plots of true values from the temporal functions.
- holford_models.R
  - `gen.lin.const()` - prduces the Z constrain matrix based off of the matrix to constrain and the columsn to orthogonalise against.
  - `my.predictMat()` - produces the model matrix for a spline model version of Holfords re-parameterisation.
  - `my.predict()` - produces predictions for a give data frame (or uses the data supplied to fit the model) using a spline model and Holfords re-parameterisation.
  - `hol.spline.fun()` - fits the spline model based on Holford re-parameterisation.
    - `data` - data frame of data including age, period, cohort, y (response) and N (total) columns.
    - `mod`  - what type of model to fit, can be APC or any sub model.
    - `knots` - list of number of knots for temporal functions.
    - `fixed` - are the splines penalised or not.
    - `colDrop` - if fitting an APC model, what slope is to be dropped. The sub models are re-parameterised in the same manor but include all linear slopes as they do not have the structural link identification problem
    - `distribution` - what is the distribution used to generate the data.
  - `hol.spline.extract()` - produces a data frame and plots of the effect and curvature estimates from the spline (penalised or not) model fit.
  - `hol.factor.X()` - produces the model matrix for a factor model version of Holfords re-parameterisation.
  - `hol.factor.fun()` - fits the factor model based on Holfords re-parameterisation.
  - `hol.factor.extract()` - produces a data frame and plots of the effect and curvature estimates from the factor model fit.
  - `run.holford.sim()` - runs the full simulation study for the factor (FA), regression smoothing spline (RSS, without penalisation) and penalised smoothing spline (PSS, with penalisation) models. 
    -  `fixedParams` - parameters for the data generation that are fixed for all distributions.
    -  `modParams` - parameters for defining the models (i.e. list of knots which changes depending on M).
    -  `distParams` - parameters for each individual distribution.
- analysis.R
  - `my.theme()` - produces the same styke for all plots.
  - `get.final.plot()` - collect the results of siumulation study and the make the plots seen in the manuscript and supplementary material.  
- Example_full_study.R
  - Runs the full simulation from top to bottom.
  - Ensure the working directory is set to source file location which contains all the `.R` files and in the `run.holford.sim()` the file paths match.
  - Choices to make:
    - Fitting the models with equal or unequal data.
    - What distribution to use  
- HMD_application.R
  - `aggregate.data()` - aggregated HMS data.
  - `gradient.plot()` - replicates the gradient plot from the supplementary material.
  - `plot.terms()` - plots the smooth functions of curvature for each temporal term from the model fit.
  - Put in your username and password in order to download the data. The HMS data is not shareable but free to download upon registering.  

## Other Files

- final_supplementary.pdf
  - Supplementary material to the manuscript.
