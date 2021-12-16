# Unequal Interval APC Models
 
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
