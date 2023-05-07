# Variable selection wrapper in presence of correlated input variables for Random Forest models, Conditional Boruta
## Highlights
• A novel method for variable selection problem for highly correlated input data

• The method is developed for a very commonly used machine learning model, namely Random Forest model

• Use of Conditional Variable importance measure based on the Random Forest model

• Select all the relevant variables among the input variables that are most relevant for predicting the output variable

• Use of simulated data to validate the approach, comparing the results with standard approaches with embedded variable selection algorithms

We compared the algorithm with the Lasso and Eastic net algorithm along with the original version of the Boruta algorithm. The ratio between the variables correctly selected/ wrongly selected variables related to the increasing correlation level. As the correlation increase the Conditional version of the Boruta algorithm.


## Example 
An example on how to run and the outcome of the algorithm can be seen here.

## Summary
The repository contain:
1. 'boruta_conditional.R': function need be read with R command 'source(boruta_conditional.R)' and than call with 
2. 'example_howtouse.R': simple example on how to use the function. 

## Function Features
The function file contain the following features:
1.'Condition_Boruta(y ~ . , data = data, doTrace = 1, getImp=getImpRfCond)' algorithm funtion
2. 'TentativeRoughFix()': function that fix the variables not yet classified
3. 'getSelectedAttributes()': function that give the subset of the selected variables
4. 'plotImpHistory()': the mean of the variable importance value calculated at each step of the algorithm
5. 'plot.Boruta()': boxplo of the variables conditional importances


### Credits
This repository includes code from the following third-party sources:

[Boruta](https://cran.r-project.org/web/packages/Boruta/index.html) - a package for feature selection in R. The original code for Boruta can be found at [https://gitlab.com/mbq/Boruta/](https://gitlab.com/mbq/Boruta/).

