# Variable selection method in the presence of correlated input variables

It is a variable selection method developed as a wrapper around the random forest model. Its purpose is to address the challenges posed by highly correlated input data, extending the original version of the Boruta model. For more details, please find a full description [here](document.pdf). 

## Highlights
* The model is a variable selection method designed for highly correlated data
* It is designed as a wrapper around the Random Forest model, providing a powerful and effective solution to the challenges posed by variables selection
* The model has been developed to address the unique complexities that arise from highly correlated input data, thereby enabling the identification of the most relevant variables.
* Rigorous evaluation of the model against established methods such as Lasso, Elastic net, VSURF, and Knockoffs variable selection has demonstrated its effectiveness.

## Example 
An example of how to run the function and analyse the outcome of the algorithm can be seen here.


## Summary
The repository contains the following:
1. 'boruta_conditional.R': function needs to be read with the R command 'source(boruta_conditional.R)'.
2. 'example_howtouse.R': a simple example of how to use the function. 
3. More details can be found in the [paper](document.pdf).
4. The comparison of the model to the established methods and a suggestion on how to choose the [parameters](parameters_comparison.pdf). 


## Function Features
The function file contains the following features:
1.'Condition_Boruta(y ~ . , data = data, doTrace = 1, getImp=getImpRfCond)' algorithm function
2. 'TentativeRoughFix()': the function that fixes the variables not yet classified
3. 'getSelectedAttributes()': the function that gives the subset of the selected variables
4. 'plotImpHistory()': the mean of the variable importance value calculated at each step of the algorithm
5. 'plot.Boruta()': boxplot of the variables' conditional importance