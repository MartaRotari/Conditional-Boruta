library(MASS)
library(corrplot)
library(Matrix)


#### Simulate some data: X=[X1,..,X20] and Y=beta1*X2 + beta2*X7
n_data <- 200
correlation_matrix <- matrix(0,20,20)
correlation_matrix[1:5,1:5] <- 0.7
diag(correlation_matrix) <- 1
data = data.frame( mvrnorm( n = n_data, mu = rep(0,20), Sigma = correlation_matrix ) )
var_noise = 0.1
noise = rnorm( n_data , mean = 0, sd = sqrt(var_noise) )
X = cbind( data$X2, data$X7)
var_beta = var_noise * solve( t(X)%*%X )
sd_beta = sqrt( diag(var_beta) )
ratio = 3  # The simulated data are such that beta/sd(beta) > ratio, where ratio > 2
beta = sd_beta*ratio
data$y = beta[1]*data$X2 + beta[2]*data$X7 + noise

####  Run Boruta Conditional    
source("boruta_conditional.R")

cb <- Condition_Boruta(y~., data = data, doTrace = 1, getImp=getImpRfCond)
print(cb)

# Fix the unclassified variables
final.boruta <- TentativeRoughFix(cb)
print(final.boruta)
# List of confirmed attributes
getSelectedAttributes(final.boruta, withTentative = F)
# Historical mean of the variables importance value
plotImpHistory(final.boruta)
# Dataframe of result
boruta.df <- attStats(cb)
print(boruta.df)
# plot the Importance box plot
plot.Boruta(cb)   
plot.Boruta(final.boruta)

#### Original Boruta algorithm
library(Boruta)

b_model = Boruta(y~., data = data, doTrace = 1)
print(b_model)
