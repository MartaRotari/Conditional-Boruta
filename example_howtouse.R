library(MASS)
library(corrplot)
library(Matrix)

#### Simulate some data: 
#  X = [X1,..,X20] multivariate normal distribute N(0,I). 1 correlation group X1,..,X5 cor = 0.8
#  Y = b1*X2 + b2*X7 + e

set.seed(123456789)
n_data <- 200
correlation_matrix <- matrix(0,20,20)
correlation_matrix[1:5,1:5] <- 0.8
diag(correlation_matrix) <- 1
data = data.frame( mvrnorm( n = n_data, mu = rep(0,20), Sigma = correlation_matrix ) )
var_noise = 0.1
noise = rnorm( n_data , mean = 0, sd = sqrt(var_noise) )
X = cbind( data$X2, data$X7)

# The simulated data are such that beta/sd(beta)=k, where k=(4,3) to ensure that the coefficients are statistically significant and appropriately scaled.
var_beta = var_noise * solve( t(X)%*%X )
sd_beta = sqrt( diag(var_beta) )
ratio = c(5,4)  
beta = sd_beta*ratio
data$y = beta[1]*data$X2 + beta[2]*data$X7 + noise


####  Run Conditional Boruta     
source('boruta_conditional.R')
cb_model = Condition_Boruta(y ~ ., data = data, doTrace = 1, getImp=getImpRfCond,mtry=20, ntree = 500)
#cb_model <- Condition_Boruta(y ~ ., data = data, mtry=20, ntree = 500, doTrace = 1, getImp=getImpRfCond)
print(cb_model)

# Decide on the unclassified variables
final.boruta <- TentativeRoughFix(cb_model)
print(final.boruta)
# List of selected variables
getSelectedAttributes(final.boruta, withTentative = F)
# The mean of the variable's importance scores over the iterations
plotImpHistory(final.boruta)
# Final variables importance table
boruta.df <- attStats(cb_model)
print(boruta.df)
# boxplot of the variables importance scores
plot.Boruta(cb_model)   
plot.Boruta(final.boruta)










#### Original Boruta algorithm could be find at https://cran.r-project.org/web/packages/Boruta/index.html

