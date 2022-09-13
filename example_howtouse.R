library(MASS)
library(corrplot)
library(Matrix)
#########################################################

#########################################################

correlation <- 0.7
n_data <- 150

######################################################### 
#####################   DATA
    
correlation_matrix <- matrix(0,20,20)
correlation_matrix[1:5,1:5] <- correlation
correlation_matrix[10:14,10:14] <- correlation
diag(correlation_matrix)<-1
data = data.frame( mvrnorm( n = n_data, mu = rep(0,20), Sigma = correlation_matrix ) )
corrplot( cor(data) )
var_noise <- 0.1
noise = rnorm( n_data , mean = 0, sd = sqrt(var_noise) )
X = cbind( data$X2, data$X11, data$X19, data$X20)
var_beta = var_noise * solve( t(X)%*%X )
sd_beta = sqrt( diag(var_beta) )
ratio = c(4,3,3,5)
beta = sd_beta*ratio
data$y = beta[1]*data$X2 + beta[2]*data$X11 + beta[3]*data$X19 + beta[4]*data$X20 + noise

######################################################### 
#####################   Run Boruta Conditional    
source("boruta_conditional.R")

cb <- Condition_Boruta(y~., data = data, doTrace = 1, getImp=getImpRfCond)
print(cb)

final.boruta <- TentativeRoughFix(cb)
print(final.boruta)

#List of confirmed attributes
getSelectedAttributes(final.boruta, withTentative = F)
plotImpHistory(final.boruta)

## Dataframe of result
boruta.df <- attStats(cb)
print(boruta.df)

# plot the Importance box plot
plot.Boruta(cb)   





