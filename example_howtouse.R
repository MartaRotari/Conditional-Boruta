library(MASS)
library(corrplot)
library(Matrix)
#########################################################

#########################################################

correlation <- 0.7


######################################################### 
#####################   DATA
    
correlation_matrix <- matrix(0,20,20)
correlation_matrix[1:5,1:5] <- correlation
correlation_matrix[10:14,10:14] <- correlation
diag(correlation_matrix)<-1
data = data.frame( mvrnorm( n = 250, mu = rep(0,20), Sigma = correlation_matrix ) )
corrplot( cor(data) )
epsilon = rnorm(250, mean = 0, sd = 0.1)
var_noise <- 0.1
noise = rnorm( 250 , mean = 0, sd = sqrt(var_noise) )
X = cbind( data$X2, data$X11, data$X19, data$X20)
var_beta = var_noise * solve( t(X)%*%X )
sd_beta = sqrt( diag(var_beta) )
ratio = c(4,3,3,5)
beta = sd_beta*ratio
data$y = beta[1]*data$X2 + beta[2]*data$X11 + beta[3]*data$X19 + beta[4]*data$X20 + noise

######################################################### 
#####################   Run Boruta Conditional    
source("cond_boruta/Boruta_cond.R")

cb <- Boruta_Conditional(y~., data = data, doTrace = 1)
final.boruta <- TentativeRoughFix(cb)
print(final.boruta)

#List of confirmed attributes
getSelectedAttributes(final.boruta, withTentative = F)
plotImpHistory(final.boruta)

## Dataframe of result
boruta.df <- attStats(cb)
print(boruta.df)

# plot the Importance box plot
plot(cb, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(cb$ImpHistory),function(i) cb$ImpHistory[is.finite(cb$ImpHistory[,i]),i])
names(lz) <- colnames(cb$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),at = 1:ncol(cb$ImpHistory), cex.axis = 0.7)

    


######################################################### 
#####################   Run Backward Boruta    
source("back_boruta/Boruta_Backward.R")

bb <- Boruta_Backward(y~., data = data, doTrace = 1)
final.boruta <- TentativeRoughFix(bb)
print(final.boruta)

#List of confirmed attributes
getSelectedAttributes(final.boruta, withTentative = F)
plotImpHistory(final.boruta)

## Dataframe of result
boruta.df <- attStats(bb)
class(boruta.df)
print(boruta.df)

# plot the Importance box plot
plot(bb, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(bb$ImpHistory),function(i) bb$ImpHistory[is.finite(bb$ImpHistory[,i]),i])
names(lz) <- colnames(bb$ImpHistory)
Labels <- sort(sapply(lz,median))
axis(side = 1,las=2,labels = names(Labels),at = 1:ncol(bb$ImpHistory), cex.axis = 0.7)









