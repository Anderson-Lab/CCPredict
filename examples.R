library(doParallel)
library(foreach)
cl<-makeCluster(8)
registerDoParallel(cl)

X <- read.csv(system.file("extdata", "X.csv", package="CCPredict"),header=FALSE)
X <- t(X)
X = scale(X,center=T,scale=T) # Scale the X data so it has a mean of 0 and a stdev of 1. Pretty standard
y <- read.csv(system.file("extdata", "y.csv", package="CCPredict"),header=FALSE)
L <- read.csv(system.file("extdata", "L.csv", package="CCPredict"),header=FALSE)
y <- as.matrix(y)
y <- factor(y[,1])
L <- as.matrix(L)
ytr=y
CRange=c(2^-8,2^-4,2^-2,2^0,2^2,2^4,2^8)
LambdaRange=c(1e-8,1e-4,1e-2,1,1e+2,1e+4,1e+8)
optimize.ccSVM(X,y,L,CRange,LambdaRange)
