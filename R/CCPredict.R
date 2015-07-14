library(kopls)
library(kernlab)
library(AUC)
#library(modeest)
library(permute)
library(foreach)
library(doParallel)

# TODO:
# 1. There are now much better and cleaner functions for running opls and svm.
#    These should be used in the optimization routines
# 2. Documentation
# 3. Examples
# 4. Cleaning code
# 5. What else is missing from the package?


#' Optimize cckOPLS parameters
#'
#' Given a range of nox and lambda values, optimize them by maximizing
#' the mean AUC (of the ROC) across all folds.
#'
#' @param X - Input data matrix (NxD)
#' @param ytr - Label matrix (Nx2)
#' @param L - Side information kernel matrix
#' @param noxRange - Range over which to optimize nox parameter (kopls)
#' @param LambdaRange - Range over which to optimize lambda parameter (cc)
#' @param kfold - Number of folds for CV
#' @return lambda, nox - optimized lambda and nox values
#'
#' @examples
#' X <- read.csv(system.file("extdata", "X.csv", package="CCPredict"),header=FALSE)
#'
#' @export
#'
optimize.cckopls <- function(X,ytr,L,noxRange,LambdaRange,kfold=2,cluster.size=8){ #optimize cckopls/kopls params
  cl<-makeCluster(cluster.size)
  registerDoParallel(cl)
  
  #kcauc <- matrix(0, nrow=length(noxRange),ncol=kfold)
  test.inxs <- generate.test.inxs(nrow(X),kfold)

  print('optimizing nox...')
  kcauc <- foreach(i=1:length(noxRange),.packages=c('kernlab','AUC','kopls','CCPredict'),.combine=rbind) %dopar% {
  #for (i in 1:length(noxRange)){
    n <- noxRange[i]
    kcauc.values <- c()
    for (j in 1:kfold){
      kcauc.values[j] <- 0
      test <- test.inxs[[j]]
      #     K <- as.kernelMatrix(crossprod(t(X[-test,])))
      K <- as.kernelMatrix(crossprod(t(X)))
      #modelCV <- koplsCV(K,ytr,1,10,nrcv=7,cvType='nfold',preProcK='mc',preProcY='mc',modelType='da')
      #first is for microarray, second for tb, something about y being different, need to go back and fix
      modelOrg <- koplsModel(K[-test,-test],ytr[-test,],1,n,preProcK='mc',preProcY='mc')
      #modelOrg <- koplsModel(K[-test,-test],ytr[-test],1,n,preProcK='mc',preProcY='mc')
      #     modelOrg <- koplsModel(K,ytr,1,n,'mc','mc')
      modelOrgPred<-koplsPredict(K[test,-test],K[test,test],K[-test,-test],modelOrg,n,rescaleY=TRUE)
      #     modelOrgPred<-koplsPredict(K,K,K,modelOrg,rescaleY=TRUE)
      labels <- factor(ytr[test,2])
      kcauc.values[j] <- auc(roc(modelOrgPred$Yhat[,2],labels))
    }
  return(kcauc.values)
  }

  b <- which.max(rowMeans(kcauc))
  
  nox <- noxRange[b[1]]
  print('finished')

  #kcauc <- matrix(0, nrow=length(LambdaRange),ncol=kfold)
  test.inxs <- generate.test.inxs(nrow(X),kfold)
  
  print('optimizing lambda...')
  kcauc <- foreach(i=1:length(LambdaRange),.packages=c('kernlab','AUC','CCPredict'),.combine=rbind) %dopar% {
  #for (i in 1:length(LambdaRange)){
    lambda <- LambdaRange[i]
    kcauc.values <- c()
    for (j in 1:kfold){
      kcauc.values[j] <- 0
      rescaled <- rescaling(X,L,lambda)
      X.new <- rescaled[[1]]
      K.new <- rescaled[[2]]
      # n.list <- rescaled[[3]]
      test <- test.inxs[[j]]
      #modelCV <- koplsCV(K.new,ytr,1,10,nrcv=7,cvType='nfold',preProcK='mc',preProcY='mc',modelType='da')
      modelOrg <- koplsModel(K.new[-test,-test],ytr[-test,],1,nox,'mc','mc')
      modelOrgPred<-koplsPredict(K.new[test,-test],K.new[test,test],K.new[-test,-test],modelOrg,rescaleY=TRUE)
      labels <- factor(ytr[test,2])
      kcauc.values[j] <- auc(roc(modelOrgPred$Yhat[,2],labels))
    }
  return (kcauc.values)
  }
  
  b <- which.max(rowMeans(kcauc))
  
  lambda <- LambdaRange[b[1]]
  
  print('finished')
  
  return(c(lambda,nox))
  

} #end of cckopls opt

#' Optimize ccSVM parameters
#' 
#' @param X - Input data matrix (NxD)
#' @param ytr - Label vector (Nx1)
#' @param L - Side information kernel matrix
#' @param CRange - Range over which to optimize C parameter (svm)
#' @param LambdaRange - Range over which to optimize lambda parameter (cc)
#' @param kfold - Number of folds for CV
#' @return lambda, C - the optimized lambda and C values
#'
#' @examples
#' X <- read.csv(system.file("extdata", "X.csv", package="CCPredict"),header=FALSE)
#'
#' @export
#'
optimize.ccSVM <- function(X,ytr,L,CRange,LambdaRange,kfold=2,cluster.size=8){ #optimize ccSVM params
  cl<-makeCluster(cluster.size)
  registerDoParallel(cl)
  
  #kcauc <- matrix(0,nrow=length(CRange),ncol=kfold) #optimize C
  
  test.inxs <- generate.test.inxs(nrow(X),kfold)
  print('optimizing C...')
  kcauc <- foreach(i=1:length(CRange),.packages=c('kernlab','AUC','CCPredict'),.combine=rbind) %dopar% {
  #for (i in 1:length(CRange)) {
    C <- CRange[i]
    #print(c)
    #results = foreach(j=1:kfold,.packages=c('kernlab','AUC','CCPredict')) %dopar% {
    kcauc.values = c()
    for (j in 1:kfold){
      kcauc.values[j] = 0
      test <- test.inxs[[j]]
      K <- as.kernelMatrix(crossprod(t(X[-test,])))
      tryCatch({
        ksvm.obj <- ksvm(K,ytr[-test],C=C,kernel='matrix')#,prob.model=T)#,type='nu-svc')
        Ktest <- as.kernelMatrix(crossprod(t(X[test,]),t(X[SVindex(ksvm.obj), ])))  
        predictions <- predict(ksvm.obj,Ktest,type='decision')
        #print(predictions)
        labels = factor(ytr[test])
        kcauc.values[j] <- auc(roc(predictions,labels))
      },
      error = function(e) {
      })
    }
    return(kcauc.values)
  }
  
  b <- which.max(rowMeans(kcauc))
  
  C <- CRange[b[1]]
  print('finished')
  
  #kcauc <- matrix(0,nrow=length(LambdaRange),ncol=kfold) #optimize lambda
 
  test.inxs <- generate.test.inxs(nrow(X),kfold)
  print('optimizing lambda...')
  if (length(LambdaRange) > 1) {
  kcauc = foreach(i=1:length(LambdaRange),.packages=c('kernlab','AUC','CCPredict'),.combine=rbind) %dopar% {
    #for (i in 1:length(LambdaRange)){
    lam <- LambdaRange[i]
    #print(lam)
    kcauc.values = c()
    for (j in 1:kfold){
      kcauc.values[j] = 0
      test <- test.inxs[[j]]
      rescaled <- rescaling(X,L,lam)
      X.new <- rescaled[[1]]
      K.new <- rescaled[[2]]
      l <- rescaled[[3]]
      tryCatch({
        ksvm.obj <- ksvm(K.new[-test,-test],y[-test],C=C,kernel='matrix')#,prob.model=T)#,type='nu-svc')
        Ktest.new <- as.kernelMatrix(crossprod(t(X.new[test,]),t(X.new[SVindex(ksvm.obj), ])))  
        predictions <- predict(ksvm.obj,Ktest.new,type='decision')
        labels <- y[test]
        kcauc.values[j] = auc(roc(predictions,labels))
      },
      error = function(e) {
      })
    }
    return(kcauc.values)
  }
  
    b <- which.max(rowMeans(kcauc))
    lambda <- LambdaRange[b[1]]
  } else {
    lambda = LambdaRange[1]
  }
  print(lambda)
  print('finished')
  stopCluster(cl)
  
  return(c(lambda,C))
  

}


#' Kernel matrix rescaling
#'
#' Uses an optimized lambda to rescale K (the kernel matrix) by way of X
#' (the data set)
#'
#' @param X - input data matrix
#' @param L - side information matrix
#' @param lambda - optimized lambda value (see Li et al.)
#'
#' @return the rescaled matrices X and K
#'
#' @export
rescaling <- function(X,L,lambda){
  
  #instead of lambda, nox
  n <- dim(X)[1]
  m <- dim(X)[2]
  
  H <- diag(n,n)-1/n*matrix(1,n,n)
  L <- H%*%L%*%H/((m-1)^2)
  
  l <- c()
  if (lambda > 0){
    for (i in 1:m){
      xi <- X[,i]
      l[i] <- sqrt(lambda*t(xi)%*%L%*%xi+1)
      X[,i] <- xi/l[i]
    }
    l = t(l)
  } else{
    l = matrix(1,m,1)
  }
  
  X.new = X
  K.new = X.new%*%t(X.new)
  
  return (list(X.new,K.new,l))
}

#' Generate Test Indices
#'
#' Generates the test indices for each fold for kfold CV
#'
#' @param n - number of samples
#' @param kfold - number of folds
#' @return test indices
#'
#' @examples
#' X <- read.csv(system.file("extdata", "X.csv", package="CCPredict"),header=FALSE)
#'
#' @export
#'
generate.test.inxs <- function(n,kfold) {
  t.inxs <- shuffle(n)
  size <- round(n/kfold)
  test.inxs <- list()
  for(i in 1:kfold){
    start <- 1 + size*(i-1)
    end <- min(nrow(X),size + size*(i-1))
    test.inxs[[i]] <- t.inxs[start:end]
  }
  return(test.inxs)
}

#' Predict with ccSVM
#'
#' @param X - Input data matrix (NxD)
#' @param y - Labels (Nx1)
#' @param L - Side information kernel matrix
#' @param test.inxs - test incides for CV
#' @param lambda - optimized lambda value (cc)
#' @param C - optimized C value (SVM)
#' @return ROC curve, labels, predicted labels, and AUC
#'
#' @examples
#' X <- read.csv(system.file("extdata", "X.csv", package="CCPredict"),header=FALSE)
#'
#' @export
#'
predict.ccsvm <- function(X,y,L,test.inxs,lambda,C) {
  res = rescaling(X,L,lambda) #isn't this the same as just using the rescaling function? or is it hiding that function?
  X.new = res[[1]]
  K.new = res[[2]]

  ksvm.obj <- ksvm(K.new[-test.inxs,-test.inxs],y[-test.inxs],C=C,kernel='matrix')
      
  #Ktest.new1 = K.new[test.inxs,test.inxs]
  Ktest.new2 <- as.kernelMatrix(crossprod(t(X.new[test.inxs,]),t(X.new[SVindex(ksvm.obj), ])))  
  # TODO: Compare the above - looks like Ktest.new1 dim = 75,75 Ktest.new2 dim = 75, 42...using Ktest.new2
  predictions <- predict(ksvm.obj,Ktest.new2,type='decision')
  roc.curve <- roc(predictions,y[test.inxs])
  m <- predictions
  
  kcauc <- auc(roc.curve)
  labels <- y[test.inxs]
  r <- roc.curve
  return(list(roc.curve=roc.curve,labels=labels,predicted.labels=m,auc=kcauc))
}
#' Predict with cckOPLS
#'
#' @param X - Input data matrix (NxD)
#' @param y - Labels (Nx1)
#' @param L - Side information kernel matrix
#' @param test.inxs - test incides for CV
#' @param lambda - optimized lambda value (cc)
#' @param nox - optimized nox value (kOPLS)
#' @return ROC curve, labels, predicted labels, and AUC
#'
#' @examples
#' X <- read.csv(system.file("extdata", "X.csv", package="CCPredict"),header=FALSE)
#'
#' @export
predict.cckopls <- function(X,y,L,test.inxs,lambda,nox) {
  # Make ytr
  values = sort(unique(y))
  ytr <- matrix(0,nrow=length(y),length(values))
  for (i in 1:length(values)) {
    ytr[y==values[i],i] <- 1
  }
  #print(ytr)

  res = rescaling(X,L,lambda)
  X.new = res[[1]]
  K.new = res[[2]]

  modelOrg <- koplsModel(K.new[-test.inxs,-test.inxs],ytr[-test.inxs,],length(values)-1,nox,'mc','mc')
  modelOrgPred<-koplsPredict(K.new[test.inxs,-test.inxs],K.new[test.inxs,test.inxs],K.new[-test.inxs,-test.inxs],modelOrg,rescaleY=TRUE)
  roc.curve <- roc(modelOrgPred$Yhat[,2],y[test.inxs])

  m <- modelOrgPred$Yhat[,2]
  kcauc <- auc(roc.curve)
  labels <- y[test.inxs]
  r <- roc.curve
  return(list(roc.curve=roc.curve,labels=labels,predicted.labels=m,auc=kcauc))
}

predict.helper <- function(X,L,lambda) {
  rescaled <- rescaling(X,L,lambda)
  X.new <- rescaled[[1]]
  K.new <- rescaled[[2]]
  return(list(X.new,K.new))
}


