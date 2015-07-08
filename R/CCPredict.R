library(kopls)
library(kernlab)
library(AUC)
#library(modeest)
#library(permute)

#' blah blah blah
#'
#' @param X blah blah
#' @param ytr blah blah
#' @param L blah blah
#' @param noxRange blah blah
#' @param LambdaRange blah blah
#' @param kfold blah blah
#' @return blah blah
#'
#' @examples
#' X <- read.csv(system.file("extdata", "X.csv", package="CCPredict"),header=FALSE)
#' X <- t(X)
#' X = scale(X,center=T,scale=T) # Scale the X data so it has a mean of 0 and a stdev of 1. Pretty standard
#' y <- read.csv(system.file("extdata", "y.csv", package="CCPredict"),header=FALSE)
#' L <- read.csv(system.file("extdata", "L.csv", package="CCPredict"),header=FALSE)
#' y <- as.matrix(y)
#' y <- factor(y[,1])
#' L <- as.matrix(L)
#' optimize.cckopls(X,y,L,c(0:3),c(1e-8,1e-4,1e-2,1,1e+2,1e+4,1e+8),5)
#'
#' @export
#'
optimize.cckopls <- function(X,ytr,L,noxRange,LambdaRange,kfold){ #optimize cckopls/kopls params
  
  kcauc <- matrix(0, nrow=length(noxRange),ncol=kfold)
  size <- round(nrow(X)/kfold)
  test.inxs <- list()
  for(i in 1:kfold){
    start <- 1 + size*(i-1)
    end <- min(nrow(X),size + size*(i-1))
    test.inxs[[i]] <- start:end
  }
  print('optimizing nox...')
  for (i in 1:length(noxRange)){
    n <- noxRange[i]
    for (j in 1:kfold){
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
      kcauc[i,j] <- auc(roc(modelOrgPred$Yhat[,2],factor(ytr[test,2])))
    }
  }
  
  a <- max(rowMeans(kcauc))
  b <- which(rowMeans(kcauc) == a)
  
  nox <- noxRange[b[1]]
  print('finished')
  
  kcauc <- matrix(0, nrow=length(LambdaRange),ncol=kfold)
  size <- round(nrow(X)/kfold)
  test.inxs <- list()
  for(i in 1:kfold){
    start <- 1 + size*(i-1)
    end <- min(nrow(X),size + size*(i-1))
    test.inxs[[i]] <- start:end
  }
  
  print('optimizing lambda...')
  for (i in 1:length(LambdaRange)){
    lambda <- LambdaRange[i]
    for (j in 1:kfold){
      rescaled <- Rescaling(X,L,lambda)
      X.new <- rescaled[[1]]
      K.new <- rescaled[[2]]
      # n.list <- rescaled[[3]]
      test <- test.inxs[[j]]
      #modelCV <- koplsCV(K.new,ytr,1,10,nrcv=7,cvType='nfold',preProcK='mc',preProcY='mc',modelType='da')
      modelOrg <- koplsModel(K.new[-test,-test],ytr[-test,],1,nox,'mc','mc')
      modelOrgPred<-koplsPredict(K.new[test,-test],K.new[test,test],K.new[-test,-test],modelOrg,rescaleY=TRUE)
      kcauc[i,j] <- auc(roc(modelOrgPred$Yhat[,2],y[test]))
    }
  }
  
  a <- max(rowMeans(kcauc))
  b <- which(rowMeans(kcauc) == a)
  
  lambda <- LambdaRange[b[1]]
  
  print('finished')
  
  return(c(lambda,nox))
  
  #   c = 1
  #   iz = matrix(0,nrow=length(LambdaRange)*length(noxRange),ncol=2)
  #   for (z in 1:length(noxRange)) {
  #     for (i in 1:length(LambdaRange)){
  #       lambda <- LambdaRange[i]
  #       kcauc_total = 0
  #       for (j in 1:kfold){
  #         rescaled <- Rescaling(X,L,lambda)
  #         X.new <- rescaled[[1]]
  #         K.new <- rescaled[[2]]
  #         test <- test.inxs[[j]]
  #         modelOrg <- koplsModel(K.new[-test,-test],ytr[-test,],1,noxRange[z],'mc','mc')
  #         modelOrgPred<-koplsPredict(K.new[test,-test],K.new[test,test],K.new[-test,-test],modelOrg,rescaleY=TRUE)
  #         kcauc_total <- kcauc_total + auc(roc(modelOrgPred$Yhat[,2],factor(ytr[test,2])))
  #         #print(auc(roc(modelOrgPred$Yhat[,2],factor(ytr[test,2]))))
  #       }
  #       kcauc_total = kcauc_total/kfold
  #       kcauc[c] = kcauc_total
  #       iz[c,1] = i
  #       iz[c,2] = z
  #       c = c + 1
  #     }
  #   }
  #   
  #   ix = which.max(kcauc)    
  #   i = iz[ix,1]
  #   z = iz[ix,2]
  #   lambda <- LambdaRange[i]
  #   nox <- noxRange[z]
  
} #end of cckopls opt

#' blah blah blah
#'
#' @param X blah blah
#' @param ytr blah blah
#' @param L blah blah
#' @param CRange blah blah
#' @param LambdaRange blah blah
#' @param kfold blah blah
#' @return blah blah
#'
#' @examples
#' X <- read.csv(system.file("extdata", "X.csv", package="CCPredict"),header=FALSE)
#' X <- t(X)
#' X = scale(X,center=T,scale=T) # Scale the X data so it has a mean of 0 and a stdev of 1. Pretty standard
#' y <- read.csv(system.file("extdata", "y.csv", package="CCPredict"),header=FALSE)
#' L <- read.csv(system.file("extdata", "L.csv", package="CCPredict"),header=FALSE)
#' y <- as.matrix(y)
#' y <- factor(y[,1])
#' L <- as.matrix(L)
#' optimize.ccSVM(X,y,L,c(2^-8,2^-4,2^-2,2^0,2^2,2^4,2^8),c(1e-8,1e-4,1e-2,1,1e+2,1e+4,1e+8))
#'
#' @export
#'
optimize.ccSVM <- function(X,ytr,L,CRange,LambdaRange,kfold=2){ #optimize ccSVM params
  cl<-makeCluster(8)
  registerDoParallel(cl)
  
  #kcauc <- matrix(0,nrow=length(CRange),ncol=kfold) #optimize C
  
  size <- round(nrow(X)/kfold)
  test.inxs <- list()
  for(i in 1:kfold){
    start <- 1 + size*(i-1)
    end <- min(nrow(X),size + size*(i-1))
    test.inxs[[i]] <- start:end
  }
  print('optimizing C...')
  kcauc = foreach(i=1:length(CRange),.packages=c('kernlab','AUC','CCPredict'),.combine=rbind) %dopar% {
  #for (i in 1:length(CRange)) {
    C <- CRange[i]
    #print(c)
    #results = foreach(j=1:kfold,.packages=c('kernlab','AUC','CCPredict')) %dopar% {
    kcauc.values = c()
    for (j in 1:kfold){
      test <- test.inxs[[j]]
      K <- as.kernelMatrix(crossprod(t(X[-test,])))
      tryCatch({
        ksvm.obj <- ksvm(K,ytr[-test],C=C,kernel='matrix',prob.model=T,type='nu-svc')
        Ktest <- as.kernelMatrix(crossprod(t(X[test,]),t(X[SVindex(ksvm.obj), ])))  
        predictions <- predict(ksvm.obj,Ktest,type='probabilities')[,2]
        labels = ytr[test]
        kcauc.values[j] = auc(roc(predictions,labels))
      },
      error = function(e) {
        kcauc.values[j] = 0
      })
    }
    return(kcauc.values)
  }
  
  a <- max(rowMeans(kcauc))
  b <- which(rowMeans(kcauc) == a)
  
  C <- CRange[b[1]]
  print(C)
  print('finished')
  
  #kcauc <- matrix(0,nrow=length(LambdaRange),ncol=kfold) #optimize lambda
  
  size <- round(nrow(X)/kfold)
  test.inxs <- list()
  for(i in 1:kfold){
    start <- 1 + size*(i-1)
    end <- min(nrow(X),size + size*(i-1))
    test.inxs[[i]] <- start:end
  } 
  
  print('optimizing lambda...')
  #foreach(i=1:length(LambdaRange),.packages=c('kernlab','AUC','CCPredict')) %dopar% { 
  kcauc = foreach(i=1:length(LambdaRange),.packages=c('kernlab','AUC','CCPredict'),.combine=rbind) %dopar% {
    #for (i in 1:length(LambdaRange)){
    lam <- LambdaRange[i]
    #print(lam)
    kcauc.values = c()
    for (j in 1:kfold){
      test <- test.inxs[[j]]
      rescaled <- Rescaling(X,L,lam)
      X.new <- rescaled[[1]]
      K.new <- rescaled[[2]]
      l <- rescaled[[3]]
      tryCatch({
        ksvm.obj <- ksvm(K.new[-test,-test],ytr[-test],C=C,kernel='matrix',prob.model=T,type='nu-svc')
        Ktest.new <- as.kernelMatrix(crossprod(t(X.new[test,]),t(X.new[SVindex(ksvm.obj), ])))  
        predictions <- predict(ksvm.obj,Ktest.new,type='probabilities')[,2]
        labels <- ytr[test]
        kcauc.values[j] = auc(roc(predictions,labels))
      },
      error = function(e) {
        kcauc.values[j] = 0
      })
    }
    return(kcauc.values)
  }
  
  a <- max(rowMeans(kcauc))
  b <- which(rowMeans(kcauc) == a)
  
  lambda <- LambdaRange[b[1]]
  print('finished')
  
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
Rescaling <- function(X,L,lambda){
  
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
