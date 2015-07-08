#optimizes nox and lambda for cckopls, and lambda and C for ccSVM, other parameters to be added as methods are added
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
#' a='blah'
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
#' @param noxRange blah blah
#' @param LambdaRange blah blah
#' @param kfold blah blah
#' @return blah blah
#'
#' @examples
#' a='blah'
#'
#' @export
#'
optimize.ccSVM <- function(X,ytr,L,CRange,LambdaRange,kfold=2){ #optimize ccSVM params
  library(kernlab)
  library(AUC)
  
  kcauc <- matrix(0,nrow=length(CRange),ncol=kfold) #optimize C
  
  size <- round(nrow(X)/kfold)
  test.inxs <- list()
  for(i in 1:kfold){
    start <- 1 + size*(i-1)
    end <- min(nrow(X),size + size*(i-1))
    test.inxs[[i]] <- start:end
  }
  library(AUC)
  print('optimizing C...')
  for (i in 1:length(CRange)){
    c <- CRange[i]
    #print(c)
    for (j in 1:kfold){
      test <- test.inxs[[j]]
      K <- as.kernelMatrix(crossprod(t(X[-test,])))
      ok = F
      while(ok == F) {
        tryCatch({
          ksvm.obj <- ksvm(K,ytr[-test],C=c,kernel='matrix',prob.model=T,type='nu-svc')
          Ktest <- as.kernelMatrix(crossprod(t(X[test,]),t(X[SVindex(ksvm.obj), ])))  
          predictions <- predict(ksvm.obj,Ktest,type='probabilities')[,2]
          labels = ytr[test]
          kcauc[i,j] <- auc(roc(predictions,labels))
          ok = T
        },
        error = function(e) {
          print('retrying ksvm')
          print('param')
          print(e)
          ok = F
        })
      }
    }
  }
  
  a <- max(rowMeans(kcauc))
  b <- which(rowMeans(kcauc) == a)
  
  C <- CRange[b[1]]
  print('finished')
  
  kcauc <- matrix(0,nrow=length(LambdaRange),ncol=kfold) #optimize lambda
  
  size <- round(nrow(X)/kfold)
  test.inxs <- list()
  for(i in 1:kfold){
    start <- 1 + size*(i-1)
    end <- min(nrow(X),size + size*(i-1))
    test.inxs[[i]] <- start:end
  }
  
  print('optimizing lambda...')
  for (i in 1:length(LambdaRange)){
    lam <- LambdaRange[i]
    #print(lam)
    for (j in 1:kfold){
      test <- test.inxs[[j]]
      rescaled <- Rescaling(X,L,lam)
      X.new <- rescaled[[1]]
      K.new <- rescaled[[2]]
      #l <- rescaled[[3]]
      ok <- F
      while(ok == F) {
        tryCatch({
          ksvm.obj <- ksvm(K.new[-test,-test],ytr[-test],C=c,kernel='matrix',prob.model=T)#,type='nu-svc')
          Ktest.new <- as.kernelMatrix(crossprod(t(X.new[test,]),t(X.new[SVindex(ksvm.obj), ])))  
          predictions <- predict(ksvm.obj,Ktest.new,type='probabilities')[,2]
          labels <- ytr[test]
          kcauc[i,j] <- auc(roc(predictions,labels))
          ok = T
        },
        error = function(e) {
          print('retrying ksvm')
          print('lambda')
          print(e)
          ok = F
        })
      }
    }   
  }
  
  a <- max(rowMeans(kcauc))
  b <- which(rowMeans(kcauc) == a)
  
  lambda <- LambdaRange[b[1]]
  print('finished')
  
  return(c(lambda,C))
  
#   kcauc <- vector(length=length(LambdaRange)*length(CRange))
#   size <- round(nrow(X)/kfold)
#   test.inxs <- list()
#   for(i in 1:kfold){
#     start <- 1 + size*(i-1)
#     end <- min(nrow(X),size + size*(i-1))
#     test.inxs[[i]] <- start:end
#   }
#   
#   c = 1
#   iz = matrix(0,nrow=length(LambdaRange)*length(CRange),ncol=2)
#   for (z in 1:length(CRange)) {
#     for (i in 1:length(LambdaRange)){
#       lambda <- LambdaRange[i]
#       kcauc_total = 0
#       for (j in 1:kfold){
#         rescaled <- Rescaling(X,L,lambda)
#         X.new <- rescaled[[1]]
#         K.new <- rescaled[[2]]
#         ok = F
#         while(ok == F) {
#           tryCatch({
#             ksvm.obj <- ksvm(K.new[-test.inxs[[j]],-test.inxs[[j]]],y[-test.inxs[[j]]],C=c,kernel='matrix',prob.model=T,type='nu-svc')
#             Ktest.new <- as.kernelMatrix(crossprod(t(X.new[test.inxs[[j]],]),t(X.new[SVindex(ksvm.obj), ])))
#             predictions <- predict(ksvm.obj,Ktest.new,type='probabilities')[,2]
#             labels = y[test.inxs[[j]]]
#             kcauc[i,j] <- auc(roc(predictions,labels))
#             ok = T
#           },
#           error = function(e) {
#             print('retrying ksvm')
#             print('param')
#             print(e)
#             ok = F
#           })
#         }
#       }
#     }
#       kcauc_total = kcauc_total/kfold
#       kcauc[c] = kcauc_total
#       iz[c,1] = i
#       iz[c,2] = z
#       c = c + 1
#     }
#   
#   ix = which.max(kcauc)    
#   i = iz[ix,1]
#   z = iz[ix,2]
#   lambda <- LambdaRange[i]
#   C <- CRange[z]
#   
  
}
