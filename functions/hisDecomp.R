require(quantmod)
require(MASS)
VARhd <- function(Estimation, #output from uhlig.reject
                  DF          #matrix of endogenous variables
                  ){
  
  ## make X and Y
  nlag <- (length(Estimation$BDraws[1,,1]))/length(Estimation$BDraws[1,1,]) 
  DATA <- DF 
  QQ   <- VARmakexy(DATA=DATA,lags=nlag,c_case=0)

  ## Retrieve and initialize variables 
  invA    <- apply(Estimation$SDraws,MARGIN=c(2,3),mean) #posterior mean of covariance matrix
  invA    <- t(chol(as.matrix(invA)))                       #inverse of the A matrix
  Fcomp   <- companionmatrix(Estimation)                    #Companion matrix
  
  #up until here
  F1      <- t(QQ$Ft)                                         # make comparable to notes
  eps     <- ginv(invA) %*% t(getResidUhlig(DF=DF,Estimation=Estimation))   # structural errors
  eps     <- eps[,-c(1:nlag)]                                 # trim lag NAs
  nvar    <- dim(Estimation$BDraws)[3]                        # number of endogenous variables
  nvarXeq <- nvar * nlag                                      # number of lagged endogenous per equation
  nvar_ex <- 0                                                # number of exogenous (excluding constant and trend)
  Y       <- QQ$Y                                             # left-hand side
  #X       <- QQ$X[,(1+det):(nvarXeq+det)]                    # right-hand side (no exogenous)
  nobs    <- nrow(Y)                                          # number of observations
  
  ## Compute historical decompositions
  # Contribution of each shock
  invA_big <- matrix(0,nvarXeq,nvar)
  invA_big[1:nvar,] <- invA
  Icomp <- cbind(diag(nvar), matrix(0,nvar,(nlag-1)*nvar))
  HDshock_big <- array(0, dim=c(nlag*nvar,nobs+1,nvar))
  HDshock <- array(0, dim=c(nvar,(nobs+1),nvar))
  
  for (j in 1:nvar){  # for each variable
    eps_big <- matrix(0,nvar,(nobs+1)) # matrix of shocks conformable with companion
    eps_big[j,2:ncol(eps_big)] <- eps[j,]
    for (i in 2:(nobs+1)){
      HDshock_big[,i,j] <- invA_big %*% eps_big[,i] + Fcomp %*% HDshock_big[,(i-1),j]
      HDshock[,i,j] <-  Icomp %*% HDshock_big[,i,j]
    } 
    
  } 
  
  HD.shock <- array(0, dim=c((nobs+nlag),nvar,nvar))   # [nobs x shock x var]
  
  for (i in 1:nvar){
    
    for (j in 1:nvar){
      HD.shock[,j,i] <- c(rep(NA,nlag), HDshock[i,(2:dim(HDshock)[2]),j])
    }
  }
  
  return(HD.shock)
  
}


VARmakexy <- function(DATA,lags,c_case){
  
  nobs <- nrow(DATA)
  
  #Y matrix 
  Y <- DATA[(lags+1):nrow(DATA),]
  Y <- DATA[-c(1:lags),]
  
  #X-matrix 
  if (c_case==0){
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    } 
  } else if(c_case==1){ #constant
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    }
    X <- cbind(matrix(1,(nobs-lags),1), X) 
  } else if(c_case==2){ # time trend and constant
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    }
    trend <- c(1:nrow(X))
    X <-cbind(matrix(1,(nobs-lags),1), t(trend))
  }
  A <- (t(X) %*% as.matrix(X)) 
  B <- (as.matrix(t(X)) %*% as.matrix(Y))
  
  Ft <- ginv(A) %*% B
  
  retu <- list(X=X,Y=Y, Ft=Ft)
  return(retu)
}

companionmatrix <- function (x) 
{
  
  K <- length(x$BDraws[1,1,])  #dimension of VAR
  p <- (length(x$BDraws[1,,1]))/length(x$BDraws[1,1,]) #lag order of VAR
  A <- t(apply(x$BDraws,MARGIN=c(2,3),mean)) #posterior mean of the lagged endogenous variables
  # A <- A[,-nrow(A)]
  companion <- matrix(0, nrow = K * p, ncol = K * p)
  companion[1:K, 1:(K * p)] <- A
  if (p > 1) {
    j <- 0
    for (i in (K + 1):(K * p)) {
      j <- j + 1
      companion[i, j] <- 1
    }
  }
  return(companion)
}

getResidUhlig <- function(DF,Estimation){
  
  #compute model fit
  xMat <- DF
  nlag <- dim(Estimation$BDraws)[3]
  temp <- list()
  for (i in 1:nlag){
    temp[[i]] <- apply(xMat,MARGIN=2,FUN=Lag,k=i)
  }
  xMat <- do.call(cbind,temp)
  xMat <- cbind(xMat)
  #coefficients
  beta <- apply(Estimation$BDraws,MARGIN=c(2,3),mean)
  #get residuals for each endogenous variables
  nEndog <- dim(Estimation$BDraws)[3]
  resid <- list()
  for (i in 1:nEndog){
  resid[[i]] <- xMat%*%beta[,i] - DF[,i]
  }
  
  resid <- do.call(cbind,resid)
  return(resid)
}
