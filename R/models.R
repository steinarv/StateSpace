LocalLevel <- function(x){
  n <- length(x)
  
  Zt=matrix(1); Tt=matrix(1); Rt=matrix(1);
  x=matrix(x, nrow=1); Ht=matrix(0); Qt=matrix(0);
  
  ssObj <- new(StateSpace, x, Zt, Tt, Rt, Ht, Qt)
  
  optFunc <- function(param){
    param <- exp(2*param)
    -ssObj$optimFunc(matrix(param[1]), matrix(param[2]))
  }
  
  optRes <- optim(c(1,1), optFunc)
  
  ssObj$Ht <- matrix(exp(2*optRes$par[1]))
  ssObj$Qt <- matrix(exp(2*optRes$par[2]))
  ssObj$kalmanFil()
  ssObj$kalmanSmooth()
  
  ssObj
}  

seasonal <- function(x, n.seas, det.seas=F){
  m <- n.seas
  Qparam <- 1+!det.seas
  
  Zt <- matrix(0, ncol=m, nrow=1)
  Tt <- matrix(0, ncol=m, nrow=m)
  Rt <- matrix(0, ncol=Qparam, nrow=m)
  
  Zt[ , 1:2] <- 1
  Tt[1,1] <- 1
  for(i in 2:m)Tt[2,i] <- -1
  for(i in 3:m)Tt[i, i-1] <- 1
  Rt[1:Qparam, ] <- diag(1,Qparam)
  
  x=matrix(x, nrow=1); Ht=matrix(0); Qt=diag(0, Qparam);
  
  ssObj <- new(StateSpace, x, Zt, Tt, Rt, Ht, Qt)
  
  optFunc <- function(param){
    param <- exp(2*param)
    -ssObj$optimFunc(matrix(param[1]), diag(param[2:(Qparam+1)], Qparam))
  }
  
  initGuess <- rep(-1, 2+!det.seas)
  optRes <- optim(initGuess, optFunc)
  
  ssObj$Ht <- matrix( exp(2*optRes$par[1]) )
  ssObj$Qt <- diag(exp(2*optRes$par[-1]), Qparam) 

  ssObj$kalmanFil()
  ssObj$kalmanSmooth()
  
  ssObj
  
}

seas_exvar <- function(x, ex.mat, n.seas, det.lev=F, det.seas=F, det.ex=rep(T, ncol(ex.mat))){
  #ex.mat must be matrix with exp.var in each column
  
  x <- matrix(x, nrow=1)
  n <-nrow(ex.mat)  
  
  m <- n.seas+ncol(ex.mat)
    
  Zt <- matrix(0, ncol=m, nrow=nrow(ex.mat)) #(for simplicity of code) will be transformed into a row vector  
  for(i in 1:n) Zt[i, ] <- c(1, 1, rep(0, n.seas-2), ex.mat[i, ])
  Zt <- matrix(t(Zt), nrow=1)
  
  Tt <- matrix(0, ncol=m, nrow=m)
  Tt[1,1] <- 1; 
  for(i in 2:n.seas)Tt[2,i] <- -1
  for(i in 3:n.seas)Tt[i, i-1] <- 1
  for(i in (n.seas+1):m)Tt[i,i] <- 1
  
  Qparam <- sum(!det.lev, !det.seas, !det.ex)
  
  Rt <- matrix(0, nrow=m, ncol=Qparam)
  if(!det.lev)Rt[1,1] <- 1
  if(!det.seas)Rt[2,2] <- 1
  if(any(!det.ex)){
  jc <- ncol(Rt); jr <- nrow(Rt)
  for(j in rev(det.ex)){if(!j){Rt[jr,jc]<-1; jc <- jc-1}; jr <- jr-1} 
  }
  
  
  Ht=matrix(0); Qt=matrix(0, ncol=Qparam, nrow=Qparam);
  ssObj <- new(StateSpace, x, Zt, Tt, Rt, Ht, Qt)
    
  
  optFunc <- function(param){
    param <- exp(2*param)
    Ht <- matrix(param[1])
    Qt <- diag(param[-1], Qparam)


    
    -ssObj$optimFunc(Ht, Qt)
  }
  
  optRes <- optim(rep(-1, Qparam+1), optFunc)
  
  ssObj$Ht <- matrix( exp(2*optRes$par[1]) )
  ssObj$Qt <- diag(exp(2*optRes$par[-1]), Qparam)
  ssObj$kalmanFil()
  ssObj$kalmanSmooth()
  
  ssObj
  
}

slope_seas_exvar <- function(x, ex.mat, n.seas, det.lev=F, det.slope=F, det.seas=F, det.ex=rep(T, ncol(ex.mat))){
  #ex.mat must be matrix with exp.var in each column
  
  x <- matrix(x, nrow=1)
  n <-nrow(ex.mat)  
  
  m <- n.seas+ncol(ex.mat)+1
  
  Zt <- matrix(0, ncol=m, nrow=nrow(ex.mat)) #(for simplicity of code) will be transformed into a row vector  
  for(i in 1:n) Zt[i, ] <- c(1, 0, 1, rep(0, n.seas-2), ex.mat[i, ])
  Zt <- matrix(t(Zt), nrow=1)
  
  Tt <- matrix(0, ncol=m, nrow=m)
  Tt[1,1:2] <- 1;
  Tt[2,2] <- 1;
  for(i in 3:(n.seas+1))Tt[3,i] <- -1
  for(i in 4:(n.seas+1))Tt[i, i-1] <- 1
  for(i in (n.seas+2):m)Tt[i,i] <- 1
  
  Qparam <- sum(!det.lev, !det.slope, !det.seas, !det.ex)
  
  Rt <- matrix(0, nrow=m, ncol=Qparam)
  if(!det.lev)Rt[1,1] <- 1
  if(!det.slope)Rt[2,2] <- 1
  if(!det.seas)Rt[3,3] <- 1
  if(any(!det.ex)){
    jc <- ncol(Rt); jr <- nrow(Rt)
    for(j in rev(det.ex)){if(!j){Rt[jr,jc]<-1; jc <- jc-1}; jr <- jr-1} 
  }
  
  
  Ht=matrix(0); Qt=matrix(0, ncol=Qparam, nrow=Qparam);
  ssObj <- new(StateSpace, x, Zt, Tt, Rt, Ht, Qt)
  
  
  optFunc <- function(param){
    param <- exp(2*param)
    Ht <- matrix(param[1])
    Qt <- diag(param[-1], Qparam)
    
    
    
    -ssObj$optimFunc(Ht, Qt)
  }
  
  optRes <- optim(rep(-1, Qparam+1), optFunc)
  
  ssObj$Ht <- matrix( exp(2*optRes$par[1]) )
  ssObj$Qt <- diag(exp(2*optRes$par[-1]), Qparam)
  ssObj$kalmanFil()
  ssObj$kalmanSmooth()
  
  ssObj
  
}