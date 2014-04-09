ssOpt <- function(x, Zt, Tt, Rt){
  Qparam <- sum(Rt)
  
  Ht <- matrix(1)
  Qt <- diag(1, Qparam)
  x=matrix(x, nrow=1);
  
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