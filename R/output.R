ssDiagnostic <- function(ssObj){
  y <- ssObj$yt[1,]   #for univariate models
  t <- ncol(ssObj$yt)
  d <- nrow(ssObj$at)
  q <- nrow(ssObj$Ht)+nrow(ssObj$Qt)
  
  stable <- min(which(abs(diff(ssObj$Ft[1, ]))<0.001*abs(mean(y))))
  #Stable when changes in Ft is less than 0.1% of mean observation value
  
  st.res <- as.vector(ssObj$vt/ssObj$Ft)
  n <- length(st.res)	
  
  AR1m <- lm(y[-1]~y[-t])
  AR1pred <- matrix(AR1m$coefficients, ncol=2)%*%rbind(1, y[-t])
  
  MEANssq <- sum((y[-(1:stable)]-mean(y[-(1:stable)]))^2)
  AR1ssq <- sum((AR1pred[1, -(1:(stable-1))]-y[-(1:stable)])^2)
  SSssq <- sum(ssObj$vt[1, -(1:stable)]^2)
  
  
  lOut <- list(Rsquared = 1-SSssq/MEANssq, AR1measure = 1-SSssq/AR1ssq, MSE = SSssq/(t-stable))
  
  
  c(list(AIC = 1/n*(-2*n*ssObj$LogL+2*(q+d))),
    #Not sure if this is right, diffuse elements are not estimated
    diagStat(st.res[-(1:d)]),
    lOut
  )	  
}

ssPlot <- function(ssObj){ 
  d <- nrow(ssObj$at)
  t <- ncol(ssObj$yt)
  st.res <- as.vector(ssObj$vt/ssObj$Ft)[-(1:d)]
  n <- length(st.res)
  Zt <- matrix(ssObj$Zt_, nrow=t, ncol=d, byrow=T)
  alfa <- Zt*t(matrix(ssObj$at[, 1:t], ncol=t))

  
  par(mfrow=c(2,2))
  plot(as.vector(ssObj$yt), type="l")
  lines(as.vector( apply(alfa,1,sum) ), type="l", col=2)
  
  
  plot(st.res, type="l")
  
  hist(st.res, breaks=floor(n*0.1))
  
  barplot(AutoCor( st.res, 1:min(round(length(st.res)*0.1, 0), 50) ), main=
            paste0("Autocorrelation from 1 to ", 
                   min(round(length(st.res)*0.1, 0), 50) )
  )
}