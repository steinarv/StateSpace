SK_gm <- function(x){
  (mean(x)-median(x))/mean(abs(x-median(x)))
}

KT_cs <- function(x, a=0.025, b=0.25){
  (quantile(x,probs=1-a,type=1)-quantile(x,probs=a,type=1))/
    (quantile(x,probs=1-b,type=1)-quantile(x,probs=b,type=1))-2.91
}

AutoCor <- function(x, lags){
  n <- length(x)
  cOut <- c()
  for(i in lags)cOut <- c( cOut, cor(x[(i+1):n],x[1:(n-i)]) )

  names(cOut) <- paste0(lags)
  cOut
}

Auto.port <- function (y, lags = 10) #Modified from package vrtest
{
  data <- y - mean(y)
  T <- length(data)
  ac <- matrix(acf(data, lag.max = lags, plot = F)$acf[, , 
                                                       1])
  ac <- ac[2:(lags + 1), 1]
  ac1 <- matrix(acf(data, lag.max = lags, type = "covariance", 
                    plot = F)$acf[, , 1])
  ac2 <- matrix(NA, nrow = lags)
  for (i in 1:lags) {
    y <- data[(i + 1):T]^2
    x <- data[1:(T - i)]^2
    t <- length(y)
    ac2[i] <- crossprod(x, y)/t
  }
  ac3 <- (ac1[2:(lags + 1), 1]^2)/ac2
  BP <- T * cumsum(ac3)
  aux <- matrix(1:lags)
  q <- 2.4
  maxro <- sqrt(T) * max(sqrt(ac3))
  pin <- 2
  if (maxro <= sqrt(q * log(T))) 
    pin <- log(T)
  Lp <- BP - aux * pin
  phat <- which.max(Lp)
  Tn <- BP[phat]
  pvalue <- 1 - pchisq(Tn, 1)
  return(list(Stat = Tn, Pvalue = pvalue, phat=phat))
}

GPH<-function(x,K){
  N=length(x)
  w=2*pi*1:N/N
  y <-log(1/length(x)*abs(fft(x))^2)
  x <- log(4*sin(w/2)^2)*-1
  fit <- lm(y[1:K]~x[1:K])
  coefficients(fit)[2]+0.5
}

diagStat <- function(x){
  auto.port.out=Auto.port(x, lags=floor(length(x)/2))
  cOut <-  c(mean(x),
             sd(x),
             max(x),
             min(x),
             quantile(x, 0.25),
             median(x, 0.5),
             quantile(x, 0.75),
             skewness(x),
             SK_gm(x),
             kurtosis(x),
             KT_cs(x),
             ad.test(x)$p.value,
             coefficients(lm(x[2:length(x)]~x[1:(length(x)-1)]))[2],
             auto.port.out$Pvalue,
             auto.port.out$phat,
             GPH(x, K=floor(length(x)^0.5))
			 )
  
  names(cOut) <- c("Mean", "SD", "Max", "Min", "Q1", "Q2", "Q3",
                   "Skewness", "skew_GM_1984", "Kurtosis", 
                   "kurt_CS_1967", "Anderson_Darling_pvalue", "AR_1","Auto.port_pvalue",
                   "Auto.port_lags", "d_GPH_1983")
  
  as.list(round(cOut,4))
}

