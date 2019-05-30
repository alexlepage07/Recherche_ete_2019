lam <- 2
a <- 2
b <- 1/50

alpha0 <- 0.1
alpha1 <- 0.2


TLS_theta <- function(t,al){
  q <- 1-al
  (exp(-t)* q / (1-(1-q)*exp(-t)))
}

TLS_theta_inv <- function(s,al){
  q <- 1-al
  log(q/s + (1-q) )
}


FFT_cop <- function(F_X){
  #algo 15
  #1 fix theta0
  theta0 <- 1
  #2 fix theta1
  theta0_1 <- 1
  #3 
  F_Xi001 <- exp(-theta0_1 * TLS_theta_inv(F_X,alpha1))
  #4
  f_Xi001 <- diff(c(0,F_Xi001))
  
  
}


fx <- c(0.2,0.4,0.4)

aa <- 2^12
nx <- length(fx)
ftx <- fft(c(fx, rep(0, aa - nx)))
fts<-exp(2 * (ftx - 1))
fs <- Re(fft(fts, T))/aa
sum(fs)
fs


