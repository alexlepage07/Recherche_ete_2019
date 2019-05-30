set.seed(2^30)
LTS_M <- function(t,alpha0){
  (-1/alpha0) * log(1-(1-exp(-alpha0))*exp(-t))
}
LTS_B <- function(t,alpha1){
  (1/(1+t))^(1/alpha1)
}

LTS_M_inv <- function(s,alpha0){
  - log( (1-exp(-alpha0 * s))/(1-exp(-alpha0))  )
}
LTS_B_inv <- function(s,alpha1){
  s^(-alpha1)-1
}

a0_1 <- -log(0.5)
a0_2 <- -log(0.1)

a1_1 <- 5
a1_2 <- 25
m <- 10^7

LTS_M(0.5,1)
LTS_M_inv(0.4835356,1)

LTS_B(0.5,1)
LTS_B_inv(2/3,1)

##algo 13



library(actuar)
sim_algo14 <- function(mm,alpha_0,alpha_1){
  M <- rep(0,mm)
  R <- rep(0,mm)
  U_0 <- rep(0,mm)
  X <- list()
  N <- rep(0,mm)
  theta_0_1 <- rep(0,mm)
  for (i in 1:mm){
    #sim theta0 aka M
    M[i] <- qlogarithmic(runif(1),1-exp(-alpha_0))
    #sim R
    R[i] <- qexp(runif(1),1)
    #sim U0
    U_0[i] <- LTS_M(R[i]/M[i],alpha_0)
    #sim N
    N[i] <- qpois(U_0[i],2)
    
    if(N[i]==0){
      X[[i]] <- 0
    }else{
      theta_0_1[i] <- sum(qgamma(runif(M[i]),1/alpha_1,1))
      R_0_1 <- qexp(runif(N[i]),1)
      U_0_1 <- LTS_M(-log(LTS_B(R_0_1/theta_0_1[i],alpha_1)),alpha_0)
      X[[i]] <- qpareto(U_0_1,3,100) 
    }
  }
  return(X)
}
XX <- sim_algo14(10000,a0_1,a1_1)
length(XX)

LOSS <- unlist(XX)
FREQ_LOSS <- unlist(lapply(XX,function(i) length(i)))
###hypothèse indép
neg_log_pois <- function(lam){
  -sum(log(dpois(FREQ_LOSS,lam)))
}
optimise(neg_log_pois,c(0,5))

neg_log_par <- function(para){
  -sum(log(dpareto(LOSS,para[1],para[2])))
}
constrOptim(c(3, 100), neg_log_par, grad = NULL, 
            ui = diag(2), ci = c(0, 0),outer.eps = .Machine$double.eps)

