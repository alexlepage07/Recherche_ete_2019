library(actuar)
set.seed(2^31)
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




sim_algo14 <- function(mm,alpha_0,alpha_1){
  M <- rep(0,mm)
  R <- rep(0,mm)
  U_0 <- rep(0,mm)
  S <- rep(0,mm)
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
      S[i] <- 0
    }else{
      theta_0_1[i] <- sum(qgamma(runif(M[i]),1/alpha_1,1))
      R_0_1 <- qexp(runif(N[i]),1)
      U_0_1 <- LTS_M(-log(LTS_B(R_0_1/theta_0_1[i],alpha_1)),alpha_0)
      S[i] <- sum(  qpareto(U_0_1,3,100) )
    }
  }
  return(S)
}
S_1_1 <- sim_algo14(m,a0_1,a1_1)

S_1_2 <- sim_algo14(m,a0_1,a1_2)
S_2_1 <- sim_algo14(m,a0_2,a1_1)
S_2_2 <- sim_algo14(m,a0_2,a1_2)
DATA <- list(S_1_1,S_1_2,S_2_1,S_2_2)

Output <- function(data){
  c(mean(data),var(data),sort(data)[length(data)*0.9],sort(data)[length(data)*0.99],
    sort(data)[length(data)*0.999],
    mean(data[data>sort(data)[length(data)*0.9]]),
    mean(data[data>sort(data)[length(data)*0.99]]),
    mean(data[data>sort(data)[length(data)*0.999]]))
}
result <- lapply(DATA,Output)


cbind(result[[1]],result[[2]],result[[3]],result[[4]])
