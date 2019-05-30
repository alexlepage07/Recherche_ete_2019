d <- 2
n1 <- 40
n2 <- 40
set.seed(2^31)


LTS_B <- function(t,q){
 (exp(-t)* q / (1-(1-q)*exp(-t)))
}
LTS_M <- function(t,q){
  log(1-q*exp(-t))/log(1-q)
}

system.time({
  m = 1000000
  M = rlogarithmic(m, 0.5)
  V = t(sapply(M, function(x){
    c(rnbinom(1, x, 0.8),rnbinom(1, x, 0.9)) + x
  }))
  U = t(sapply(1:m, function(i){
    t1 = rexp(n1) / V[i,1]
    t2 = rexp(n2) / V[i,2]
    c(
      LTS_M(
        -log(LTS_B(t1, 0.8)),
        0.5),
      LTS_M(
        -log(LTS_B(t2, 0.9)),
        0.5)
    )
  }))
})
qq <- c(0.05+0.005*1:40,0.1+0.005*1:40)
X <- t(sapply(1:10^6,function(i) qbinom(U[i,],10,qq) ))

S <- sapply(1:10^6,function(i) sum(X[i,]))
mean(S)
var(S)
S <- sort(S)

VaR_S <- function(k) S[ceiling(k*length(S))]
TVaR_S <- function(k) mean(S[S>VaR_S(k)])


kappa <- c(0.9, 0.99, 0.999, 0.9999)
sapply(kappa, VaR_S)
sapply(kappa, TVaR_S)

###pacakge
library(copula)
library(nCopula)
structure_c <- LOG(0.5, NULL, list(GEO(0.8, 1:40, NULL),
                                   GEO(0.9,1:40, NULL)))
U <- rCompCop(10^6, structure_c)
LOG
mean(S)
var(S)

S <- sort(S)

VaR_S <- function(k) S[ceiling(k*length(S))]
TVaR_S <- function(k) mean(S[S>VaR_S(k)])
t <- new("Log_Child", parameter = as.character(0.5), 
         arg = 1, dimension = length(1), name = "Logarithmic distribution", 
         type = "Child", obj = "Log")

kappa <- c(0.9, 0.99, 0.999, 0.9999)
sapply(kappa, VaR_S)
sapply(kappa, TVaR_S)

###function

rCopula_single <- function(m,M_loi,para_M,B_loi,para_B){
  if(M_loi == 'GEO'){
    sim_M <- function(mm){rgeom(mm, para_M)+1}
    LTS_M <- function(t,q){
      (exp(-t)* q / (1-(1-q)*exp(-t)))
    }
  }
  if(M_loi == 'LOG'){
    sim_M <- function(mm){rlogarithmic(mm, para_M)}
    LTS_M <- function(t,q){
      log(1-q*exp(-t))/log(1-q)
    }
  }
  if(M_loi == 'GAMMA'){
    sim_M <- function(mm){rgamma(mm, para_M,1)}
    LTS_M <- function(t,a){
      (1/(1-t))^a
    }
  }
  if(B_loi == 'GEO'){
    sim_B <- function(mm){rgeom(mm, para_M)+1}
    LTS_B <- function(t,q){
      (exp(-t)* q / (1-(1-q)*exp(-t)))
    }
  }
  if(B_loi == 'LOG'){
    sim_B <- function(mm){rlogarithmic(mm, para_M)}
    LTS_B <- function(t,q){
      log(1-q*exp(-t))/log(1-q)
    }
  }
  if(B_loi == 'GAMMA'){
    sim_B <- function(mm){rgamma(mm, para_M,1)}
    LTS_B <- function(t,a){
      (1/(1-t))^a
    }
  }
  
  M = sim_M(m)
  V = t(sapply(M, function(x){
    c(sum(sim_B(x)),sum(sim_B(x)))
  }))
  U = t(sapply(1:m, function(i){
    t1 = rexp(n1) / V[i,1]
    t2 = rexp(n2) / V[i,2]
    c(
      LTS_M(
        -log(LTS_B(t1, 0.8)),
        0.5),
      LTS_M(
        -log(LTS_B(t2, 0.9)),
        0.5)
    )
  }))
}