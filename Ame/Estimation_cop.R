set.seed(2^30)
LTS_M <- function(t,alpha0){
  (-1/alpha0) * log(1-(1-exp(-alpha0))*exp(-t))
}

TLS_geo <- function(t, q){
  q * exp(-t) / (1 - (1 - q) * exp(-t))
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

# Simulation --------------------------------------------------------------
lam <- 1

##algo 13
N_MAX <- qpois(0.99999,lam)


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
    M[i] <- qgeom(runif(1),1-exp(-alpha_0))+1
    #sim R
    R[i] <- qexp(runif(1),1)
    #sim U0
    U_0[i] <- TLS_geo(R[i]/M[i],1-exp(-alpha_0))
    #sim N
    N[i] <- qpois(U_0[i],lam)
    
    if(N[i]==0){
      X[[i]] <- rep(0,N_MAX)
    }else{
      theta_0_1[i] <- sum(qgamma(runif(M[i]),1/alpha_1,1))
      R_0_1 <- qexp(runif(N[i]),1)
      U_0_1 <- LTS_M(-log(LTS_B(R_0_1/theta_0_1[i],alpha_1)),alpha_0)
      X[[i]] <- c(qpareto(U_0_1,3,100),rep(0,N_MAX-length(U_0_1))) 
    }
  }
  return(X)
}
m <- 10000
XX <- sim_algo14(m,a0_1,a1_1)
length(XX)

LOSS <- unlist(XX)
LOSS <- LOSS[LOSS>0]
FREQ_LOSS <- unlist(lapply(XX,function(i) sum(i>0)))
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


# Estimation --------------------------------------------------------------


###estimation copule archimedienne

TLS_geo <- function(t, q){
  q * exp(-t) / (1 - (1 - q) * exp(-t))
}

TLS_geo_deriv <- function(t, q){
  (-q * exp(t)) / ((q + exp(t) - 1)^2)
}

TLS_geo_deriv_deriv <- function(t, q){
  (q * exp(t) * (-q + exp(t) + 1)) / ((q + exp(t) - 1)^3)
}

TLS_geo_inv <- function(y, q){
  ans <- -log(y / (q + y * (1 - q)))
  return(ans)
}

TLS_geo_inv_deriv <- function(y, q){
  q / (y * (q * (y - 1) - y))
}

P_geo <- function(t, q){
  q * t / (1 - (1 - q) * t)
}

P_deriv <- function(t, q){
  q / ((q - 1) * t + 1)^2
}

P_deriv_deriv <- function(t, q){
  -2 * (q - 1) * q / (((q - 1) * t + 1)^3)
}

TLS_gam <- function(t, al, be=1){
  ans <- (be / (be + t))^al
  return(ans)
}

density_cop_ROOT <- function(u, param,TLS_M_inv,TLS_M_inv_deriv,TLS_M_deriv_deriv){
  x <- sapply(u, function(i) TLS_M_inv(i, param))
  TLS_M_deriv_deriv(sum(x), param) * TLS_M_inv_deriv(u[1], param) * TLS_M_inv_deriv(u[2], param)
}

ll=function(u,x){
  s <- 0
  u1 <- u[1]
  u2 <- u[2:(N_MAX+1)]
  
  u2 <- u2[u2>0]
   if(length(u2)==0){
     s <- s+log(u1)
   }else{
     W <- expand.grid(u1,u2)#,u3)
     l <- dim(W)[1]
     for(i in 1:l) s <- s + log(density_cop_ROOT(unlist(W[i,]),x,TLS_geo_inv,TLS_geo_inv_deriv,TLS_geo_deriv_deriv))
   }
  return(s)
}



##### Estimation du M (root) --------------------------------------------------

##calculs de U empiriques

#frequences
F_N_emp <- Vectorize(function(x){
  sum(FREQ_LOSS <= x)/length(FREQ_LOSS)
})
U0 <- F_N_emp(FREQ_LOSS)

#U0 <- U0[FREQ_LOSS>0]
n <- length(U0)
#Montants

X <- matrix(rep(0,N_MAX*m),m,N_MAX)
U <- matrix(rep(0,(N_MAX+1)*n),n,(N_MAX+1))
for (i in 1:N_MAX){
  X[,i] <- unlist(lapply(XX, function(x) x[i]))
}

F_X_emp <- Vectorize(function(x){
  sum(LOSS <= x)/length(LOSS)
})

U[,1] <- U0
for (i in 2:(N_MAX+1)){
#  F_X_emp <- Vectorize(function(x){
#    sum(X[,i-1] <= x)/length(X[,i-1])
#  })
#  U[,i] <- F_X_emp(X[FREQ_LOSS>0,i-1])
  U[,i] <- F_X_emp(X[,i-1])
}
U[3:8,]

X[3:8,]

neg_logvrais <- function(x,UU) -sum(sapply(1:n,function(j) ll(UU[j,],x)))

(gam_M <- optimize(neg_logvrais, c(0.4, 0.9), UU = U)$minimum)

# THETA (0, 1)

#library(copula)
#library(nCopula)
#structure <- GEO(1-exp(-a0_1), 1, list(GAMMA(1 / a1_1, 1:12, NULL)))

PGF_M <- "(gamma)*(z) / (1 - (1-(gamma))*(z))"
PGFInv_M <- "1 / (((gamma)/(z)) + (1 - (gamma)))"
Lap_M <- "(gamma)*exp(-(z)) / (1 - (1 - (gamma)) * exp(-(z)))"
LapInv_M <- "-log(1 / (((gamma)/(z)) + (1 - (gamma))))"

Lap_B <- "(1 / (1 + (z)))^(alpha)"
LapInv_B <- "((z)^(-1/(alpha)) - 1)"

PGF_M
C01 <- stringr::str_replace_all(PGF_M, "z", as.character(Lap_B))
somme <- "z1 + z2 + z3 + z4 + z5+ z6 + z7 + z8 "#+ z9 + z10 + z11 + z12 "

C01 <- stringr::str_replace_all(C01, "z", as.character(somme))



PGFInv_M_U <- list()
LapInv_B_U <- list()
for(i in 1:N_MAX){
  PGFInv_M_U[[i]] <- stringr::str_replace_all(PGFInv_M, "z", paste("u",as.character(i),sep = ""))
  LapInv_B_U[[i]] <- stringr::str_replace_all(LapInv_B, "z", PGFInv_M_U[[i]])
}
C01

#


for(i in 1:N_MAX){
  C01 <- stringr::str_replace(C01,paste("z",as.character(i),sep = ""), LapInv_B_U[[i]])
}
for(i in 1:N_MAX){
  C01 <- stringr::str_replace(C01,paste("z",as.character(i),sep = ""), LapInv_B_U[[i]])
}

C01

dC01 <- Deriv::Deriv(C01, "u1", cache.exp = F)


for ( i in 2:N_MAX){
  time <- system.time(dC01 <- Deriv::Deriv(dC01, paste("u",as.character(i),sep = ""), cache.exp = F))
  print(i)
  print(time)
}

C01 <- stringr::str_replace_all(C01, "gamma", as.character(gam_M))
dC01 <- stringr::str_replace_all(dC01, "gamma", as.character(gam_M))

dC01 <- parse(text = dC01)

dC01

UU[,9]

neg_logv_C01 <- function(alpha, data){
  u1 <- data[,2]
  u2 <- data[,3]
  u3 <- data[,4]
  u4 <- data[,5]
  u5 <- data[,6]
  u6 <- data[,7]
  u7 <- data[,8]
  u8 <- data[,9]
  -sum(log(eval(dC01)))
}

(alpha1<- optimize(neg_logv_C01, c(0, 1), data = UU)$minimum)


