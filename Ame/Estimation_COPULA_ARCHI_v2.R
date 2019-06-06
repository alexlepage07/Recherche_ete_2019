set.seed(2^30)
##algo 13



library(actuar)
sim_algo <- function(mm,alpha_0,alpha_1,q_M_gen,q_B_gen,
                     TLS_M_gen,TLS_B_gen,q_freq,q_sev){
  
  q_M <- q_M_gen(alpha_0)
  q_B <- q_B_gen(alpha_1)
  TLS_M <- TLS_M_gen(alpha_0)
  TLS_B <- TLS_B_gen(alpha_1)
  M <- rep(0,mm)
  R <- rep(0,mm)
  U_0 <- rep(0,mm)
  X <- list()
  N <- rep(0,mm)
  theta_0_1 <- rep(0,mm)
  for (i in 1:mm){
    #sim theta0 aka M
    M[i] <- q_M(runif(1))
    #sim R
    R[i] <- qexp(runif(1),1)
    #sim U0
    U_0[i] <- TLS_M(R[i]/M[i])
    #sim N
    N[i] <- q_freq(U_0[i])
    
    if(N[i]==0){
      X[[i]] <- 0
    }else{
      theta_0_1[i] <- sum(q_B(runif(M[i])))
      R_0_1 <- qexp(runif(N[i]),1)
      U_0_1 <- TLS_M(-log(TLS_B(R_0_1/theta_0_1[i])))
      X[[i]] <- q_sev(U_0_1) 
    }
  }
  return(X)
}

#fonction LTS
LTS_M <- function(t,alpha0){
  (-1/alpha0) * log(1-(1-exp(-alpha0))*exp(-t))
}

LTS_M_gen <- function(para){
  function(t){
    LTS_M(t,para)
  }
}

qu_M_gen <- function(para){
  function(p){
    qlogarithmic(p,1-exp(-para))
  }
}

LTS_B <- function(t,alpha1){
  (1/(1+t))^(1/alpha1)
}

LTS_B_gen <- function(para){
  function(t){
    LTS_B(t,para)
  }
}

qu_B_gen <- function(para){
  function(p){
    qgamma(p,1/para,1)
  }
}

a0_1 <- -log(0.5)
a0_2 <- -log(0.1)

a1_1 <- 5
a1_2 <- 25


# Simulation --------------------------------------------------------------

lam <- 2

m <- 10^4
# XX <- sim_algo(m,a0_1,a1_1,qu_M_gen,qu_B_gen,LTS_M_gen,LTS_B_gen,
#                function(p) qpois(p,lam),function(p) qpareto(p,3,100))
XX <- sim_algo(m,a0_1,a1_1,qu_M_gen,qu_B_gen,LTS_M_gen,LTS_B_gen,
                function(p) qbinom(p,5,0.4),function(p) qexp(p,1/100))
length(XX)

LOSS <- unlist(XX)

LOSS <- LOSS[LOSS>0]
mean(LOSS)
FREQ_LOSS <- unlist(lapply(XX,function(i) sum(i>0)))
mean(FREQ_LOSS)

S <- unlist(lapply(XX,function(i) sum(i)))
mean(S)
###hypothèse indép pois et pareto
# neg_log_pois <- function(lam){
#   -sum(log(dpois(FREQ_LOSS,lam)))
# }
# optimise(neg_log_pois,c(0,5))
# 
# neg_log_par <- function(para){
#   -sum(log(dpareto(LOSS,para[1],para[2])))
# }
# constrOptim(c(3, 100), neg_log_par, grad = NULL, 
#             ui = diag(2), ci = c(0, 0),outer.eps = .Machine$double.eps)


# Estimation --------------------------------------------------------------
N_MAX <- 5
# N_MAX <- qpois(0.99999,lam)
##### Estimation du M (root) et B en même temps --------------------------------------------------


X <- matrix(rep(0,N_MAX*m),m,N_MAX)

LOSS <- list()

for (i in 1:N_MAX){
  X[,i] <- unlist(lapply(XX, function(x) x[i]))
}

X[3:10,]


#library(copula)
#library(nCopula)
#structure <- GEO(1-exp(-a0_1), 1, list(GAMMA(1 / a1_1, 1:12, NULL)))

# PGF_M <- "(gamma)*(z) / (1 - (1-(gamma))*(z))"
# PGFInv_M <- "1 / (((gamma)/(z)) + (1 - (gamma)))"
# 
# Lap_M <- "(gamma)*exp(-(z)) / (1 - (1 - (gamma)) * exp(-(z)))"
# LapInv_M <- "-log(1 / (((gamma)/(z)) + (1 - (gamma))))"

F_x <- "(1-exp(-pa3 * X))"
##M est logarithmique
PGF_M <- "log(1 - (pa1)*(z)) / log(1 - (pa1))"
PGFInv_M <- "((1 - (1 - (pa1))^(z))/(pa1))"

Lap_M <- "log(1 - (pa1) * exp(-(z))) / log(1 - (pa1))"
LapInv_M <- "-log((1 - (1 - (pa1))^(z)) / (pa1))"
##B est gamma
Lap_B <- "(1 / (1 + (z)))^(1/pa2)"
LapInv_B <- "((z)^(-(pa2)) - 1)"

Lap_M
C01 <- stringr::str_replace_all(Lap_M, "z", "z0 - log(z1)")
C01

LAP_INV_M_u0 <- stringr::str_replace_all(LapInv_M, "z", "u0")
LAP_INV_M_u0
C01 <- stringr::str_replace_all(C01, "z0", as.character(LAP_INV_M_u0))
C01

dC01<- list()

for( j in 1:N_MAX){
  somme <- "z1"
  if(j>1){
    for(i in 2:j){
      somme <- paste(paste(somme, " + z",sep=""),as.character(i),sep="")
    }
  }

  #somme <- "z1 + z2 + z3 + z4 + z5 "#+ z6 + z7 + z8 + z9 + z10 + z11 + z12 "
  som_Lap <- stringr::str_replace_all(Lap_B, "z", as.character(somme))
  
  C01_j <- stringr::str_replace_all(C01, "z1", as.character(som_Lap))
  
  PGFInv_M_U <- list()
  LapInv_B_U <- list()
  for(i in 1:j){
    PGFInv_M_U[[i]] <- stringr::str_replace_all(PGFInv_M, "z", F_x)
    PGFInv_M_U[[i]] <- stringr::str_replace_all(PGFInv_M_U[[i]], "X", paste("x",as.character(i),sep = ""))
    LapInv_B_U[[i]] <- stringr::str_replace_all(LapInv_B, "z", PGFInv_M_U[[i]])
  }
  
  
  for(i in 1:j){
    C01_j <- stringr::str_replace(C01_j,paste("z",as.character(i),sep = ""), LapInv_B_U[[i]])
  }
  
  # for(i in 1:N_MAX){
  #   C01 <- stringr::str_replace(C01,paste("z",as.character(i),sep = ""), LapInv_B_U[[i]])
  # }
  
  dC01_x <- Deriv::Deriv(C01_j, "x1", cache.exp = F)
  if(j>1){
    for ( i in 2:j){
      time <- system.time(dC01_x <- Deriv::Deriv(dC01_x, paste("x",as.character(i),sep = ""), cache.exp = F))
      print(i)
      print(time)
    }
  }
  
  dC01[[j]] <- dC01_x

}


dC01[[5]]


dC01_parse <- lapply(dC01,function(l) parse(text = l))

calc_deriv_dC01 <- function(para,U,derive){
  for(i in 1:length(para)) { 
    nam <- paste("pa", i, sep = "")
    assign(nam, para[i])
  }
  u0 <- U[1]
  for(i in 1:length(U)) { 
    nam <- paste("x", i, sep = "")
    assign(nam, U[i+1])
  }
  return(eval(dC01_parse[[derive]]))
}

# calc_deriv_dC01(rep(0.5,2),rep(0.9,6),1)
# 
# calc_deriv_dC01(rep(0.5,2),c(U_N[9,2],U_X[9,1:2]),2)

neg_log_vrais_C01 <- function(para,data_N=FREQ_LOSS, data_X=X){
  
  nn <- length(data_N)
  F_N <- pbinom(data_N,5,para[4])
  F_N_1 <- pbinom(data_N-1,5,para[4])
  d <- rep(0,nn)
  
  for( i in 1:nn){
    N <- data_N[i]
    if(N==0){
      d[i] <- F_N[i]
    }else{
      d[i] <-  calc_deriv_dC01(para[1:3],c(F_N[i],data_X[i,1:N]),N) -  calc_deriv_dC01(para[1:3],c(F_N_1[i],data_X[i,1:N]),N)
    }
  }
  #print(para)
  return(-sum(log(d)))
}

val_depart <- c(0.5,5, 1/100,0.4)
system.time(mle <- constrOptim(val_depart, neg_log_vrais_C01, grad = NULL, 
                   ui = diag(4), ci = c(0, 0,0,0),outer.eps = 0.01 )
)
##.Machine$double.eps

mle

hist(FREQ_LOSS)

XX[[1]]
XX[[2]]
XX[[3]]
XX[[4]]
XX[[5]]

