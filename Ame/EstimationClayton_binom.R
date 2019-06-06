Cclayton_generateur <- function(delta){
  function(u){
    (sum(u^(-delta) )- length(u)+1)^(-1/delta)
  }
}


pCopula_generateur <- function(delta,FFN,FFX){
  CC <- Cclayton_generateur(delta)
  function(X){
    end <- length(X)
    vec <- c(FFN(X[1]),FFX(X[2:end]))
    CC(vec)
  }
}

C1_2_cond <- Vectorize(function(u,u2,alpha){
  ((u * u2^(alpha+1))^(- alpha / (alpha+1)) - u2^(-alpha)+1)^(-1/alpha)
})

# Simulation --------------------------------------------------------------
set.seed(2^30)
sim_clayton <- function(m,dim,alpha){
  U <- matrix(rep(0,m*dim),m,dim)
  Y_cond <- matrix(rep(0,m*dim),m,dim)
  V <- matrix(runif(m*dim),m,dim)
  
  for (i in 1:m){
    theta <- rgamma(1,1/alpha,1)
    Y_cond[i,] <- (-1/theta) * log(1-V[i,])
    U[i,] <- (1/(1+Y_cond[i,]))^(1/alpha)
  }
  return(U)
}

# UU <- sim_clayton(1000000,5,5)
# UU_N_X <- sim_clayton(1000000,2,10)
# 
# plot(UU[1:10000,1],UU[1:10000,5])
# cor(UU[1:10000,1],UU[1:10000,5],method = "kendall")

# Estimation --------------------------------------------------------------

q <- 0.4

al <- 1000



alpha_0 <- 10


mm <- 10^4
UU <- sim_clayton(mm,6,alpha_0)

FREQ <- qbinom(UU[,1],5,q)

#CHECK!
FREQ[FREQ>=5] #pour dérivé
#CHECK!
n <- 5

LOSS <- list()

X <- matrix(rep(0,mm*n),mm,5)

X_t <- matrix(rep(0,mm*n),mm,5)

for (i in 1:mm){
  X_t[i,] <- qexp(UU[i,2:6],1/al)
  if(FREQ[i]==0){
    LOSS[[i]] <- 0
  }else{
    LOSS[[i]] <- qexp(UU[i,2:(FREQ[i]+1)],1/al)
    X[i,1:FREQ[i]] <- qexp(UU[i,2:(FREQ[i]+1)],1/al)
  }
}


summary(X_t)

LOSS_unlist <- unlist(LOSS)

mean(LOSS_unlist)

LOSS_unlist <- LOSS_unlist[LOSS_unlist>0]




mean(LOSS_unlist)
mean(FREQ)

plot(ecdf(LOSS_unlist))



calcule_deriv <- function(F_x,clay_str,dim){
  C_clay <- list()
  
  dC01 <- list()
  for (j in 1:dim){
    somme <- "u0^(-pa1) + u1^(-pa1)"
    if(j>1){
      for(i in 2:j){
        somme <- paste(paste(paste(somme, " + u",sep=""),as.character(i),sep=""),"^(-pa1)",sep="")
      }
    }
    somme <- paste(somme,as.character(j),sep=" - ")
 
    C_clay[[j]] <- stringr::str_replace_all(clay_str,"U",somme)
    #print(C_clay[[j]])
    for ( i in 1:j){
      C_clay[[j]] <- stringr::str_replace(C_clay[[j]], paste("u",as.character(i),sep=""), F_x)
    }
    
    for(i in 1:j){
      C_clay[[j]] <- stringr::str_replace(C_clay[[j]], "X", paste("x",as.character(i),sep=""))
    }
    
    dC01_x <- Deriv::Deriv(C_clay[[j]], "x1", cache.exp = F)
    if(j>1){
      for ( i in 2:j){
        time <- system.time(dC01_x <- Deriv::Deriv(dC01_x, paste("x",as.character(i),sep = ""), cache.exp = F))
        print(i)
        print(time)
      }
    }
    
    dC01[[j]] <- dC01_x
  }
  
  dC01_parse <- list()
  dC01_parse <- lapply(dC01,function(l) parse(text = l))
  return(dC01_parse)
}
  
generateur_evalue_deriv <- function(dC01_v){
  function(para,U,derive){
    
    for(i in 1:length(para)) { 
      nam <- paste("pa", i, sep = "")
      assign(nam, para[i])
    }
    u0 <- U[1]
    for(i in 1:length(U)) { 
      nam <- paste("x", i, sep = "")
      assign(nam, U[i+1])
    }
    return(eval(dC01_v[[derive]]))
  }
}
F_x <- "(1-exp(-pa2 * X))"
clay_str <- "( U )^(-1/pa1)"
dCC01 <- calcule_deriv(F_x = F_x,clay_str = clay_str,5)

calc_deriv_dC01 <- generateur_evalue_deriv(dCC01)

 
neg_log_vrais_C01 <- function(para,data_N=FREQ, data_X=X){
  
  nn <- length(data_N)
  F_N <- pbinom(data_N,5,para[3])
  F_N_1 <- pbinom(data_N-1,5,para[3])
  d <- rep(0,nn)
  for( i in 1:nn){
    N <- data_N[i]
    if(N==0){
      d[i] <- pbinom(N,5,para[3])
    }else{
      d[i] <-  calc_deriv_dC01(para[1:2],c(F_N[i],data_X[i,1:N]),N)-  calc_deriv_dC01(para[1:2],c(F_N_1[i],data_X[i,1:N]),N)
    }
  }
  #print(para)
  #print(d)
  return(-sum(log(d)))
}



pbinom(0,5,0.2)

val_depart <- c(alpha_0,1/al, q)

system.time(mle <- constrOptim(val_depart, neg_log_vrais_C01, grad = NULL, 
                               ui = diag(3), ci = c(0, 0,0),outer.eps = .Machine$double.eps )
)

##.Machine$double.eps

mle


