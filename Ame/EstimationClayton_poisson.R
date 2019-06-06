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

lam <- 1
qpois(0.999999,lam)
al <- 100



alpha_0 <- 5

n <- 8
mm <- 10^4
UU <- sim_clayton(mm,n,alpha_0)

FREQ <- qpois(UU[,1],lam)

#CHECK!
FREQ[FREQ>=10] #pour dérivé
#CHECK!



LOSS <- list()

X <- matrix(rep(0,mm*(n-1)),mm,n-1)

X_t <- matrix(rep(0,mm*(n-1)),mm,n-1)

for (i in 1:mm){
  X_t[i,] <- qexp(UU[i,2:n],1/al)
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
        somme <- paste0(somme, " + u",i,"^(-pa1)")
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
dCC01 <- calcule_deriv(F_x = F_x,clay_str = clay_str,max(FREQ))

calc_deriv_dC01 <- generateur_evalue_deriv(dCC01)


neg_log_vrais_C01 <- function(para,data_N=FREQ, data_X=X,pFreq= ppois){
  
  nn <- length(data_N)
  F_N <- pFreq(data_N,para[3])
  F_N_1 <- pFreq(data_N-1,para[3])
  d <- rep(0,nn)
  for( i in 1:nn){
    N <- data_N[i]
    if(N==0){
      d[i] <- F_N[i]
    }else{
      d[i] <-  calc_deriv_dC01(para[1:2],c(F_N[i],data_X[i,1:N]),N)-  calc_deriv_dC01(para[1:2],c(F_N_1[i],data_X[i,1:N]),N)
    }
  }
  #print(para)
  #print(d)
  return(-sum(log(d)))
}



val_depart <- c(alpha_0,1/al, lam)

system.time(mle <- constrOptim(val_depart, neg_log_vrais_C01, grad = NULL, 
                               ui = diag(3), ci = c(0, 0,0),outer.eps = .Machine$double.eps )
)

##.Machine$double.eps

mle
LOSS[[1]]
LOSS[[2]]
LOSS[[3]]
LOSS[[4]]
LOSS[[5]]

hist(FREQ)

max(FREQ)
##hypothèse indép 
neg_log_pois <- function(lam){
  -sum(log(dpois(FREQ,lam)))
}
optimise(neg_log_pois,c(0,5))

neg_log_par <- function(para){
  -sum(log(dexp(LOSS_unlist,1/para)))
}
optimise(neg_log_par,c(0,1000))

