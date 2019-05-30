lam <- 5
l_exp <- 1
theta <- c(1/3,2/3)

cGumbel_generateur <- function(alpha){
  function(u){
    exp(-(sum(  (-log(u))^alpha))^(1/alpha))
  }
}

cFrank_generateur <- function(alpha){
  function(u){
    (-1/alpha) * log(1+ prod(exp(-alpha*u)-1)/(exp(-alpha)-1)^(length(u)-1)  )
  }
}

pCopula_generateur <- function(alpha,FFN,FFV,cop_generateur){
  CC <- cop_generateur(alpha)
  function(n,v){
    vec <- c(FFN[n+1],FFV[v])
    CC(vec)
  }
}

dCopula <- function(n,v,alpha,FFN,FFV,cop_generateur){
  pCopula <- pCopula_generateur(alpha,FFN,FFV,cop_generateur)
  d <- pCopula(n,v)
  if(n > 0){
    d <- d- pCopula(n-1,v)
  }
  if(v-1 > 0){
    d <- d- pCopula(n,v-1)
  }
  if( (n>0) &&(v-1 >0)){
    d <- d+ pCopula(n-1,v-1)
  }
  return(d)
}



generateur_f_nv <- function(alpha,theta,pois_para,exp_para,cop_generateur,fin_n=300,fin_v=1900){
  n <- 0:fin_n
  v <- 0:fin_v
  F_N <- ppois(n,pois_para)
  F_V <- pgeom(v,1-theta)
  v <- v+1
  ff<- matrix(c(rep(n,each=fin_v+1),rep(v,fin_n+1),rep(0,(fin_n+1)*(fin_v+1))),(fin_n+1)*(fin_v+1),3)
  for (i in 1:((fin_n+1)*(fin_v+1))) {
    ff[i,3]<- dCopula(ff[i,1],ff[i,2],alpha,F_N,F_V,cop_generateur)
  }
  return(ff)
}



# Gumbel ------------------------------------------------------------------


f_1_1 <- generateur_f_nv(1,1/3,lam,1,cGumbel_generateur)
sum(f_1_1[,3])
f_1_2 <- generateur_f_nv(1,2/3,lam,1,cGumbel_generateur)
sum(f_1_2[,3])
f_2_1 <- generateur_f_nv(1.5,1/3,lam,1,cGumbel_generateur)
sum(f_2_1[,3])
f_2_2 <- generateur_f_nv(1.5,2/3,lam,1,cGumbel_generateur)
sum(f_2_3[,3])
f_3_1 <- generateur_f_nv(4,1/3,lam,1,cGumbel_generateur)
sum(f_3_1[,3])
f_3_2 <- generateur_f_nv(4,2/3,lam,1,cGumbel_generateur)
sum(f_3_1[,3])


nn <- 0:300
vv <- 1:1901
beta <- 1/(1-theta)

(E_S_1_1 <- sum(f_1_1[,3] *f_1_1[,1]*f_1_1[,2]/beta[1] ))
(E_S_1_2 <- sum(f_1_2[,3] *f_1_2[,1]*f_1_2[,2]/beta[2] ))
(E_S_2_1 <- sum(f_2_1[,3] *f_2_1[,1]*f_2_1[,2]/beta[1] ))
(E_S_2_2 <- sum(f_2_2[,3] *f_2_2[,1]*f_2_2[,2]/beta[2] ))
(E_S_3_1 <- sum(f_3_1[,3] *f_3_1[,1]*f_3_1[,2]/beta[1] ))
(E_S_3_2 <- sum(f_3_2[,3] *f_3_2[,1]*f_3_2[,2]/beta[2] ))

(E_S_3_2 <- sum(f_3_2[,3] *f_3_2[,1]*f_3_2[,2]/beta[2] ))

F_s <- function(s,ff,alpha,theta,pois_para=5,fin_n=300,fin_v=1900){
  n <- 0:fin_n
  v <- 0:fin_v
  F_N <- ppois(n,pois_para)
  F_V <- pgeom(v,1-theta)
  v <- v+1
  beta <- 1/(1-theta)
  sum(ff[ff[,1]==0,3]) + sum ( sapply(1:300,function(i){
    sapply(1:1901,function(j){
      dCopula(i,j,alpha,F_N,F_V) * pgamma(s,i*j,beta) }  )
  }
  ))
}

VaR_S <- function(kapp,ff,theta){
  beta <- 1/(1-theta)
  FF_s <- function(s){
    sum(ff[ff[,1]==0,3]) + sum ( sapply(302:(301*1901),function(i){
     ff[i,3]* pgamma(s,ff[i,1]*ff[i,2],beta) 
    }  
    ))
  }
  optimise(function(x) abs(FF_s(x)-kapp),c(0,100))$minimum
}

kappa <- c(0.9,0.99,0.999)
#system.time(VAR_1_1_9 <- VaR_S(0.9,f_1_1,1/3))

VaR_1_1 <- sapply(kappa,function(k) VaR_S(k,f_1_1,1/3))
VaR_1_2 <- sapply(kappa,function(k) VaR_S(k,f_1_2,2/3))
VaR_2_1 <- sapply(kappa,function(k) VaR_S(k,f_2_1,1/3))
VaR_2_2 <- sapply(kappa,function(k) VaR_S(k,f_2_2,2/3))
VaR_3_1 <- sapply(kappa,function(k) VaR_S(k,f_3_1,1/3))
VaR_3_2 <- sapply(kappa,function(k) VaR_S(k,f_3_2,2/3))

#system.time(VAR_1_2_9 <- VaR_S(0.9,f_1_2,1,2/3))
#system.time(VAR_1_2_99 <- VaR_S(0.99,f_1_2,1,2/3))
#system.time(VAR_1_2_999 <- VaR_S(0.999,f_1_2,1,2/3))

#system.time(VAR_2_1_9 <- VaR_S(0.9,f_2_1,1.5,1/3))
#system.time(VAR_2_1_99 <- VaR_S(0.99,f_2_1,1.5,1/3))
#system.time(VAR_2_1_999 <- VaR_S(0.999,f_2_1,1.5,1/3))

#system.time(VAR_2_2_9 <- VaR_S(0.9,f_2_2,1.5,2/3))
#system.time(VAR_2_2_99 <- VaR_S(0.99,f_2_2,1.5,2/3))
#system.time(VAR_2_2_999 <- VaR_S(0.999,f_2_2,1.5,2/3))

#system.time(VAR_3_1_9 <- VaR_S(0.9,f_3_1,4,1/3))
#system.time(VAR_3_1_99 <- VaR_S(0.99,f_3_1,4,1/3))
#system.time(VAR_3_1_999 <- VaR_S(0.999,f_3_1,4,1/3))

#system.time(VAR_3_2_9 <- VaR_S(0.9,f_3_2,4,2/3))
#system.time(VAR_3_2_99 <- VaR_S(0.99,f_3_2,4,2/3))
#system.time(VAR_3_2_999 <- VaR_S(0.999,f_3_2,4,2/3))

cbind(VaR_1_1,VaR_2_1,VaR_3_1)
cbind(VaR_1_2,VaR_2_2,VaR_3_2)

#f_1_1<- matrix(c(rep(n,each=1901),rep(v,301),rep(0,301*1901)),301*1901,3)

TVaRS <- function(kapp,VaR,ff,beta){
  sum ( sapply(1:(301*1901),function(i){
      ff[i,1] *ff[i,2]*ff[i,3] /beta * (1-pgamma(VaR,ff[i,1]*ff[i,2]+1,beta)) }  )
  )/(1-kapp)
}


TVAR_1_1 <- mapply(function(kap,V) TVaRS(kap,V,f_1_1,beta[1]),kappa,VaR_1_1)
TVAR_1_2 <- mapply(function(kap,V) TVaRS(kap,V,f_1_2,beta[2]),kappa,VaR_1_2)
TVAR_2_1 <- mapply(function(kap,V) TVaRS(kap,V,f_2_1,beta[1]),kappa,VaR_2_1)
TVAR_2_2 <- mapply(function(kap,V) TVaRS(kap,V,f_2_2,beta[2]),kappa,VaR_2_2)
TVAR_3_1 <- mapply(function(kap,V) TVaRS(kap,V,f_3_1,beta[1]),kappa,VaR_3_1)
TVAR_3_2 <- mapply(function(kap,V) TVaRS(kap,V,f_3_2,beta[2]),kappa,VaR_3_2)

cbind(c(VaR_1_1,TVAR_1_1),c(VaR_2_1,TVAR_2_1),c(VaR_3_1,TVAR_3_1))
cbind(c(VaR_1_2,TVAR_1_2),c(VaR_2_2,TVAR_2_2),c(VaR_3_2,TVAR_3_2))


# Frank -------------------------------------------------------------------



f_1_1 <- generateur_f_nv(-10,1/3,lam,1,cFrank_generateur)
sum(f_1_1[,3])
f_1_2 <- generateur_f_nv(-5,1/3,lam,1,cFrank_generateur)
sum(f_1_2[,3])
f_1_3 <- generateur_f_nv(5,1/3,lam,1,cFrank_generateur)
sum(f_1_3[,3])
f_1_4 <- generateur_f_nv(10,1/3,lam,1,cFrank_generateur)
sum(f_1_4[,3])


nn <- 0:300
vv <- 1:1901
beta <- 1/(1-theta)

(E_S_1_1 <- sum(f_1_1[,3] *f_1_1[,1]*f_1_1[,2]/beta[1] ))
(E_S_1_2 <- sum(f_1_2[,3] *f_1_2[,1]*f_1_2[,2]/beta[1] ))
(E_S_1_3 <- sum(f_1_3[,3] *f_1_3[,1]*f_1_3[,2]/beta[1] ))
(E_S_1_4 <- sum(f_1_4[,3] *f_1_4[,1]*f_1_4[,2]/beta[1] ))



kappa <- c(0.9,0.99,0.999)
#system.time(VAR_1_1_9 <- VaR_S(0.9,f_1_1,1/3))

VaR_1_1 <- sapply(kappa,function(k) VaR_S(k,f_1_1,1/3))
VaR_1_2 <- sapply(kappa,function(k) VaR_S(k,f_1_2,1/3))
VaR_1_3 <- sapply(kappa,function(k) VaR_S(k,f_1_3,1/3))
VaR_1_4 <- sapply(kappa,function(k) VaR_S(k,f_1_4,1/3))



cbind(VaR_1_1,VaR_1_2,VaR_1_3,VaR_1_4)


#f_1_1<- matrix(c(rep(n,each=1901),rep(v,301),rep(0,301*1901)),301*1901,3)

TVaRS <- function(kapp,VaR,ff,beta){
  sum ( sapply(1:(301*1901),function(i){
    ff[i,1] *ff[i,2]*ff[i,3] /beta * (1-pgamma(VaR,ff[i,1]*ff[i,2]+1,beta)) }  )
  )/(1-kapp)
}


TVAR_1_1 <- mapply(function(kap,V) TVaRS(kap,V,f_1_1,beta[1]),kappa,VaR_1_1)
TVAR_1_2 <- mapply(function(kap,V) TVaRS(kap,V,f_1_2,beta[1]),kappa,VaR_1_2)
TVAR_1_3 <- mapply(function(kap,V) TVaRS(kap,V,f_1_3,beta[1]),kappa,VaR_1_3)
TVAR_1_4 <- mapply(function(kap,V) TVaRS(kap,V,f_1_4,beta[1]),kappa,VaR_1_4)


cbind(c(VaR_1_1,TVAR_1_1),c(VaR_1_2,TVAR_1_2),c(VaR_1_3,TVAR_1_3),c(VaR_1_4,TVAR_1_4))
