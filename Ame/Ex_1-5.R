###1.5

# 5 -----------------------------------------------------------------------


n <- c(0,1,2)
X <- seq(1,20)
fX <- ( (4/(3+X))^3 - (4/(4+X))^3 )/(1-(4/24)^3)
fN <- c(0.4,0.5,0.1)
sum(fX)
FN <- cumsum(fN)
FX <- cumsum(fX)

Cclayton_generateur <- function(delta){
  function(u){
    (sum(u^(-delta) )- length(u)+1)^(-1/delta)
  }
}

pCopula_generateur <- function(delta,FFN,FFX){
  CC <- Cclayton_generateur(delta)
  function(n,k1,k2){
    vec <- c(FFN[n+1],FFX[k1],FFX[k2])
    CC(vec)
  }
}

dCopula <- function(n,x1,x2,delta,FFN,FFX){
  pCopula <- pCopula_generateur(delta,FFN,FFX)
  d <- pCopula(n,x1,x2)
  if(n > 0){
    d <- d- pCopula(n-1,x1,x2)
  }
  if(x1-1 > 0){
    d <- d- pCopula(n,x1-1,x2)
  }
  if(x2-1 > 0){
    d <- d- pCopula(n,x1,x2-1) 
  }
  if( (n > 0) && (x1-1 >0)){
    d <- d+ pCopula(n-1,x1-1,x2)
  }
  if( (x1-1 >0) && (x2 -1 > 0)){
    d <- d+ pCopula(n,x1-1,x2-1 )
  }
  if( (n>0) && (x2 -1 > 0)){
    d <- d+ pCopula(n-1,x1,x2-1) 
  }
  if( (n>0) &&(x1-1 >0) && (x2 -1 > 0)){
    d <- d- pCopula(n-1,x1-1,x2-1)
  }
  return(d)
}


f_2<- matrix(c(rep(n,each=400),rep(rep(X,each=20),3), rep(X,60),rep(0,1200)),1200,4)
f_5<- matrix(c(rep(n,each=400),rep(rep(X,each=20),3), rep(X,60),rep(0,1200)),1200,4)
f_10<- matrix(c(rep(n,each=400),rep(rep(X,each=20),3), rep(X,60),rep(0,1200)),1200,4)
for (i in 1:1200) {
  f_2[i,4]<- dCopula(f_2[i,1],f_2[i,2],f_2[i,3],2,FN,FX)
  f_5[i,4]<- dCopula(f_5[i,1],f_5[i,2],f_5[i,3],5,FN,FX)
  f_10[i,4]<- dCopula(f_10[i,1],f_10[i,2],f_10[i,3],10,FN,FX)
}

sum(f_2[,4])
sum(f_5[,4])
sum(f_10[,4])

f_S <- function(x,ff){
  s <- 0
  if(x==0){
    s <- sum(ff[ff[,1]==0,4])
  }else{
    s <- sum(ff[(ff[,1]==1)&(ff[,2]==x),4]) + sum(ff[(ff[,1]==2)&(ff[,2]+ff[,3]==x),4])
  }
  return(s)
}

fS_2 <- sapply(0:40,function(k) f_S(k,f_2))
fS_5 <- sapply(0:40,function(k) f_S(k,f_5))
fS_10 <- sapply(0:40,function(k) f_S(k,f_10))

sum(fS_2)
sum(fS_5)
sum(fS_10)

FS_2 <- cumsum(fS_2)
FS_5 <- cumsum(fS_5)
FS_10 <- cumsum(fS_10)

df <- data.frame(Delta = rep(c(2,5,10),each=41), F_S = c(FS_2,FS_5,FS_10))


plot(0:40,FS_2,type="S",ylab = "F_S(x)",xlab="x")
lines(FS_5,col="blue",type="S")
lines(FS_10,col="red",type="S")

(E_S_2 <- sum( (0:40) * fS_2))
(E_S_5 <- sum( (0:40) * fS_5))
(E_S_10 <- sum( (0:40) * fS_10))


(Stoploss_S_2 <- sapply(0:40,function(k) sum( pmax((0:40)-k,0) * fS_2)))
(Stoploss_S_5 <- sapply(0:40,function(k) sum( pmax((0:40)-k,0) * fS_5)))
(Stoploss_S_10 <- sapply(0:40,function(k) sum( pmax((0:40)-k,0) * fS_10)))

(TVaR_S_2 <- sapply(0:40,function(k) sum( pmax((0:40)-k,0) * fS_2)/(1-FS_2[k+1]) + k  ))
(TVaR_S_5 <- sapply(0:40,function(k) sum( pmax((0:40)-k,0) * fS_5)/(1-FS_5[k+1]) + k  ))
(TVaR_S_10 <- sapply(0:40,function(k) sum( pmax((0:40)-k,0) * fS_10)/(1-FS_10[k+1]) + k  ))



# 6 -----------------------------------------------------------------------


n <- c(0,1,2)
X <- seq(1,20)
fX <- ( (4/(3+X))^3 - (4/(4+X))^3 )/(1-(4/24)^3)
fN <- c(0.4,0.5,0.1)
sum(fX)
FN <- cumsum(fN)
FX <- cumsum(fX)

Cclayton_generateur <- function(delta){
  function(u){
    (sum(u^(-delta) )- length(u)+1)^(-1/delta)
  }
}

pCopula_generateur <- function(delta,FFN,FFX){
  CC_1 <- Cclayton_generateur(delta[1])
  CC_2 <- Cclayton_generateur(delta[2])
  function(n,k1,k2){
    CC_1(c(FFN[n+1],CC_2(c(FFX[k1],FFX[k2]))))
  }
}

dCopula <- function(n,x1,x2,delta,FFN,FFX){
  pCopula <- pCopula_generateur(delta,FFN,FFX)
  d <- pCopula(n,x1,x2)
  if(n > 0){
    d <- d- pCopula(n-1,x1,x2)
  }
  if(x1-1 > 0){
    d <- d- pCopula(n,x1-1,x2)
  }
  if(x2-1 > 0){
    d <- d- pCopula(n,x1,x2-1) 
  }
  if( (n > 0) && (x1-1 >0)){
    d <- d+ pCopula(n-1,x1-1,x2)
  }
  if( (x1-1 >0) && (x2 -1 > 0)){
    d <- d+ pCopula(n,x1-1,x2-1 )
  }
  if( (n>0) && (x2 -1 > 0)){
    d <- d+ pCopula(n-1,x1,x2-1) 
  }
  if( (n>0) &&(x1-1 >0) && (x2 -1 > 0)){
    d <- d- pCopula(n-1,x1-1,x2-1)
  }
  return(d)
}


f_2_5<- matrix(c(rep(n,each=400),rep(rep(X,each=20),3), rep(X,60),rep(0,1200)),1200,4)

for (i in 1:1200) {
  f_2_5[i,4]<- dCopula(f_2_5[i,1],f_2_5[i,2],f_2_5[i,3],c(2,5),FN,FX)
}

sum(f_2_5[,4])


f_S <- function(x,ff){
  s <- 0
  if(x==0){
    s <- sum(ff[ff[,1]==0,4])
  }else{
    s <- sum(ff[(ff[,1]==1)&(ff[,2]==x),4]) + sum(ff[(ff[,1]==2)&(ff[,2]+ff[,3]==x),4])
  }
  return(s)
}

fS_2_5 <- sapply(0:40,function(k) f_S(k,f_2_5))


sum(fS_2_5)

FS_2_5 <- cumsum(fS_2_5)

df <- data.frame(Delta = rep(c(2,5,10,"imbriquée"),each=41), F_S = c(FS_2,FS_5,FS_10,FS_2_5))


plot(0:40,FS_2,type="S",ylab = "F_S(x)",xlab="x")
lines(0:40,FS_5,col="blue",type="S")
lines(0:40,FS_10,col="green",type="S")
lines(0:40,FS_2_5,col="red",type="S")

(E_S_2_5 <- sum( (0:40) * fS_2_5))


(Stoploss_S_2_5 <- sapply(0:40,function(k) sum( pmax((0:40)-k,0) * fS_2_5)))


(TVaR_S_2_5 <- sapply(0:40,function(k) sum( pmax((0:40)-k,0) * fS_2_5)/(1-FS_2_5[k+1]) + k  ))

FS_2_5
FS_2

CC_t1 <- function(u1,u2,u3){
  ( u1^(-2) + ((u2^(-5) + u3^(-5)-1)^(-1/5))^(-2)  -1  )^(-1/2)
}
CC_t2 <- function(u1,u2,u3){
  ( u1^(-2) + u2^(-2) + u3^(-2)  -2  )^(-1/2)
}

CC_t1(0.5,0.5,0.5)
CC_t2(0.5,0.5,0.5)
