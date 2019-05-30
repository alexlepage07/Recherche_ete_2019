alpha <- 0.5

Cclayton <- function(u){
  (sum( u^(-alpha) ) - length(u) + 1)^(-1/alpha)
}

pCopula <- function(U,copul){
  copul(U)
}

FrN <- c(0.562,0.918,1.000)
FrX <- c(0.29,0.71,1.00)
pCopula(c(FrN[1],FrX[1],FrX[1]),Cclayton)

dCopula <- function(n,x1,x2,FrN,FrX,copul){
  d <- pCopula(c(FrN[n+1],FrX[x1],FrX[x2]),copul)
  if(n > 0){
    d <- d- pCopula(c(FrN[n+1-1],FrX[x1],FrX[x2]),copul)
  }
  if(x1-1 > 0){
    d <- d- pCopula(c(FrN[n+1],FrX[x1-1],FrX[x2]),copul)
  }
  if(x2-1 > 0){
    d <- d- pCopula(c(FrN[n+1],FrX[x1],FrX[x2-1]),copul)
  }
  if( (n > 0) && (x1-1 >0)){
    d <- d+ pCopula(c(FrN[n+1-1],FrX[x1-1],FrX[x2]),copul)
  }
  if( (x1-1 >0) && (x2 -1 > 0)){
    d <- d+ pCopula(c(FrN[n+1],FrX[x1-1],FrX[x2-1]),copul)
  }
  if( (n>0) && (x2 -1 > 0)){
    d <- d+ pCopula(c(FrN[n+1-1],FrX[x1],FrX[x2-1]),copul)
  }
  if( (n>0) &&(x1-1 >0) && (x2 -1 > 0)){
    d <- d- pCopula(c(FrN[n+1-1],FrX[x1-1],FrX[x2-1]),copul)
  }
return(d)
}

p3<- matrix(c(rep(0:2,each=9),rep(rep(1:3,each=3),3), rep(1:3,9),rep(0,27)),27,4)
for (i in 1:27) {
  p3[i,4]<- dCopula(p3[i,1],p3[i,2],p3[i,3],FrN,FrX,Cclayton)
}


sum(p3[,4])
  
  