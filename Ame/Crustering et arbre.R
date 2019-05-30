##exemple 2

K <- matrix(c(1,0.407,0.192,0.192,0.192,0.192,
              0.407,1,0.192,0.192,0.192,0.192,
              0.192,0.192,1,0.665, 0.441, 0.441,
              0.192,0.192,0.665,1,0.441, 0.441,
              0.192,0.192,0.441,0.441,1,0.545,
              0.192,0.192,0.441,0.441,0.545,1),6,6)

distance <- function(KK){
  col <- dim(KK)[2]
  ligne <- dim(KK)[1]
  res <- matrix(rep(0,length(KK)),ligne,col)
  
  sapply(1:ligne,function(i){
    sapply(1:col,function(j){
      if(i==j){
        res[i,j] <- 100000
      }else{
        res[i,j] <- max(abs(KK[i,]-KK[j,]))
      }
      
      })})
}

D1 <- distance(K)

D <- 1-K

ind <- which(D1==min(D1), arr.ind=TRUE)

ind

ave <-( D[ind[1,1],]+D[ind[1,2],])/2
ave
D
cluster <- function(KK){
  DD <- distance(KK)
  

  level <- 0
  m <- 0
  
  while (condition) {
    col <- dim(DD)[2]-1
    ligne <- dim(DD)[1]-1
    min_dist_coord <- which(DD==min(DD), arr.ind=TRUE)[1]
    DD <- 1-KK
    average <- ( DD[min_dist_coord[1]]+DD[min_dist_coord[2]])/2
    m <- m+1
    res <- matrix(rep(0,ligne*col),ligne,col)
    
    sapply(1:ligne,function(i){
      res[i,] <- DD[i,]
      })
  }
  
}

