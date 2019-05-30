# Exemple pour Diamilatou et Alexandre
# 
# Date : 2018-05-08
#
#
nlignes<-27
ncolonnes<-4
pNX1X2<-matrix(0,nlignes,ncolonnes)
pNX1X2
p1<-0.8
p2<-1-p1
a1<-0.2
a2<-0.5
b1<-0.4
b2<-0.9

j<-0
for (n in 0:2)
{
  for (k1 in 1:3)
  {
    for (k2 in 1:3)
    {
      j<-j+1
            pNX1X2[j,]<-c(n,k1,k2,p1*dbinom(n,2,a1)*dbinom(k1-1,2,b1)*dbinom(k2-1,2,b1)+p2*dbinom(n,2,a2)*dbinom(k1-1,2,b2)*dbinom(k2-1,2,b2))
    }
  }
}

pNX1X2
sum(pNX1X2[,4])

n <- 0:2

ff <- Vectorize(function(s){
  sum <- 0
  sum <- sum + sum(pNX1X2[(pNX1X2[,1]==1) &((pNX1X2[,2])==s) ,4])
    sum <- sum + sum(pNX1X2[(pNX1X2[,1]==2) &((pNX1X2[,2]+pNX1X2[,3])==s) ,4])
  return(sum)
})

pNX1X2[pNX1X2[,1]==0,4]

f_S <- c(sum(pNX1X2[pNX1X2[,1]==0,4]),ff(1:6))
sum(f_S)
f_S
(E_S <- sum((0:6) * f_S))
(V_S <- sum((0:6)^2 * f_S) - E_S^2)


f_N <- c(sum(pNX1X2[pNX1X2[,1]==0,4]),sum(pNX1X2[pNX1X2[,1]==1,4]),sum(pNX1X2[pNX1X2[,1]==2,4]))
f_X1 <- c(sum(pNX1X2[pNX1X2[,2]==1,4]),sum(pNX1X2[pNX1X2[,2]==2,4]),sum(pNX1X2[pNX1X2[,2]==3,4]))
f_X2 <- c(sum(pNX1X2[pNX1X2[,3]==1,4]),sum(pNX1X2[pNX1X2[,3]==2,4]),sum(pNX1X2[pNX1X2[,3]==3,4]))
#2

ff_ind <- Vectorize(function(s){
  sum <- 0
  if(s<=3){
    sum <- f_N[2]* f_X1[s]
  }
  sum <- sum+ f_N[3]*sum(sapply(1:3,function(i) sapply(1:3,function(j)  sum((i+j==s)*f_X1[i]*f_X2[j]) )))
})

f_S_ind <- c(f_N[1],ff_ind(1:6))
sum(f_S_ind)

f_S_ind
f_S

(E_S_ind <- sum((0:6) * f_S_ind))
(V_S_ind <- sum((0:6)^2 * f_S_ind) - E_S_ind^2)



