m <- 5
q <- 0.3
beta <- 1/10

F_S <- Vectorize(function(x){
  dbinom(0,m,q) + sum ( dbinom(1:5,m,q) * pexp(x/(1:5),beta))
})

F_S(c(0,10,20,30,40))

m*q/beta
