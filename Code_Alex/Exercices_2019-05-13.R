# semaine du 13 mai 2019
#============================ exercice 3 =================================

F_S <- function(x){
    sapply(x, function(k)
        if (k == 0) 
            return(dbinom(0, 5, 0.3))
        else
            dbinom(0, 5, 0.3) + sum(sapply(1:5, function(n)
                dbinom(n, 5, 0.3) * pgamma(k, 1, scale = n*10)))
    )
}

F_S((0:3)*10)
VaR_S <- function(k) uniroot(function(x) F_S(x) - k, c(0, 1e+10))$root
VaR_0.99 <- VaR_S(0.99)

sum(1-F_S(0:VaR_0.99)) ; 5*0.3 * 10

#============================ exercice 4 =================================
f_N <- c(0.4, 0.5, 0.1)
F_N <- c(0.4, 0.9, 1)


f_X <- function(k) {
    if (k %in% 1:20)
        return(216/215 * ((4 / (3+k))^3 - (4 / (4+k))^3))
    return(0)
    }

F_X <- function(x){
    if (x == 0)
        return(0)
    sum(sapply(1:x, function(k)f_X(k)))
}


C <- function(n, x_1, x_2, delta){
    if (n %in% 0:2 & x_1 %in% 1:20 & x_2 %in% 1:20)
        return((F_N[n+1]^(-delta) + F_X(x_1)^(-delta) + F_X(x_2)^(-delta) - 2)^(-1/delta))
    return(0)
}

f_nX1X2 <- function(n, x_1, x_2, delta){
    C(n, x_1, x_2, delta) -
        C(n-1, x_1, x_2, delta) - C(n, x_1-1, x_2, delta) - C(n, x_1, x_2-1, delta) +
        C(n-1, x_1-1, x_2, delta) + C(n, x_1-1, x_2-1, delta) + C(n-1, x_1, x_2-1, delta) -
        C(n-1, x_1-1, x_2-1, delta)
}


construire_tableau <- function(N, X, delta){
    nrow__ = length(N) * length(X)^2
    
    n <- floor(seq(N[1], N[length(N)]+0.99999999, 1/length(X)^2))
    
    X_1 <- rep(
        floor(seq(X[1], X[length(X)]+0.999999999, 1/length(X))),
        length(N))
    
    X_2 <- rep(X, nrow__ / length(X))
    
    table <- as.matrix(cbind(n, X_1, X_2))
    
    P <- numeric(nrow(table))
    for (i in 1:nrow(table)){
        P[i] <- as.numeric( f_nX1X2(table[i,1], table[i,2], table[i,3], delta))
    }
    sum(P)
    table <- cbind(table, P)
    return(table)
}


Tableau_2 <- construire_tableau(N=0:2, X=1:20, delta = 2)
Tableau_5 <- construire_tableau(N=0:2, X=1:20, delta = 5)
Tableau_10 <- construire_tableau(N=0:2, X=1:20, delta = 10)

head(Tableau_2)
head(Tableau_5)
head(Tableau_10)

source("Exemple_2019-05-08.R")
Modele_dependance(Tableau_2)
Modele_dependance(Tableau_5)
Modele_dependance(Tableau_10)

Modele_dependance(Tableau_2)$TVaR_S
Modele_dependance(Tableau_5)$TVaR_S
Modele_dependance(Tableau_10)$TVaR_S

#============================ exercice 5 =================================
f_N <- c(0.4, 0.5, 0.1)
F_N <- c(0.4, 0.9, 1)


f_X <- function(k) {
    if (k %in% 1:20)
        return(216/215 * ((4 / (3+k))^3 - (4 / (4+k))^3))
    return(0)
}

F_X <- function(x){
    if (x == 0)
        return(0)
    sum(sapply(1:x, function(k)f_X(k)))
}


C <- function(n, x_1, x_2, delta){
    if (n %in% 0:2 & x_1 %in% 1:20 & x_2 %in% 1:20)
        return((F_N[n+1]^(-2) + (F_X(x_1)^(-5) + F_X(x_2)^(-5) -1)^(2/5) - 1)^(-1/2))
    return(0)
}

f_nX1X2 <- function(n, x_1, x_2, delta){
    C(n, x_1, x_2, delta) -
        C(n-1, x_1, x_2, delta) - C(n, x_1-1, x_2, delta) - C(n, x_1, x_2-1, delta) +
        C(n-1, x_1-1, x_2, delta) + C(n, x_1-1, x_2-1, delta) + C(n-1, x_1, x_2-1, delta) -
        C(n-1, x_1-1, x_2-1, delta)
}

Tableau <- construire_tableau(N=0:2, X=1:20, delta = 2)

source("Exemple_2019-05-08.R")
Modele_dependance(Tableau)
