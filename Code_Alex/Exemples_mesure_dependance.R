

tau_kendall <- function(X,Y){
    n <- length(X)
    concord <- outer(1:n,1:n, function(i,j) sign((X[i] - X[j]) * (Y[i] - Y[j])))
    E_n <- (length(concord[concord == 0]) - n) / 2 # On ajoute la probabilité d'avoir des valeurs égales.
    P_n <- sum(concord[concord > 0]) / 2
    tau <- (4 * P_n + 2 * E_n) / (n * (n - 1)) - 1
    return(tau)
}
 
round(tau_kendall(X,Y),8)
cor(X,Y,method = "kendall")

library(copula)
tau(fgmCopula(-1)) # Permet de trouver le tau de kendall en fonction du paramètre de dépendance
plot(alpha <- (-100:100)/100, sapply(alpha, function(a) tau(fgmCopula(a))), type="l")

simul <- rCopula(1e+4, fgmCopula(1))
simul <- qpois(simul,1)
cor(simul, method = "kendall")

plot(N <- (1:30), 
     sapply(N, function(n)
         tau_kendall(rbinom(100, n, 0.5), rbinom(100, n, 0.5))),
     ylim = c(-0.2222,0.2222)
     )

{
# Scénario 1 - tiré de l'article de Christian Genest : Everything
X <- c(-2.224, -1.538, -0.807, 0.024, 0.052, 1.324)
Y <- c(0.431, 1.035, 0.586, 1.465, 1.115, -0.847)

#Scénario 2 - Simulation perso - Avec une distribution discrète, 
# les résultats obtenus donnent des résultats différents de la fonction
# préprogrammé dans R. Cela implique que cette dernière utilise une 
# méthode de traitement des valeurs répétées.
X <- rpois(20,10)
Y <- rpois(20,10)
rbind(X,Y)

#Scénario 3 - Simulation perso
X <- rexp(20,50)
Y <- rexp(20,50)
} # Exemples


rho_spearman <- function(X,Y){
    n <- length(X)
    R_i <- rank(X, ties.method = "average")
    S_i <- rank(Y, ties.method = "average")
    rho <- 12 / (n * (n + 1) * (n - 1)) * sum(R_i*S_i) - 3 * (n + 1) / (n - 1)
    return(rho)
}

rho_spearman(X,Y)
cor(X,Y,method = "spearman")
