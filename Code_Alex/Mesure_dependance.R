library(xtable)
library(copula)
library(nCopula)
library(rlist)


#======================================= Tau de Kendall ===========================
#---------------------------------- Calcul du tau de Kendall empirique ------------

tau_kendall_empirique <- function(X,Y){
    # Fonction qui calcule le tau de Kendall empirique lorsqu'au moins une des 
    # deux v.a. est discète.
    n <- length(X)
    concord <- outer(1:n,1:n, function(i,j) sign((X[i] - X[j]) * (Y[i] - Y[j])))
    E_n <- (length(concord[concord == 0]) - n) # On ajoute la probabilité d'avoir des valeurs égales.
    P_n <- sum(concord[concord > 0])
    tau <- (2 * P_n + E_n) / (n * (n - 1)) - 1
    return(tau)
}


bootstrap_tau <- function(data, nb_sous_intervalles = nrow(data) / 10,
                          nb_samples = 100) {
    # Fonction pour utiliser la méthode de ré-échantillonage (bootstrap) sur les
    # données
    tau <- numeric(nb_samples)
    for (i in 1:nb_samples){
        indices <- sample(1:nrow(data), nb_sous_intervalles)
        tau[i] <- tau_kendall_empirique(data[indices, 1],
                                        data[indices, 2])
    }
    return(tau)
}


{
    # Scénario 1 - tiré de l'article de Christian Genest : Everything
    X <- c(-2.224, -1.538, -0.807, 0.024, 0.052, 1.324)
    Y <- c(0.431, 1.035, 0.586, 1.465, 1.115, -0.847)
    
    #Scénario 2 - Simulation perso - Avec une distribution discrète, 
    # les résultats obtenus donnent des résultats différents de la fonction
    # préprogrammé dans R. Cela implique que cette dernière utilise une 
    # méthode de traitement des valeurs répétées.
    set.seed(20190618)
    UU <- rCopula(20, fgmCopula(0.5))
    X <- qpois(UU[,1],10)
    Y <- qpois(UU[,2],10)
    table_XY <- rbind(X,Y)
    xtable(table_XY, digits = 0)
    
    #Scénario 3 - Simulation perso
    X <- rexp(20,50)
    Y <- rexp(20,50)
} # Scénarios
#---------------------------------- Calcul du tau de Kendall théorique ------------

tau_kendall_theorique <- function(F_N, F_X,
                                  para,
                                  Copule,
                                  alpha,
                                  n_max, x_max){
    # Fonction qui permet de calculer le tau de kendall avec une variable aléatoire
    # discrète et une autre qui est continue.
    #
    # F_N peut être une fonction pré-programmée de R, mais F_X doit être une fonction
    # dérivable (pas de boucle, ni de condition, ni d'appel à d'autres fonctions).
    # 
    # Para doit être une liste nommée de la forme 
    # para = list(N=list(size=5, prob=0.3), X=list(par1=1, par2=2, ...))
    #
    # alpha est le paramètre de dépendance de la copule.
            
    F_N. <- function(n) {
        # Paramétriser la fonction de répartition de N.
        do.call(F_N, list.flatten(list(q = n, para$N)))
    }
    F_X. <- function(x) {
        # Paramétriser la fonction de répartition de N.
        do.call(F_X, list.flatten(list(q = x, para$X)))
    }
    
    f_N <- function(n) F_N.(n) - F_N.(n-1)
    
    F_NX <- function(n, x) {
        # Fonction de répartition conjointe de N et X.
        # f <- pCompCop(struct_copule, T , F)
        f <- function(U) pCopula(U, Copule(alpha, use.indepC = "FALSE"))
        return(f(c(F_N.(n), F_X.(x))))
    }
    
    f_NX <- function(n, x) {
        eps <- .Machine$double.eps^0.25
        f <- function(n,x) (F_NX(n, x + eps) - F_NX(n, x)) / eps
        return(f(n, x) - f(n - 1, x))
    }
    
    return(sum(sapply(0:n_max, function(n)
        4 * integrate(
                Vectorize(function(x)
                    F_NX(n - 1, x) * f_NX(n, x)),
                0, x_max,
                subdivisions = 100L,
                rel.tol = .Machine$double.eps ^ 0.25
            )$value +
            f_N(n) ^ 2)) - 1)
}

{
## Scénario 1 : Avec une copule archimédienne hiérarchique
# alpha <- 5
# struct_copule <- LOG(1-exp(-alpha), 1, list(LOG(1-exp(-alpha), 1:2, NULL)))
# F_N <- ppois
# F_X <- function(q, par1) 1-exp(-par1 * q)
# n_max <- qpois(0.99999, 2)
# x_max <- qexp(0.99999, 1/100)
# para <- list(N=list(lambda=2), X=list(par1=1/100))
# tau_kendall_theorique(F_N, F_X, struct_copule, para, alpha, n_max, x_max)

# Scénario 2
alpha <- 5
Copule <- frankCopula
F_N <- ppois
F_X <- function(q, rate) 1-exp(-rate * q)
n_max <- qpois(0.99999, 2)
x_max <- qexp(0.99999, 1/100)
para <- list(N=list(lambda=2), X=list(rate=1/100))

tau_theorique <- tau_kendall_theorique(F_N, F_X,
                                        para,
                                        Copule,
                                        alpha,
                                        n_max, x_max)

# Simulation avec calcul empirique
set.seed(20190618)
UU <- rCopula(1e+4, Copule(alpha))
# UU <- rCompCop(1e+4, struct_copule)
X <- qpois(UU[,1],2)
Y <- qexp(UU[,2], 1/100)

tau_empirique <- tau_kendall_empirique(X, Y)

(tbl_resultats <- rbind(
    "tau empirique (mixte)" = round(tau_empirique, 4),
    "tau théorique (mixte)" = round(tau_theorique, 4),
    "Tau théorique (continu)" = tau(Copule(alpha)),
    "Fonction cor de R" = cor(X, Y, method = "kendall")
))
xtable(tbl_resultats, digits = 4)
} # Scénarios
#---------------------------------- Inversion du tau de Kendall -------------------

inversion_tau_kendall <- function(F_N, F_X,
                         para,
                         Copule,
                         bornes,
                         Tau_Kendall,
                         n_max, x_max){
    # Fonction qui permet de trouver le paramètre de dépendance à l'aide du 
    # tau de Kendall empirique lorsque le modèle comprend une v.a. discrète et
    # une v.a. continue.
    #
    #=========== Bornes: ========
    # Frank: (-infty, infty) / {0}
    # Gumbel : [1, infty]
    # AMH : [-1,1]
    # Joe : [1, infty)
    # Clayton : (0, infty)
    # ===========================
    temps <- system.time(
        alpha <- uniroot(function(alpha)
            tau_kendall_theorique(F_N, F_X,
                                  para,
                                  Copule,
                                  alpha,
                                  n_max, x_max) - Tau_Kendall,
        bornes)$root)
    return(list("alpha"=alpha, "temps"=temps[[3]]))
}

#---------------------------------------------------------------------------------
# Scénario 1 : Frank(5), N suit Pois(2), X suit Exp(1/100) -----------

# Les paramètres
alpha <- 5
Copule <- frankCopula
F_N <- ppois
F_X <- pexp
para <- list(N=list(lambda=2), X=list(rate=1/100))
n_max <- qpois(0.99999, para$N$lambda)
x_max <- qexp(0.99999, para$X$rate)

# Calcul du alpha à l'aide du tau de Kendal empirique
tau_empirique <- numeric(nsim <- 1e+3)
set.seed(20190618)
for (i in 1:nsim) {
    UU <- rCopula(1e+3, Copule(alpha))
    X <- qpois(UU[,1], para$N$lambda)
    Y <- qexp(UU[,2], para$X$rate)
    tau_empirique[i] <- tau_kendall_empirique(X, Y)
}
hist(tau_empirique, freq =T, xlab = "tau_n", ylab = "Fréquence", breaks = 20, main = "")
(var_tau <- var(tau_empirique))
(mean_tau <- mean(tau_empirique))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))
tau_kendall_theorique(F_N, F_X, para, Copule, alpha, n_max, x_max)

# Analyse graphique pour connaître les bornes d'optimisation
Tau <- sapply(c(-5:-1, 1:10), function(a)
    tau_kendall_theorique(F_N, F_X, para, Copule, a, n_max, x_max))

plot(c(-5:-1, 1:10), Tau,
     ylab="tau de kendall",
     xlab="alpha",
     type="l"
)
axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
axis(1,at=-5:10, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
abline(a=tau_empirique, b=0, col="green")

# Bornes trouvées graphiquement.
bornes <- c(4.5, 5.5)

# Optimisation à partir du tau empirique
alpha_trouve <- inversion_tau_kendall(F_N, F_X,
                             para,
                             Copule,
                             bornes,
                             mean_tau,
                             n_max, x_max)
(tbl <- rbind(
    "Vrai paramètre" = alpha,
    "Paramètre trouvé" = alpha_trouve$alpha,
    "Temps d'optimisation" = alpha_trouve$temps
))
xtable(tbl, digits=4)

# Optimisation à partir du tau théorique
# tau_theorique <- tau_kendall_theorique(F_N, F_X,
#                                        para,
#                                        Copule,
#                                        alpha,
#                                        n_max, x_max)
# 
# alpha_trouve <- inversion_tau_kendall(F_N, F_X,
#                              para,
#                              Copule,
#                              bornes,
#                              tau_theorique,
#                              n_max, x_max)
# (tbl <- rbind(
#     "Vrai paramètre" = alpha,
#     "Paramètre trouvé" = alpha_trouve$alpha,
#     "Temps d'optimisation" = alpha_trouve$temps
# ))
# xtable(tbl, digits=8)

# Scénario 2 : AMH(0.5), N suit Pois(2), X suit Exp(1/100) -----------

# Les paramètres
alpha <- 0.5
Copule <- amhCopula
F_N <- ppois
F_X <- pexp
para <- list(N=list(lambda=2), X=list(rate=1/100))
n_max <- qpois(0.99999, para$N$lambda)
x_max <- qexp(0.99999, para$X$rate)

# Calcul du alpha à l'aide du tau de Kendal empirique
tau_empirique <- numeric(nsim <- 1e+3)
set.seed(20190618)
for (i in 1:nsim) {
    UU <- rCopula(1e+3, Copule(alpha))
    X <- qpois(UU[,1], para$N$lambda)
    Y <- qexp(UU[,2], para$X$rate)
    tau_empirique[i] <- tau_kendall_empirique(X, Y)
}
(var_tau <- var(tau_empirique))
(mean_tau <- mean(tau_empirique))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))
tau_kendall_theorique(F_N, F_X, para, Copule, alpha, n_max, x_max)

# Analyse graphique pour connaître les bornes d'optimisation
Tau <- sapply((-9:9)/10, function(a)
    tau_kendall_theorique(F_N, F_X, para, Copule, a, n_max, x_max))

plot((-9:9)/10, Tau,
     ylab="tau de kendall",
     xlab="alpha",
     type="l",
)
axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
axis(1, at=(-9:9)/10, tck = 1, lty = 2, col = "grey",) # L'axe des abscisses
abline(a=tau_empirique, b=0, col="green")

# Bornes trouvées graphiquement.
bornes <- c(0.45, 0.55)

# Optimisation à partir du tau empirique
alpha_trouve <- inversion_tau_kendall(F_N, F_X,
                             para,
                             Copule,
                             bornes,
                             mean_tau,
                             n_max, x_max)
(tbl <- rbind(
    "Vrai paramètre" = alpha,
    "Paramètre trouvé" = alpha_trouve$alpha,
    "Temps d'optimisation" = alpha_trouve$temps
))
xtable(tbl, digits=4)

# # Optimisation à partir du tau théorique
# tau_theorique <- tau_kendall_theorique(F_N, F_X,
#                                        para,
#                                        Copule,
#                                        alpha,
#                                        n_max, x_max)
# 
# alpha_trouve <- inversion_tau_kendall(F_N, F_X,
#                              para,
#                              Copule,
#                              bornes,
#                              tau_theorique,
#                              n_max, x_max)
# (tbl <- rbind(
#     "Vrai paramètre" = alpha,
#     "Paramètre trouvé" = alpha_trouve$alpha,
#     "Temps d'optimisation" = alpha_trouve$temps
# ))
# xtable(tbl, digits=8)

# Scénario 3 : gumbel(0.5), N suit Pois(2), X suit Exp(1/100) -----------

# Les paramètres
alpha <- 5
Copule <- gumbelCopula
F_N <- ppois
F_X <- function(q, rate) 1-exp(-rate * q)
para <- list(N=list(lambda=2), X=list(rate=1/100))
n_max <- qpois(0.99999, para$N$lambda)
x_max <- qexp(0.99999, para$X$rate)

# Calcul du alpha à l'aide du tau de Kendal empirique
tau_empirique <- numeric(nsim <- 1e+3)
set.seed(20190618)
for (i in 1:nsim) {
    UU <- rCopula(1e+3, Copule(alpha))
    X <- qpois(UU[,1], para$N$lambda)
    Y <- qexp(UU[,2], para$X$rate)
    tau_empirique[i] <- tau_kendall_empirique(X, Y)
}
(var_tau <- var(tau_empirique))
(mean_tau <- mean(tau_empirique))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))
tau_kendall_theorique(F_N, F_X, para, Copule, alpha, n_max, x_max)

# Analyse graphique pour connaître les bornes d'optimisation
Tau <- sapply(1:10, function(a)
    tau_kendall_theorique(F_N, F_X, para, Copule, a, n_max, x_max))

plot(1:10, Tau,
     ylab="tau de kendall",
     xlab="alpha",
     type="l",
     xlim=c(1,10),
)
axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
axis(1, at= 1:10, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
abline(a=tau_empirique, b=0, col="green")

# Bornes trouvées graphiquement.
bornes <- c(4.5, 5.5)

# Optimisation à partir du tau empirique
alpha_trouve <- inversion_tau_kendall(F_N, F_X,
                             para,
                             Copule,
                             bornes,
                             mean_tau,
                             n_max, x_max)
(tbl <- rbind(
    "Vrai paramètre" = alpha,
    "Paramètre trouvé" = alpha_trouve$alpha,
    "Temps d'optimisation" = alpha_trouve$temps
))
xtable(tbl, digits=4)

# # Optimisation à partir du tau théorique
# tau_theorique <- tau_kendall_theorique(F_N, F_X,
#                                        para,
#                                        Copule,
#                                        alpha,
#                                        n_max, x_max)
# 
# alpha_trouve <- inversion_tau_kendall(F_N, F_X,
#                              para,
#                              Copule,
#                              bornes,
#                              tau_theorique,
#                              n_max, x_max)
# (tbl <- rbind(
#     "Vrai paramètre" = alpha,
#     "Paramètre trouvé" = alpha_trouve$alpha,
#     "Temps d'optimisation" = alpha_trouve$temps
# ))
# xtable(tbl, digits=8)


# Scénario 4 : Joe(0.5), N suit Pois(2), X suit Exp(1/100) -----------

# Les paramètres
alpha <- 5
Copule <- joeCopula
F_N <- ppois
F_X <- function(q, rate) 1-exp(-rate * q)
para <- list(N=list(lambda=2), X=list(rate=1/100))
n_max <- qpois(0.99999, para$N$lambda)
x_max <- qexp(0.99999, para$X$rate)

# Calcul du alpha à l'aide du tau de Kendal empirique
tau_empirique <- numeric(nsim <- 1e+3)
set.seed(20190618)
for (i in 1:nsim) {
    UU <- rCopula(1e+3, Copule(alpha))
    X <- qpois(UU[,1], para$N$lambda)
    Y <- qexp(UU[,2], para$X$rate)
    tau_empirique[i] <- tau_kendall_empirique(X, Y)
}
(var_tau <- var(tau_empirique))
(mean_tau <- mean(tau_empirique))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))
tau_kendall_theorique(F_N, F_X, para, Copule, alpha, n_max, x_max)

# Analyse graphique pour connaître les bornes d'optimisation
Tau <- sapply(1:10, function(a)
    tau_kendall_theorique(F_N, F_X, para, Copule, a, n_max, x_max))

plot(1:10, Tau,
     ylab="tau de kendall",
     xlab="alpha",
     type="l",
     xlim=c(1,10),
)
axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
axis(1, at= 1:10, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
abline(a=tau_empirique, b=0, col="green")

# Bornes trouvées graphiquement.
bornes <- c(4.5, 6.5)

# Optimisation à partir du tau empirique
alpha_trouve <- inversion_tau_kendall(F_N, F_X,
                                      para,
                                      Copule,
                                      bornes,
                                      mean_tau,
                                      n_max, x_max)
(tbl <- rbind(
    "Vrai paramètre" = alpha,
    "Paramètre trouvé" = alpha_trouve$alpha,
    "Temps d'optimisation" = alpha_trouve$temps
))
xtable(tbl, digits=4)

# # Optimisation à partir du tau théorique
# tau_theorique <- tau_kendall_theorique(F_N, F_X,
#                                        para,
#                                        Copule,
#                                        alpha,
#                                        n_max, x_max)
# 
# alpha_trouve <- inversion_tau_kendall(F_N, F_X,
#                                       para,
#                                       Copule,
#                                       bornes,
#                                       tau_theorique,
#                                       n_max, x_max)
# (tbl <- rbind(
#     "Vrai paramètre" = alpha,
#     "Paramètre trouvé" = alpha_trouve$alpha,
#     "Temps d'optimisation" = alpha_trouve$temps
# ))
# xtable(tbl, digits=8)


#---------------------------------------------------------------------------------
# Scénario 5 : Modèle collectif du risque - Logarithmique-Gamma ------------------

{
alpha0 <- 5
alpha1 <- 1/4
Copule0 <- frankCopula
Copule1 <- claytonCopula
F_N <- pbinom
F_X <- pexp
para <- list(N=list(size=7, prob=0.4), X=list(rate=1/100))
n_max <- qbinom(0.99999, para$N$size, para$N$prob)
x_max <- qexp(0.99999, para$X$rate)
nsim <- 1e+4
}# Les paramètres

{
set.seed(20190618)
DATA_train <- rCompCop(nsim, LOG(1-exp(-alpha0), 1, list(GAMMA(alpha1, 1:n_max, NULL))))
for (i in 1:nsim){
    DATA_train[i, 1] <- N <- qbinom(DATA_train[i,1], para$N$size, para$N$prob)
    if (N==0){ 
        DATA_train[i,-1] <- rep(NaN, n_max)
    }else{
        DATA_train[i,-1] <- c(qexp(DATA_train[i, 2:(N + 1)], para$X$rate), rep(NaN, n_max - N))}
}
colnames(DATA_train) <- c("N", sapply(1:n_max, function(i) paste0("X", i)))

# head(DATA_train) ; tail(DATA_train)
# summary(DATA_train)
    
}# Simulation d'une base de données

{
    lst_tau <- list()
    couples_NX <- list()
    for (i in 1:n_max) {
        couples_NX <- append(couples_NX,
                             list(DATA_train[which(DATA_train[, 1] >= i), c(1, i+1)]))
        
        if (nrow(couples_NX[[i]]) > 1e+3){
            lst_tau <- append(lst_tau, list(bootstrap_tau(couples_NX[[i]])))
        }else{
            if (nrow(couples_NX[[i]]) > 1e+2){
                lst_tau <- append(lst_tau,
                          tau_kendall_empirique(couples_NX[[i]][, 1],
                                                couples_NX[[i]][, 2]))
            }else break
        }
    }
} # Calcul des taus de Kendall entre N et chacun des X_i en utilisant du 
# rééchantillonage si le nombre de données est significatif.

tau <- unlist(lst_tau[1:1])
(var_tau <- var(tau))
(mean_tau <- mean(tau))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))

tau_kendall_theorique(F_N, F_X, para, Copule0, alpha0, n_max, x_max)

{
Tau_graph <- sapply(c(-5:-1, 1:10), function(a)
    tau_kendall_theorique(F_N, F_X, para, Copule0, a, n_max, x_max))

plot(c(-5:-1, 1:10), Tau_graph,
     ylab="tau de kendall",
     xlab="alpha",
     type="l"
)
axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
axis(1,at=-5:10, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
abline(a=mean_tau, b=0, col="green")
}# Analyse graphique

# Bornes trouvées graphiquement.
bornes <- c(4, 5.5)

alpha_trouve <- inversion_tau_kendall(F_N, F_X,
                                      para,
                                      Copule0,
                                      bornes,
                                      mean_tau,
                                      n_max, x_max)
(tbl <- rbind(
    "Vrai paramètre" = alpha0,
    "Paramètre trouvé" = alpha_trouve$alpha,
    "Temps d'optimisation" = alpha_trouve$temps
))


# Scénario 6 : Modèle collectif du risque - Géométrique-Gamma ------------------

{
    alpha0 <- 0.5
    alpha1 <- 1/4
    Copule0 <- amhCopula
    Copule1 <- claytonCopula
    F_N <- pbinom
    F_X <- pexp
    para <- list(N=list(size=7, prob=0.4), X=list(rate=1/100))
    n_max <- qbinom(0.99999, para$N$size, para$N$prob)
    x_max <- qexp(0.99999, para$X$rate)
    nsim <- 1e+4
}# Les paramètres

{
    set.seed(20190618)
    DATA_train <- rCompCop(nsim, GEO(1-alpha0, 1, list(GAMMA(alpha1, 1:n_max, NULL))))
    for (i in 1:nsim){
        DATA_train[i, 1] <- N <- qbinom(DATA_train[i,1], para$N$size, para$N$prob)
        if (N==0){ 
            DATA_train[i,-1] <- rep(NaN, n_max)
        }else{
            DATA_train[i,-1] <- c(qexp(DATA_train[i, 2:(N + 1)], para$X$rate), rep(NaN, n_max - N))}
    }
    colnames(DATA_train) <- c("N", sapply(1:n_max, function(i) paste0("X", i)))
    
    # head(DATA_train) ; tail(DATA_train)
    # summary(DATA_train)
    
}# Simulation d'une base de données

{
    lst_tau <- list()
    couples_NX <- list()
    for (i in 1:n_max) {
        couples_NX <- append(couples_NX,
                             list(DATA_train[which(DATA_train[, 1] >= i), c(1, i+1)]))
        
        if (nrow(couples_NX[[i]]) > 1e+3){
            lst_tau <- append(lst_tau, list(bootstrap_tau(couples_NX[[i]])))
        }else{
            if (nrow(couples_NX[[i]]) > 1e+2){
                lst_tau <- append(lst_tau,
                                  tau_kendall_empirique(couples_NX[[i]][, 1],
                                                        couples_NX[[i]][, 2]))
            }else break
        }
    }
} # Calcul des taus de Kendall entre N et chacun des X_i en utilisant du 
# rééchantillonage si le nombre de données est significatif.

tau <- unlist(lst_tau[1:1])
(var_tau <- var(tau))
(mean_tau <- mean(tau))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))

tau_kendall_theorique(F_N, F_X, para, Copule0, alpha0, n_max, x_max)

{
    Tau <- sapply((-9:9)/10, function(a)
        tau_kendall_theorique(F_N, F_X, para, Copule0, a, n_max, x_max))
    
    plot((-9:9)/10, Tau,
         ylab="tau de kendall",
         xlab="alpha",
         type="l",
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=(-9:9)/10, tck = 1, lty = 2, col = "grey",) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
}# Analyse graphique

# Bornes trouvées graphiquement.
bornes <- c(0.4, 0.55)

alpha_trouve <- inversion_tau_kendall(F_N, F_X,
                                      para,
                                      Copule0,
                                      bornes,
                                      mean_tau,
                                      n_max, x_max)
(tbl <- rbind(
    "Vrai paramètre" = alpha0,
    "Paramètre trouvé" = alpha_trouve$alpha,
    "Temps d'optimisation" = alpha_trouve$temps
))


#---------------------------------------------------------------------------------
#======================================= Rho de Spearman ==========================
rho_spearman <- function(X,Y){
    n <- length(X)
    R_i <- rank(X, ties.method = "average")
    S_i <- rank(Y, ties.method = "average")
    rho <- 12 / (n * (n + 1) * (n - 1)) * sum(R_i*S_i) - 3 * (n + 1) / (n - 1)
    return(rho)
}

rho_spearman(X,Y)
cor(X,Y,method = "spearman")

#####