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
F_X <- function(q, beta) 1-exp(-beta * q)
n_max <- qpois(0.99999, 2)
x_max <- qexp(0.99999, 1/100)
para <- list(N=list(lambda=2), X=list(beta=1/100))

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

depend_optim <- function(F_N, F_X,
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
set.seed(20190618)
UU <- rCopula(1e+3, Copule(alpha))
X <- qpois(UU[,1], para$N$lambda)
Y <- qexp(UU[,2], para$X$rate)
tau_empirique <- tau_kendall_empirique(X, Y)
cor(X,Y, method = "kendall")

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
alpha_trouve <- depend_optim(F_N, F_X,
                             para,
                             Copule,
                             bornes,
                             tau_empirique,
                             n_max, x_max)
(tbl <- rbind(
    "Vrai paramètre" = alpha,
    "paramètre trouvé" = alpha_trouve$alpha,
    "Temps d'optimisation" = alpha_trouve$temps
))
xtable(tbl, digits=4)

# Optimisation à partir du tau théorique
tau_theorique <- tau_kendall_theorique(F_N, F_X,
                                       para,
                                       Copule,
                                       alpha,
                                       n_max, x_max)

alpha_trouve <- depend_optim(F_N, F_X,
                             para,
                             Copule,
                             bornes,
                             tau_theorique,
                             n_max, x_max)
(tbl <- rbind(
    "Vrai paramètre" = alpha,
    "paramètre trouvé" = alpha_trouve$alpha,
    "Temps d'optimisation" = alpha_trouve$temps
))
xtable(tbl, digits=8)

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
set.seed(20190618)
UU <- rCopula(1e+3, Copule(alpha))
X <- qpois(UU[,1], para$N$lambda)
Y <- qexp(UU[,2], para$X$rate)
tau_empirique <- tau_kendall_empirique(X, Y)
rbind(
    "tau empirique" = tau_empirique,
    "Fonction cor" = cor(X,Y, method = "kendall"))

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
alpha_trouve <- depend_optim(F_N, F_X,
                             para,
                             Copule,
                             bornes,
                             tau_empirique,
                             n_max, x_max)
(tbl <- rbind(
    "Vrai paramètre" = alpha,
    "paramètre trouvé" = alpha_trouve$alpha,
    "Temps d'optimisation" = alpha_trouve$temps
))
xtable(tbl, digits=4)

# Optimisation à partir du tau théorique
tau_theorique <- tau_kendall_theorique(F_N, F_X,
                                       para,
                                       Copule,
                                       alpha,
                                       n_max, x_max)

alpha_trouve <- depend_optim(F_N, F_X,
                             para,
                             Copule,
                             bornes,
                             tau_theorique,
                             n_max, x_max)
(tbl <- rbind(
    "Vrai paramètre" = alpha,
    "paramètre trouvé" = alpha_trouve$alpha,
    "Temps d'optimisation" = alpha_trouve$temps
))
xtable(tbl, digits=8)

# Scénario 3 : gumbel(0.5), N suit Pois(2), X suit Exp(1/100) -----------

# Les paramètres
alpha <- 5
Copule <- gumbelCopula
F_N <- ppois
F_X <- function(q, beta) 1-exp(-beta * q)
para <- list(N=list(lambda=2), X=list(beta=1/100))
n_max <- qpois(0.99999, para$N$lambda)
x_max <- qexp(0.99999, para$X$beta)

# Calcul du alpha à l'aide du tau de Kendal empirique
set.seed(20190618)
UU <- rCopula(1e+3, Copule(alpha))
X <- qpois(UU[,1], para$N$lambda)
Y <- qexp(UU[,2], para$X$beta)
tau_empirique <- tau_kendall_empirique(X, Y)
rbind(
    "tau empirique" = tau_empirique,
    "Fonction cor" = cor(X,Y, method = "kendall"))

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
alpha_trouve <- depend_optim(F_N, F_X,
                             para,
                             Copule,
                             bornes,
                             tau_empirique,
                             n_max, x_max)
(tbl <- rbind(
    "Vrai paramètre" = alpha,
    "paramètre trouvé" = alpha_trouve$alpha,
    "Temps d'optimisation" = alpha_trouve$temps
))
xtable(tbl, digits=4)

# Optimisation à partir du tau théorique
tau_theorique <- tau_kendall_theorique(F_N, F_X,
                                       para,
                                       Copule,
                                       alpha,
                                       n_max, x_max)

alpha_trouve <- depend_optim(F_N, F_X,
                             para,
                             Copule,
                             bornes,
                             tau_theorique,
                             n_max, x_max)
(tbl <- rbind(
    "Vrai paramètre" = alpha,
    "paramètre trouvé" = alpha_trouve$alpha,
    "Temps d'optimisation" = alpha_trouve$temps
))
xtable(tbl, digits=8)

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