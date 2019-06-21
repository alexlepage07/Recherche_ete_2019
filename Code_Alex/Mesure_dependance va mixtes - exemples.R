source("Mesure_dependance va mixtes.R")

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
alpha0 <- 6
alpha1 <- 4
Copule0 <- frankCopula

struct_copule <- function(alpha0, alpha1){
    LOG(1-exp(-alpha0), NULL, list(GAMMA(1/alpha1, 1:2, NULL)))
}

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
    nb_couple <- sum(!is.na(DATA_train[,-1]))
    couples_NX <- matrix(ncol=2, nrow=nb_couple)
    k <- 0
    for (i in 1:nrow(DATA_train)) {
        if (DATA_train[i, 1] == 0)
            next
        for (j in 1:n_max){
            if (is.na(DATA_train[i, j + 1]))
                next
            else
                couples_NX[k <- k + 1,] <- DATA_train[i, c(1, j+1)]
        }
    }
    
    tau_0 <- bootstrap_tau_mixte(couples_NX)
    
} # Calcul des taus de Kendall entre N et chacun des X_i en utilisant du 
# rééchantillonage.
summary(couples_NX)

(var_tau <- var(tau_0))
(mean_tau <- mean(tau_0))
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

{
bornes <- c(5, 6.5)
    
temps_0 <- system.time(
    alpha0_n <- inversion_tau_kendall(F_N, F_X,
                                      para,
                                      Copule0,
                                      bornes,
                                      mean_tau,
                                      n_max, x_max)
)

(tbl_0 <- rbind(
    "Vrai paramètre" = alpha0,
    "Paramètre trouvé" = alpha0_n,
    "Temps d'optimisation" = temps_0[[3]]
))
} # Alpha0

{
    lst_tau <- list()
    for (i in 1:(n_max - 1)) {
        for (j in (i + 1):n_max) {
            lst_tau <- append(lst_tau,
                              list(bootstrap_tau_continues(DATA_train[, c(i+1, j+1)])))
        }
    }
}# Calcul des taus entre les X_i

tau_1 <- unlist(lst_tau[1:2])
(var_tau <- var(tau_1, na.rm = T))
(mean_tau <- mean(tau_1, na.rm = T))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))
tau_kendall_theorique_continues(struct_copule,
                                alpha0_n, alpha1)

{
    Tau_graph <- sapply(c(0.1,1:3), function(a)
        tau_kendall_theorique_continues(struct_copule,
                                        alpha0_n, alpha1=a))
    
    plot(0:3, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1,at=0:3, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
}# Analyse graphique

bornes <- c(7, 9)

temps_1 <- system.time(
    alpha1_n <- inversion_tau_kendall_XX(struct_copule, alpha0_n, mean_tau, bornes)
)

(tbl_1 <- rbind(
    "Vrai paramètre" = alpha1,
    "Paramètre trouvé" = alpha1_n,
    "Temps d'optimisation" = temps_1[[3]]
))


#####
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
            lst_tau <- append(lst_tau, list(bootstrap_tau_mixte(couples_NX[[i]])))
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
