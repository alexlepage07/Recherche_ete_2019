library(pracma)


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
    N <- qpois(UU[,1], para$N$lambda)
    X <- qexp(UU[,2], para$X$rate)
    tau_empirique[i] <- tau_kendall_empirique(N, X)
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
DATA_train <- rCompCop(nsim, LOG(1-exp(-alpha0), 1, list(GAMMA(1/alpha1, 1:n_max, NULL))))
for (i in 1:nsim){
    DATA_train[i, 1] <- N <- qbinom(DATA_train[i,1], para$N$size, para$N$prob)
    if (N==0){ 
        DATA_train[i,-1] <- rep(NaN, n_max)
    }else{
        DATA_train[i,-1] <- c(qexp(DATA_train[i, 2:(N + 1)], para$X$rate), rep(NaN, n_max - N))}
}
colnames(DATA_train) <- c("N", sapply(1:n_max, function(i) paste0("X", i)))

}# Simulation d'une base de données


set.seed(20190618)
tau_0 <- tau_NX(DATA_train, nsim=10, silent=F)

var_tau <- var(tau_0)
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
bornes <- c(6, 7)
    
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

set.seed(20190618)
tau_1 <- tau_XX(DATA_train, nsim=30, silent=F)


(var_tau <- var(tau_1, na.rm = T))
(mean_tau <- mean(tau_1, na.rm = T))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))
tau_kendall_theorique_continues(struct_copule,
                                alpha0_n, alpha1)

{
    Tau_graph <- sapply(c(0.001,1:5), function(a)
        tau_kendall_theorique_continues(struct_copule,
                                        alpha0_n, alpha1=a))
    
    plot(0:5, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=0:5, tck=1, lty = 2, col = "grey", labels = F) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
    
    plot(0.001, Tau_graph[1],
         ylab="tau de kendall",
         xlab="alpha",
         type="p"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=(0:2)/10, tck=1, lty = 2, col = "grey", labels = F) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
}# Analyse graphique

bornes <- c(1e-12, 0.001)

temps_1 <- system.time(
    alpha1_n <- inversion_tau_kendall_XX(struct_copule, alpha0_n, mean_tau, bornes)
)

(tbl_1 <- rbind(
    "Vrai paramètre" = alpha1,
    "Paramètre trouvé" = alpha1_n,
    "Temps d'optimisation" = temps_1[[3]]
))

tbl_tot <- cbind(tbl_0, tbl_1)
colnames(tbl_tot) <- c("alpha0", "alpha1")
tbl_tot

# Scénario 6 : Modèle collectif du risque - Logarithmique-Logarithmique ----------

{
    alpha0 <- 6
    alpha1 <- 4
    Copule0 <- frankCopula
    
    struct_copule <- function(alpha0, alpha1){
        LOG(1-exp(-alpha0), NULL, list(LOG(1-exp(-alpha1), 1:2, NULL)))
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
    DATA_train <- rCompCop(nsim, LOG(1-exp(-alpha0), 1, list(LOG(1-exp(-alpha1), 1:n_max, NULL))))
    for (i in 1:nsim){
        DATA_train[i, 1] <- N <- qbinom(DATA_train[i,1], para$N$size, para$N$prob)
        if (N==0){ 
            DATA_train[i,-1] <- rep(NaN, n_max)
        }else{
            DATA_train[i,-1] <- c(qexp(DATA_train[i, 2:(N + 1)], para$X$rate), rep(NaN, n_max - N))}
    }
    colnames(DATA_train) <- c("N", sapply(1:n_max, function(i) paste0("X", i)))
    
}# Simulation d'une base de données


set.seed(20190618)
tau_0 <- tau_NX(DATA_train, nsim=10, silent=F)

var_tau <- var(tau_0)
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
    bornes <- c(6, 7)
    
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

set.seed(20190618)
tau_1 <- tau_XX(DATA_train, nsim=30, silent=F)


(var_tau <- var(tau_1, na.rm = T))
(mean_tau <- mean(tau_1, na.rm = T))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))
tau_kendall_theorique_continues(struct_copule,
                                alpha0_n, alpha1)

{
    Tau_graph <- sapply(c(0.001,1:5), function(a)
        tau_kendall_theorique_continues(struct_copule,
                                        alpha0_n, alpha1=a))
    
    plot(0:5, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=0:5, tck=1, lty = 2, col = "grey", labels = F) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
    
    plot(0.001, Tau_graph[1],
         ylab="tau de kendall",
         xlab="alpha",
         type="p"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=(0:2)/10, tck=1, lty = 2, col = "grey", labels = F) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
}# Analyse graphique

bornes <- c(1e-16, 1)

temps_1 <- system.time(
    alpha1_n <- inversion_tau_kendall_XX(struct_copule, alpha0_n, mean_tau, bornes)
)

(tbl_1 <- rbind(
    "Vrai paramètre" = alpha1,
    "Paramètre trouvé" = alpha1_n,
    "Temps d'optimisation" = temps_1[[3]]
))

tbl_tot <- cbind(tbl_0, tbl_1)
colnames(tbl_tot) <- c("alpha0", "alpha1")
tbl_tot

# Scénario 7 : Modèle collectif du risque - Logarithmique-géométrique ------------

{
    alpha0 <- 6
    alpha1 <- 0.5
    Copule0 <- frankCopula
    
    struct_copule <- function(alpha0, alpha1){
        LOG(1-exp(-alpha0), NULL, list(GEO(1-alpha1, 1:2, NULL)))
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
    DATA_train <- rCompCop(nsim, LOG(1-exp(-alpha0), 1, list(GEO(1-alpha1, 1:n_max, NULL))))
    for (i in 1:nsim){
        DATA_train[i, 1] <- N <- qbinom(DATA_train[i,1], para$N$size, para$N$prob)
        if (N==0){ 
            DATA_train[i,-1] <- rep(NaN, n_max)
        }else{
            DATA_train[i,-1] <- c(qexp(DATA_train[i, 2:(N + 1)], para$X$rate), rep(NaN, n_max - N))}
    }
    colnames(DATA_train) <- c("N", sapply(1:n_max, function(i) paste0("X", i)))
    
}# Simulation d'une base de données


set.seed(20190618)
tau_0 <- tau_NX(DATA_train, nsim=10, silent=F)

var_tau <- var(tau_0)
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
    bornes <- c(6, 7)
    
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

set.seed(20190618)
tau_1 <- tau_XX(DATA_train, nsim=30, silent=F)


(var_tau <- var(tau_1, na.rm = T))
(mean_tau <- mean(tau_1, na.rm = T))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))
tau_kendall_theorique_continues(struct_copule,
                                alpha0_n, alpha1)

{
    Tau_graph <- sapply((1:9)/10, function(a)
        tau_kendall_theorique_continues(struct_copule,
                                        alpha0_n, alpha1=a))
    
    plot((1:9)/10, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=(1:9)/10, tck=1, lty = 2, col = "grey", labels = F) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
    
}# Analyse graphique

bornes <- c(1e-16, 0.00001)

temps_1 <- system.time(
    alpha1_n <- inversion_tau_kendall_XX(struct_copule, alpha0_n, mean_tau, bornes)
)

(tbl_1 <- rbind(
    "Vrai paramètre" = alpha1,
    "Paramètre trouvé" = alpha1_n,
    "Temps d'optimisation" = temps_1[[3]]
))

tbl_tot <- cbind(tbl_0, tbl_1)
colnames(tbl_tot) <- c("alpha0", "alpha1")
tbl_tot

#---------------------------------------------------------------------------------
# Scénario 8 : Modèle collectif du risque - Géométrique-Gamma ------------------

{
    alpha0 <- 0.7
    alpha1 <- 4
    Copule0 <- amhCopula
    struct_copule <- function(alpha0, alpha1){
        GEO(1-alpha0, NULL, list(GAMMA(1/alpha1, 1:2, NULL)))
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
    DATA_train <- rCompCop(nsim, GEO(1-alpha0, 1, list(GAMMA(1/alpha1, 1:n_max, NULL))))
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


set.seed(20190618)
tau_0 <- tau_NX(DATA_train, nsim=10, silent=F)

var_tau <- var(tau_0)
(mean_tau <- mean(tau_0))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))

tau_kendall_theorique(F_N, F_X, para, Copule0, alpha0, n_max, x_max)

{
    Tau_graph <- sapply((-9:9)/10, function(a)
        tau_kendall_theorique(F_N, F_X, para, Copule0, a, n_max, x_max))
    
    plot((-9:9)/10, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=(-9:9)/10, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
}# Analyse graphique

{
    bornes <- c(0.8, 0.9)
    
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


set.seed(20190618)
tau_1 <- tau_XX(DATA_train, nsim=10, silent=T)

(var_tau <- var(tau_1, na.rm = T))
(mean_tau <- mean(tau_1, na.rm = T))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))
tau_kendall_theorique_continues(struct_copule,
                                alpha0_n, alpha1)

{
    Tau_graph <- sapply(c(2:6), function(a)
        tau_kendall_theorique_continues(struct_copule,
                                        alpha0_n, alpha1=a))
    
    plot(2:6, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=2:6, tck=1, lty = 2, col = "grey", labels = F) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
    
}# Analyse graphique

bornes <- c(4.5, 5)

temps_1 <- system.time(
    alpha1_n <- inversion_tau_kendall_XX(struct_copule, alpha0_n, mean_tau, bornes)
)

(tbl_1 <- rbind(
    "Vrai paramètre" = alpha1,
    "Paramètre trouvé" = alpha1_n,
    "Temps d'optimisation" = temps_1[[3]]
))

tbl_tot <- cbind(tbl_0, tbl_1)
colnames(tbl_tot) <- c("alpha0", "alpha1")
tbl_tot


# Scénario 9 : Modèle collectif du risque - Géométrique-Géométrique ------------------

{
    alpha0 <- 0.7
    alpha1 <- 0.8
    Copule0 <- amhCopula
    struct_copule <- function(alpha0, alpha1){
        GEO(1-alpha0, NULL, list(GEO(1-alpha1, 1:2, NULL)))
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
    DATA_train <- rCompCop(nsim, GEO(1-alpha0, 1, list(GEO(1-alpha1, 1:n_max, NULL))))
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


set.seed(20190618)
tau_0 <- tau_NX(DATA_train, nsim=10, silent=F)

var_tau <- var(tau_0)
(mean_tau <- mean(tau_0))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))

tau_kendall_theorique(F_N, F_X, para, Copule0, alpha0, n_max, x_max)

{
    Tau_graph <- sapply((1:9)/10, function(a)
        tau_kendall_theorique(F_N, F_X, para, Copule0, a, n_max, x_max))
    
    plot((1:9)/10, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=(1:9)/10, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
}# Analyse graphique

{
    bornes <- c(0.8, 0.9)
    
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


set.seed(20190618)
tau_1 <- tau_XX(DATA_train, nsim=10, silent=T)

(var_tau <- var(tau_1, na.rm = T))
(mean_tau <- mean(tau_1, na.rm = T))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))
tau_kendall_theorique_continues(struct_copule,
                                alpha0_n, alpha1)

{
    Tau_graph <- sapply((1:9)/10, function(a)
        tau_kendall_theorique_continues(struct_copule,
                                        alpha0_n, alpha1=a))
    
    plot((1:9)/10, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=(1:9)/10, tck=1, lty = 2, col = "grey", labels = F) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
    
}# Analyse graphique

bornes <- c(1e-16, 1e-1)

temps_1 <- system.time(
    alpha1_n <- inversion_tau_kendall_XX(struct_copule, alpha0_n, mean_tau, bornes)
)

(tbl_1 <- rbind(
    "Vrai paramètre" = alpha1,
    "Paramètre trouvé" = alpha1_n,
    "Temps d'optimisation" = temps_1[[3]]
))

tbl_tot <- cbind(tbl_0, tbl_1)
colnames(tbl_tot) <- c("alpha0", "alpha1")
tbl_tot

# Scénario 10 : Modèle collectif du risque - Géométrique-Logarithmique ------------

{
    alpha0 <- 0.7
    alpha1 <- 4
    Copule0 <- amhCopula
    struct_copule <- function(alpha0, alpha1){
        GEO(1-alpha0, NULL, list(LOG(1-exp(-alpha1), 1:2, NULL)))
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
    DATA_train <- rCompCop(nsim, GEO(1-alpha0, 1, list(LOG(1-exp(-alpha1), 1:n_max, NULL))))
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


set.seed(20190618)
tau_0 <- tau_NX(DATA_train, nsim=10, silent=F)

var_tau <- var(tau_0)
(mean_tau <- mean(tau_0))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))

tau_kendall_theorique(F_N, F_X, para, Copule0, alpha0, n_max, x_max)

{
    Tau_graph <- sapply((1:9)/10, function(a)
        tau_kendall_theorique(F_N, F_X, para, Copule0, a, n_max, x_max))
    
    plot((1:9)/10, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=(1:9)/10, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
}# Analyse graphique

{
    bornes <- c(0.8, 0.9)
    
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


set.seed(20190618)
tau_1 <- tau_XX(DATA_train, nsim=10, silent=T)

(var_tau <- var(tau_1, na.rm = T))
(mean_tau <- mean(tau_1, na.rm = T))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))
tau_kendall_theorique_continues(struct_copule,
                                alpha0_n, alpha1)

{
    Tau_graph <- sapply((1:9), function(a)
        tau_kendall_theorique_continues(struct_copule,
                                        alpha0_n, alpha1=a))
    
    plot((1:9), Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=(1:9), tck=1, lty = 2, col = "grey", labels = F) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
    
}# Analyse graphique

bornes <- c(3.5, 4)

temps_1 <- system.time(
    alpha1_n <- inversion_tau_kendall_XX(struct_copule, alpha0_n, mean_tau, bornes)
)

(tbl_1 <- rbind(
    "Vrai paramètre" = alpha1,
    "Paramètre trouvé" = alpha1_n,
    "Temps d'optimisation" = temps_1[[3]]
))

tbl_tot <- cbind(tbl_0, tbl_1)
colnames(tbl_tot) <- c("alpha0", "alpha1")
tbl_tot

#---------------------------------------------------------------------------------
