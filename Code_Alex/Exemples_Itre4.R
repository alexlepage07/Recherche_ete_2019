# Article de Itre : Hierarchical Archimedian copulas...
#============================Exemple 1 ========================
source("Code_simul_copule.R")
nsim <- 1e+4
# 1)
colMeans(actrisk.rcompcop(nsim,"geo", 0.5, 2, "geo", 0.1))
# 2)
colMeans(actrisk.rcompcop(nsim,"geo", 0.5, 2, "log", 0.1))
# 3)
colMeans(actrisk.rcompcop(nsim,"geo", 0.5, 2, "gamma", 3))

#============================Exemple 2 ========================
# 1)
colMeans(actrisk.rcompcop(nsim,"log", 0.5, 2, "geo", 0.1))
# 2)
colMeans(actrisk.rcompcop(nsim,"log", 0.5, 2, "log", 0.1))
# 3)
colMeans(actrisk.rcompcop(nsim,"log", 0.5, 2, "gamma", 3))

# Avec une gamma comme mère
# 1)
colMeans(actrisk.rcompcop(nsim,"gamma", 2, 2, "geo", 0.1))
# 2)
colMeans(actrisk.rcompcop(nsim,"gamma", 2, 2, "log", 0.1))
# 3)
colMeans(actrisk.rcompcop(nsim,"gamma", 2, 2, "gamma", 3))

# Simulation de S
# U <- actrisk.rcompcop(nsim,"geo", 0.5, 2, "gamma", 3)
# N <- qpois(U[,1],1)
# S <- c()
# for (i in 1:nrow(U)){
#     if (N[i] == 0)
#         S <- append(S, 0)
#     else
#         S <- append(S, sum(qgamma(U[i, 2:N[i]+1], 1, 0.001)))
# }
# mean(S)

#============================Exemple 3 ========================
colMeans(actrisk.rcompcop(nsim,
                          "geo", 0.5,
                          c(2, 2),
                          c("log", "gamma"), c(0.1, 3)))

# =========================== Exemple 4 =======================
colMeans(actrisk.rcompcop(nsim,
                          "geo", 0.1,
                          c(2, 2),
                          c("gamma", "gamma"), c(0.04, 0.2)))
# =========================== Exemple 5 =======================
# Pour une copule hiérarchique à deux niveaux
p1 <- 0.5; p2 <- 0.1; alpha <- 1/30
n <- 2

L_theta_01i <- function(t) {
    LST_geom(p1, -log(
        LST_geom(p2, -log(
            LST_gamma(alpha, t)
        ))))
}

L_theta_0i <- function(t) {
    LST_geom(p1, -log(
        LST_gamma(alpha, t)
    ))
}

exemple5 <- function(p, alpha, n){
    p_1 <- p[1]; p_2 <- p[2]
    M_0 <- rgeom(1, p1) + 1
    M_01 <- rnbinom(1, M_0, p2) + M_0
    U <- c()
    for (i in 1:2) {
        theta_01i <- rgamma(1, M_01 * alpha, 1)
        U_i <- L_theta_01i(rexp(n) / theta_01i)
        U <- append(U, U_i)
    }
    theta_02 <- rgamma(1, M_0 * alpha, 1)
    U <- append(U, L_theta_0i(rexp(n) / theta_02))
    return(U)
}


colMeans(actrisk.simul(nsim, exemple5, c(p1,p2), alpha, n))

# =========================== Exemple 6 =======================
source("Code_simul_copule.R")
nsim <- 1e+6
U <- actrisk.rcompcop(nsim,
                      "log", 0.5, 
                      c(40, 40), 
                      c("geo", "geo"), 
                      c(0.8, 0.9)
                      )[,-1]

# library(nCopula)
# U <- rCompCop(nsim, LOG(0.5, NULL, list(GEO(0.8, 1:40, NULL),
#                                         GEO(0.9, 1:40, NULL))))


simul_S <- function(U){
    U_1 <- U[,1:(ncol(U)/2)]
    U_2 <- U[,-(1:(ncol(U)/2))]

    simul_X <- function(U_si, s){
        sapply(1:length(U_si), function(i)
            qbinom(U_si[i], 10, 0.05 * s + 0.005 * i))
    }

    S <- numeric(nrow(U))
    for(i in 1:nrow(U)){
        X_1 <- simul_X(U_1[i,], 1)
        X_2 <- simul_X(U_2[i,], 2)
        S[i] <- sum(X_1 + X_2)
    }
    return(S)
}
S <- simul_S(U)

mean(S)
var(S)

S <- sort(S)
summary(S)

VaR_S <- function(k) S[ceiling(k*length(S))]
TVaR_S <- function(k) mean(S[S>VaR_S(k)])


kappa <- c(0.9, 0.99, 0.999, 0.9999)
sapply(kappa, VaR_S)
sapply(kappa, TVaR_S)
