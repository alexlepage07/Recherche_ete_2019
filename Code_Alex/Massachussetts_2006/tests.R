library(copula)


# avec z = 1 : Equicorrelation matrix 
Sig_rho2.z1 <- function(rho, k) {
    # Fonction qui construit la matrice de corrélation de rho2 avec z=1
    matrix(c(rep(c(1, rep(rho[2], k)), k-1),1), nrow = k, byrow = T)
}
Sig_rho2.z1.inv <- function(rho, k) {
    (diag(k) - rho[2] / (1 + (k-1) * rho[2]) * matrix(1, k, k)) / (1 - rho[2])
}
Sig_rho1.2.z1 <- function(rho, k) {
    # Fonction qui construit la matrice de corrélation de rho1 et rho2 avec z=1
    mat <- rbind(rep(rho[1], k), Sig_rho2.z1(rho, k))
    mat <- cbind(c(1, rep(rho[1], k)), mat)
    return(mat)
}
Sig_rho1.2.z1.inv <- function(rho, k) {
    denom <- 1 + (k - 1) * rho[2] - k * rho[1] ^ 2
    mat <- matrix((diag(k) - (rho[2] - rho[1] ^ 2) / denom * matrix(1, k, k)) / (1 - rho[2]), k)
    mat <- rbind(-rho[1] / denom * t(rep(1, k)) , mat)
    mat <- cbind(c(1 + rho[1] ^ 2 * k / denom,-rho[1] / denom * rep(1, k)), mat)
    return(mat)
}

# Avec z = 2 : Autoregressive correlation matrix
Sig_rho2.z2 <- function(rho, k) {
    # Fonction qui construit la matrice de corrélation de rho2 avec z=2
    outer(0:(k-1), 0:(k-1), function(i,j) rho[2]^abs(i-j))
}
Sig_rho2.z2.inv <- function(rho, k) {
    # Fonction qui construit la matrice de corrélation inverse de rho2 avec z=2
    solve(Sig_rho2.z2(rho, k))
}
Sig_rho1.2.z2 <- function(rho, k) {
    # Fonction qui construit la matrice de corrélation de rho1 et rho2 avec z=2
    mat <- rbind(rep(rho[1], k), Sig_rho2.z2(rho, k))
    mat <- cbind(c(1, rep(rho[1], k)), mat)
    return(mat)
}
Sig_rho1.2.z2.inv <- function(rho, k) {
    # Fonction qui construit la matrice de corrélation inverse de rho1 et rho2 avec z=2
    solve(Sig_rho1.2.z2(rho, k))
}

mu.norm <- function(U, rho, z=1) {
    # Fonction qui calcule le paramètre de moyenne de la copule gaussienne
    k <- length(U)
    Sig_rho.inv <- ifelse(z == 1, Sig_rho2.z1.inv, Sig_rho2.z2.inv)
    
    rho[1] * t(rep(1, k)) %*% Sig_rho.inv(rho, k) %*% (qnorm(U))
}
sig.norm <- function(U, rho, z=1) {
    # Fonction qui calcule le paramètre d'écart-type de la copule gaussienne
    k <- length(U)
    Sig_rho.inv <- ifelse(z == 1, Sig_rho2.z1.inv, Sig_rho2.z2.inv)
    
    sqrt(1 - (rho[1] * t(rep(1, k))) %*% Sig_rho.inv(rho, k) %*% (rho[1] * rep(1, k)))
}


parameter.settings <- matrix(
    # Paramètres initiaux présentés dans le tableau 1 de l'article.
    c(rep(-2.5,12), rep(0.5,12), rep(1, 6), rep(1.5, 6), 
      rep(8, 12), rep(-0.1, 12), rep(0.3, 12), rep(0.7, 12),
      rep(c(rep(-0.05, 2), rep(0.05, 2), rep(0.1, 2)), 2),
      rep(c(0.1,0.05),6)),
    nrow = 12
)

scenario <- 1
nsim <- 1e+6
rho <- parameter.settings[scenario, 8:9]
U <- rCopula(nsim, normalCopula(P2p(Sig_rho1.2.z1(rho, 2)), 3, dispstr = "un"))

# rho <- numeric(2)
# rho[1] <- cor(U[,1], U[,2], method = "spearman")
# rho[2] <- cor(U[,2], U[,3], method = "spearman")
# 2 * sin(pi * rho / 6)
# 
# rho1 <- optimize(function(para)
#     - sum(log(dCopula(U, normalCopula(P2p(Sig_rho1.2.z1(c(para, rho[2]), 2)), 3, dispstr = "un")
#     ))),
#     interval =  c(-1, sqrt(((2 - 1) * rho[2] + 1) / 2)))$minimum


N <- qpois(U[,1], 0.1)
Y <- qexp(U[,2:3], 0.01)

rho <- numeric(2)
rho1 <- sapply(1:10, function(i)
    cor(rank(N, ties.method = "random") / nsim,
        rank(Y[,1]) / nsim,
        method = "spearman"))
rho[1] <- mean(rho1)
rho[1] - 1.645 * sqrt(var(rho1)) 
rho[1] + 1.645 * sqrt(var(rho1))
rho[2] <- cor(Y[,1], Y[,2], method = "spearman")
rho <- 2 * sin(pi * rho / 6)


bounds <- matrix(c(-1, 1,
                   -1, 1), ncol = 2, byrow = T)
n_par <- nrow(bounds)
ui <- rbind( diag(n_par), -diag(n_par) )
ci <- c( bounds[,1], - bounds[,2] )
i <- as.vector(is.finite(bounds))
ui <- ui[i,]
ci <- ci[i]

rho <- constrOptim(rho,
                    function(para)
                        -sum(log(dCopula(cbind(ppois(N, 0.1), pexp(Y, 0.01)), normalCopula(P2p(Sig_rho1.2.z1(rho, 2)), 3, dispstr = "un")) )),
                    grad = NULL, ui=ui, ci=ci)$par
rho

rho1 <- optimize(function(para)
    - sum(log(dCopula(cbind(ppois(N, 0.1), pexp(Y, 0.01)),
                      normalCopula(P2p(Sig_rho1.2.z1(c(para, rho[2]), 2)), 3, dispstr = "un")
    ))),
    interval =  c(-1, sqrt(((2 - 1) * rho[2] + 1) / 2)))$minimum

