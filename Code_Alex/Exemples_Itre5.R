# Exemples de l'article Itre5
#========================= Expemle 9 =================================
#------------------------- Copule de Gumbel --------------------------
library(copula)
C_NV <- function(n, v, alpha, theta) {
    pCopula(c(ppois(n, 5), pgeom(v - 1, 1 - theta)),
            gumbelCopula(
                param = alpha,
                dim = 2,
                use.indepC = "FALSE"
            ))
}


f_NV <- function(n, v, alpha, theta){
    C_NV(n, v, alpha, theta) - 
        C_NV(n-1, v, alpha, theta) - C_NV(n, v-1, alpha, theta) +
        C_NV(n-1, v-1, alpha, theta)
}


k_v <- ceiling(uniroot(function(v) C_NV(1000, v, 4, 2/3) - 0.99999999,
                      c(0, 50))$root)

k_n <- ceiling(uniroot(function(n) C_NV(n, k_v, 1, 2/3) - 0.99999999,
                       c(0, 50))$root)


GEN_f_NV <- function(alpha, theta) {
    n <- 0:k_n
    v <- 1:k_v
    
    F_N <- ppois(n, 5)
    F_V <- pgeom(1 - v, 1 - theta)
    
    mat_ <- matrix(c(rep(n, each  = k_v),
                     rep(v, times = k_n + 1),
                     rep(0, times = (k_n + 1) * k_v)), (k_n + 1) * k_v, 3)
    
    for (i in 1:((k_n + 1) * k_v)) {
        mat_[i, 3] <- f_NV(mat_[i, 1], mat_[i, 2], alpha, theta)
    }
    return(mat_)
}
matrices_0.33 <- lapply(c(1, 1.5, 4), function(a) GEN_f_NV(a, 1/3))
matrices_0.66 <- lapply(c(1, 1.5, 4), function(a) GEN_f_NV(a, 2/3))


F_S <- function(s, matrice, theta) {
    if (s == 0)
        return(dpois(0, 5))
    sum(sapply(1:nrow(matrice), function(i)
        matrice[i, 3] * pgamma(s, prod(matrice[i, 1:2]), 1 / (1 - theta))))
}
F_S(1, GEN_f_NV(1.5, 1/3), 1/3)

E_S <- function(matrice, theta){
    sum(sapply(1:nrow(matrice), function(i)
        matrice[i, 3] * prod(matrice[i, 1:2]) * (1 - theta)))
}
unlist(lapply(matrices_0.33, function(mat) round(E_S(mat, 1/3), 2)))
unlist(lapply(matrices_0.66, function(mat) round(E_S(mat, 2/3), 2)))


V_S <- function(matrice, theta) {
    E_L2 <- sum(sapply(1:nrow(matrice), function(i)
        matrice[i, 3] * (prod(matrice[i, 1:2]))^2))
    E_L <- sum(sapply(1:nrow(matrice), function(i)
        matrice[i, 3] * (prod(matrice[i, 1:2]))))
    return((E_L + E_L2 - E_L^2) * (1 - theta) ^ 2)
}
unlist(lapply(matrices_0.33, function(mat) round(V_S(mat, 1/3), 2)))
unlist(lapply(matrices_0.66, function(mat) round(V_S(mat, 2/3), 2)))


kappa <- c(0.9, 0.99, 0.999)


VaR_S <- function(kappa, matrice, theta) {
    sapply(kappa, function(k) {
        uniroot(function(s)
            F_S(s, matrice, theta) - k, c(0, 100))$root
    })
}
lapply(matrices_0.33, function(mat) round(VaR_S(kappa, mat, 1/3), 2))
lapply(matrices_0.66, function(mat) round(VaR_S(kappa, mat, 2/3), 2))


TVaR_S <- function(kappa, matrice, theta) {
    rep <- c()
    beta <- 1 / (1 - theta)
    for (k in kappa){
        VaR_ <- VaR_S(k, matrice, theta)
        TVaR_ <- sum(sapply(1:nrow(matrice), function(i)
            matrice[i, 3] * prod(matrice[i, 1:2]) / beta *
                (1 - pgamma(VaR_, prod(matrice[i, 1:2]) + 1, beta))
            )) / (1 - k)
        rep <- append(rep, TVaR_)
    }
    return(rep)
}
lapply(matrices_0.33, function(mat) round(TVaR_S(kappa, mat, 1/3), 2))
lapply(matrices_0.66, function(mat) round(TVaR_S(kappa, mat, 2/3), 2))

#------------------------- Copule de Frank ---------------------------
C_NV <- function(n, v, alpha, theta) {
    pCopula(c(ppois(n, 5), pgeom(v - 1, 1 - theta)),
            frankCopula(
                param = alpha,
                dim = 2,
                use.indepC = "FALSE"
            ))
}


f_NV <- function(n, v, alpha, theta){
    C_NV(n, v, alpha, theta) - 
        C_NV(n-1, v, alpha, theta) - C_NV(n, v-1, alpha, theta) +
        C_NV(n-1, v-1, alpha, theta)
}


k_v <- ceiling(uniroot(function(v) C_NV(1000, v, 4, 2/3) - 0.99999999,
                       c(0, 50))$root)

k_n <- ceiling(uniroot(function(n) C_NV(n, k_v, 1, 2/3) - 0.99999999,
                       c(0, 50))$root)


GEN_f_NV <- function(alpha, theta) {
    n <- 0:k_n
    v <- 1:k_v
    
    F_N <- ppois(n, 5)
    F_V <- pgeom(1 - v, 1 - theta)
    
    mat_ <- matrix(c(rep(n, each  = k_v),
                     rep(v, times = k_n + 1),
                     rep(0, times = (k_n + 1) * k_v)), (k_n + 1) * k_v, 3)
    
    for (i in 1:((k_n + 1) * k_v)) {
        mat_[i, 3] <- f_NV(mat_[i, 1], mat_[i, 2], alpha, theta)
    }
    return(mat_)
}
matrices_0.33 <- lapply(c(-10, -5, 5, 10), function(a) GEN_f_NV(a, 1/3))


F_S <- function(s, matrice, theta) {
    if (s == 0)
        return(dpois(0, 5))
    sum(sapply(1:nrow(matrice), function(i)
        matrice[i, 3] * pgamma(s, prod(matrice[i, 1:2]), 1 / (1 - theta))))
}


E_S <- function(matrice, theta){
    sum(sapply(1:nrow(matrice), function(i)
        matrice[i, 3] * prod(matrice[i, 1:2]) * (1 - theta)))
}
unlist(lapply(matrices_0.33, function(mat) round(E_S(mat, 1/3), 2)))


V_S <- function(matrice, theta) {
    E_L2 <- sum(sapply(1:nrow(matrice), function(i)
        matrice[i, 3] * (prod(matrice[i, 1:2]))^2))
    E_L <- sum(sapply(1:nrow(matrice), function(i)
        matrice[i, 3] * (prod(matrice[i, 1:2]))))
    return((E_L + E_L2 - E_L^2) * (1 - theta) ^ 2)
}
unlist(lapply(matrices_0.33, function(mat) round(V_S(mat, 1/3), 2)))


kappa <- c(0.9, 0.99, 0.999)


VaR_S <- function(kappa, matrice, theta) {
    sapply(kappa, function(k) {
        uniroot(function(s)
            F_S(s, matrice, theta) - k, c(0, 100))$root
    })
}
lapply(matrices_0.33, function(mat) round(VaR_S(kappa, mat, 1/3), 2))


TVaR_S <- function(kappa, matrice, theta) {
    rep <- c()
    beta <- 1 / (1 - theta)
    for (k in kappa){
        VaR_ <- VaR_S(k, matrice, theta)
        TVaR_ <- sum(sapply(1:nrow(matrice), function(i)
            matrice[i, 3] * prod(matrice[i, 1:2]) / beta *
                (1 - pgamma(VaR_, prod(matrice[i, 1:2]) + 1, beta))
        )) / (1 - k)
        rep <- append(rep, TVaR_)
    }
    return(rep)
}
lapply(matrices_0.33, function(mat) round(TVaR_S(kappa, mat, 1/3), 2))

#========================= Expemle 14 ================================
library(actuar)


# Loi logistique
LST_logarithmic <- function(g, t) {
    log(1 - g * exp(-t)) / log(1 - g)
}

# Loi gamma
LST_gamma <- function(alpha, t) {
    (1 + t) ^ (-alpha)
}


VaR_S <- function(kappa, S) {
    sapply(kappa, function(k)
        sort(S)[floor(length(S) * k)])
}

TVaR_S <- function(kappa, S) {
    sapply(kappa, function(k)
        mean(S[S > VaR_S(k, S)]))
}

{
# library(nCopula)
#     simul_NXX <- function(alphas, nsim = 1e+6) {
#     n_max <- qpois(0.999999999, 2)
#     
#     structure_C <-
#         LOG(1 - exp(-alphas[1]), 1, list(GAMMA(1 / alphas[2], 1:n_max, NULL)))
#     
#     UU <- rCompCop(nsim, structure_C)
#     NXX <- cbind(qpois(UU[, 1], 2), qpareto(UU[,-1], 3, 100))
#     
#     S <- numeric(nrow(NXX))
#     for (i in 1:nrow(NXX)) {
#         n <- NXX[i, 1]
#         if (n == 0)
#             S[i] <- 0
#         else
#             S[i] <- sum(sapply(2:(n+1), function(j)
#                 NXX[i, j]))
#     }
#     
#     E_S <- mean(S)
#     V_S <- var(S)
#     kappa <- c(0.9, 0.99, 0.999)
#     VaR_S <- VaR_S(kappa, S)
#     TVaR_S <- TVaR_S(kappa, S)
#     
#     return(list("E_S"=E_S, "V_S"=V_S, "kappa"=kappa,
#                 "VaR_S"=VaR_S, "TVaR_S"=TVaR_S))
# }
} # Code avec le package à Simon-Pierre

simul_S <- function(alphas, nsim = 1e+6) {
    g <- 1 - exp(-alphas[1])
    S <- numeric(nsim)
    k_n <- 2
    uu <- runif(nsim * (k_n + 3))
    j <- 0
    for (i in 1:nsim)   {
        theta_0 <- qlogarithmic(uu[j <- j + 1], g)
        U_0 <- LST_logarithmic(g, qexp(uu[j <- j + 1]) / theta_0)
        N <- qpois(U_0, 2)
        if (N == 0){
            S[i] <- 0
        }else{
            theta_01 <- qgamma(uu[j <- j + 1], theta_0 / alphas[2], 1)
            j <- j + 1
            U_01 <- LST_logarithmic(g, -log(LST_gamma(
                    1 / alphas[2], 
                    qexp(uu[j:(j + N - 1)]) / theta_01
                )))
            j <- j + N - 1
            XX <- qpareto(U_01, 3, 100)
            S[i] <- sum(XX)
        }
    }
    
    E_S <- mean(S)
    V_S <- var(S)
    kappa <- c(0.9, 0.99, 0.999)
    VaR_S <- VaR_S(kappa, S)
    TVaR_S <- TVaR_S(kappa, S)
    
    return(list("E_S"=E_S, "V_S"=V_S, "kappa"=kappa,
                "VaR_S"=VaR_S, "TVaR_S"=TVaR_S))
}

alphas <- list(c(0.7, 5),
               c(2.3, 5),
               c(2.3, 25)
               )
lapply(alphas, simul_S, 1e+7)

#========================= Exemple 16 ================================
library(actuar)

LST_geom_inv <- function(p, u) {
    log(p * (1 / u - 1) + 1)
}

LST_geom <- function(p, t) {
    p / (exp(t) - (1 - p))
}

# LST_NegBinom_inv <- function(r, p, u) {
#     log(p * (u ^ (-1 / r) - 1) + 1)
# }

fgp_N <- function(t, f_n) {
    sapply(t, function(tt)
        sum(sapply(1:length(f_n), function(n)
            tt ^ (n-1) * f_n[n])))
}

# alpha_0 <- 0.1 ; alpha_1 <- 0.2 ; h <- 1
# theta_0 <- 1 ; theta_01 <- 1

algorithme_15 <- function(alpha_0, alpha_1, h=1, TOL=1e-10) {
    
    p_0 <- 1 - alpha_0
    p_1 <- 1 - alpha_1
    p_01 <- p_1 / p_0
    
    F_Xi_theta <- function(u, theta_01) {
        exp(-theta_01 * LST_geom_inv(p_1, u))
    }
    F_N_theta <- function(u, theta_0) {
        exp(-theta_0 * LST_geom_inv(p_0, u))
    }
    
    k_0 <- qztgeom(1 - TOL, p_0)
    k_01 <- qztnbinom(1 - TOL, alpha_0, p_01)
    k_x <- qgamma(1 - TOL, 2, 1 / 50)
    k_n <- qpois(1 - TOL, 2)
    
    F_X <- pgamma((0:(k_x / h)) * h, 2, 1 / 50)
    F_N <- ppois(0:k_n, 2)
    
    nbr_S <- length(F_X) + 1
    f_S <- rep(0, nbr_S)
    f_S_theta_0 <- rep(0, nbr_S)
    
    for (theta_0 in 1:k_0) {
        for (theta_01 in 1:k_01) {
            
            F_X_array <- F_Xi_theta(F_X, theta_01)
            f_X_array <- c(F_X_array, 1) - c(0, F_X_array)

            fft_f_X <- fft(f_X_array)
            
            F_N_array <- F_N_theta(F_N, theta_0)
            f_N_array <- c(F_N_array, 1) - c(0, F_N_array)

            fft_f_S <- fgp_N(fft_f_X, f_N_array)
            
            f_S_theta_0 <- f_S_theta_0 + 
                Re(fft(fft_f_S, inverse = T)) / nbr_S  * 
                dztnbinom(theta_01, alpha_0, p_01)
        }
        f_S <- f_S + f_S_theta_0 * dztgeom(theta_0, p_0)
    }
    return(f_S)
}


h <- 1

f_S <- algorithme_15(0.1,0.2, h)
sum(f_S)

nb_S <- length(F_S) - 1
S <- (0:nb_S) * h

F_S <- cumsum(f_S)

(E_S <- sum(1 - F_S))
sum(S * f_S)

(V_S <- sum((S - E_S) ^ 2 * f_S))

VaR_S <- function(kappa, F_S, h = 1) {
    VaR_k <- sapply(kappa, function(k) {
        ifelse(k <= F_S[1],
               S[1],
               uniroot(function(m)
                   F_S[m + 1] - k, c(0, nb_S))$root * h)
    })
    return(VaR_k)
}
TVaR_S <- function(kappa, f_S, h = 1) {
    TVaR_ <- c()
    F_S <- cumsum(f_S)
    for (k in kappa) {
        VaR_k <- VaR_S(k, F_S, h)
        index <- which(S > VaR_k)
        
        TVaR_k <-
            (sum(S[index] * f_S[index]) +
                 VaR_k * (F_S[VaR_k / h + 1] - k)
             ) / (1 - k)
        
        TVaR_ <- append(TVaR_, TVaR_k)
    }
    return(TVaR_)
}

kappa <- c(0.9, 0.99, 0.999)
VaR_S(kappa, F_S, h)
TVaR_S(kappa, f_S, h)

#========================= Exemple 19 ================================
library(combinat)


# Loi logistique
L_M <- function(t, alpha) {
    - log(1 - (1 - exp(-alpha)) * exp(-t)) / alpha
}
L_M.inv <- function(u, alpha) {
    - log((1 - exp(-alpha * u)) / (1 - exp(-alpha)))
}

# Loi gamma
L_B <- function(t, alpha) {
    (1 + t) ^ (-1 / alpha)
}
L_B.inv <- function(u, alpha) {
    u ^ (-alpha) - 1
}

f_J <- function(j){
    if(j==2)
        return(0.2281946)
    if(j==9)
        return(0.5429590)
    if(j==14)
        return(0.2288464)
    else
        return(0)
    }
F_J <- function(j){
    if(j>=14)
        return(1)
    if(j>=9)
        return(0.7711536)
    if(j>=2)
        return(0.2281946)
    else
        return(0)
}

pCop <- function(n, jj, alphas) {
    a0 <- alphas[1]
    a1 <- alphas[2]
    
    L_M(L_M.inv(pbinom(n, 3, 0.1), a0) -
            log(L_B(sum(
                sapply(jj, function(j)
                    L_B.inv(exp(-L_M.inv(F_J(j), a0)), a1)
                    )), a1)),
        a0)
}

dCop <- function(n, jj, alphas) {
    k <- length(c(n, jj))
    mat <- numeric(k)
    for (i in 1:(k - 1)) {
        perm <- unique(permn(c(rep(1, i), rep(0, k - i))))
        nperm <- length(perm)
        
        mat <- rbind(mat, matrix(unlist(perm),
                                 nperm, k, byrow = T))
    }
    mat <- rbind(mat, numeric(k) + 1)
    
    f_NJJ <- sum(sapply(1:nrow(mat), function(i)
        (-1) ^ sum(mat[i,]) * pCop(n - mat[i, 1], jj - mat[i,-1], alphas)))
    
    return(f_NJJ)
}

f_L <- function(l, alphas) {
    if(l == 0)
        return(pbinom(0, 3, 0.1))
    if(l == 1)
        return(dCop(1, 1, alphas))
    if(l==2)
        return(dCop(1, 2, alphas) + 
                   dCop(2, c(1, 1), alphas))
    else
        dCop(1, l, alphas) +
        sum(sapply(1:l, function(j)
            dCop(2, c(j, l - j), alphas))) +
        sum(sapply(1:l, function(j)
            sum(sapply(1:(l - j), function(i)
                dCop(3, c(i, j, l - j - i), alphas)))))
}

exemple19 <- function(alphas, beta_=1.960312){
    L <- 0:42
    kappa <- c(0.9, 0.99, 0.999)
    
    f_L_array <- sapply(L, f_L, alphas)
    E_L <- sum(L * f_L_array)
    V_L <- sum(L ^ 2 * f_L_array) - E_L ^ 2
    
    F_S <- function(s, f_L_array) {
        f_L_array[1] + sum(f_L_array[-1] * pgamma(s, L[-1], beta_))
    }
    
    VaR_S <- function(kappa, f_L_array) {
        sapply(kappa, function(k)
            ifelse(k == f_L_array[1],
                   L[1],
                   uniroot(function(s)
                           F_S(s, f_L_array) - k, c(0, 60)
                   )$root))
    }
    VaR_S_k <- VaR_S(kappa, f_L_array)
    
    TVaR_S <- function(VaR_S_k, kappa, f_L_array) {
        sapply(1:length(kappa), function(i)
            sum(f_L_array * L / beta_ * 
                    (1 - pgamma(VaR_S_k[i], L + 1 , beta_))
                ) / (1 - kappa[i]))
    }
    TVaR_S_k <- TVaR_S(VaR_S_k, kappa, f_L_array)
    
    return(
        list(
            "E_S" = E_L/beta_,
            "V_S" = (E_L + V_L) / beta_^2,
            "kappa" = kappa,
            "VaR_S" = VaR_S_k,
            "TVaR_S" = TVaR_S_k
        )
    )
}

alphas <- list(c(-log(0.5), 5),
               c(-log(0.1), 5),
               c(-log(0.1), 25))

lapply(alphas, exemple19)
