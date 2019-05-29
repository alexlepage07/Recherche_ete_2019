#============================= Code pour simuler des copules ======================
library(actuar)
# Loi géométrique
LST_geom <- function(p, t) {
    p / (exp(t) - (1 - p))
}

# LST_geom_inv <- function(p, u) {
#     log(p * (1 - u) / u + 1)
# }


#Loi logistique
LST_logarithmic <- function(g, t) {
    log(1 - g * exp(-t)) / log(1 - g)
}

# LST_logarithmic_inv <- function(g, u) {
#     -log((1 - (1 - g) ^ u) / g)
# }


# Loi gamma
LST_gamma <- function(alpha, t) {
    (1 + t) ^ (-alpha)
}

# LST_gamma_inv <- function(alpha, t) {
#     u ^ (-1 / alpha) - 1
# }


# Fonctions génératrices de probs
actrisk.compcop <- function(mother, par_mother, n, child, par_child) {
        d <- length(child)
        if (d != length(par_child) | d != length(n)){
            return("Le nombre de paramètres doit
                   être égal au nombre de dimensions (enfants)")}
        
        if (mother == "geo") {
            M <- rgeom(1, par_mother) + 1
            
            L_M <- function(t)
                LST_geom(par_mother, t)
        }
        if (mother == "log") {
            M <- rlogarithmic(1, par_mother)
            
            L_M <- function(t)
                LST_logarithmic(par_mother, t)
        }
        if (mother == "gamma") {
            M <- rgamma(1, par_mother) + 1
            
            L_M <- function(t)
                LST_gamma(par_mother, t)
        }
        
        U <- c(L_M(rexp(1) / M))
        
        for (i in 1:d) {
            if (child[i] == "geo") {
                theta_i <- rnbinom(1, M, par_child[i]) + M
                
                L_theta <- function(t)
                    L_M(-log(LST_geom(par_child[i], t)))
            }
            if (child[i] == "log") {
                theta_i <- sum(rlogarithmic(M, par_child[i]))
                
                L_theta <- function(t)
                    L_M(-log(LST_logarithmic(par_child[i], t)))
            }
            if (child[i] == "gamma") {
                theta_i <- rgamma(1, M * par_child[i], 1)
                
                L_theta <- function(t)
                    L_M(-log(LST_gamma(par_child[i], t)))
            }
            
            U <- append(U, L_theta(rexp(n[i]) / theta_i))
        }
        return(U)
    }

actrisk.simul <- function(nsim, FUNC, ...) {
    return(t(sapply(1:nsim, function(i)
        FUNC(...))))
}

actrisk.rcompcop <- function(nsim, mother, par_mother, n, child, par_child) {
    actrisk.simul(nsim,
                  actrisk.compcop,
                  mother,
                  par_mother,
                  n,
                  child,
                  par_child)
}