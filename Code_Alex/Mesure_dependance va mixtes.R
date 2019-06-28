library(xtable)
library(copula)
library(nCopula)
library(rlist)
library(pracma)


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

# À partir d'une matrice de données représentant le modèle collectif du risque
tau_NX <- function(DATA, nsim = 10, silent = T) {
    # À partir d'une matrice contenant des observations issues du modèle collectif 
    # du risque, calcule le tau de Kendall entre N et les X_i, de façon
    # itérative, en prenant aléatoirement un X_i.
    tau <- numeric(nsim)
    temps <- 0
    nb_couples <- sum(DATA[, 1] > 0)
    for (realisation in 1:nsim) {
        temps <- temps + system.time({
            couples_NX <- matrix(ncol = 2, nrow = nb_couples)
            k <- 0
            for (i in 1:nrow(DATA)) {
                if (DATA[i, 1] == 0)
                    next
                couples_NX[k <- k + 1, ] <- c(N <- DATA[i, 1],
                                              sample(DATA[i, 2:(N + 1)], 1))
            }
            tau[realisation] <-
                tau_kendall_empirique(couples_NX[, 1], couples_NX[, 2])
        })[3]
        
        if (silent == F) {
            tbl_iter <- cbind(
                "tau" = round(tau[realisation], 6),
                "temps ecoule" = paste(round(temps), "secondes")
            )
            rownames(tbl_iter) <- paste0(realisation, "/", nsim)
            print(tbl_iter)
        }
    }
    return(tau)
} 

tau_XX <- function(DATA, nsim = 50, silent = T) {
    # À partir d'une matrice contenant des observations issues du modèle collectif 
    # du risque, calcule le tau de Kendall entre les couples (X_i, X_j), de façon
    # itérative, en prenant les couples de façon aléatoire.
    tau <- numeric(nsim)
    temps <- 0
    nb_couples <- sum(DATA[, 1] > 1)
    for (realisation in 1:nsim) {
        temps <- temps + system.time({
            couples_XX <- matrix(ncol = 2, nrow = nb_couples)
            k <- 0
            for (i in 1:nrow(DATA)) {
                if (DATA[i, 1] < 2)
                    next
                N <- DATA[i, 1]
                couples_XX[k <- k + 1,] <- sample(DATA[i, 2:(N + 1)], 2)
            }
            tau[realisation] <- corKendall(couples_XX)[1, 2]
        })[3]
        
        if (silent == F) {
            tbl_iter <- cbind(
                "tau" = round(tau[realisation], 6),
                "temps ecoule" = paste(round(temps), "secondes"))
            rownames(tbl_iter) <- paste0(realisation, "/", nsim)
            print(tbl_iter)
        }
    }
    return(tau)
} 


bootstrap_tau_mixte <- function(data, 
                                nb_sous_intervalles = min(nrow(data) / 10, 1e+4),
                                nb_samples = min(100, nb_sous_intervalles)) {
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

bootstrap_tau_continues <- function(data,
                                    nb_sous_intervalles = min(nrow(data) / 10, 1e+4),
                                    nb_samples = min(100, nb_sous_intervalles)) {
    # Fonction pour utiliser la méthode de ré-échantillonage (bootstrap) sur les
    # données
    tau <- numeric(nb_samples)
    for (i in 1:nb_samples){
        indices <- sample(1:nrow(data), nb_sous_intervalles)
        tau[i] <- corKendall(data[indices, ])[1,2]
    }
    return(tau)
}

#---------------------------------- Calcul du tau de Kendall théorique ------------

tau_kendall_theorique <- function(F_N, F_X,
                                  para,
                                  Copule,
                                  alpha,
                                  n_max, x_max){
    # Fonction qui permet de calculer le tau de kendall avec une variable aléatoire
    # discrète et une autre qui est continue.
    #
    # F_N et F_X peuvent être des fonctions pré-programmées de R.
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
        # Paramétriser la fonction de répartition de X.
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

tau_kendall_theorique_continues <- function(struct_copule,
                                            alpha0, alpha1){
    # Fonction qui admet n'importe quelle copule archimédienne hiérarchique
    # et qui retourne le tau de Kendall théorique pour les cas continus.
    #
    # struct_copule correspond à la structure de la copule archimédienne 
    # hiérarchique. Elle doit sous forme de fonction de alpha.
    #
    # alpha est le paramètre de dépendance de la copule.
    
    C_ <- function(u_1, u_2) {
        # Fonction de répartition conjointe de N et X.
        f <- pCompCop(struct_copule(alpha0, alpha1), T , F)
        # f <- function(U) pCopula(U, Copule(alpha, use.indepC = "FALSE"))
        return(f(c(u_1, u_2)))
    }
    
    c_ <- function(u_1, u_2) {
        eps <- .Machine$double.eps^0.25
        f <- function(u_1, u_2) (C_(u_1 + eps, u_2) - C_(u_1, u_2)) / eps
        f. <- function(u_1, u_2) (f(u_1, u_2 + eps) - f(u_1, u_2)) / eps
        return(f.(u_1, u_2))
    }
    
    return(
        4 * integral2(
            Vectorize(function(u_1, u_2)
                C_(u_1, u_2) * c_(u_1, u_2)),
            0, 1,
            0, 1, 
            reltol = .Machine$double.eps^0.25,
            maxlist = 1e+3
        )$Q - 1
    )
}


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
    alpha <- uniroot(
        function(alpha)
            tau_kendall_theorique(F_N, F_X,
                                  para,
                                  Copule,
                                  alpha,
                                  n_max, x_max) - Tau_Kendall,
        bornes
    )$root
    return(alpha)
}


inversion_tau_kendall_XX <- function(struct_copule,
                                     alpha0,
                                     Tau,
                                     bornes) {
    # Fonction qui permet de trouver le paramètre de dépendance à l'aide du
    # tau de Kendall empirique lorsque le modèle comprend des v.a. strictement continues
    alpha <- uniroot(
        function(alpha1)
            tau_kendall_theorique_continues(struct_copule,
                                            alpha0, alpha1) - Tau,
        bornes)$root
    return(alpha)
}

#======================================= Rho de Spearman ==========================
rho_spearman <- function(X,Y){
    n <- length(X)
    R_i <- rank(X, ties.method = "average")
    S_i <- rank(Y, ties.method = "average")
    rho <- 12 / (n * (n + 1) * (n - 1)) * sum(R_i*S_i) - 3 * (n + 1) / (n - 1)
    return(rho)
}
