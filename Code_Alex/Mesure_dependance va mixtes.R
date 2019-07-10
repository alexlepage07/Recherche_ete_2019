library(xtable)
library(copula)
library(nCopula)
library(rlist)
library(pracma)


#======================================= Tau de Kendall ===========================
#---------------------------------- Calcul du tau de Kendall empirique ------------
tau_kendall_empirique <- function(Data) {
    # Fonction qui calcule le tau de Kendall empirique lorsqu'au moins une des
    # deux v.a. est discète.
    
    # X <- Data[,1] ; Y <- Data[,2]
    # n <- length(X)
    # concord <- outer(1:n,1:n, function(i,j) sign((X[i] - X[j]) * (Y[i] - Y[j])))
    # E_n <- (length(concord[concord == 0]) - n) # On ajoute la probabilité d'avoir des valeurs égales.
    # P_n <- sum(concord[concord > 0])
    # tau <- (2 * P_n + E_n) / (n * (n - 1)) - 1
    
    n <- nrow(na.omit(Data))
    E_n <- P_n <- 0
    for (i in 1:(n - 2)) {
        d <- Data[-1,]
        d[, 1] <- d[, 1] - Data[1, 1]
        d[, 2] <- d[, 2] - Data[1, 2]
        E_n <- E_n + sum(d[, 1] == 0)
        P_n <- P_n + sum(sign(d[, 1]) == sign(d[, 2]))
        Data <- Data[-1,]
    }
    tau <- (4 * P_n + 2 * E_n) / (n * (n - 1)) - 1
    
    return(tau)
}

# À partir d'une matrice de données représentant le modèle collectif du risque
constr_mat_NX <- function(Data){
    # Fonction qui construit une matrice dont qui crée
    # des paires de (N, X_i) de façon aléatoire.
    Data <- Data[(Data[,1] > 0),]
    mat_NX <- matrix(ncol = 2, nrow = nrow(Data))
    
    for (i in 1:nrow(Data)) {
        N <- Data[i, 1]
        X_i <- sample.int(N, 1)
        mat_NX[i,] <- c(N, Data[i, X_i + 1])
    }
    return(mat_NX)
}

tau_NX <- function(DATA, nsim = 100) {
    # À partir d'une matrice contenant des observations issues du modèle collectif 
    # du risque, calcule le tau de Kendall entre N et les X_i, de façon
    # itérative, en prenant aléatoirement un X_i.
    tau <- sapply(1:nsim, function(k)
        tau_kendall_empirique(constr_mat_NX(DATA)))
    
    mean_tau <- mean(tau)
    var_tau <- var(tau)
    IC_tau <- c(mean_tau - qnorm(0.975) * var_tau,
                mean_tau + qnorm(0.975) * var_tau)
    
    output <- list("moyenne"=mean_tau, "Intervalle de confiance a 95%"=IC_tau,
                   "variance"=var_tau, "liste"=tau)
    return(output)
} 


constr_mat_XX <- function(Data){
    # Fonction qui construit une matrice dont qui crée
    # des paires de (N, X_i) de façon aléatoire.
    Data <- Data[(Data[,1] > 1),]
    mat_XX <- matrix(ncol = 2, nrow = nrow(Data))
    
    for (i in 1:nrow(Data)) {
        N <- Data[i, 1]
        X_i <- sample.int(N, 2)
        mat_XX[i,] <- Data[i, X_i + 1]
    }
    return(mat_XX)
}

tau_XX <- function(DATA, nsim = 100) {
    # À partir d'une matrice contenant des observations issues du modèle collectif 
    # du risque, calcule le tau de Kendall entre les couples (X_i, X_j), de façon
    # itérative, en prenant les couples de façon aléatoire.
    tau <- sapply(1:nsim, function(k)
        corKendall(constr_mat_XX(DATA))[1,2])
    
    mean_tau <- mean(tau)
    var_tau <- var(tau)
    IC_tau <- c(mean_tau - qnorm(0.975) * var_tau,
                mean_tau + qnorm(0.975) * var_tau)
    
    output <- list("moyenne"=mean_tau, "Intervalle de confiance a 95%"=IC_tau,
                   "variance"=var_tau, "liste"=tau)
    return(output)
} 

#---------------------------------- Inversion du tau de Kendall -------------------
inversion_tau_kendall <- function(F_N, para, Copule, bornes,
                                  Tau_Kendall, n_max){
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
            tau_kendall_theorique(F_N,
                                  para,
                                  Copule,
                                  alpha,
                                  n_max) - Tau_Kendall,
        bornes
    )$root
    return(alpha)
}


#---------------------------------- Calcul du tau de Kendall théorique ------------
tau_kendall_theorique <- function(F_N, para, Copule,
                                  alpha, n_max){
    # Fonction qui permet de calculer le tau de kendall avec une variable aléatoire
    # discrète et une autre qui est continue.
    #
    # F_N peut être une fonction pré-programmée de R.
    # 
    # Para doit être une liste nommée de la forme 
    # para = list(N=list(size=5, prob=0.3))
    #
    # alpha est le paramètre de dépendance de la copule.
            
    F_N. <- function(n) {
        # Paramétriser la fonction de répartition de N.
        do.call(F_N, list.flatten(list(n, para$N)))
    }
    
    f_N <- function(n) {
        if (n == 0)
            F_N.(n)
        else
            F_N.(n) - F_N.(n-1)}
    
    F_NX <- function(n, u) {
        # Fonction de répartition conjointe de N et X.
        # f <- pCompCop(struct_copule, T , F)
        f <- function(U) pCopula(U, Copule(alpha, use.indepC = "FALSE"))
        return(f(c(F_N.(n), u)))
    }
    
    f_NX <- function(n, u) {
        eps <- .Machine$double.eps^0.25
        f <- function(n, u) (F_NX(n, u + eps) - F_NX(n, u)) / eps
        if (n==0)
            return(f(n, u))
        return(f(n, u) - f(n - 1, u))
    }
    
    return(sum(sapply(0:n_max, function(n)
        4 * integrate(
                Vectorize(function(u)
                    F_NX(n - 1, u) * f_NX(n, u)),
                0, 1)$value +
            f_N(n) ^ 2)) - 1)
}


# Fonctions obtenues grâce à Diamilatou pour calculer le tau de Kendall théorique
# de copules bivariées archiméridennes quelconques dont les v.a. sont continues.
ConstrCopule <- function(mere = c("geo", "log"),
                         enfant = c("geo", "log", "gamma")) {
    # Fonction qui permet de trouver la fonction de densité bivariée des
    # principales copules archimédiennes hiérarchiques.
        mere <- match.arg(mere)
        enfant <- match.arg(enfant)
        
        if (mere == "geo")
            strm <-
            c("(1-alpha0)/(exp(u)-alpha0)",
              "-log(1/(((1-alpha0)/(u))+alpha0))")
        else
            strm <-
            c(" -log(1-(1-exp(-alpha0))*exp(-(u)))/alpha0",
                "-log((1-exp(-alpha0*(u)))/(1-exp(-alpha0)))")
        
        if (enfant == "geo")
            stre <-
            c(" (1-alpha1)/(exp(u)-alpha1)",
              "-log(1/(((1-alpha1)/(u))+alpha1))")
        else if (enfant == "log")
            stre <- c(
                "-log(1-(1-exp(-alpha1))*exp(-(u)))/alpha1",
                "-log((1-exp(-alpha1*(u)))/(1-exp(-alpha1)))")
        else
            stre <- c(" 1/((1+(u))^(1/alpha1))", "(u)^(-alpha1)-1")
        
        Laptheta <-
            str_replace_all(strm[1], "u", paste0("-log(", as.character(stre[1]), ")"))
        Lapthetainv <-
            str_replace_all(stre[2], "u", paste0("exp(-", as.character(strm[2]), ")"))
        
        FrCop <- "chr"
        FrCop <- str_replace_all(Laptheta, "u", "u1+u2")
        
        Lapthetainv <- str_replace_all(Lapthetainv, "u", "u1")
        FrCop <- str_replace_all(FrCop, "u1", as.character(Lapthetainv))
        
        Lapthetainv <- str_replace_all(Lapthetainv, "u1", "u2")
        FrCop <- str_replace_all(FrCop, "u2", as.character(Lapthetainv))
        
        fCop <- Deriv(Deriv(as.character(FrCop), "u1"), "u2")
        expres <-
            paste0("(", as.character(FrCop), ")*(", as.character(fCop), ")")
        
        expr <- parse(text = expres)
        
        return(expr)
}

tau_knd.continu <- function(alpha0, alpha1, expr) {
    # Fonction qui calcule le tau théorique d'une copule bivariéev 
    # archimédienne quelconque dont les v.a. sont continues.
    PEval <- function(u1, u2, alpha0, alpha1, expr) eval(expr)
    
    probb <- integral2(function(u1, u2)
        PEval(u1, u2, alpha0, alpha1, expr), 0, 1, 0, 1)$`Q`
    return(4 * probb - 1)
}

Invert_tau_knd.continu <- function(tau, alpha0, bornes, expr){
    # Fonction qui inverse le tau de Kendall pour une copule bivariée 
    # archimédienne quelconque dont les v.a. sont continues.
    return(uniroot(function(alpha1)
        tau_knd.continu(alpha0, alpha1, expr) - tau, bornes)$root)
}

# Anciennes versions ----
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

bootstrap_tau_mixte <- function(data, 
                                nb_sous_intervalles = min(nrow(data) / 10, 1e+3),
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

#======================================= Rho de Spearman ==========================
rho_spearman <- function(X,Y){
    n <- length(X)
    R_i <- rank(X, ties.method = "average")
    S_i <- rank(Y, ties.method = "average")
    rho <- 12 / (n * (n + 1) * (n - 1)) * sum(R_i*S_i) - 3 * (n + 1) / (n - 1)
    return(rho)
}
