library(copula)
library(nCopula)
library(Deriv)
library(stringr)

# source("../../Mesure_dependance va mixtes.R")

simul_modele.collectif <- function(n, q, beta, alpha0, alpha1, nsim){
    # Fonction qui permet de simuler le modèle collectif du risque avec une loi 
    # binomiale comme loi de fréquence, une loi exponentielle comme loi de sévérité
    # et une copule archimédienne hiérarchique LOG-LOG
    DATA_train <- rCompCop(nsim, LOG(1-exp(-alpha0), 1, list(GAMMA(1/alpha1, 1:nb_xi, NULL))))
    
    for (i in 1:nsim){
        DATA_train[i, 1] <- N <- qbinom(DATA_train[i,1], nb_xi, q)
        if (N==0){ 
            DATA_train[i,-1] <- rep(NaN, nb_xi)
        }else{
            DATA_train[i,-1] <- c(qexp(DATA_train[i, 2:(N + 1)], beta), rep(NaN, nb_xi - N))}
    }
    
    colnames(DATA_train) <- c("N", sapply(1:(ncol(DATA_train)-1),function(i) paste0("X", i)))
    return(DATA_train)
}

nb_xi <- 5
q <- 1/2
beta <- 1/100
alpha0 <- 4
alpha1 <- 6
nsim <- n_obs <- 1e+3

set.seed(2^12)
DATA_train <- simul_modele.collectif(nb_xi, q, beta, alpha0, alpha1, nsim)

# head(DATA_train)
# summary(DATA_train)
# barplot(table(DATA_train[,1]))

# Méthode 1: Maximum de vraisemblance complète
TLS <- function(mere=c("geo", "log"),
                enfant=c("geo", "log", "gamma")) {
    # Fonction qui sort une chaîne de charactère pour
    # la TLS de theta et la TLS inverse de theta.
    mere <- match.arg(mere)
    enfant <- match.arg(enfant)
    
    TLS_log <- "-1 / ai * log(1 - (1 - exp(-ai)) * exp(-T))"
    TLS_geo <- "(1 - ai) / (exp(T) - ai)"
    TLS_gamma <- "(1 + T)^(-1 / ai)"
    TLS_log.inv <- "-log((1 - exp(-ai * ui)) / (1 - exp(-ai)))"
    TLS_geo.inv <- "log((1 - ai) / ui + ai)"
    TLS_gamma.inv <- "(ui^{-ai} - 1)"
    
    TLS_mere <- get(paste0("TLS_", mere))
    TLS_mere.inv <- get(paste0("TLS_", mere, ".inv"))
    TLS_enfant <- get(paste0("TLS_", enfant))
    TLS_enfant.inv <- get(paste0("TLS_", enfant, ".inv"))
    
    
    # TLS de theta
    TLS_Theta <- str_replace_all(TLS_mere, "ai", "a0")
    TLS_Theta <- str_replace_all(TLS_Theta, "T",
                           paste0("-log(", TLS_enfant,")"))
    TLS_Theta <- str_replace_all(TLS_Theta, "ai", "a1")
    # TLS_Theta <- Simplify(TLS_Theta)
    
    # TLS inverse de theta
    TLS_Theta.inv <- str_replace_all(TLS_enfant.inv, "ai", "a1")
    TLS_Theta.inv <- str_replace_all(TLS_Theta.inv, "ui",
                           paste0("exp(-",TLS_mere.inv,")"))
    TLS_Theta.inv <- str_replace_all(TLS_Theta.inv, "ai", "a0")
    # TLS_Theta.inv <- Simplify(TLS_Theta.inv)
    
    return(list("TLS_Theta" = TLS_Theta, 
                "TLS_Theta.inv" = TLS_Theta.inv))
}


derivee <- function(TLS_, d=2) {
    # Fonction qui trouve la densité d'une copule
    # archimédienne de d dimensions.
    Psi <- TLS_$TLS_Theta
    Psi.inv <- TLS_$TLS_Theta.inv
    
    deriv_1 <- Deriv(Psi, "T")
    for (i in 1:(d - 1)) {
        deriv_1 <- Deriv(deriv_1, "T")
    }
    deriv_1 <- str_replace(deriv_1, "T", 
                             paste0("sum(sapply(U, function(ui)",
                                    Psi.inv, "))"))
    # deriv_1 <- paste0("(-1)^", d, " * (", deriv_1, ")")
    
    
    deriv_2 <- Deriv(Psi, "T")
    deriv_2 <- str_replace(deriv_1, "T", Psi.inv)
    deriv_2 <- paste0("prod(sapply(U, function(ui)",
                      deriv_2, "))")
    
    deriv_tot <- paste("(",deriv_1, ") / (", deriv_2,")")
    
    return(parse(text = deriv_tot))
}


list_derivees <- function(mere=c("geo", "log"),
                          enfant=c("geo", "log", "gamma"), d_max) {
    # Fonction qui sort la liste des dérivées nécessaires pour calculer
    # la vraisemblance du modèle collectif du risque.
    mere <- match.arg(mere)
    enfant <- match.arg(enfant)
    TLS_ <- TLS(mere, enfant)
    
    lapply(2:d_max, function(d)derivee(TLS_, d))
}


dCopule <- function(U, a0, a1, liste_derivees) {
    # fonction qui évalue la densité d'une copule archimédienne
    # à partir de la liste des dérivées fournie par list_derivees.
    U <- U[!is.na(U)]
    
    dCopule <- function(U, a0, a1)
        eval(liste_derivees[[length(U) - 1]])
    
    print(dCopule(U, a0, a1))
    return(dCopule(U, a0, a1))
}


fct_score <- function(Data, a0, a1, liste_derivees, F_X, ...) {
    # Fonction a minimiser afin de maxisimer la vraisemblance.
    # ... correspond aux arguments de la fonction F_X
    Data <- Data[Data[,1] > 1,]
    Data <- Data[,-1]
    Data <- F_X(Data, ...)
    
    -sum(log(apply(Data, 1, dCopule, a0, a1, liste_derivees)))
}


# Tests et scénario ----
temps_deriv <- system.time(
    # Max d: - mere = log : 13 ,
    #        - mere = geo : plus de 6 = très long.
    liste_derivees <- list_derivees("log", "gamma", nb_xi)
)[[3]]

fct_Score_XX <- function(a1){
    fct_score(DATA_train, alpha0, a1, liste_derivees, pexp, beta)}

Data <- DATA_train
Data <- Data[Data[,1] > 1,]
Data <- pexp(Data[,-1], beta)
summary(apply(Data, 1, dCopule, alpha0, alpha1, liste_derivees))

temps_optimize <- system.time({
    plot(domaine <- 1:10, sapply(domaine, fct_Score_XX), type="l")
    axis(2, tck = 1, lty = 2, col = "grey")
    axis(1, tck=1, lty = 2, col = "grey",)
})[[3]]

temps_optimize <- temps_optimize + system.time(
    mle_alpha1 <- optimize(fct_Score_XX, interval = c(1e-6, 10))$minimum
)

