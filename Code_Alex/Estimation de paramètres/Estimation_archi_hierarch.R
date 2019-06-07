# On commence par tester avec une gentille petite copule de Clayton
library(copula)
library(nCopula)
library(Deriv)
library(stringr)
library(xtable)


# ================================== Simulations des données d'entraînement ================
source("../Code_simul_copule.R")
n <- 5
q <- 2/5
beta <- 1/100
alpha0 <- 0.5
alpha1 <- 5

nsim <- 1e+3

DATA_train <- actrisk.rcompcop(nsim,"log", 1-exp(-alpha0), 5, "gamma", 1/alpha1)
DATA_train <- cbind(qbinom(DATA_train[,1], n, q),
                    qexp(DATA_train[,-1], beta)
                    )

(sommaire_train <- summary(DATA_train))
xtable(sommaire_train)

par(mfrow=c(1,2))
hist(DATA_train[,1], breaks = 0:5, probability = T)
hist(dbinom(0:5, n, q), breaks = 0:5, probability = T)
plot(ecdf(DATA_train[,-1]))
plot(0:800, pexp(0:800, beta), type="l")
par(mfrow=c(1,1))

# ================================== Estimation des paramètres d'entraînement ======================
para <- c(q=2/5, beta=1/100, alpha0=0.5, alpha1=5)

F_N <- function(x0, pa1) pbinom(x0, 5, pa1)
F_X <- "1 - exp(-x_i * pa2)"

LST.Log <- "-1 / pa3 * log(1 - (1 - exp(-pa3)) * exp(-T)) "
LST.Log.inv <- "- log((1 - exp(-pa3 * U)) / (1 - exp(-pa3)))"
LST.Gamma <- "(1 + T)^(- 1 / pa4)"
LST.Gamma.inv <- "(U^(-pa4) - 1)"

str_copule_ext <- str_replace(LST.Log, "T",
                              paste0(
                                  str_replace(LST.Log.inv, "U", "F_N(x0, pa1)"),
                                  " + log(", LST.Gamma,")"
                                  )
                              )
str_copule_int <- str_replace(LST.Gamma.inv, "U",
                              paste0("exp(-",LST.Log.inv,")"))


func_str_copule <- function(str_copule_ext, str_copule_int, nb_xi) {
    # Génère la chaîne de caractères qui permettra d'effectuer les dérivées.
    str_tot <- str_copule_ext
    str_int <- "0"

    for (i in 1:nb_xi) {
        str_int <- paste(str_int, "+", str_copule_int)
        str_int <- str_replace(str_int,"U", paste0("(",F_X,")"))
        str_int <- str_replace(str_int,"x_i", paste0("x",i))
    }
    str_tot <- str_replace(str_tot, "T", str_int)
    return(str_tot)
}

chain_derivative <- function(str_copule, nb_xi){
    # Fonction qui effectue les dérivations en chaîne sur le texte généré précédemment
    X_i <- sapply(1:nb_xi, function(i) paste0("x",i))
    derivees <- str_copule
    
    for (xi in X_i){
        derivees <- Deriv(derivees, xi, cache.exp = T)
    }
    return(derivees)
}

temps_deriv <- system.time(
    # Calcul des dérivées
    derivees <- lapply(1:5, function(n)
        parse(text = chain_derivative(func_str_copule(str_copule_ext, str_copule_int, n), n)))
)

generateur_evalue_deriv <- function(derivee){
    # Générateur permettant d'utiliser la dérivée à titre de fonction évaluable.
    # à partir de la liste de dérivées saisie en argument.
    function(n, xx, para){
        # Fonction générée à l'aide des dérivées.
        for (i in 1:length(para)) {
            assign(paste0("pa", i), para[i])
        }
        x0 <- n
        for (i in 1:length(xx)) {
            assign(paste0("x", i), xx[i])
        }
        return(eval(derivee[[length(xx)]]))
    }
}
densite <- generateur_evalue_deriv(derivees)

# densite(1,c(1), para)
# densite(2,c(100, 100), para)
# densite(3,c(200, 200, 200), para)
# densite(4,c(500, 500, 500, 500), para)
# densite(5,c(900, 900, 900, 900, 900), para)

dCopule <- function(n, xx=NULL, para){
    # Fonction de densité de la copule
    if (n == 0)
        d <- F_N(0, para[1])
    else
        d <- densite(n, xx, para) - densite(n-1, xx, para)
    return(unname(d))
}
# dCopule(0,NULL, para)
# dCopule(1,c(1), para)
# dCopule(2,c(100, 100), para)
# dCopule(3,c(200, 200, 200), para)
# dCopule(4,c(500, 500, 500, 500), para)
# dCopule(5,c(900, 900, 900, 900, 900), para)

fct_Score <- function(para, Data = DATA_train){
    # Fonction de score: On cherchera à minimiser la log-vraisemblance négative
    neg_log_vrais <- 0
    for (i in 1:nrow(Data)){
        N <- Data[i, 1]
        XX <- Data[i, 2:(N + 1)]
        neg_log_vrais <- neg_log_vrais - log(dCopule(N, XX, para))
    }
    return(neg_log_vrais)
}

val_depart <- c("q"=mean(DATA_train[,1] / 5),
                "beta"=1/mean(DATA_train[,-1]),
                "alpha0"=0.1,
                "alpha1"=0.1)


temps_solv <- system.time(
    mle <- constrOptim(val_depart, 
                       fct_Score,
                       grad = NULL, 
                       ui = diag(4),
                       ci = c(0, 0, 0, 0),
                       outer.eps = 1e-5 )
)
(resultats <- rbind("Estimateurs" = round(mle$par, 4), "Vrais paramètres" = para))
(rbind(temps_deriv, temps_solv))
xtable(resultats)
