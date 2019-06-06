# On commence par tester avec une gentille petite copule de Clayton
library(copula)
library(Deriv)
library(stringr)
library(xtable)


# ================================== Simulations des données d'entraînement ================
lambda <- 2
beta <- 1/100
alpha <- 6
nsim <- 1e+4
nb_xi <- ceiling(qpois(0.99999, lambda))

DATA_train <- rCopula(nsim, claytonCopula(alpha, dim = nb_xi))
DATA_train <- cbind(qpois(DATA_train[,1], lambda),
                    qexp(DATA_train[,-1], beta)
                    )
summary(DATA_train)
nb_xi <- max(DATA_train[,1])

# ================================== Estimation des paramètres d'entraînement ======================
para <- c(lambda=1, beta=1/100, alpha=6)


F_N <- function(x0, pa1) ppois(x0, pa1)
F_X <- "1 - exp(-x_i * pa2)"
str_copule_ext <- "(U - n + 1) ^ (-1/pa3)"
str_copule_int <- " u ^ (-pa3)"


func_str_copule <- function(str_copule_ext, str_copule_int, nb_xi) {
    # Génère la chaîne de caractères qui permettra d'effectuer les dérivées.
    str_tot <- str_copule_ext
    str_int <- str_replace(str_copule_int, "u", "F_N(x0, pa1)")

    for (i in 1:nb_xi) {
        str_int <- paste(str_int, "+", str_copule_int)
        str_int <- str_replace(str_int,"u", paste0("(",F_X,")"))
        str_int <- str_replace(str_int,"x_i", paste0("x",i))
    }
    str_tot <- str_replace(str_tot, "U", str_int)
    str_tot <- str_replace(str_tot, "n", as.character(nb_xi + 1))
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

time <- system.time(
    # Évaluation des dérivées
    derivees <- lapply(1:nb_xi, function(n)
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


dCopule <- function(n, xx=NULL, para){
    # Fonction de densité de la copule
    if (n == 0)
        d <- F_N(0, para[1])
    else
        d <- densite(n, xx, para) - densite(n-1, xx, para)
    return(unname(d))
}


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

val_depart <- para

time <- time + system.time(
    mle <- constrOptim(val_depart, 
                       fct_Score,
                       grad = NULL, 
                       ui = diag(3),
                       ci = c(0, 0, 0),
                       outer.eps = 1e-5 )
)
(resultats <- rbind("Estimateurs" = round(mle$par, 2), "Vrais paramètres" = para))
time
xtable(resultats)
load("Clayton_Poisson.RData")
