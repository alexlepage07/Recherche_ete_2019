# Estimation de paramètres du modèle collectif du risque avec structure de dépendance
# sur des données simulées.
#
# Le modèle est construit avec une copule archimédienne hiérachique dont la v.a. mère
# suit une loi logarithmique de même que la v.a. enfant.
library(copula)
library(nCopula)
library(Deriv)
library(stringr)
library(xtable)
library(ggplot2)
library(psych)


# ================================== Simulations des données d'entraînement ================
# Paramètres de simulation
n <- 4
q <- 2/5
beta <- 1/100
alpha0 <- 5
alpha1 <- 7
nsim <- 1e+4

# Simulation:
#
# source("../Code_simul_copule.R")
# DATA_train <- actrisk.rcompcop(nsim,"log", 1-exp(-alpha0), n, "log", 1-exp(-alpha1))
DATA_train <- rCompCop(nsim, LOG(1-exp(-alpha0), 1, list(LOG(1-exp(-alpha1), 1:n, NULL))))
DATA_train <- cbind(qbinom(DATA_train[,1], n, q),
                    qexp(DATA_train[,-1], beta)
                    )
colnames(DATA_train) <- c("N", sapply(1:(ncol(DATA_train)-1),function(i) paste0("X", i)))
nb_xi <- n

# Sommaire des résultats de simulation
(sommaire_train <- summary(DATA_train))
xtable(sommaire_train)  # Permet de  convertir un tableau de R à LaTeX.


#Graphiques de goodness of fit pour les données simulées.
datas <- data.frame(c(DATA_train[,1], qbinom((0:100)/100, n, q)),
                    Source <- c(rep("Empirique", length(DATA_train[,1])),rep("Théorique", 101)))

ggplot(datas,  aes(datas[, 1], fill = Source)) + 
    geom_histogram(alpha = 0.3, aes(y = ..density..), position = 'identity', binwidth = 1)+
    xlab("n") + ylab("Probabilité") +
    theme(legend.title = element_blank())


datas <- data.frame(c(DATA_train[,-1], qexp((0:100)/100, beta)),
                    Source <- c(rep("Simul", length(DATA_train[,-1])), rep("théorique", 101)))

ggplot() + 
    geom_histogram(alpha = 0.3, aes(x= DATA_train[,-1], y = ..density.., fill = "Empirique"), position = 'identity')+
    geom_density(alpha = 0.3, aes(x= qexp((0:100)/100, beta), y = ..density.., fill = "Théorique")) + 
    xlab("x") + ylab("Densité") +
    theme(legend.title = element_blank())


# Scatterplot
stats_ordre <- sapply(1:ncol(DATA_train), function(j)
    rank(DATA_train[,j], ties.method = "first") / (nsim + 1))
colnames(stats_ordre) <- c("N", sapply(1:(ncol(DATA_train)-1),function(i) paste0("X", i)))
pairs.panels(stats_ordre, density = F, ellipses = F, method = "kendall", pch=".")

# ================================== Estimation des paramètres ============================
para <- c("q"=q, "beta"=beta, "alpha0"=alpha0, "alpha1"=alpha1)
nb_para <- length(para)

F_N <- function(x0, pa1) pbinom(x0, n, pa1)
F_X <- "1 - exp(-x_i * pa2)"

LST.Log_M <- "-1 / pa3 * log(1 - (1 - exp(-pa3)) * exp(-T)) "
LST.Log_M.inv <- "- log((1 - exp(-pa3 * U)) / (1 - exp(-pa3)))"
LST.Log_B <- "-1 / pa4 * log(1 - (1 - exp(-pa4)) * exp(-T)) "
LST.Log_B.inv <- "- log((1 - exp(-pa4 * U)) / (1 - exp(-pa4)))"

str_copule_ext <- str_replace(LST.Log_M, "T",
                              paste0(
                                  str_replace(LST.Log_M.inv, "U", "F_N(x0, pa1)"),
                                  " + log(", LST.Log_B,")"
                                  )
                              )
str_copule_int <- str_replace(LST.Log_B.inv, "U",
                              paste0("exp(-",LST.Log_M.inv,")"))


func_str_copule <- function(str_copule_ext, str_copule_int, nb_xi) {
    # Fonction qui génère la chaîne de caractères qui permettra d'effectuer les dérivées.
    str_tot <- str_copule_ext
    str_int <- "0"

    for (i in 1:nb_xi) {
        str_int <- paste(str_int, "-", str_copule_int)
        str_int <- str_replace(str_int,"U", paste0("(",F_X,")"))
        str_int <- str_replace(str_int,"x_i", paste0("x",i))
    }
    str_tot <- str_replace(str_tot, "T", str_int)
    return(str_tot)
}

chain_derivative <- function(str_copule, nb_xi){
    # Fonction qui effectue les dérivations en chaîne sur le texte donné en argument
    X_i <- sapply(1:nb_xi, function(i) paste0("x",i))
    derivees <- str_copule
    
    for (xi in X_i){
        derivees <- Deriv(derivees, xi, cache.exp = T)
    }
    return(derivees)
}

derivees <- list()
temps_deriv <- numeric(nb_xi)
for (n in 1:(nb_xi)){
    # Calcul des dérivées
    temps_deriv[n] <- system.time(
        derivees <- append(derivees, parse(text = 
                              chain_derivative(func_str_copule(str_copule_ext, str_copule_int, n), n)))
    )[[3]]
    print(c("Dérivée"=n, "temps"=temps_deriv[n]))
}

# Graphique du temps de dérivation
plot(temps_deriv[1:nb_xi], type="l",
     xlab="nb de dérivées partielles",
     ylab="temps de dérivation")

temps_deriv <- sum(temps_deriv)


generateur_evalue_deriv <- function(derivee){
    # Générateur permettant d'utiliser la dérivée à titre de fonction évaluable
    # à partir de la liste de dérivées saisie en argument.
    function(n, xx, para){
        # Fonction générée à l'aide des dérivées.
        for (i in 1:nb_para) {
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

val_depart <- c("q" = mean(DATA_train[,1]) / nb_xi,
                "beta" = 1/mean(DATA_train[,-1]),
                "alpha0" = 1,
                "alpha1" = 1)

temps_solv <- system.time(
    # Estimation des paramètres
    mle <- constrOptim(val_depart, fct_Score, grad = NULL, 
                       ui = diag(nb_para),
                       ci = rep(0, nb_para),
                       outer.eps = 1e-2 )
)

# Tableaux permettant de présenter les résultats
(resultats <- rbind("Valeurs de départ"=round(val_depart,4),
                    "Estimateurs" = round(mle$par, 4),
                    "Vrais paramètres" = round(para,4)))
(temps_tot <- rbind("Temps de dérivation"=temps_deriv,
                    "Temps d'estimation"=temps_solv[[3]]))
# Conversion en LaTeX
xtable(resultats, digits = 4)
xtable(temps_tot)
