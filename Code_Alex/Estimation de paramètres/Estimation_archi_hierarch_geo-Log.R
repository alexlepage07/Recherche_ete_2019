# On commence par tester avec une gentille petite copule de Clayton
library(copula)
library(nCopula)
library(Deriv)
library(stringr)
library(xtable)


# ================================== Simulations des données d'entraînement ================
source("../Code_simul_copule.R")
n <- 4
q <- 1/2
beta <- 1/100
alpha0 <- 0.5
alpha1 <- 0.8
nsim <- 1e+4

DATA_train <- actrisk.rcompcop(nsim,"log", 1-exp(-alpha0), 5, "gamma", 1/alpha1)
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
pairs.panels(DATA_train, density = F, ellipses = F, method = "spearman", pch=".")

# ================================== Estimation des paramètres d'entraînement ======================
para <- c("q"=q, "beta"=beta, "alpha0"=alpha0, "alpha1"=alpha1)

F_N <- function(x0, pa1) pbinom(x0, 5, pa1)
F_X <- "1 - exp(-x_i * pa2)"

LST.M <- "(1 - pa3)/(exp(T) - pa3)"
LST.M.inv <- "log((1 - pa3) / U + pa3)"
LST.B <- "-1 / pa4 * log(1 - (1 - exp(-pa4)) * exp(-T)) "
LST.B.inv <- "- log((1 - exp(-pa4 * U)) / (1 - exp(-pa4)))"

str_copule_ext <- str_replace(LST.M, "T",
                              paste0(
                                  str_replace(LST.M.inv, "U", "F_N(x0, pa1)"),
                                  " - log(", LST.B,")"
                                  )
                              )
str_copule_int <- str_replace(LST.B.inv, "U",
                              paste0("exp(-",LST.M.inv,")"))


func_str_copule <- function(str_copule_ext, str_copule_int, nb_xi) {
    # Génère la chaîne de caractères qui permettra d'effectuer les dérivées.
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
    # Fonction qui effectue les dérivations en chaîne sur le texte généré précédemment
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
plot(temps_deriv[1:5], type="l",
     xlab="nb de dérivées partielles",
     ylab="temps de dérivation")

temps_deriv <- sum(temps_deriv)


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

# densite(1,c(40), para)
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
                "alpha0"=0.5,
                "alpha1"=0.5)


temps_solv <- system.time(
    mle <- constrOptim(val_depart, 
                       fct_Score,
                       grad = NULL, 
                       ui = diag(4),
                       ci = c(0, 0, 0, 0),
                       outer.eps = 1e-5 )
)
(resultats <- rbind("Estimateurs" = round(mle$par, 4), "Vrais paramètres" = round(para,4)))
(temps_tot <- rbind("temps de dérivation"=temps_deriv, "temps d'estimation"=temps_solv[[3]]))
xtable(resultats)
xtable(temps_tot)
