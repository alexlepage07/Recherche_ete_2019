# On commence par tester avec une gentille petite copule de Clayton
library(copula)
library(Deriv)
library(stringr)
library(xtable)
library(ggplot2)
library(psych)


# ================================== Simulations des données d'entraînement ================
n <- 5
q <- 2/5
alpha <- 1.5
beta <- 1/100
nsim <- 1e+4
xtable(data.frame(n, q, alpha, beta))


DATA_train <- rCopula(nsim, claytonCopula(alpha, dim = 6))
DATA_train <- cbind(qbinom(DATA_train[,1], n, q),
                    qexp(DATA_train[,-1], beta)
                    )
colnames(DATA_train) <- c("N", sapply(1:(ncol(DATA_train)-1),function(i) paste0("X", i)))

(sommaire_train <- summary(DATA_train))
xtable(sommaire_train)

datas <- data.frame(c(DATA_train[,1], qbinom((0:100)/100, n, q)),
                    Source <- c(rep("Empirique", length(DATA_train[,1])),rep("Théorique", 101)))

ggplot(datas, aes(datas[,1], fill = Source)) + 
    geom_histogram(alpha = 0.3, aes(y = ..density..), position = 'identity', binwidth = 1)+
    xlab("n") + ylab("Probabilité") +
    theme(legend.title = element_blank())


datas <- data.frame(c(DATA_train[,-1], qexp((0:100)/100, beta)),
                    Source <- c(rep("Simul", length(DATA_train[,-1])),rep("théorique", 101)))

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


# ================================== Estimation des paramètres d'entraînement ======================
para <- c("q"=q, "beta"=beta, "alpha"=alpha)


F_N <- function(x0, pa1) pbinom(x0, n, pa1)
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

temps_deriv <- system.time(
    # Calcul des dérivées
    derivees <- lapply(1:5, function(n)
        parse(text = chain_derivative(func_str_copule(str_copule_ext, str_copule_int, n), n)))
)[[3]]

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

val_depart <- c("q"=mean(DATA_train[,1] / 5),
                "beta"=1/mean(DATA_train[,-1]),
                "alpha"=1)


temps_solv <- system.time(
    mle <- constrOptim(val_depart, 
                       fct_Score,
                       grad = NULL, 
                       ui = diag(3),
                       ci = c(0, 0, 0),
                       outer.eps = 1e-2 )
)

(resultats <- rbind("Valeurs de départ"=round(val_depart,4),
                    "Estimateurs" = round(mle$par, 4),
                    "Vrais paramètres" = round(para,4)))
(temps_tot <- rbind("Temps de dérivation"=temps_deriv,
                    "Temps d'estimation"=temps_solv[[3]]))
xtable(resultats, digits = 4)
xtable(temps_tot)

# load("Clayton_Binomial.RData")

