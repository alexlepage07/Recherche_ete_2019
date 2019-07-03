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
library(combinat)

source("../../Mesure_dependance va mixtes.R")

set.seed(20190702)
# ================================== Simulations des données d'entraînement ================
# Paramètres de simulation
nb_xi <- 5
q <- 2/5
beta <- 1/100
alpha0 <- 8
alpha1 <- 7
nsim <- n_obs <- 1e+4

DATA_train <- rCompCop(nsim, LOG(1-exp(-alpha0), 1, list(LOG(1-exp(-alpha1), 1:5, NULL))))

for (i in 1:nsim){
    DATA_train[i, 1] <- N <- qbinom(DATA_train[i,1], nb_xi, q)
    if (N==0){ 
        DATA_train[i,-1] <- rep(NaN, nb_xi)
    }else{
        DATA_train[i,-1] <- c(qexp(DATA_train[i, 2:(N + 1)], beta), rep(NaN, nb_xi - N))}
}

colnames(DATA_train) <- c("N", sapply(1:(ncol(DATA_train)-1),function(i) paste0("X", i)))
head(DATA_train)

# Sommaire des résultats de simulation
(sommaire_train <- summary(DATA_train))
# xtable(sommaire_train)  # Permet de  convertir un tableau de R à LaTeX.

# ====================== Estimation des paramètres ds v.a. N et X ======================
# Estimation des paramètres de N
(mle_binom <-  optimize(function(par) -sum(log(dbinom(DATA_train[, 1], nb_xi, par))),
                       lower = 1e-6, upper = 1-1e-6)$minimum)
q

# Estimation des paramètres de X
lst_XX <- array(DATA_train[, -1])
lst_XX <- lst_XX[!is.na(lst_XX)]
summary(lst_XX)

# plot(domaine <- (1:100)/1000, 
#      sapply(domaine, function(para) -sum(log(dexp(lst_XX, para)))),
#      type="l")
# axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
# axis(1, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
# abline(a=tau_0, b=0, col="green")

(mle_exp <-  optimize(function(par) -sum(log(dexp(lst_XX, par))), 
                      interval = c(1/200, 1/55))$minimum)
beta


# Graphiques de goodness of fit pour les données simulées.
datas <- data.frame(c(DATA_train[,1], qbinom((0:100)/100, nb_xi, mle_binom)),
                    Source <- c(rep("Empirique", length(DATA_train[,1])),rep("Théorique", 101)))

ggplot(datas,  aes(datas[, 1], fill = Source)) + 
    geom_histogram(alpha = 0.3, aes(y = ..density..), position = 'identity', binwidth = 1)+
    xlab("n") + ylab("Probabilité") +
    theme(legend.title = element_blank())

Q <- 0
for (n in unique(DATA_train[,1])){
    Q <- Q + (sum(DATA_train[,1] == n) - n_obs * dbinom(n, nb_xi, mle_binom)) ^ 2 /
        sum(DATA_train[,1] == n)
}
chi2 <- rbind(
              "Binomiale" = c(
                  Q,
                  qchisq(0.95, length(unique(DATA_train[,1])) - 2, lower.tail = T),
                  pchisq(Q, length(unique(DATA_train[,1])) - 2, lower.tail = F)
              ))
colnames(chi2) <- c("Statistique", "valeur critique" , "P-value")
chi2


datas <- data.frame(c(lst_XX, qexp((0:100)/100, mle_exp)),
                    Source <- c(rep("Simul", length(lst_XX)),
                                rep("théorique", 101)))

ggplot() + 
    geom_histogram(alpha = 0.3, aes(x= lst_XX, y = ..density.., fill = "Empirique"), position = 'identity')+
    geom_density(alpha = 0.3, aes(x= qexp((0:100)/100, beta), y = ..density.., fill = "Théorique")) + 
    xlab("x") + ylab("Densité") +
    theme(legend.title = element_blank())

ks.test(lst_XX, pexp, mle_exp)

# ====================== Estimation du paramètre de la v.a. M. ============================
Copule0 <- frankCopula
para <- c(q, beta, alpha0)

F_N <- function(x0) pbinom(x0, nb_xi, mle_binom)
F_X <- "(1 - exp(-x * mle_exp))"

str_copule <- "-1 / pa1 * log( 1 + 1 / (exp(-pa1) - 1) * (exp(-pa1 * u1) - 1) * (exp(-pa1 * u2) - 1) )"
str_copule <- str_replace(str_copule, "u1", "F_N(x0)")
str_copule <- str_replace(str_copule, "u2", F_X)

derivees <- Deriv(str_copule, "x", cache.exp = T)
derivees <-  parse(text = derivees)


generateur_evalue_deriv <- function(derivee){
    # Générateur permettant d'utiliser la dérivée à titre de fonction évaluable.
    # à partir de la liste de dérivées saisie en argument.
    function(n, x, para){
        # Fonction générée à l'aide des dérivées.
        pa1 <- para
        x0 <- n
        return(eval(derivee))
    }
}
densite <- generateur_evalue_deriv(derivees)


f_NX <- function(n, x, para){
    # Fonction de densité de la copule
    if (n == 0)
        d <- F_N(0)
    else
        d <- densite(n, x, para) - densite(n-1, x, para)
    return(unname(d))
}


fct_Score_NX <- function(para, Data = DATA_train) {
    # Fonction de score: On cherchera à minimiser la log-vraisemblance négative
    neg_log_vrais <- 0
    for (i in 1:nrow(Data)) {
        N <- Data[i, 1]
        if (N == 0) {
            neg_log_vrais <- neg_log_vrais - log(F_N(0))
        } else{
            for (j in 1:N) {
                X <- Data[i, j + 1]
                neg_log_vrais <-
                    neg_log_vrais - log(f_NX(N, X, para))
                # print(c(N, X, f_NX(N,X,para)))}
            }
        }
    }
    return(neg_log_vrais)
}


plot(domaine <- 1:10, sapply(domaine, fct_Score_NX), type="l")
axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
axis(1, at=domaine, tck=1, lty = 2, col = "grey",) # L'axe des abscisses

temps_solv <- system.time(
    mle_M <- optimize(fct_Score_NX, interval = c(5, 9))
)

# Tableaux permettant de présenter les résultats
(resultats <- rbind("Estimateurs" = round(mle_M$par, 4),
                    "Vrais paramètres" = round(alpha0,4)))
(temps_tot <- temps_solv[[3]])
# Conversion en LaTeX
xtable(resultats, digits = 4)
xtable(temps_tot)

alpha0 <- mle_M$par


# ====================== Estimation des paramètres des v.a. M, N et X. ==========================
Copule0 <- frankCopula
para <- c(q, beta, alpha0)

F_N <- function(x0, pa1) pbinom(x0, nb_xi, pa1)
F_X <- "(1 - exp(-x * pa2))"

str_copule <- "-1 / pa3 * log( 1 + 1 / (exp(-pa3) - 1) * (exp(-pa3 * u1) - 1) * (exp(-pa3 * u2) - 1) )"
str_copule <- str_replace(str_copule, "u1", "F_N(x0, pa1)")
str_copule <- str_replace(str_copule, "u2", F_X)

derivees <- Deriv(str_copule, "x", cache.exp = T)
derivees <-  parse(text = derivees)


generateur_evalue_deriv <- function(derivee){
    # Générateur permettant d'utiliser la dérivée à titre de fonction évaluable.
    # à partir de la liste de dérivées saisie en argument.
    function(n, x, para){
        # Fonction générée à l'aide des dérivées.
        for (i in 1:length(para)) {
            assign(paste0("pa", i), para[i])
        }
        x0 <- n
        return(eval(derivee))
    }
}
densite <- generateur_evalue_deriv(derivees)


f_NX <- function(n, x, para){
    # Fonction de densité de la copule
    if (n == 0)
        d <- F_N(0, pa1)
    else
        d <- densite(n, x, para) - densite(n-1, x, para)
    return(unname(d))
}


fct_Score_NX <- function(para, Data = DATA_train) {
    # Fonction de score: On cherchera à minimiser la log-vraisemblance négative
    neg_log_vrais <- 0
    for (i in 1:nrow(Data)) {
        N <- Data[i, 1]
        if (N == 0) {
            neg_log_vrais <- neg_log_vrais - log(F_N(0, para[1]))
        } else{
            for (j in 1:N) {
                X <- Data[i, j + 1]
                neg_log_vrais <-
                    neg_log_vrais - log(f_NX(N, X, para))
                # print(c(N, X, f_NX(N,X,para)))}
            }
        }
    }
    return(neg_log_vrais)
}

{
    tau_0 <- mean(tau_NX(DATA_train, nsim=5, silent=F))
    
    
    Tau_graph <- sapply(domaine <- 5:10, function(a)
        tau_kendall_theorique(pbinom,
                              list(N=list(size=nb_xi, prob=mle_binom)),
                              Copule0, a, nb_xi))
    
    plot(domaine, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=domaine, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
    abline(a=tau_0, b=0, col="green")
    
    
    bornes <- c(8, 9)
    
    temps_0 <- system.time(
        alpha0_n <- inversion_tau_kendall(pbinom,
                                          list(N=list(size=nb_xi, prob=mle_binom)),
                                          Copule0,
                                          bornes,
                                          tau_0,
                                          nb_xi)
    )
    
    (tbl_0 <- rbind(
        "Vrai paramètre" = alpha0,
        "Paramètre trouvé" = alpha0_n,
        "Temps d'optimisation" = temps_0[[3]]
    ))
} # Estimation alpha0 par la méthode des moments.

# Valeurs de départ de l'optimisation
(val_depart <- c(
    "prob" = mean(DATA_train[, 1]) / nb_xi,
    "beta" = 1 / mean(DATA_train[,-1], na.rm = T),
    "alpha0" = alpha0_n
))

# Bornes de l'optimisation
bounds <- matrix(c(1e-6, 1 - 1e-6,
                   1e-6, Inf,
                   1e-6, Inf), ncol = 2, byrow = T)
# Convertir les constraintes en matrices ui et ci.
n_par <- nrow(bounds)
ui <- rbind( diag(n_par), -diag(n_par) )
ci <- c( bounds[,1], - bounds[,2] )
# Retirer les valeurs infinies
i <- as.vector(is.finite(bounds))
ui <- ui[i,]
ci <- ci[i]

temps_solv <- system.time(
    # Estimation des paramètres
    mle_NX <- constrOptim(val_depart, fct_Score_NX, grad=NULL,
                         ui=ui, ci=ci)
)

# Tableaux permettant de présenter les résultats
(resultats <- rbind("Valeurs de départ"=round(val_depart,4),
                    "Estimateurs" = round(mle_NX$par, 4),
                    "Vrais paramètres" = round(para,4)))
(temps_tot <- temps_solv[[3]])
# Conversion en LaTeX
xtable(resultats, digits = 4)
xtable(temps_tot)

alpha0 <- mle_NX$par[3]

# ====================== Estimation des paramètres des v.a. theta et X. ==========================
F_X. <- function(x) {
    # Fonction de répartition marginale de X.
    eval(parse(text = F_X), list(pa2 = beta))
}

Copule1 <- function(xx, alpha1) {
    # Définission de la structure de la copule à utiliser.
    k <- length(xx)
    LOG(1 - exp(-alpha0), NULL, list(LOG(1 - exp(-alpha1), 1:k, NULL)))
}

F_XX <- function(xx, alpha1, struct_Copule){
    # Fonction de répartition conjointe de X_1,...,X_k.
    f <- pCompCop(struct_Copule(xx, alpha1), T, F)
    return(f(c(F_X.(xx))))
}

f_XX <- function(xx, alphas, struct_Copule) {
    # Fonction de densité conjointe de X_1,...,X_k.
    eps <- .Machine$double.eps^0.25
    k <- length(xx)
    mat <- numeric(k)
    for (i in 1:(k - 1)) {
        perm <- unique(permn(c(rep(1, i), rep(0, k - i))))
        nperm <- length(perm)
        
        mat <- rbind(mat, matrix(unlist(perm),
                                 nperm, k, byrow = T))
    }
    mat <- rbind(mat, numeric(k) + 1)
    
    f_XX <- sum(sapply(1:nrow(mat), function(i)
        (-1) ^ sum(mat[i,]) * F_XX(xx - mat[i,] * eps, alphas, struct_Copule))) / eps
    
    return(f_XX)
}

fct_Score_XX <- function(para, struct_Copule=Copule1, Data = DATA_train) {
    # Fonction de score: On cherchera à minimiser la log-vraisemblance négative
    neg_log_vrais <- 0
    for (i in 1:nrow(Data)) {
        N <- Data[i, 1]
        if (N > 2)
            next
        neg_log_vrais <- neg_log_vrais - log(f_XX(Data[i, 2:(N+1)], para, struct_Copule))
    }
    return(neg_log_vrais)
}

plot(domaine <- (5:10), sapply(domaine, fct_Score_NX), type="l")
axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
axis(1, at=domaine, tck=1, lty = 2, col = "grey",) # L'axe des abscisses

temps_solv <- system.time(
    mle_theta <- optimize(fct_Score_XX, interval = c(6.5, 7.5))
)

# Tableaux permettant de présenter les résultats
(resultats <- rbind("Estimateurs" = round(mle_theta$par, 4),
                    "Vrais paramètres" = round(alpha1,4)))
(temps_tot <- temps_solv[[3]])
# Conversion en LaTeX
xtable(resultats, digits = 4)
xtable(temps_tot)
