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
library(graphics)

source("../../Mesure_dependance va mixtes.R")
#

# ====================== Simulations des données d'entraînement ================
simul_modele.collectif <- function(n, q, beta, alpha0, alpha1, nsim){
    # Fonction qui permet de simuler le modèle collectif du risque avec une loi 
    # binomiale comme loi de fréquence, une loi exponentielle comme loi de sévérité
    # et une copule archimédienne hiérarchique LOG-LOG
    DATA_train <- rCompCop(nsim, LOG(1-exp(-alpha0), 1, list(LOG(1-exp(-alpha1), 1:nb_xi, NULL))))
    
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
# Paramètres de simulation
nb_xi <- 5
q <- 3/5
beta <- 1/100
alpha0 <- 8
alpha1 <- 6
nsim <- n_obs <- 1e+3

set.seed(2^12)
DATA_train <- simul_modele.collectif(nb_xi, q, beta, alpha0, alpha1, nsim)

head(DATA_train)
(sommaire_train <- summary(DATA_train))
# xtable(sommaire_train)  # Permet de  convertir un tableau de R à LaTeX.

# data_exemple1 <- simul_modele.collectif(nb_xi, q, beta, 10, 10, 1e+4)
# sommaire <- cbind(t(summary(data_exemple1[,1])), "NA's"=0)
# for (i in 2:ncol(data_exemple1[,-1])){
#     if (length(summary(data_exemple1[,i])) < 7)
#         sommaire <- rbind(
#             sommaire,
#             cbind(t(summary(data_exemple1[,i])), "NA's"=0)
#             )
#     else
#         sommaire <- rbind(sommaire, summary(data_exemple1[,i]))
# }
# sommaire
# rownames(sommaire) <- c("N", paste0("X",1:4))
# xtable(head(data_exemple1))
# xtable(sommaire)
# mean(data_exemple1[,-1],na.rm = T)

# ====================== Estimation des paramètres des v.a. N et X ======================
# Estimation des paramètres de N
temps_N <- system.time(
    mle_binom <-  optimize(function(par) -sum(log(dbinom(DATA_train[, 1], nb_xi, par))),
                           lower = 1e-6, upper = 1-1e-6)$minimum
)[[3]]
mle_binom; q
temps_N

# Estimation des paramètres de X

# arr_XX <- array(DATA_train[,-1])
# arr_XX <- arr_XX[!is.na(arr_XX)]
# summary_X <- summary(arr_XX)
# (mle_exp <- 1/mean(arr_XX))

mat_X <- matrix(nrow = sum(DATA_train[,1] > 0), ncol=100)
for (j in 1:100){
    k <- 0
    for (i in 1:nrow(DATA_train)){
        if ((N <- DATA_train[i, 1]) == 0)
            next
        mat_X[k <- k + 1, j] <- sample(DATA_train[i, 2:(N + 1)], 1)
    }
}

arr_X <- array(mat_X)
(summary_X <- summary(arr_X))


temps_X <- system.time(
(mle_exp <-  optimize(function(par) -sum(log(dexp(arr_X, par))), 
                      interval = c(1/summary_X[[5]], 1/summary_X[[2]]))$minimum)
)[[3]]
mle_exp; beta
temps_X

# Graphiques de goodness of fit pour les données simulées. ----
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

temps_M <- system.time({
    
    plot(domaine <- 1:10, sapply(domaine, fct_Score_NX), type="l")
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=domaine, tck=1, lty = 2, col = "grey",) # L'axe des abscisses

    mle_M <- optimize(fct_Score_NX, interval = c(6, 9))$minimum
})[[3]]


# ====================== Estimation du paramètre de la v.a. theta ==========================
# Méthode de dérivation numérique ----
# F_X. <- function(x) {
#     # Fonction de répartition marginale de X.
#     eval(parse(text = F_X), list(pa2 = mle_exp))
# }
# 
# Copule1 <- function(xx, alpha1) {
#     LOG(1 - exp(-alpha0), NULL, list(LOG(1 - exp(-alpha1), 1:length(xx), NULL)))
# }
# 
# F_XX <- function(xx, alpha1, struct_Copule){
#     # Fonction de répartition conjointe de X_1,...,X_k.
#     f <- pCompCop(struct_Copule(xx, alpha1), T, F)
#     return(unname(f(c(F_X.(xx)))))
# }
# 
# 
# mat_deriv <- function(d){
#     # Fonction qui permet de créer la matrice de 1 et de 0 qui servira à
#     # produire la dérivation numérique.
#     mat <- numeric(d)
#     for (i in 1:(d - 1)) {
#         perm <- unique(permn(c(rep(1, i), rep(0, d - i))))
#         nperm <- length(perm)
#         
#         mat <- rbind(mat, matrix(unlist(perm),
#                                  nperm, d, byrow = T))
#     }
#     mat <- rbind(mat, numeric(d) + 1)
#     return(mat)
# }
# 
# temps_mat_deriv <- system.time(
#     # Création des matrices de dérivation
#     lst_mat_deriv <- lapply(2:nb_xi, mat_deriv)
# )
# 
# f_XX <- function(xx, alpha, struct_Copule) {
#     # Fonction de densité conjointe de X_1,...,X_k.
#     eps <- .Machine$double.eps^0.25
#     mat <- lst_mat_deriv[[length(xx) - 1]]
#     
#     f_XX <- sum(sapply(1:nrow(mat), function(i)
#         (-1) ^ (sum(mat[i,])+1) * F_XX(xx - mat[i,] * eps, alpha, struct_Copule))) / eps
# 
#     return(abs(f_XX))
# }
# 
# fct_Score_XX <- function(alpha, struct_Copule=Copule1, Data = DATA_train) {
#     # Fonction de score: On cherchera à minimiser la log-vraisemblance négative
#     neg_log_vrais <- 0
#     for (i in 1:nrow(Data)) {
#         N <- Data[i, 1]
#         if (N < 2)
#             next
#         neg_log_vrais <- neg_log_vrais - log(f_XX(Data[i, 2:(N+1)], alpha, struct_Copule))
#     }
#     return(neg_log_vrais)
# }
# 
# temps_theta_deriv_num <- system.time({
#     
#     plot(domaine <- (5:10), sapply(domaine, fct_Score_XX), type="l")
#     axis(2, tck = 1, lty = 2, col = "grey")
#     axis(1, at=domaine, tck=1, lty = 2, col = "grey")
# 
# 
#     mle_theta_deriv_num <- optimize(fct_Score_XX, interval = c(6, 8))
# })
# 
# # Tableaux permettant de présenter les résultats
# (resultats <- rbind("Estimateurs" = round(mle_theta$par, 4),
#                     "Vrais paramètres" = round(alpha1,4)))
# (temps_tot <- temps_solv[[3]])
# # Conversion en LaTeX
# xtable(resultats, digits = 4)
# xtable(temps_tot)


# Méthode de dérivation symbolique ----
mle_alpha0 <- mle_M

F_X <- "1 - exp(-x_i * mle_exp)"

LST.Log_M <- "-1 / mle_alpha0 * log(1 - (1 - exp(-mle_alpha0)) * exp(-T)) "
LST.Log_M.inv <- "- log((1 - exp(-mle_alpha0 * U)) / (1 - exp(-mle_alpha0)))"
LST.Log_B <- "-1 / pa1 * log(1 - (1 - exp(-pa1)) * exp(-T)) "
LST.Log_B.inv <- "log((1 - exp(-pa1 * U)) / (1 - exp(-pa1)))"

str_copule_ext <- str_replace(LST.Log_M, "T",
                              paste0("-log(", LST.Log_B,")"))
str_copule_int <- str_replace(LST.Log_B.inv, "U",
                              paste0("exp(-",LST.Log_M.inv,")"))


func_str_copule <- function(str_copule_ext, str_copule_int, nb_xi) {
    # Fonction qui génère la chaîne de caractères qui permettra d'effectuer les dérivées.
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
    # Fonction qui effectue les dérivations en chaîne sur le texte donné en argument
    X_i <- sapply(1:nb_xi, function(i) paste0("x",i))
    derivees <- str_copule
    
    for (xi in X_i){
        derivees <- Deriv(derivees, xi, cache.exp = T)
    }
    return(derivees)
}

derivees <- list()
temps_deriv <- numeric(2)
for (n in 1:(2)){
    # Calcul des dérivées
    temps_deriv[n] <- system.time(
        derivees <- append(derivees, parse(text = 
                                               chain_derivative(func_str_copule(str_copule_ext, str_copule_int, n), n)))
    )[[3]]
    print(c("Dérivée"=n, "temps"=temps_deriv[n]))
}

generateur_evalue_deriv <- function(derivee){
    # Générateur permettant d'utiliser la dérivée à titre de fonction évaluable
    # à partir de la liste de dérivées saisie en argument.
    function(xx, para){
        # Fonction générée à l'aide des dérivées.
        for (i in 1:length(para)) {
            assign(paste0("pa", i), para[i])
        }
        for (i in 1:length(xx)) {
            assign(paste0("x", i), xx[i])
        }
        return(eval(derivee[[length(xx)]]))
    }
}
densite <- generateur_evalue_deriv(derivees)

dCopule <- function(xx, para){
    # Fonction de densité de la copule
    return(densite(xx, para))
}



fct_Score_XX <- function(para, Data=mat_XX) {
    # Fonction de score: On cherchera à minimiser la log-vraisemblance négative
    return(-sum(log(apply(Data, 1, dCopule, para))))
}


temps_theta <- system.time({
    mat_XX <- constr_mat_XX(DATA_train)
    plot(domaine <- 6:9, sapply(domaine, fct_Score_XX), type="l")
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
})[[3]]

arr_mle_theta <- numeric(nsim <- 30)
for (iter in 1:nsim) {
  temps_theta <- temps_theta + system.time({
    
    mat_XX <- constr_mat_XX(DATA_train)
    arr_mle_theta[iter] <- optimize(fct_Score_XX, interval = c(6, 8))$minimum
    
  })[[3]]
  
  print(c(paste0(iter, "/", nsim),
          paste0(round(temps_theta %/% 60),":", round(temps_theta %% 60)),
          arr_mle_theta[iter]))
}

# sapply(1:3, function(i) hist(arr_mle_M.N.X[, i]))
(mle_theta <- mean(arr_mle_theta))
var_theta <- var(arr_mle_theta)
IC_theta <- matrix(
  c(mle_theta - qnorm(0.975) * var_theta,
    mle_theta + qnorm(0.975) * var_theta),
  ncol=2)


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


mat_NX <- constr_mat_NX(DATA_train)
fct_Score_NX <- function(para, Data=mat_NX) {
    # Fonction de score: On cherchera à minimiser la log-vraisemblance négative
    neg_log_vrais <- 0
    for (i in 1:nrow(Data)){
        neg_log_vrais <- neg_log_vrais - log(f_NX(Data[i, 1], Data[i, 2], para))
    }
    return(neg_log_vrais)
}


# Estimation alpha0 par la méthode des moments.
set.seed(20190618)
(tau_0  <- tau_NX(DATA_train, 30)$moyenne)

temps_M_moments <- system.time({
    
  Tau_graph <- sapply(domaine <- 5:10, function(a)
      tau_kendall_theorique(pbinom,
                            list(N = list(size=nb_xi, prob=mle_binom)),
                            Copule0, a, nb_xi))
  
  plot(domaine, Tau_graph,
       ylab="tau de kendall",
       xlab="alpha",
       type="l")
  axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
  axis(1, at=domaine, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
  abline(a=tau_0, b=0, col="green")

  bornes <- c(7.5, 8.5)
  
  alpha0_n <- inversion_tau_kendall(pbinom,
                                    list(N=list(size=nb_xi, prob=mle_binom)),
                                    Copule0, bornes, tau_0, nb_xi)
})[[3]]


# Valeurs de départ de l'optimisation
(val_depart <- c(
    "prob" = mean(DATA_train[, 1]) / nb_xi,
    "beta" = mle_exp,
    "alpha0" = alpha0_n
))

# Bornes de l'optimisation
bounds <- matrix(c(1e-6, 1 - 1e-6,
                   1/summary_X[[5]], 1/summary_X[[2]],
                   1e-6, 10), ncol = 2, byrow = T)
# Convertir les constraintes en matrices ui et ci.
n_par <- nrow(bounds)
ui <- rbind( diag(n_par), -diag(n_par) )
ci <- c( bounds[,1], - bounds[,2] )
# Retirer les valeurs infinies
i <- as.vector(is.finite(bounds))
ui <- ui[i,]
ci <- ci[i]

temps_M.N.X <- 0
arr_mle_M.N.X <- matrix(ncol=3, nrow=nsim <- 30)
for (iter in 1:nsim) {
    temps_M.N.X <- temps_M.N.X + system.time({
        
        mat_NX <- constr_mat_NX(DATA_train)
        arr_mle_M.N.X[iter,] <- constrOptim(val_depart, fct_Score_NX, grad=NULL,
                                            ui=ui, ci=ci, outer.eps = 1e-3)$par
        
    })[[3]]
    
    print(c(paste0(iter, "/", nsim),
            paste0(round(temps_M.N.X %/% 60),":", round(temps_M.N.X %% 60)),
            arr_mle_M.N.X[iter,]))
}

# sapply(1:3, function(i) hist(arr_mle_M.N.X[, i]))
mle_M.N.X <- apply(arr_mle_M.N.X, 2, mean)
var_M.N.X <- apply(arr_mle_M.N.X, 2, var)
IC_M.N.X <- matrix(
    c(mle_M.N.X - qnorm(0.975) * var_M.N.X,
      mle_M.N.X + qnorm(0.975) * var_M.N.X),
    ncol=2)


# ====================== Estimation du paramètre des v.a. X et theta ==========================
mle_alpha0 <- mle_M.N.X[3]

F_X <- "1 - exp(-x_i * pa1)"
LST.Log_M <- "-1 / mle_alpha0 * log(1 - (1 - exp(-mle_alpha0)) * exp(-T)) "
LST.Log_M.inv <- "- log((1 - exp(-mle_alpha0 * U)) / (1 - exp(-mle_alpha0)))"
LST.Log_B <- "-1 / pa2 * log(1 - (1 - exp(-pa2)) * exp(-T)) "
LST.Log_B.inv <- "log((1 - exp(-pa2 * U)) / (1 - exp(-pa2)))"

str_copule_ext <- str_replace(LST.Log_M, "T",
                              paste0("-log(", LST.Log_B,")"))
str_copule_int <- str_replace(LST.Log_B.inv, "U",
                              paste0("exp(-",LST.Log_M.inv,")"))


func_str_copule <- function(str_copule_ext, str_copule_int, nb_xi) {
    # Fonction qui génère la chaîne de caractères qui permettra d'effectuer les dérivées.
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
    # Fonction qui effectue les dérivations en chaîne sur le texte donné en argument
    X_i <- sapply(1:nb_xi, function(i) paste0("x",i))
    derivees <- str_copule
    
    for (xi in X_i){
        derivees <- Deriv(derivees, xi, cache.exp = T)
    }
    return(derivees)
}

derivees <- list()
temps_deriv <- numeric(2)
for (n in 1:(2)){
    # Calcul des dérivées
    temps_deriv[n] <- system.time(
        derivees <- append(derivees, parse(text = 
            chain_derivative(func_str_copule(str_copule_ext, str_copule_int, n), n)))
    )[[3]]
    print(c("Dérivée"=n, "temps"=temps_deriv[n]))
}

generateur_evalue_deriv <- function(derivee){
    # Générateur permettant d'utiliser la dérivée à titre de fonction évaluable
    # à partir de la liste de dérivées saisie en argument.
    function(xx, para) {
        # Fonction générée à l'aide des dérivées.
        for (i in 1:2) {
            assign(paste0("pa", i), para[i])
            assign(paste0("x", i), xx[i])
        }
        return(eval(derivee[[2]]))
    }
}
dCopule <- generateur_evalue_deriv(derivees)


mat_XX <- constr_mat_XX(DATA_train)
fct_Score_XX <- function(para, Data=mat_XX) {
    # Fonction de score: On cherchera à minimiser la log-vraisemblance négative
    return(-sum(log(apply(Data, 1, dCopule, para))))
}


# Estimation de alpha1 par la méthode des moments.
Expr <- ConstrCopule("log", "log")
set.seed(20190618)
tau_1 <- tau_XX(DATA_train, nsim=100)$moyenne


# list_Tau <- sapply(domaine <- 1:10, function(a0)
#   Invert_tau_knd.continu(
#     tau_XX(simul_modele.collectif(nb_xi, mle_binom, beta, a0, 6, 1e+4), nsim=100)$moyenne,
#     a0, c(1e-6,7), Expr)
#   )
# 
# plot(domaine, list_Tau,
#      ylab="alpha1",
#      xlab="alpha0",
#      type="l"
# )
# axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
# axis(1, at=domaine, tck=1, lty = 2, col = "grey", labels = F) # L'axe des abscisses
# abline(a=6, b=0, col="green")


temps_theta_moments <- system.time({
    
    Tau_graph <- sapply(domaine <- (1:10), function(a)
      tau_knd.continu(alpha0_n, a, Expr))
    
    plot(domaine, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=domaine, tck=1, lty = 2, col = "grey", labels = F) # L'axe des abscisses
    abline(a=tau_1, b=0, col="green")
    
    bornes <- c(3, 5)
    
    alpha1_n <- Invert_tau_knd.continu(tau_1, alpha0_n, bornes, Expr)
})[[3]]

    
# Valeurs de départ de l'optimisation
(val_depart <- c(
    mle_M.N.X[2],
    "alpha1" = alpha1_n
))

# Bornes de l'optimisation
bounds <- matrix(c(1 / summary_X[[5]], 1 / summary_X[[2]],
                   1, 10),
                 ncol = 2,
                 byrow = T)
# Convertir les constraintes en matrices ui et ci.
n_par <- nrow(bounds)
ui <- rbind(diag(n_par), -diag(n_par))
ci <- c(bounds[,1], - bounds[,2])
# Retirer les valeurs infinies
i <- as.vector(is.finite(bounds))
ui <- ui[i,]
ci <- ci[i]


temps_X.theta <- 0
arr_mle_X.theta <- matrix(ncol=2, nrow=nsim <- 30)
for (iter in 1:nsim) {
  temps_X.theta <- temps_X.theta + system.time({
    
    mat_XX <- constr_mat_XX(DATA_train)
    arr_mle_X.theta[iter,] <- constrOptim(val_depart, fct_Score_XX, grad=NULL,
                                          ui=ui, ci=ci, outer.eps = 1e-03)$par
    
  })[[3]]
  
  print(c(paste0(iter, "/", nsim),
          paste0(round(temps_X.theta %/% 60),":", round(temps_X.theta %% 60)),
          arr_mle_X.theta[iter,]))
}

mle_X.theta <- apply(arr_mle_X.theta, 2, mean)
var_X.theta <- apply(arr_mle_X.theta, 2, var)
IC_X.theta <- matrix(
  c(mle_X.theta - qnorm(0.975) * var_X.theta,
    mle_X.theta + qnorm(0.975) * var_X.theta),
  ncol=2)


# Résultats --------------------------------------------------------------------------

tbl_resultats <- rbind(
    "Vrais paramètres" = c(q, beta, alpha0, alpha1),
    "Méthode IFM" = c(mle_binom, mle_exp, mle_M, mle_theta),
    "Méthode composite" = c(mle_M.N.X[1], mle_M.N.X[2], mle_M.N.X[3], mle_X.theta[2]),
    "Méthode des moments" = c(
        mean(DATA_train[, 1]) / nb_xi,
        1 / mean(DATA_train[, -1], na.rm = T),
        alpha0_n, alpha1_n
    )
)
tbl_resultats <- cbind(tbl_resultats,
      c(
          NaN,
          temps_N + temps_X + temps_M + temps_theta,
          temps_M.N.X + temps_X.theta,
          temps_M_moments + temps_theta_moments
      )
)
colnames(tbl_resultats) <- c("q", "beta", "alpha0", "alpha1", "temps")

tbl_resultats

xtable(tbl_resultats, digits=4)

