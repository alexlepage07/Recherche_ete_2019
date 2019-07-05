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


# ====================== Simulations des données d'entraînement ================
# Paramètres de simulation
set.seed(20190702)
nb_xi <- 10
q <- 1/5
beta <- 1/100
alpha0 <- 0.6
alpha1 <- 0.4
nsim <- n_obs <- 1e+4

DATA_train <- rCompCop(nsim, GEO(1-alpha0, 1, list(GEO(1-alpha1, 1:nb_xi, NULL))))

for (i in 1:nsim){
    DATA_train[i, 1] <- N <- qbinom(DATA_train[i,1], nb_xi, q)
    if (N==0){ 
        DATA_train[i,-1] <- rep(NaN, nb_xi)
    }else{
        DATA_train[i,-1] <- c(qexp(DATA_train[i, 2:(N + 1)], beta), rep(NaN, nb_xi - N))}
}

colnames(DATA_train) <- c("N", sapply(1:(ncol(DATA_train)-1),function(i) paste0("X", i)))
head(DATA_train)
# write.csv(DATA_train, file = "Data_LOG-LOG")

# Sommaire des résultats de simulation
(sommaire_train <- summary(DATA_train))
# xtable(sommaire_train)  # Permet de  convertir un tableau de R à LaTeX.

# ====================== Estimation des paramètres des v.a. N et X ======================
# Estimation des paramètres de N
temps_N <- system.time(
    mle_binom <-  optimize(function(par) -sum(log(dbinom(DATA_train[, 1], nb_xi, par))),
                           lower = 1e-6, upper = 1-1e-6)$minimum
)[[3]]
mle_binom; q
temps_N

# Estimation des paramètres de X
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
Copule0 <- amhCopula
para <- c(q, beta, alpha0)

F_N <- function(x0) pbinom(x0, nb_xi, mle_binom)
F_X <- "(1 - exp(-x * mle_exp))"

str_copule <- "u1 * u2 / (1 - pa1 * (1-u1) * (1-u2))"
str_copule <- str_replace_all(str_copule, "u1", "F_N(x0)")
str_copule <- str_replace_all(str_copule, "u2", F_X)

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
    
    plot(domaine <- (1:10)/10, sapply(domaine, fct_Score_NX), type="l")
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=domaine, tck=1, lty = 2, col = "grey",) # L'axe des abscisses


    mle_M <- optimize(fct_Score_NX, interval = c(0.55, 0.75))
})[[3]]

# ====================== Estimation du paramètre de la v.a. theta ==========================
mle_alpha0 <- mle_M$minimum

F_X <- "1 - exp(-x_i * mle_exp)"

LST.M <- "(1 - mle_alpha0) / (exp(T) - mle_alpha0)"
LST.M.inv <- "log((1 - mle_alpha0) / U + mle_alpha0)"
LST.B <- "(1 - pa1) / (exp(T) - pa1)"
LST.B.inv <- "log((1 - pa1) / U + pa1)"

str_copule_ext <- str_replace(LST.M, "T",
                              paste0("-log(", LST.B,")"))
str_copule_int <- str_replace(LST.B.inv, "U",
                              paste0("exp(-",LST.M.inv,")"))


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
func_str_copule(str_copule_ext, str_copule_int, 2)

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


mat_XX <- combn2(na.omit(DATA_train[1,-1]))
for (i in 2:nrow(DATA_train)) {
    N <- DATA_train[i, 1]
    if (N < 2)
        next
    mat_XX <- rbind(mat_XX, combn2(DATA_train[i, 2:(N + 1)]))
}


fct_Score_XX <- function(para, Data=mat_XX) {
    # Fonction de score: On cherchera à minimiser la log-vraisemblance négative
    return(-sum(log(apply(Data, 1, dCopule, para))))
}

temps_theta <- system.time({
    
    plot(domaine <- (4:8)/10, sapply(domaine, fct_Score_XX), type="l")
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=domaine, tck=1, lty = 2, col = "grey",) # L'axe des abscisses


    mle_theta <- optimize(fct_Score_XX, interval = c(0.55, 0.75))
})


# ====================== Estimation des paramètres des v.a. M, N et X. ==========================
Copule0 <- amhCopula
para <- c(q, beta, alpha0)

F_N <- function(x0, pa1) pbinom(x0, nb_xi, pa1)
F_X <- "(1 - exp(-x * pa2))"

str_copule <- "u1 * u2 / (1 - pa1 * (1-u1) * (1-u2))"
str_copule <- str_replace_all(str_copule, "u1", "F_N(x0, pa1)")
str_copule <- str_replace_all(str_copule, "u2", F_X)

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
    
    temps_M_moments <- system.time({
        
        Tau_graph <- sapply(domaine <- (1:5)/5, function(a)
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
        
        
        bornes <- c(8, 8.5)
        
        alpha0_n <- inversion_tau_kendall(pbinom,
                                          list(N=list(size=nb_xi, prob=mle_binom)),
                                          Copule0, bornes, tau_0, nb_xi)
    })[[3]]
} # Estimation alpha0 par la méthode des moments.

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

temps_M.N.X <- system.time(
    # Estimation des paramètres
    mle_M.N.X <- constrOptim(val_depart, fct_Score_NX, grad=NULL,
                             ui=ui, ci=ci)
)[[3]]


# ====================== Estimation du paramètre des v.a. X et theta ==========================
# Méthode d'estimation par décomposition de la densité conjointe.
mle_alpha0 <- mle_M.N.X$par[3]

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


mat_XX <- combn2(na.omit(DATA_train[1,-1]))
for (i in 2:nrow(DATA_train)) {
    N <- DATA_train[i, 1]
    if (N < 2)
        next
    mat_XX <- rbind(mat_XX, combn2(DATA_train[i, 2:(N + 1)]))
}


fct_Score_XX <- function(para, Data=mat_XX) {
    # Fonction de score: On cherchera à minimiser la log-vraisemblance négative
    return(-sum(log(apply(Data, 1, dCopule, para))))
}


{
set.seed(20190618)
mean_tau <- mean(tau_XX(DATA_train, nsim=30, silent=T), na.rm=T)

struct_copule <- function(alpha0, alpha1) {
    LOG(1 - exp(-alpha0), NULL, list(LOG(1 - exp(-alpha1), 1:2, NULL)))}

temps_theta_moments <- system.time({
    
    Tau_graph <- sapply(domaine <- (1:5)*2, function(a)
        tau_kendall_theorique_continues(struct_copule,
                                        mle_alpha0, alpha1 = a))
    
    plot(domaine, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=domaine, tck=1, lty = 2, col = "grey", labels = F) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
    
    bornes <- c(5, 6)
    
    alpha1_n <- inversion_tau_kendall_XX(struct_copule, mle_alpha0, mean_tau, bornes)
})

} # Estimation de alpha1 par la méthode des moments.

# Valeurs de départ de l'optimisation
(val_depart <- c(
    mle_M.N.X$par[2],
    "alpha1" = alpha1_n
))

# Bornes de l'optimisation
bounds <- matrix(c(1e-06, 1 - 1e-06,
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

temps_X.theta <- system.time(
    # Estimation des paramètres
    mle_X.theta <- constrOptim(val_depart, fct_Score_XX, grad=NULL,
                          ui=ui, ci=ci, outer.eps = 1e-03)
)
mle_M.N.X$par[2]; mle_X.theta[1]

# Résultats --------------------------------------------------------------------------
{
    tbl_resultats <- rbind(
        "Vrais paramètres" = c(q, beta, alpha0, alpha1),
        "Méthode IFM" = c(mle_binom, mle_exp, mle_M$minimum, mle_theta$minimum),
        "méthode composite" = c(mle_M.N.X$par[1], mle_M.N.X$par[2], mle_M.N.X$par[3], mle_),
        "méthode des moments" = c(
            mean(DATA_train[, 1]) / nb_xi,
            1 / mean(DATA_train[, -1], na.rm = T),
            alpha0_n,
            alpha1_n
        )
    )
    cbind(tbl_resultats,
          c(
            "s.o.",
            temps_N + temps_X + temps_M + temps_theta,
            temps_M.N.X + temps_X.theta,
            temps_M_moments + temps_theta_moments
          )
    )
    colnames(tbl_resultats) <- c(q, beta, alpha0, alpha1)
} # Production du tableau des résultats

tbl_resultats

xtable(tbl_resultats, digits=4)

