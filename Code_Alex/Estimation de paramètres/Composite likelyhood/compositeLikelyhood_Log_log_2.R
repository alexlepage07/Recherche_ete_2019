n<- 15000
alpha0ll<-4
alpha1ll<-7
q <- 0.4
m <- 7
beta<- 1/100

ULOGLOG<- rCompCop(n,LOG(1-exp(-alpha0ll), 1, list(LOG(1-exp(-alpha1ll), 1:m, NULL))))

mat<- matrix(NA,n,m+1)
mat[,1]<- qbinom(ULOGLOG[,1],m,q)
for (i in 1:n) {
  count <- mat[i,1]
  if (count !=0) mat[i,2:(count+1)]<- qexp(ULOGLOG[i,2:(count+1)],beta)
}
mat

library(copula)
library(nCopula)
library(Deriv)
library(stringr)
library(xtable)
library(ggplot2)
library(psych)
library(combinat)

source("../../Mesure_dependance va mixtes.R")

DATA_train <- mat
nb_xi <- max(DATA_train[,1])
n_obs <- nrow(DATA_train)
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

# lst_XX <- array(DATA_train[,-1])
# lst_XX <- lst_XX[!is.na(lst_XX)]
# summary_X <- summary(lst_XX)
# temps_X <- system.time(
#     (mle_exp <-  optimize(function(par) -sum(log(dexp(DATA_train[,-1], par))), 
#                           interval = c(1/summary_X[[5]], 1/summary_X[[2]]))$minimum)
# )[[3]]
# mle_exp; beta
# temps_X

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


datas <- data.frame(c(arr_X, qexp((0:100)/100, mle_exp)),
                    Source <- c(rep("Simul", length(arr_X)),
                                rep("théorique", 101)))

ggplot() + 
    geom_histogram(alpha = 0.3, aes(x= arr_X, y = ..density.., fill = "Empirique"), position = 'identity')+
    geom_density(alpha = 0.3, aes(x= qexp((0:100)/100, beta), y = ..density.., fill = "Théorique")) + 
    xlab("x") + ylab("Densité") +
    theme(legend.title = element_blank())

ks.test(unique(arr_X), pexp, mle_exp)

# ====================== Estimation du paramètre de la v.a. M. ============================
Copule0 <- frankCopula
para <- c(q, beta, alpha0 <- alpha0ll)

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
    
    
    mle_M <- optimize(fct_Score_NX, interval = c(3, 5))
})[[3]]


# ====================== Estimation du paramètre de la v.a. theta ==========================
mle_alpha0 <- mle_M$minimum

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
    
    plot(domaine <- (1:5)*2, sapply(domaine, fct_Score_XX), type="l")
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
    
    
    mle_theta <- optimize(fct_Score_XX, interval = c(6, 10))
})
