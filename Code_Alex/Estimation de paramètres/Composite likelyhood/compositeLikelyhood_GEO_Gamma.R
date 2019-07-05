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

lambda <- 2
nb_xi <- qpois(1 - 1e-06, lambda)

alpha <- 0.8
beta <- 1/100

alpha0 <- 0.8
alpha1 <- 7

nsim <- n_obs <- 1e+4

DATA_train <- rCompCop(nsim, GEO(1 - alpha0, 1, list(GAMMA(1 / alpha1, 1:nb_xi, NULL))))

for (i in 1:nsim){
    DATA_train[i, 1] <- N <- qpois(DATA_train[i,1], lambda)
    if (N==0){ 
        DATA_train[i,-1] <- rep(NaN, nb_xi)
    }else{
        DATA_train[i,-1] <- c(qgamma(DATA_train[i, 2:(N + 1)], alpha, beta), rep(NaN, nb_xi - N))}
}

colnames(DATA_train) <- c("N", sapply(1:(ncol(DATA_train)-1),function(i) paste0("X", i)))
head(DATA_train)
summary(DATA_train)
DATA_train <- DATA_train[,1:11]
summary(DATA_train)


nb_xi <- max(DATA_train[,1])
n_obs <- nrow(DATA_train)


# ====================== Estimation des paramètres des v.a. N et X ======================
# Estimation des paramètres de N
summary_N <- summary(DATA_train[,1])
temps_N <- system.time(
    mle_N <- optimize(function(par) -sum(log(dpois(DATA_train[, 1], par))),
                           lower = summary_N[[2]], upper = summary_N[[5]])$minimum
)[[3]]
mle_N; lambda
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
)rbind(mle_X$par,
      c(alpha, beta))
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