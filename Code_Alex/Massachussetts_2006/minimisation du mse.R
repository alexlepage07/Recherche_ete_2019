library(pscl)
library(statmod)
library(tweedie)
library(copula)


#==== Fonctions pour évaluer les paramètres de dépendance de la copule gaussienne. ====
# avec z = 1 : Equicorrelation matrix 
Sig_rho2.z1 <- function(rho, k) {
    # Fonction qui construit la matrice de corrélation de rho2 avec z=1
    matrix(c(rep(c(1, rep(rho[2], k)), k-1),1), nrow = k, byrow = T)
}
Sig_rho2.z1.inv <- function(rho, k) {
    (diag(k) - rho[2] / (1 + (k-1) * rho[2]) * matrix(1, k, k)) / (1 - rho[2])
}
Sig_rho1.2.z1 <- function(rho, k) {
    # Fonction qui construit la matrice de corrélation de rho1 et rho2 avec z=1
    mat <- rbind(rep(rho[1], k), Sig_rho2.z1(rho, k))
    mat <- cbind(c(1, rep(rho[1], k)), mat)
    return(mat)
}
Sig_rho1.2.z1.inv <- function(rho, k) {
    denom <- 1 + (k - 1) * rho[2] - k * rho[1] ^ 2
    mat <- matrix((diag(k) - (rho[2] - rho[1] ^ 2) / denom * matrix(1, k, k)) / (1 - rho[2]), k)
    mat <- rbind(-rho[1] / denom * t(rep(1, k)) , mat)
    mat <- cbind(c(1 + rho[1] ^ 2 * k / denom,-rho[1] / denom * rep(1, k)), mat)
    return(mat)
}

# Avec z = 2 : Autoregressive correlation matrix
Sig_rho2.z2 <- function(rho, k) {
    # Fonction qui construit la matrice de corrélation de rho2 avec z=2
    outer(0:(k-1), 0:(k-1), function(i,j) rho[2]^abs(i-j))
}
Sig_rho2.z2.inv <- function(rho, k) {
    # Fonction qui construit la matrice de corrélation inverse de rho2 avec z=2
    solve(Sig_rho2.z2(rho, k))
}
Sig_rho1.2.z2 <- function(rho, k) {
    # Fonction qui construit la matrice de corrélation de rho1 et rho2 avec z=2
    mat <- rbind(rep(rho[1], k), Sig_rho2.z2(rho, k))
    mat <- cbind(c(1, rep(rho[1], k)), mat)
    return(mat)
}
Sig_rho1.2.z2.inv <- function(rho, k) {
    # Fonction qui construit la matrice de corrélation inverse de rho1 et rho2 avec z=2
    solve(Sig_rho1.2.z2(rho, k))
}

# Paramètre de la copule gaussienne
mu.norm <- function(U, rho, z=1) {
    # Fonction qui calcule le paramètre de moyenne de la copule gaussienne
    k <- length(U)
    Sig_rho.inv <- ifelse(z == 1, Sig_rho2.z1.inv, Sig_rho2.z2.inv)
    
    rho[1] * t(rep(1, k)) %*% Sig_rho.inv(rho, k) %*% (qnorm(U))
}
sig.norm <- function(U, rho, z=1) {
    # Fonction qui calcule le paramètre d'écart-type de la copule gaussienne
    k <- length(U)
    Sig_rho.inv <- ifelse(z == 1, Sig_rho2.z1.inv, Sig_rho2.z2.inv)
    
    sqrt(1 - (rho[1] * t(rep(1, k))) %*% Sig_rho.inv(rho, k) %*% (rho[1] * rep(1, k)))
}

dGaussian.copula <- function(n, Y, x, rho, z=1) {
    # Fonction qui calcule la densité d'une copule gaussienne
    if(typeof(x) == "list") x <- as.integer(unlist(x))
    if(typeof(Y) == "list") Y <- as.double(unlist(Y))
    
    p_R.1 <- prob(x)
    if (n == 0) return(1 - p_R.1)
    if (n > 7) return(1)
    
    U <- F2(Y,x)
    k <- length(Y)
    Sig <- sig.norm(U, rho, z)
    mu <- mu.norm(U, rho, z)
    
    Sig_rho <- ifelse(z == 1, Sig_rho2.z1, Sig_rho2.z2)
    Sig_rho.inv <- ifelse(z == 1, Sig_rho2.z1.inv, Sig_rho2.z2.inv)
    det.Sig_rho <- ifelse (z == 1,
                           (1 + (k - 1) * rho[2]) * (1 - rho[2]) ^ (k - 1),
                           (1 - rho[2] ^ 2) ^ k)
    
    c_Y <- det.Sig_rho ^ (-0.5) * exp(- t(qnorm(U)) %*%
                                          (Sig_rho.inv(rho, k) - diag(1, k)) %*%
                                          qnorm(U) / 2)
    
    c_Y <- c_Y * prod(f2(Y, x))
    
    d <- p_R.1 * c_Y * (pnorm((qnorm(F1(n, x)) - mu) / Sig) -
                            pnorm((qnorm(F1(n-1, x)) - mu) / Sig))
    
    # print(paste(n, d))
    return(d)
}
rGaussian.copula <- function(n_sim, rho, k = max(FREQ$N), z = 1) {
    # Fonction qui permet de simuler des réalisations d'une copule gaussienne.
    Sig_rho <- ifelse(z == 1, Sig_rho1.2.z1, Sig_rho1.2.z2)
    rCopula(n_sim, normalCopula(P2p(Sig_rho(rho, k)), k + 1, dispstr = "un"))
}

VaR <- function(x){
    # Valeur au risque empirique avec un seuil de 0.005
    sort(x)[ceiling(length(x) * 0.995)]
}
TVaR <- function(x) {
    # TVaR empirique avec un seuil de 0.005
    mean(x[x > VaR(x)])

}

simul.S <- function(n_sim, x, rho, k = max(FREQ.train$N), z = 1) {
    # Fonction qui permet de faire de la simulation Monte-Carlo
    # pour estimer les mesures de risque théorique du modèle.
    U <- rGaussian.copula(n_sim, rho, k, z)
    
    Couples <- cbind(
        "N" = sapply(1:n_sim, function(i)
            F1.inv(U[i, 1], unlist(x))),
        "Y" = matrix(
            sapply(1:n_sim, function(i)
                qgamma(U[i, -1], 1 / Y.disp, 1 / exp(unlist(x) %*% Y.coef) / Y.disp)),
            nrow = n_sim, byrow = T
        ))
    
    S <- numeric(n_sim)
    for (i in 1:n_sim) {
        N <- Couples[i,1]
        if (N == 0)
            S[i] <- 0
        else
            S[i] <- sum(Couples[i, 2:(N + 1)])
    }

    return(list(
        "E[S]" = mean(S),
        "VaR_0.995" = VaR(S),
        "TVaR_0.995" = TVaR(S)
    ))
}

randomize.Y <- function(data, n) {
    # Fonction qui prend n montants de sinistres aléatoirement
    # parmi ceux disponibles pour chacune des polices.
    data <- data[data[, 1] >= n, ]
    
    Y <- matrix(sapply(1:nrow(data), function(i)
        sample(unlist(data[, -1][[i]]), n)),
        ncol = n,
        byrow = T)
    
    cbind(unlist(data[, 1]), Y)
}
spearmans.Rho <- function(data, nsim=30) {
    # Fonction qui retourne le rho de spearman pour une base de données d'assurance.
    cor.spearman <- matrix(NA, nsim, 2)
    for (i in 1:nsim) {
        couples.YY <- randomize.Y(data, 2)
        cor.spearman[i,1] <- cor(couples.YY[,1], couples.YY[,2], method = "spearman")
        cor.spearman[i,2] <- cor(couples.YY[,2], couples.YY[,3], method = "spearman")
        print(paste0(i, "/", nsim))
    }
    res <- list()
    res$mean <- apply(cor.spearman, 2, mean)
    res$variance <- apply(cor.spearman, 2, var)
    res$IC <- matrix(c(res$mean - qnorm(0.975) * sqrt(res$variance),
                res$mean + qnorm(0.975) * sqrt(res$variance)), 2)
    return(res)
}


#==== Importation du jeu de données ====
setwd("C:/Users/Alex/Desktop/Recherche_ete_2019/Code_Alex/Massachussetts_2006")
# Réarrangement des données en facteurs
AGG <- read.csv("COVARIATES.txt")
AGG$class4 <- AGG$class4 %% 10
AGG$class4[AGG$class4 == 1] <- "A"
AGG$class4[AGG$class4 == 2] <- "S"
AGG$class4[AGG$class4 %in% c(3, 4)] <- "M"
AGG$class4[AGG$class4 == 5] <- "B"
AGG$class4[AGG$class4 %in% (6:9)] <- "I"
AGG$class4 <- factor(AGG$class4)

AGG$tgroup <- factor(AGG$tgroup, 1:6)
AGG$pol_id <- factor(AGG$pol_id)
AGG$clm_id <- factor(AGG$clm_id)
AGG$losspaid[is.na(AGG$losspaid)] <- 0
summary(AGG)
head(AGG)
# Réorganiser les données pour faciliter le travail.
train.id <- 1:1e+6

AGG.train <- AGG[train.id,]
AGG.test <- AGG[-train.id,]

# Données de sinistres
SEV.train <- AGG.train[AGG.train$losspaid > 0,]
SEV.train <- aggregate(SEV.train$losspaid,
                       by = list(SEV.train$pol_id,
                                 SEV.train$class4,
                                 SEV.train$tgroup),
                       list)
colnames(SEV.train) <- c("pol_id", "class4", "tgroup", "losspaid")

SEV.test <- AGG.test[AGG.test$losspaid > 0,]
SEV.test <- aggregate(AGG.test$losspaid,
                       by = list(AGG.test$pol_id,
                                 AGG.test$class4,
                                 AGG.test$tgroup),
                       list)
colnames(SEV.test) <- c("pol_id", "class4", "tgroup", "losspaid")


FREQ.train <- aggregate(AGG.train$losspaid,
                  by = list(AGG.train$pol_id,
                            AGG.train$class4,
                            AGG.train$tgroup),
                  function(x) sum(x>0))
colnames(FREQ.train) <- c("pol_id", "class4", "tgroup", "N")

FREQ.test <- aggregate(AGG.test$losspaid,
                        by = list(AGG.test$pol_id,
                                  AGG.test$class4,
                                  AGG.test$tgroup),
                        function(x) sum(x>0))
colnames(FREQ.test) <- c("pol_id", "class4", "tgroup", "N")


SUM.train <- aggregate(AGG.train$losspaid,
                       by = list(AGG.train$pol_id, AGG.train$tgroup, AGG.train$class4),
                       sum)
colnames(SUM.train) <- c("pol_id", "class4", "tgroup", "S")

SUM.test <- aggregate(AGG.test$losspaid,
                       by = list(AGG.test$pol_id, AGG.test$tgroup, AGG.test$class4),
                       sum)
colnames(SUM.test) <- c("pol_id", "class4", "tgroup", "S")

#==== Paramétrisation des marginales ====
#---- Fréquence
N.mod <- hurdle(N ~ class4 + tgroup, data = FREQ.train)
N.coef <- N.mod$coefficients$count
R.coef <- N.mod$coefficients$zero

prob <- function(x) {
    # Fonction qui calcule la probabilité d'avoir un accident.
    p <- exp(x %*% R.coef)
    return(p / (1 + p))
}
F1 <- function(n, x) {
    # Fonction qui calcule la c.d.f de la v.a. de fréquence.
    p_N.0 <- dpois(0, exp(x %*% N.coef))
    if (n < 0) return(0)
    p_R.1 <- prob(x)
    (1 - p_R.1) + p_R.1 * (ppois(n, exp((x) %*% (N.coef))) - p_N.0) / (1 - p_N.0)
}
F1.inv <- function(p, x) {
    p.0 <- 1 - prob(x)
    if(p <= p.0) return(0)
    floor(uniroot(function(n) F1(n, x) - p, c(0, 50))$root)
}

#---- Sévérité
EstimParamSeverite <- function(nsim) {
    # Fonction qui estime les paramètres de la v.a. des
    # montants de sinistre en sélectionnant un montant
    # par police d'assurance aléatoirement et en applicant
    # du rééchantillonnage.
    
    Y.disp_ <- numeric(nsim)
    Y.coef_ <- matrix(NA, nsim, 10)
    
    for (sim in 1:nsim){
        print(paste0(sim, "/", nsim))
        
        rand.Y <- randomize.Y(cbind(FREQ.train$N[FREQ.train$N > 0], SEV.train$losspaid), 1)
        gamma.data <- cbind(SEV.train[, 1:3], rand.Y[, 2])
        colnames(gamma.data) <- colnames(SEV.train)
        
        gamma.model <- glm(losspaid ~ class4 + tgroup, data = gamma.data,
                           family = Gamma(link = "log"))
        
        Y.disp_[sim] <- MASS::gamma.dispersion(gamma.model)
        Y.coef_[sim,] <- gamma.model$coefficients
    }
    
    Y.disp <- mean(Y.disp_, na.rm=T)
    Y.coef <- apply(Y.coef_, 2, mean, na.rm=T)
    return(list("dispersion"=Y.disp, "coefficients"=Y.coef))
}
Y.param <- EstimParamSeverite(30)

Y.disp <- Y.param$dispersion
Y.coef <- Y.param$coefficients

F2 <- function(Y, x) {
    pgamma(Y, 1 / Y.disp, exp(-x %*% Y.coef) / Y.disp)
}
f2 <- function(Y, x) {
    dgamma(Y, 1 / Y.disp, exp(-x %*% Y.coef) / Y.disp)
}
F2.inv <- function(p, x) {
    qgamma(p, 1 / Y.disp, exp(-x %*% Y.coef) / Y.disp)
}


#==== Paramétrisation de la copule ====
AGG.SUM.train_mean <- aggregate(SUM.train$S, list(SUM.train$class4, SUM.train$tgroup), mean)
AGG.SUM.train_VaR <- aggregate(SUM.train$S, list(SUM.train$class4, SUM.train$tgroup), VaR)
AGG.SUM.train_TVaR <- aggregate(SUM.train$S, list(SUM.train$class4, SUM.train$tgroup), TVaR)
AGG.SUM.train <- cbind(AGG.SUM.train_mean, AGG.SUM.train_VaR$x, AGG.SUM.train_TVaR$x)
colnames(AGG.SUM.train) <- c("tgroup", "class4", "mean.S", "VaR", "TVaR")

(rho_ <- spearmans.Rho(cbind(FREQ.train$N[FREQ.train$N > 0], SEV.train$losspaid), 100))
rho.IC <- rho_$IC
rho <- 2 * sin(pi * rho_$mean / 6)


#---- Par minimisation du MSE ----
x.train <- model.matrix(mean.S ~ class4 + tgroup, data = AGG.SUM.train)
fct_score <- function(rho) {
    # Fonction à minimiser pour estimer les paramètres de dépendance.
    MSE <- mean((sapply(1:30, function(i) 
        simul.S(1e+4, x.train[i,], rho)$'E[S]') - AGG.SUM.train$mean.S)^2)
    print(paste0("Rho = ",rho[1], ", ", rho[2], ", MSE = ", MSE))
    return(MSE)
}

EstimParamDependance <- function(){
    # Fonction qui estime les paramètres de dépendance en minimisant
    # L'erreur quadratique sur les valeurs agrégées.
    k <- max(FREQ.train$N)
    
    bounds <- matrix(c(-1, rho.IC[2], sqrt(((k - 1) * rho.IC[2] + 1) / k), rho.IC[4]), 2)
    ui <- rbind( diag(2), -diag(2) )
    ci <- c( bounds[,1], - bounds[,2] )
    
    constrOptim(rho, fct_score, grad = NULL, ui, ci,
                outer.eps = 1e-05, mu = 1e-03, outer.iterations = 100)
}

MinimizedSquareError.rho <- EstimParamDependance()
rho.optim <- MinimizedSquareError.rho$par
rho.optim <- c(0.116208594597012, 0.177315720293825) # MSE = 25439.5143346203
rho.optim <- rho


#---- Par maximisation de la vraisemblance ----
EstimParamDependance <- function(nsim=100){
    # Fonction qui permet de trouver les paramètres de dépendance
    # d'une copule gaussienne par optimisation numérique en utilisant un bootstrap.
    #
    # nsim correspond au nombre de simulations du bootstrap
    k <- max(FREQ.train$N)
    x.train <- model.matrix(N ~ class4 + tgroup, data = FREQ.train[FREQ.train$N > 0,])
    
    bounds <- matrix(c(-1, 0, sqrt(((k - 1) * rho[2] + 1) / k), 1), 2)
    ui <- rbind( diag(2), -diag(2) )
    ci <- c( bounds[,1], - bounds[,2] )
    
    fct_score <- function(rho) {
        -sum(log(sapply(sample(nrow(SEV.train), 1e+4), function(i)
            dGaussian.copula(
                n = length(SEV.train$losspaid[[i]]),
                Y = unlist(SEV.train$losspaid[[i]]),
                x = x.train[i,],
                rho = rho
            ))))
    }
    
    values <- sapply(1:nsim, function(sim) {
        constrOptim(rho, fct_score, grad = NULL,
                ui = ui, ci = ci, outer.eps = 1e-2)$par
        })
    
    Resultat <- list()
    Resultat$mean <- apply(values, 1, mean)
    Resultat$variance <- apply(values, 1, var)
    Resultat$IC <- matrix(c(Resultat$mean - qnorm(0.975) * sqrt(Resultat$variance), 
                     Resultat$mean + qnorm(0.975) * sqrt(Resultat$variance)), 2)
    return(Resultat)
}
(MLE.rho <- EstimParamDependance(100))
rho.optim <- MLE.rho$mean


#=== Modèle de tweedie ====
Tweedie.mod <- glm(S ~ class4 + tgroup, data = SUM.train,
                   family = tweedie(1.1, 0))
summary.Tweedie <- summary(Tweedie.mod)

mean_Tweedie.train <- aggregate(fitted(Tweedie.mod), list(SUM.train$class4, SUM.train$tgroup), mean)

mean_Tweedie.test <- predict(Tweedie.mod, newdata = SUM.test, type = "response")
mean_Tweedie.test <- aggregate(mean_Tweedie.test, list(SUM.test$class4, SUM.test$tgroup), mean)

VaR_Tweedie.train <- sapply(1:30, function(i)
    qtweedie(0.995, 1.1, mu = exp(Tweedie.mod$coefficients %*% x.train[i,]), phi = summary.Tweedie$dispersion))
TVaR_Tweedie.train <- sapply(1:30, function(i)
    TVaR(rtweedie(1e+5, 1.1, mu = exp(Tweedie.mod$coefficients %*% x.train[i,]), phi = summary.Tweedie$dispersion)))

#==== Résultats ====
x.train <- model.matrix(mean.S ~ class4 + tgroup, data = AGG.SUM.train)
Predictions <- sapply(1:30, function(i)
    simul.S(1e+4, x.train[i,], rho.optim))


RESULT.train <- cbind(
    "tgroup" = AGG.SUM.train$tgroup,
    "class4" = AGG.SUM.train$class4,
    "mean.S" = AGG.SUM.train$mean.S,
    "VaR_n" = AGG.SUM.train$VaR,
    "TVaR_n" = AGG.SUM.train$TVaR,
    "E[S]" = unlist(Predictions[1,]),
    "VaR_0.995" = unlist(Predictions[2,]),
    "TVaR_0.995" = unlist(Predictions[3,]),
    "E[S]_Tweedie" = Predict_Tweedie.train$x,
    "VaR_Tweedie" = VaR_Tweedie.train,
    "TVaR_Tweedie" = TVaR_Tweedie.train
)

AGG.SUM.test.mean <- aggregate(SUM.test$S, list(SUM.test$class4, SUM.test$tgroup), mean)
AGG.SUM.test.VaR <- aggregate(SUM.test$S, list(SUM.test$class4, SUM.test$tgroup), VaR)
AGG.SUM.test.TVaR <- aggregate(SUM.test$S, list(SUM.test$class4, SUM.test$tgroup), TVaR)
AGG.SUM.test <- cbind(AGG.SUM.test.mean, AGG.SUM.test.VaR$x, AGG.SUM.test.TVaR$x)
colnames(AGG.SUM.test) <- c("tgroup", "class4", "mean.S", "VaR", "TVaR")
RESULT.test <- cbind(
    "tgroup" = AGG.SUM.test$tgroup,
    "class4" = AGG.SUM.test$class4,
    "mean.S" = AGG.SUM.test$mean.S,
    "VaR_n" = AGG.SUM.test$VaR,
    "TVaR_n" = AGG.SUM.test$TVaR,
    "E[S]" = unlist(Predictions[1,]),
    "VaR_0.995" = unlist(Predictions[2,]),
    "TVaR_0.995" = unlist(Predictions[3,]),
    "E[S]_Tweedie" = Predict_Tweedie.test$x,
    "VaR_Tweedie" = VaR_Tweedie.train,
    "TVaR_Tweedie" = TVaR_Tweedie.train
)

RESULT.train
RESULT.test

MSE <- matrix(c(
    mean((RESULT.train[,3] - RESULT.train[,6])^2),
    mean((RESULT.train[,3] - RESULT.train[,9])^2),
    mean((RESULT.train[,4] - RESULT.train[,7])^2),
    mean((RESULT.train[,4] - RESULT.train[,10])^2),
    mean((RESULT.train[,5] - RESULT.train[,8])^2),
    mean((RESULT.train[,5] - RESULT.train[,11])^2),
    mean((RESULT.test[,3] - RESULT.test[,6])^2),
    mean((RESULT.test[,3] - RESULT.test[,9])^2),
    mean((RESULT.test[,4] - RESULT.test[,7])^2),
    mean((RESULT.test[,4] - RESULT.test[,10])^2),
    mean((RESULT.test[,5] - RESULT.test[,8])^2),
    mean((RESULT.test[,5] - RESULT.test[,11])^2),
    2837, 2984, 21187813, 28210104, NA, NA
    ), 6
)
colnames(MSE) <- c("train", "test", "Papier")
rownames(MSE) <- c("Mean", "Mean Tweedie", "VaR", "VaR Tweedie", "TVaR", "TVaR Tweedie")
round(MSE, 0)

save.image("C:/Users/Alex/Desktop/Recherche_ete_2019/Code_Alex/Massachussetts_2006/minimisation du mse.RData")
