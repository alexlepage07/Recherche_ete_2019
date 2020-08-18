# Document R servant à essayer les exemples donnés dans l'article "On Copula-based Collective Risk Models"
# écrit par Rosy Oh, Jae Youn Ahn et Woojoo Lee.

library(pscl)
library(statmod)
library(tweedie)
library(copula)
library(utils)


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

dGaussian.copula <- function(n, Y, x, N.coef, Y.coef, Y.disp, rho, z=1) {
    # Fonction qui calcule la densité d'une copule gaussienne
    if(typeof(x) == "list") x <- as.double(unlist(x))
    if(typeof(Y) == "list") Y <- as.double(unlist(Y))

    F1 <- function(n, x) {
      ppois(n, exp(t(x) %*% N.coef))
    }
    # if (n == 0) return(F1(0, x)) 
    F2 <- function(y, x) {
      pgamma(y,
             1 / Y.disp,
             exp(-t(x) %*% Y.coef) / Y.disp)
    }
    f2 <- function(y, x) {
      dgamma(y,
             1 / Y.disp,
             exp(-t(x) %*% Y.coef) / Y.disp)
    }
    
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
    
    d <- c_Y * (pnorm((qnorm(F1(n, x)) - mu) / Sig) - 
                        pnorm((qnorm(F1(n-1, x)) - mu) / Sig))
    # print(paste(n, Y, d))
    # if (d <= .Machine$double.eps) return(1)
    return(d)
}

optim_rho <- function(data, x, N.coef, Y.coef, Y.disp,
                      bornes, nsim=30, length.ech=1000){
    # Fonction qui permet de trouver les paramètres de dépendance
    # d'une copule gaussienne par optimisation numérique en utilisant un bootstrap.
    #
    # nsim correspond au nombre de simulations du bootstrap et 
    # length.ech correspond à la longueur des échantillons du bootstrap.
    
    # Bornes d'optimisation ui et ci.
    bounds <- matrix(bornes, ncol = 2, byrow = T)
    n_par <- nrow(bounds)
    ui <- rbind( diag(n_par), -diag(n_par) )
    ci <- c( bounds[,1], - bounds[,2] )
    # Retirer les valeurs infinies
    i <- as.vector(is.finite(bounds))
    ui <- ui[i,]
    ci <- ci[i]

    Depart <- numeric(2)
    Depart[1] <- cor(as.double(unlist(data[,1])), as.double(unlist(data[,2])), method="pearson")
    Depart[2] <- cor(as.double(unlist(data[,2])), as.double(unlist(data[,3])), method="pearson")
    
    mle_rho <- matrix(ncol=2, nrow=nsim)
    for (j in 1:nsim) {
        # Bootstrap et optimisation numérique.
        sample.id <- sample(1:nrow(data), length.ech)
        
        mle_rho[j, ] <- constrOptim(Depart, function(rho) {
            -sum(log(sapply(sample.id, function(i)
              dGaussian.copula(
                n = data[i, 1],
                Y = data[i, -1],
                x = x[i,],
                N.coef = N.coef,
                Y.coef = Y.coef,
                Y.disp = Y.disp,
                rho = rho
              ))))
        },
        grad = NULL, ui = ui, ci = ci, outer.eps = 1e-2)$par
        
        Depart <- mle_rho[j,]
        
        print(paste0(j/nsim*100,"%", " : ", list(mle_rho[j, ])))
    }
    mean_rho <- apply(mle_rho, 2, mean, na.rm = T)
    var_rho <- apply(mle_rho, 2, var)
    IC_rho <- matrix(c(mean_rho - qnorm(0.975) * sqrt(var_rho),
                       mean_rho + qnorm(0.975)* sqrt(var_rho)),
                     ncol=2)
    return(list("mean"=mean_rho, "variance"=var_rho, "IC"=IC_rho))
}


#==== Section 7 : Numercial study ====

parameter.settings <- matrix(
    # Paramètres initiaux présentés dans le tableau 1 de l'article.
    c(rep(-2.5,12), rep(0.5,12), rep(1, 6), rep(1.5, 6), 
      rep(8, 12), rep(-0.1, 12), rep(0.3, 12), rep(0.7, 12),
      rep(c(rep(-0.05, 2), rep(0.05, 2), rep(0.1, 2)), 2),
        rep(c(0.1,0.05),6)),
    nrow = 12
)

section7.simul <- function(nsim, scenario, z=1) {
    # Fonction qui permet de simuler des scénarios en fonction des paramètres.
    x_i <- rbinom(nsim, 1, 0.5)
    w_i <- rbinom(nsim, 1, 0.5)
    rho <- parameter.settings[scenario, 8:9]
    N.para <-  parameter.settings[scenario, 1:3]
    Y.para <- parameter.settings[scenario, 4:6]
    v <- parameter.settings[scenario, 7]
    U <- rCopula(nsim, normalCopula(P2p(Sig_rho1.2.z1(rho, 2)),3 , dispstr = "un"))
    
    Couples <- cbind(
      "N" = sapply(1:nsim, function(i)
        qpois(U[i, 1], exp(t(c(1, x_i[i], w_i[i])) %*% N.para))),
      "Y" = matrix(
        sapply(1:nsim, function(i)
          qgamma(U[i,-1], 1 / v, exp(-t(c(1, x_i[i], w_i[i])) %*% Y.para) / v)),
        nrow = nsim, byrow = T
      ))

    cbind(Couples, x_i, w_i)
}


tbl_result <- array(dim = c(12, 9, nsim <- 100))
# tbl_result <- as.data.frame(tbl_result)
colnames(tbl_result) <- c("beta0", "beta1", "beta2",
                          "gamma0", "gamma1", "gamma2", "v",
                          "rho1", "rho2")
for (simul in 1:nsim){
  for (scenario in 1:12){
    # Boucle qui recré chacun des scénario de la section 7 de l'article.
    
    print(paste0("simul ", simul, "/", nsim,", ", "Scenario ", scenario, "/12"))
    
    DATA <- as.data.frame(section7.simul(5e+3, scenario))
    
    N.mod <- glm(N ~ x_i + w_i, data = DATA, family = poisson("log"))
    N.coef <- N.mod$coefficients
    
    Y.mod <- glm(DATA[,2] ~ x_i + w_i, data = DATA, family = Gamma("log"))
    Y.disp <- MASS::gamma.dispersion(Y.mod)
    Y.coef <- Y.mod$coefficients
    k <- 2
    
    # (Depart <- c(cor(rank(DATA[,1], ties.method = "random"), DATA[,2], method = "spearman"),
    #              cor(DATA[,2], DATA[,3], method = "spearman")))
    # 
    # (rho2 <- 2 * sin(pi * Depart[2] / 6))
    # rho1 <- Depart[1]
    # DATA <- as.matrix(DATA)
    
    # system.time({
    #   rho1 <- optim(Depart[1], function(rho1)
    #     - sum(log(sapply(1:nrow(DATA), function(i)
    #       dGaussian.copula(
    #         n = DATA[i,1],
    #         Y = DATA[i, 2:3],
    #         x = (c(1, DATA[i,4:5])),
    #         N.coef , Y.coef, Y.disp,
    #         rho = c(rho1, rho2)
    #       )))),
    #     lower = -0.2, upper = sqrt(rho2), method="L-BFGS-B")$par
    # })[[3]]
    # # Optim est plus long que Optimize
    
    # # system.time({
    #   rho1 <- optimize(function(rho1)
    #     - sum(log(sapply(1:nrow(DATA), function(i)
    #       dGaussian.copula(
    #         n = DATA[i,1],
    #         Y = DATA[i, 2:3],
    #         x = (c(1, DATA[i,4:5])),
    #         N.coef , Y.coef, Y.disp,
    #         rho = c(rho1, rho2)
    #       )))),
    #     interval =  c(-0.2, sqrt(((k - 1) * rho2 + 1) / k)))$minimum
    # # })[[3]]
    
    # system.time({
    mle_rho <- optim_rho(data=DATA[,1:3],
                         x=cbind(1, DATA[,4:5]),
                         N.coef=N.coef,
                         Y.coef=Y.coef,
                         Y.disp=Y.disp,
                         bornes = c(-0.2, 0.2, 0, 0.2),
                         nsim = 1,
                         length.ech = 5000)
    # })[[3]]
    rho1 <- mle_rho$mean[1] ; rho2 <- mle_rho$mean[2]
    # La full maximum likelyhood pour les paramètres de dépendance
    # Est beaucoup plus longue, mais n'apporte pas beaucoup plus de
    # précision par rapport à la méthode par parties.
      
    tbl_result[scenario,,simul] <-
      c(
        (round(summary(N.mod)$coefficients[, 1], 1)),
        (round(summary(Y.mod)$coefficients[, 1], 1)),
        round(summary(Y.mod)$dispersion, 1),
        round(rho1, 2),
        round(rho2, 2)
      )
    # print(tbl_result[scenario,,simul])
  }
  print(tbl_result[,,simul])
}

tbl_result.means <- array(dim = c(12, 9))
for (i in 1:12) {
  for (j in 1:9) {
    tbl_result.means[i,j] <- mean(tbl_result[i,j,], na.rm = T)
  }
}
colnames(tbl_result.means) <- c("beta0", "beta1", "beta2",
                         "gamma0", "gamma1", "gamma2", "v",
                         "rho1", "rho2")
round(tbl_result.means, 2)

# Tableau 2 : Biais
tbl_biais <- (tbl_result.means - parameter.settings) / (parameter.settings) * 100
colnames(tbl_biais) <- c("beta0", "beta1", "beta2",
                         "gamma0", "gamma1", "gamma2", "v",
                         "rho1", "rho2")
round(tbl_biais, 2)

# Tableau 3 : MSE
tbl_MSE <- diag(t(tbl_result.means - parameter.settings) %*% (tbl_result.means - parameter.settings)) / nsim
tbl_MSE <- t(tbl_MSE)
colnames(tbl_MSE) <- c("beta0", "beta1", "beta2",
                         "gamma0", "gamma1", "gamma2", "v",
                         "rho1", "rho2")
round(tbl_MSE, 4)


#==== Section 8 : Real data analysis ====

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
SEV <- AGG[train.id,]
SEV <- SEV[SEV$losspaid>0,]
LOSS <- aggregate(SEV$losspaid,
                  by = list(SEV$pol_id,
                            SEV$class4,
                            SEV$tgroup),
                  list)
colnames(LOSS) <- c("pol_id", "class4", "tgroup", "losspaid")
FREQ <- aggregate(AGG[train.id, ]$losspaid,
               by = list(AGG[train.id, ]$pol_id,
                         AGG[train.id, ]$class4,
                         AGG[train.id, ]$tgroup),
               function(x) sum(x>0))

colnames(FREQ) <- c("pol_id", "class4", "tgroup", "N")
# nrow(FREQ); nrow(SEV)
summary(unlist(FREQ$N))
# summary(unlist(SEV$losspaid))
barplot(table(FREQ$N))

# tgroup : 1 = least risky, 6 = most risky
# class4 : A = Adult, B = Business, I = <3 yrs experience, M = 3-6 yrs experience, S = Senior

#==== Analyse de la fréquence ===================
n_obs <- nrow(FREQ)

# Matrice d'aggrégation de la moyenne de la fréquence.
agg_means_AGG <- aggregate(FREQ$N, by = list(FREQ$class4, FREQ$tgroup), mean)
mat_means_AGG <- matrix(agg_means_AGG$x, ncol = 5, byrow = T)
colnames(mat_means_AGG) <- agg_means_AGG$Group.1[1:5]
rownames(mat_means_AGG) <- 1:6
round(mat_means_AGG,3)
# J'arrive à reproduire le pattern attendu par rapport à ce qui est donné
# dans le tableau 4 du papier de Oh, Ahn et Lee 2019 même si les valeurs
# ne sont pas les mêmes.

#Matrice d'aggrégation de la distribution des features
agg_means_AGG <- aggregate(FREQ$pol_id, 
                           by = list(FREQ$class4, FREQ$tgroup), 
                           function(i) length(i) / n_obs)
mat_means_AGG <- matrix(agg_means_AGG$x, ncol = 5, byrow = T)
colnames(mat_means_AGG) <- agg_means_AGG$Group.1[1:5]
rownames(mat_means_AGG) <- 1:6
round(mat_means_AGG * 100,2) # Valeurs données en pourcentages.
# J'ai une répartion très similaire à ce qui est présenté dans le tableau
# 4 du papier de Oh, Ahn et Lee.

# sort(FREQ$N)[ceiling(0.9999 * n_obs)]

hurdle.model <- hurdle(N ~ class4 + tgroup, data = FREQ)
(Sumry_1 <- summary(hurdle.model))

{
# hurdle.model_2 <-
#     hurdle(CompteDeclm_id ~ class4 + tgroup | 1,
#            data = AGG[train.id, ] ,
#            offset = log(earnexpo))
# (Sumry_2 <- summary(hurdle.model_2))

# Poisson.model <- glm(unlist(N) ~ class4 + tgroup, data = FREQ[FREQ$N < 5,], family = poisson("log"))
# (Sumry_3.1 <- summary(Poisson.model))
# AIC(hurdle.model, Poisson.model, k = log(nrow(AGG)))
# Sumry_3.1$aic - 20 + 2*Sumry_1$loglik
# qchisq(0.99, 10)

# zero.infl.model <- zeroinfl(CompteDeclm_id ~ class4 + tgroup,
#                             data = AGG[train.id, ],
#                             offset = log(earnexpo))
# (Sumry_3 <- summary(zero.infl.model))
# 
# hurdle.model_3 <-
#     hurdle(CompteDeclm_id ~ class4 + tgroup,
#            data = AGG[train.id, ],
#            offset = log(earnexpo), dist = "negbin")
# (Sumry_4 <- summary(hurdle.model))
# 
# # Test d'AIC
# AIC(hurdle.model, hurdle.model_2, hurdle.model_3, zero.infl.model)
# (Sumry_1$loglik - Sumry_2$loglik) *2
# qchisq(0.95, 9)
# # Le modèle de Hurdle 3 (celui qui utilise tous les paramètres pour sa 
# # régression binomiale et qui est jumelé à une loi binomiale négative)
# # est le plus adéquat.
# (Sumry_1$loglik - Sumry_4$loglik) *2
# # Cependant, la version avec la loi Poisson (Hurdle.model_1) est une
# # bonne simplification
} # Tests avec d'autres modèles possibles.
N.coef <- hurdle.model$coefficients$count
R.coef <- hurdle.model$coefficients$zero


AGG.pred <- predict(hurdle.model, newdata = AGG[-train.id, 3:5],
                     type = "response")

agg_pred_AGG <- aggregate(AGG.pred, by = list(AGG$class4[-train.id], AGG$tgroup[-train.id]), mean)
agg_N.test <- aggregate(AGG[-train.id,]$losspaid,
                        by = list(AGG[-train.id, ]$pol_id,
                                  AGG[-train.id, ]$class4,
                                  AGG[-train.id, ]$tgroup),
                        function(x) sum(x > 0))
agg_means_AGG <- aggregate(agg_N.test$x,
                           by = list(agg_N.test$Group.2,
                                     agg_N.test$Group.3),
                           mean)
tbl_agg_AGG <- cbind(agg_means_AGG, agg_pred_AGG$x, agg_pred_AGG$x - agg_means_AGG$x)
colnames(tbl_agg_AGG) <- c("Class", "Territory", "Empirical mean", "theorical mean", "Écarts")
tbl_agg_AGG
sum((tbl_agg_AGG$`Empirical mean` - tbl_agg_AGG$`theorical mean`)^2 / tbl_agg_AGG$`Empirical mean`)
qchisq(0.95, 30 - 20 -1)

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
  ceiling(uniroot(function(n) F1(n, x) - p, c(0, 50))$root)
}

#==== Analyse de la sévérité =================
hist(SEV$losspaid)
sum(SEV$losspaid == 25000) / length(SEV$losspaid)

# Matrice d'aggrégation de la moyenne de la sévérité
agg_means_sev <- aggregate(unlist(SEV$losspaid), by = list(SEV$class4, SEV$tgroup), mean)
mat_means_sev <- matrix(agg_means_sev$x, ncol = 5, byrow = T)
colnames(mat_means_sev) <- agg_means_sev$Group.1[1:5]
rownames(mat_means_sev) <- 1:6
mat_means_sev
# J'arrive à reproduire le pattern attendu par rapport à ce qui est donné
# dans le tableau 4 du papier de Oh, Ahn et Lee 2019 même si les valeurs
# ne sont pas les mêmes.

# Modèle de sévérité gamma construit avec un glm
randomize.Y <- function(data, n) {
  data <- data[data[,1] >= n,]
  Y <- matrix(sapply(1:nrow(data), function(i) sample(unlist(data[,-1][[i]]), n)), ncol=n, byrow = T)
  cbind(unlist(data[,1]), Y)
}
# summary(randomize.Y(cbind(FREQ$N[FREQ$N > 0], LOSS$losspaid)))

Y.disp_ <- numeric(nsim <- 30)
Y.coef_ <- matrix(NA, nsim, 10)
for (sim in 1:nsim){
  print(paste0(sim, "/", nsim))
  rand.Y <- randomize.Y(cbind(FREQ$N[FREQ$N > 0], LOSS$losspaid), 1)
  gamma.data <- cbind(LOSS[, 1:3], rand.Y[, 2])
  colnames(gamma.data) <- colnames(LOSS)
  
  gamma.model <- glm(losspaid ~ class4 + tgroup, data = gamma.data,
                     family = Gamma(link = "log"))
  # Dans l'article, Oh, Ahn et Lee ont pris le lien log.
  Y.disp_[sim] <- MASS::gamma.dispersion(gamma.model)
  Y.coef_[sim,] <- gamma.model$coefficients
  
  print(list("dispersion" = Y.disp_[sim], "coefficients" = Y.coef_[sim, ]))
  }

Y.disp <- mean(Y.disp_, na.rm=T)
Y.coef <- apply(Y.coef_, 2, mean, na.rm=T)

F2 <- function(Y, x) {
  pgamma(Y, 1 / Y.disp, exp(-x %*% Y.coef) / Y.disp)
}
f2 <- function(Y, x) {
  dgamma(Y, 1 / Y.disp, exp(-x %*% Y.coef) / Y.disp)
}
F2.inv <- function(p, x) {
  qgamma(p, 1 / Y.disp, exp(-x %*% Y.coef) / Y.disp)
}


#--- Mesures de l'adéquation du modèle gamma 
SEV.test <- AGG[-train.id,]
SEV.test <- SEV.test[SEV.test$losspaid > 0,]
x.mat_test <- model.matrix(losspaid ~ class4 + tgroup, data = SEV.test)
x.mat_train <- model.matrix(losspaid ~ class4 + tgroup, data = SEV)


# # Analyse graphique de l'adéquation 
#
# U <- F2(SEV.test$losspaid, x.mat_test)
# q <- qgamma(U, 1 / Y.disp, 1 / Y.disp / exp(x.mat_test %*% Y.coef))
# plot(SEV.test$losspaid, q, type = "l",
#      xlab = "quantiles empiriques",
#      ylab = "quantiles théoriques")
# abline(a=0, b=1, col="red")
# axis(2, tck = 1, lty = 2, col = "grey")
# axis(1, tck=1, lty = 2, col = "grey")
# # Graphiquement, le modèle semble adéquat.


#--- Tests d'adéquation
# Sur les données d'entraînement
Y.obs <- aggregate(rand.Y[,2], list(SEV[FREQ$N > 0,]$class4, SEV[FREQ$N > 0,]$tgroup), mean)
x.mat_train <- model.matrix(x ~ Group.1 + Group.2, data = Y.obs)
sev.pred <- exp(x.mat_train %*% Y.coef)

cbind("Predictions" = sev.pred, "Observations" = Y.obs$x)

mat_sev.pred <- matrix(sev.pred, ncol = 5, byrow = T)
colnames(mat_sev.pred) <- agg_means_sev$Group.1[1:5]
rownames(mat_sev.pred) <- 1:6
mat_sev.pred

# Chi2 de pearson
(Chi2_pearson <- sum((sev.pred - Y.obs$x)^2) / var(sev.pred))
qchisq(0.99, nrow(SEV) - 11)
# On trouve que le modèle est adéquat si on prend un sinistre par assuré de façon aléatoire.

# Kolmogorov-Smirnov
ks.test(sev.pred, Y.obs$x)


#==== Modèle d'aggrégation =================

x.mat_train <- model.matrix(losspaid ~ class4 + tgroup, data=SEV)
x.mat_train <- aggregate(x.mat_train,
                         list(SEV$pol_id, SEV$class4, SEV$tgroup),
                         data.table::first)[,-(2:3)]
k <- max(FREQ$N)
N <- FREQ[FREQ$N > 0, ]$N


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

spearmans.Rho <- function(data, nsim=100) {
  # Fonction qui retourne le rho de spearman pour une base de données d'assurance.
  pb <- txtProgressBar(min = 1, max = nsim, style = 3)
  cor.spearman <- numeric(nsim)
  for (i in 1:nsim) {
    couples.YY <- randomize.Y(data, 2)
    cor.spearman[i] <- cor(couples.YY[,2], couples.YY[,3])
    setTxtProgressBar(pb, i)
  }
  rep <- list()
  rep$mean <- mean(cor.spearman)
  rep$variance <- var(cor.spearman)
  rep$IC <- c(rep$mean - qnorm(0.975) * sqrt(rep$variance),
              rep$mean + qnorm(0.975) * sqrt(rep$variance))
  close(pb)
  return(rep)
}

# fct_Score <- function(rho) {
#   # Fonction à minimiser afin de maximiser la vraisemblance du
#   # paramètre de dépendance rho1.
#   
#   print(rho)
#   
#   Score <- - sum(log(sapply(1:nrow(LOSS), function(i)
#     dGaussian.copula(
#       n = N[[i]],
#       Y = unlist(LOSS$losspaid[[i]]),
#       x = x.mat_train[i, -1],
#       rho = c(rho, rho2.mean)
#     ))))
#   return(Score)
# }

optim_rho <- function(data, x, Depart, bornes, nsim=50, length.ech=1e+3){
  # Fonction qui permet de trouver les paramètres de dépendance
  # d'une copule gaussienne par optimisation numérique en utilisant un bootstrap.
  #
  # nsim correspond au nombre de simulations du bootstrap et 
  # length.ech correspond à la longueur des échantillons du bootstrap.
  
  # Bornes d'optimisation ui et ci.
  bounds <- matrix(bornes, ncol = 2, byrow = T)
  ui <- rbind( diag(2), -diag(2) )
  ci <- c( bounds[,1], - bounds[,2] )
  
  pb <- txtProgressBar(min = 1, max = nsim, style = 3)
  mle_rho <- matrix(ncol=2, nrow=nsim)
  for (j in 1:nsim) {
    # Bootstrap et optimisation numérique.
    sample.id <- sample(nrow(LOSS), length.ech)
    
    mle_rho[j,] <- constrOptim(Depart, function(rho) { 
      - sum(log(sapply(sample.id, function(i)
      dGaussian.copula(
        n = N[[i]],
        Y = unlist(LOSS$losspaid[[i]]),
        x = x.mat_train[i, -1],
        rho = rho
      ))))},
        grad = NULL, ui = ui, ci = ci, outer.eps = 1e-2)$par
    
    Depart <- mle_rho[j,]
    
    setTxtProgressBar(pb, j)
  }
  
  mean_rho <- apply(mle_rho, 2, mean, na.rm = T)
  var_rho <- apply(mle_rho, 2, var)
  IC_rho <- matrix(c(mean_rho - qnorm(0.975) * sqrt(var_rho),
                   mean_rho + qnorm(0.975)* sqrt(var_rho)),
                   ncol=2)
  names(mean_rho) <-  names(var_rho) <- c("rho1","rho2")
  dimnames(IC_rho) <- list(c("rho1","rho2"), c("inf","sup"))
  close(pb)
  return(list("mean"=mean_rho, "variance"=var_rho, "IC"=IC_rho))
}


rho1 <- cor(N, rand.Y[,2], method = "spearman")
rho1 <- 2 * sin(pi * rho1 / 6)

(rho2 <- spearmans.Rho(cbind(FREQ$N[FREQ$N > 0], LOSS$losspaid)))
rho2.mean <- 2 * sin(pi * rho2$mean / 6)
rho2.IC <- 2 * sin(pi * rho2$IC / 6)

# unique(copula::rho(normalCopula(P2p(Sig_rho1.2.z1(c(0.2, rho1$mean), k - 1)), k, dispstr = "un")))[2]
# 6/pi * asin(Depart$mean / 2)
# # On voit que le rho de Spearman empirique est un excellent point de départ pour l'optimisation
# # numérique puisque pour la copule gaussienne, celui-ci est pratiquement identique à ses
# # paramètres de dépendance.

# rho1 <- optimize(fct_Score, lower = -1,
#                     upper = sqrt(((k - 1) * rho2.mean + 1) / k))$minimum
# rho <- c(rho1, rho2.mean)


mle_rho <- optim_rho(data = cbind(N, LOSS$losspaid),
                     x = x.mat_train[,-1],
                     Depart = c(rho1, rho2.mean),
                     bornes = c(-1, sqrt(((k - 1) * rho2.mean + 1) / k),
                                rho2.IC[1], rho2.IC[2]),
                     nsim=50,
                     length.ech=1e+3)
rho <- mle_rho$mean


#==== Reproduction du tableau 6 - Résultats ----
rGaussian.copula <- function(n_sim, rho, k = max(FREQ$N), z = 1) {
    # Fonction qui permet de simuler des réalisations d'une copule gaussienne.
    Sig_rho <- ifelse(z == 1, Sig_rho1.2.z1, Sig_rho1.2.z2)
    rCopula(n_sim, normalCopula(P2p(Sig_rho(rho, k)), k+1, dispstr = "un"))
}

simul.S <- function(n_sim, x, rho, k = max(FREQ$N), z = 1) {
    # Fonction qui permet de faire de la simulation Monte-Carlo
    # pour estimer les mesures de risque théorique du modèle.
    U <- rGaussian.copula(n_sim, rho, k, z)
    
    Couples <- cbind(
      "N" = sapply(1:n_sim, function(i)
        F1.inv(U[i,1], unlist(x))),
      "Y" = matrix(
        sapply(1:n_sim, function(i)
          qgamma(U[i,-1], 1 / Y.disp, 1 / exp(unlist(x) %*% Y.coef) / Y.disp)),
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
    mean_S <- mean(S)
    VaR_0.995 <- sort(S)[ceiling(0.995 * length(S))]
    return(list(
        "E[S]" = mean_S,
        "VaR_0.995" = VaR_0.995
    ))
}

# Mesures empiriques
AGG_train <- AGG[train.id,]
S.empirique <- aggregate(AGG_train$losspaid,
                         by = list(AGG_train$pol_id, AGG_train$tgroup, AGG_train$class4),
                         sum)
mean.S.empirique <- aggregate(S.empirique$x,
                              by = list(S.empirique$Group.2, S.empirique$Group.3),
                              mean)
mean.S.empirique$x <- round(mean.S.empirique$x)
VaR <- function(x) sort(x)[ceiling(length(x) * 0.995)]
VaR.S.empirique <- aggregate(S.empirique$x,
                             by = list(S.empirique$Group.2, S.empirique$Group.3),
                             VaR)
tbl_6.empirique <- cbind(mean.S.empirique, VaR.S.empirique$x)
colnames(tbl_6.empirique) <- c("Territory", "Class", "mean(S)", "VaR_n")
tbl_6.empirique <- data.table::data.table(tbl_6.empirique)
tbl_6.empirique <- data.table::setorderv(tbl_6.empirique, c('Territory', 'Class'), c(1, 1))

# Mesures théoriques
proposed <- matrix(ncol=2, nrow=30)
risk_groups <- unique(x.mat_train[,-1])
for (i in 1:30) {
  proposed[i,] <- round(unlist(simul.S(5e+3, unlist(risk_groups[i,]), rho)))
}
Class <- as.matrix(risk_groups[1:5]) %*% 0:4 + 1
Territory <- as.matrix(risk_groups[c(1,6:10)]) %*% 0:5 + 1

tbl_6.theorique <- cbind(Territory, Class, proposed)
colnames(tbl_6.theorique) <- c("Territory", "Class", "E[S]", "VaR_0.995")
tbl_6.theorique <- data.table::data.table(tbl_6.theorique)
tbl_6.theorique <- data.table::setorderv(tbl_6.theorique, c('Territory', 'Class'), c(1, 1))
tbl_6 <- cbind(tbl_6.empirique, tbl_6.theorique[,3:4])


# Modèle de Tweedie

AGG_2 <- S.empirique
colnames(AGG_2) <- c("pol_id", "class4", "tgroup", "losspaid")
# head(AGG_2)
# summary(AGG_2)

# prof <- tweedie.profile(
#   losspaid ~ class4 + tgroup,
#   data = AGG_2,
#   p.vec = seq(1.05, 1.25, 0.05),
#   method = "saddlepoint",
#   control = list(maxit = 400)
# )
p.max <- 1
Tweedie.mod <- glm(losspaid ~ class4 + tgroup, data = AGG_2,
                  family = tweedie(p.max, 0))
summary(Tweedie.mod)
Tweedie.mod$coefficients

x.agg.train <- model.matrix(losspaid ~ class4 + tgroup,
                          data = AGG_2)

tbl_6.tweedie <- aggregate(fitted(Tweedie.mod), list(AGG_2$class4, AGG_2$tgroup), mean)
tbl_6.tweedie$x <- round(tbl_6.tweedie$x)
# VaR_tweedie.pred <- aggregate(fitted(Tweedie.mod), list(AGG_2$class4, AGG_2$tgroup), VaR) # ne marche pas.
colnames(tbl_6.tweedie) <- c("Territory", "Class", "E[S].Tweedie")
tbl_6.tweedie <- data.table::data.table(tbl_6.tweedie)
tbl_6.tweedie <- data.table::setorderv(tbl_6.tweedie, c('Territory', 'Class'), c(1, 1))
tbl_6 <- cbind(tbl_6, tbl_6.tweedie[,3])

# Le modèle surestime la somme des réclamations...
MSE.train.mean <- mean((tbl_6$`mean(S)` - tbl_6$`E[S]`)^2)
MSE.train.Tweedie <- mean((tbl_6$`mean(S)` - tbl_6$`E[S].Tweedie`)^2)
MSE.train.VaR <- mean((tbl_6$VaR_n - tbl_6$VaR_0.995)^2)


#---- Avec les données de test ----

AGG_test <- AGG[-train.id,]
S.empirique <- aggregate(AGG_test$losspaid,
                         by = list(AGG_test$pol_id, AGG_test$tgroup, AGG_test$class4),
                         sum)
mean.S.empirique <- aggregate(S.empirique$x,
                              by = list(S.empirique$Group.2, S.empirique$Group.3),
                              mean)
mean.S.empirique$x <- round(mean.S.empirique$x)
VaR <- function(x) sort(x)[ceiling(length(x) * 0.995)]
VaR.S.empirique <- aggregate(S.empirique$x,
                             by = list(S.empirique$Group.2, S.empirique$Group.3),
                             VaR)
tbl_6.empirique <- cbind(mean.S.empirique, VaR.S.empirique$x)
colnames(tbl_6.empirique) <- c("Territory", "Class", "mean(S)", "VaR_n")
tbl_6.empirique <- data.table::data.table(tbl_6.empirique)
tbl_6.empirique <- data.table::setorderv(tbl_6.empirique, c('Territory', 'Class'), c(1, 1))

# Mesures théoriques
proposed <- matrix(ncol=2, nrow=30)

x.mat_test <- model.matrix(losspaid ~ class4 + tgroup, data=SEV.test)
x.mat_test <- aggregate(x.mat_test,
                         list(SEV.test$pol_id, SEV.test$class4, SEV.test$tgroup),
                         data.table::first)[,-(2:3)]

risk_groups <- unique(x.mat_test[,-1])
for (i in 1:30) {
  proposed[i,] <- round(unlist(simul.S(1e+4, unlist(risk_groups[i,]), rho)))
}
Class <- as.matrix(risk_groups[1:5]) %*% 0:4 + 1
Territory <- as.matrix(risk_groups[c(1,6:10)]) %*% 0:5 + 1

tbl_6.theorique <- cbind(Territory, Class, proposed)
colnames(tbl_6.theorique) <- c("Territory", "Class", "E[S]", "VaR_0.995")
tbl_6.theorique <- data.table::data.table(tbl_6.theorique)
tbl_6.theorique <- data.table::setorderv(tbl_6.theorique, c('Territory', 'Class'), c(1, 1))
tbl_6_test <- cbind(tbl_6.empirique, tbl_6.theorique[,3:4])


# Modèle de Tweedie

AGG_2 <- S.empirique
colnames(AGG_2) <- c("pol_id", "class4", "tgroup", "losspaid")
tweedie.pred <- round(predict(Tweedie.mod, newdata = AGG_2, type = "response"))
mean_tweedie.pred <- aggregate(tweedie.pred, list(AGG_2$class4, AGG_2$tgroup), mean)
colnames(mean_tweedie.pred) <- c("Territory", "Class", "E[S].Tweedie")
tbl_6.tweedie <- data.table::data.table(tbl_6.tweedie)
tbl_6.tweedie <- data.table::setorderv(tbl_6.tweedie, c('Territory', 'Class'), c(1, 1))
tbl_6_test <- cbind(tbl_6_test, tbl_6.tweedie[,3])



MSE.test.mean <- mean((tbl_6_test$`mean(S)` - tbl_6_test$`E[S]`)^2)
MSE.test.VaR <- mean((tbl_6_test$VaR_n - tbl_6_test$VaR_0.995)^2)
MSE.test.tweedie <- mean((tbl_6_test$`E[S].Tweedie` - tbl_6_test$`mean(S)`)^2)

#---- Compilation des résultats ----
tbl_6
tbl_6_test

round(rbind("Mean" = c("Test" = MSE.test.mean, "train"=MSE.train.mean,
                    "Article"=2837),
      "Tweedie" = c(MSE.test.tweedie, MSE.train.Tweedie, 2984),
      "VaR" = c(MSE.test.VaR, MSE.train.VaR, 21187813)
      ))
