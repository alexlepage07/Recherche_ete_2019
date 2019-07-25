# Document R servant à essayer les exemples donnés dans l'article "On Copula-based Collective Risk Models"
# écrit par Rosy Oh, Jae Youn Ahn et Woojoo Lee.

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
    Depart[1] <- cor(as.double(unlist(data[,1])), as.double(unlist(data[,2])), method="spearman")
    Depart[2] <- cor(as.double(unlist(data[,2])), as.double(unlist(data[,3])), method="spearman")
    
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
        
        # print(paste0(j/nsim*100,"%", " : ", list(mle_rho[j, ])))
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


tbl_result <- array(dim = c(12, 9, nsim <- 10))
# tbl_result <- as.data.frame(tbl_result)
colnames(tbl_result) <- c("beta0", "beta1", "beta2",
                          "gamma0", "gamma1", "gamma2", "v",
                          "rho1", "rho2")
for (simul in 1:nsim){
  for (scenario in 1:12){
    # Boucle qui recré chacun des scénario de la section 7 de l'article.
    
    print(paste0("Scenario ", scenario, "/12"))
    
    DATA <- as.data.frame(section7.simul(5e+3, scenario))
    
    N.mod <- glm(N ~ x_i + w_i, data = DATA, family = poisson("log"))
    N.coef <- N.mod$coefficients
    
    Y.mod <- glm(DATA[,2] ~ x_i + w_i, data = DATA, family = Gamma("log"))
    Y.disp <- MASS::gamma.dispersion(Y.mod)
    Y.coef <- Y.mod$coefficients
    k <- 2
    
    (Depart <- c(cor(rank(DATA[,1], ties.method = "random"), DATA[,2], method = "spearman"),
                 cor(DATA[,2], DATA[,3], method = "spearman")))
    
    # (rho2 <- 2 * sin(pi * Depart[2] / 6))
    # DATA <- as.matrix(DATA)
    
    # system.time({
    #   mle_rho <- optim(Depart[1], function(rho1)
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
    
    # system.time({
      # mle_rho1 <- optimize(function(rho1)
      #   - sum(log(sapply(1:nrow(DATA), function(i)
      #     dGaussian.copula(
      #       n = DATA[i,1],
      #       Y = DATA[i, 2:3],
      #       x = (c(1, DATA[i,4:5])),
      #       N.coef , Y.coef, Y.disp,
      #       rho = c(rho1, rho2)
      #     )))),
      #   interval =  c(-0.1, sqrt(((k - 1) * rho2 + 1) / k)))$minimum
    # })[[3]]
    
    # system.time({
    #   mle_rho <- optim_rho(data=DATA[,1:3],
    #                        x=cbind(1, DATA[,4:5]),
    #                        N.coef=N.coef,
    #                        Y.coef=Y.coef,
    #                        Y.disp=Y.disp,
    #                        bornes = c(-0.2, 0.2, 0, 0.2),
    #                        nsim = 1,
    #                        length.ech = 5000)
    # })[[3]]
    # mle_rho1 <- mle_rho$mean[1] ; rho2 <- mle_rho$mean[2]
    # # La full maximum likelyhood pour les paramètres de dépendance
    # # Est beaucoup plus longue, mais n'apporte pas beaucoup plus de
    # # précision par rapport à la méthode par parties.
      
    tbl_result[scenario,,simul] <-
      c(
        (round(summary(N.mod)$coefficients[, 1], 1)),
        (round(summary(Y.mod)$coefficients[, 1], 1)),
        round(summary(Y.mod)$dispersion, 1),
        round(mle_rho1, 2),
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
FREQ <- read.csv("FREQ.txt")
SEV <- read.csv("SEVERITY.txt")
# tgroup : 1 = least risky, 6 = most risky
# class4 : A = Adult, B = Business, I = <3 yrs experience, M = 3-6 yrs experience, S = Senior

#==== Analyse de la fréquence ===================
# Réarrangement des données en facteurs
FREQ$class4 <- FREQ$class4 %% 10
FREQ$class4[FREQ$class4 == 1] <- "A"
FREQ$class4[FREQ$class4 == 2] <- "S"
FREQ$class4[FREQ$class4 %in% c(3, 4)] <- "M"
FREQ$class4[FREQ$class4 == 5] <- "B"
FREQ$class4[FREQ$class4 %in% (6:9)] <- "I"
FREQ$class4 <- factor(FREQ$class4)

FREQ$tgroup <- factor(FREQ$tgroup)
FREQ$pol_id <- factor(FREQ$pol_id)
summary(FREQ)
head(FREQ)
attach(FREQ)

n_obs <- nrow(FREQ)

# Matrice d'aggrégation de la moyenne de la fréquence.
agg_means_freq <- aggregate(CompteDeclm_id, by = list(class4, tgroup), mean)
mat_means_freq <- matrix(agg_means_freq$x, ncol = 5, byrow = T)
colnames(mat_means_freq) <- agg_means_freq$Group.1[1:5]
rownames(mat_means_freq) <- 1:6
round(mat_means_freq,3)
# J'arrive à reproduire le pattern attendu par rapport à ce qui est donné
# dans le tableau 4 du papier de Oh, Ahn et Lee 2019 même si les valeurs
# ne sont pas les mêmes.

#Matrice d'aggrégation de la distribution des features
agg_means_freq <- aggregate(pol_id, by = list(class4, tgroup), 
                            function(i) length(i) / n_obs)
mat_means_freq <- matrix(agg_means_freq$x, ncol = 5, byrow = T)
colnames(mat_means_freq) <- agg_means_freq$Group.1[1:5]
rownames(mat_means_freq) <- 1:6
round(mat_means_freq * 100,2) # Valeurs données en pourcentages.
# J'ai une répartion très similaire à ce qui est présenté dans le tableau
# 4 du papier de Oh, Ahn et Lee.


train.id <- 1:1e+6
test.id <- (max(train.id)+1):n_obs

hurdle.model <-
    hurdle(CompteDeclm_id ~ class4 + tgroup,
           data = FREQ[train.id, ],
           offset = log(earnexpo))
(Sumry_1 <- summary(hurdle.model))

{
# hurdle.model_2 <-
#     hurdle(CompteDeclm_id ~ class4 + tgroup | 1,
#            data = FREQ[train.id, ] ,
#            offset = log(earnexpo))
# (Sumry_2 <- summary(hurdle.model_2))

# Poisson.model <-
#     glm(CompteDeclm_id ~ class4 + tgroup,
#            data = FREQ[train.id, ] ,
#            offset = log(earnexpo),
#         family = poisson("log"))
# (Sumry_3.1 <- summary(Poisson.model))
# AIC(hurdle.model_2, Poisson.model)
# Sumry_3.1$aic - 20 + 2*Sumry_2$loglik
# qchisq(0.95,1)

# zero.infl.model <- zeroinfl(CompteDeclm_id ~ class4 + tgroup,
#                             data = FREQ[train.id, ],
#                             offset = log(earnexpo))
# (Sumry_3 <- summary(zero.infl.model))
# 
# hurdle.model_3 <-
#     hurdle(CompteDeclm_id ~ class4 + tgroup,
#            data = FREQ[train.id, ],
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


freq.pred <- predict(hurdle.model, newdata = FREQ[test.id, 3:5],
                     type = "response")

agg_pred_freq <- aggregate(freq.pred, by = list(class4[test.id], tgroup[test.id]), mean)
tbl_agg_freq <- cbind(agg_means_freq, agg_pred_freq$x, agg_pred_freq$x - agg_means_freq$x)
colnames(tbl_agg_freq) <- c("Class", "Territory", "Empirical mean", "theorical mean", "Écarts")
tbl_agg_freq
mean((tbl_agg_freq$`Empirical mean` - tbl_agg_freq$`theorical mean`)^2)


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
# Réarrangement des données en facteurs
SEV$class4 <- SEV$class4 %% 10
SEV$class4[SEV$class4 == 1] <- "A"
SEV$class4[SEV$class4 == 2] <- "S"
SEV$class4[SEV$class4 %in% c(3, 4)] <- "M"
SEV$class4[SEV$class4 == 5] <- "B"
SEV$class4[SEV$class4 %in% (6:9)] <- "I"
SEV$class4 <- factor(SEV$class4)

SEV$tgroup <- factor(SEV$tgroup, 1:6)
SEV$pol_id <- factor(SEV$pol_id)
# SEV <- SEV[SEV$losspaid>0,]
summary(SEV)
hist(SEV$losspaid)
# La distribution marginale de la sévérité peut être une gamma, 
# une exponentielle ou une lognormale.
head(SEV)
attach(SEV)

# Matrice d'aggrégation de la moyenne de la sévérité
agg_means_sev <- aggregate(SEV$losspaid, by = list(SEV$class4, SEV$tgroup), mean)
mat_means_sev <- matrix(agg_means_sev$x, ncol = 5, byrow = T)
colnames(mat_means_sev) <- agg_means_sev$Group.1[1:5]
rownames(mat_means_sev) <- 1:6
mat_means_sev
# J'arrive à reproduire le pattern attendu par rapport à ce qui est donné
# dans le tableau 4 du papier de Oh, Ahn et Lee 2019 même si les valeurs
# ne sont pas les mêmes.

# Modèle de sévérité gamma construit avec un glm
train.id <- which(SEV$pol_id %in% FREQ[train.id,]$pol_id)
x.mat_train <- model.matrix(losspaid ~ class4 + tgroup, data=SEV[train.id,])
{
# bounds <- matrix(c(1e-6, Inf,
#                    rep(c(-Inf, Inf), 10))
#                    , ncol = 2, byrow = T)
# # Convertir les constraintes en matrices ui et ci.
# n_par <- nrow(bounds)
# ui <- rbind( diag(n_par), -diag(n_par) )
# ci <- c( bounds[,1], - bounds[,2] )
# # Retirer les valeurs infinies
# i <- as.vector(is.finite(bounds))
# ui <- ui[i,]
# ci <- ci[i]
# 
# gamma.model <- constrOptim(c(5, rep(1e-6, 10)), function(para)
#     - sum(log(dgamma(
#         SEV[train.id, 2], para[1], x.mat %*% para[-1]
#     ))),
#     grad = NULL, ui = ui, ci = ci)
# round(gamma.model$par,3)
}
gamma.model <- glm(losspaid ~ class4 + tgroup, data = SEV[train.id, ],
                   family = Gamma(link = "log"))
# Dans l'article, Oh, Ahn et Lee ont pris le lien log.
(sumry_gam <- summary(gamma.model))

Y.disp <- sumry_gam$dispersion
Y.coef <- gamma.model$coefficients

F2 <- function(Y, x) {
  pgamma(Y, 1 / Y.disp, 1 / Y.disp / exp(x %*% Y.coef))
}
f2 <- function(Y, x) {
  dgamma(Y, 1 / Y.disp, 1 / Y.disp / exp(x %*% Y.coef))
}
F2.inv <- function(p, x) {
  qgamma(p, 1 / Y.disp, 1 / Y.disp / exp(x %*% Y.coef))
}

# Analyse de l'adéquation graphiquement
x.mat_test <- model.matrix(losspaid ~ class4 + tgroup, data=SEV[-train.id,])
U <- F2(SEV[-train.id, ]$losspaid, x.mat_test)
q <- qgamma(U, 1 / Y.disp, 1 / Y.disp / exp(x.mat_test %*% Y.coef))
plot(SEV[-train.id, ]$losspaid, q, type = "l",
     xlab = "quantiles empiriques",
     ylab = "quantiles théoriques")
abline(a=0, b=1, col="red")
axis(2, tck = 1, lty = 2, col = "grey")
axis(1, tck=1, lty = 2, col = "grey")
# Graphiquement, le modèle semble adéquat.


# Mesure de l'adéquation du modèle gamma avec le Chi2 de pearson
sev.pred <- predict(gamma.model, newdata = SEV[-train.id, ], type = "response")

agg_pred_sev <- aggregate(sev.pred, by = list(class4[-train.id], tgroup[-train.id]), mean)
tbl_agg_sev <- cbind(agg_means_sev, agg_pred_sev$x, agg_pred_sev$x - agg_means_sev$x)
colnames(tbl_agg_sev) <- c("Class", "Territory", "Empirical mean", "theorical mean", "Écarts")
head(tbl_agg_sev)
(MSE <- mean((tbl_agg_sev$Écarts)^2))
sqrt(MSE) # L'adéquation est satisfaisante.

(Chi2_pearson <- sum((sev.pred - SEV[-train.id, ]$losspaid)^2 / var(sev.pred)))
qchisq(0.95, nrow(SEV[-train.id, ]) - 10)
# Le Chi-2 de Pearson est pourris...

ks.test(sev.pred, SEV[-train.id, ]$losspaid)
# En revenche, selon le test de Kolmogorov-Smirnov, le modèle est adéquat.

{
# mean.y <- mean(SEV$losspaid[train.id])
# var.y <- var(SEV$losspaid[train.id])
# mle_gamma <- constrOptim(c(mean.y^2/var.y, mean.y/var.y),
#             function(para) 
#                 -sum(log(dgamma(SEV$losspaid[train.id], para[1], para[2]))),
#             grad = NULL, ui=diag(2), ci=c(0,0))
# 
# # Test du ratio de vraisemblance
# sumry_gam$aic
# 2 * (2 + mle_gamma$value)
# (sumry_gam$aic - 20 - 2*mle_gamma$value)
# qchisq(0.95, 8)
# # Le modèle avec le glm est significativement meilleur.
}# Modèle de sévérité gamma construit sans glm


#==== Modèle d'aggrégation =================

dGaussian.copula <- function(n, Y, x, rho, z=1) {
  # Fonction qui calcule la densité d'une copule gaussienne
  if(typeof(x) == "list") x <- as.integer(unlist(x))
  if(typeof(Y) == "list") Y <- as.double(unlist(Y))
  
  p_R.1 <- prob(x)
  if (n == 0) return(1 - p_R.1)
  if (n > 10) return(1)
  
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
  if (d <= .Machine$double.eps) return(1)
  return(d)
}

randomize.Y <- function(data, n) {
  data <- data[data[,1] > n,]
  Y <- matrix(sapply(1:nrow(data), function(i) sample(unlist(data[,-1][[i]]), n)), ncol=n, byrow = T)
  cbind(unlist(data[,1]), Y)
}

spearmans.Rho <- function(nsim, data) {
  # Fonction qui retourne le rho de spearman pour une base de données d'assurance.
  cor.spearman <- matrix(ncol = 2, nrow=nsim)
  for (i in 1:nsim) {
    couples.NY <- randomize.Y(data, 1)
    couples.YY <- randomize.Y(data, 2)
    cor.spearman[i, 1] <- cor(rank(couples.NY[,1], ties.method = "random"), couples.NY[,2], method = "spearman")
    cor.spearman[i, 2] <- cor(couples.YY[,2], couples.YY[,3], method = "spearman")
  }
  
  mean_rho <- apply(cor.spearman, 2, mean)
  var_rho <- apply(cor.spearman, 2, var)
  IC_rho <- matrix(c(mean_rho - qnorm(0.975) * sqrt(var_rho),
                     mean_rho + qnorm(0.975)* sqrt(var_rho)),
                   ncol=2)
  return(list("mean"=mean_rho, "variance"=var_rho, "IC"=IC_rho))
}

optim_rho <- function(data, x, Depart, bornes, nsim=100, length.ech=1000){
  # Fonction qui permet de trouver les paramètres de dépendance
  # d'une copule gaussienne par optimisation numérique en utilisant un bootstrap.
  #
  # nsim correspond au nombre de simulations du bootstrap et 
  # length.ech correspond à la longueur des échantillons du bootstrap.
  
  # Bornes d'optimisation ui et ci.
  bounds <- matrix(bornes, ncol = 2, byrow = T)
  ui <- rbind( diag(2), -diag(2) )
  ci <- c( bounds[,1], - bounds[,2] )

  mle_rho <- matrix(ncol=2, nrow=nsim)
  for (j in 1:nsim) {
    # Bootstrap et optimisation numérique.
    sample.id <- sample(1:nrow(data), length.ech)
    
    mle_rho[j,] <- constrOptim(Depart, function(rho) {
      -sum(log(sapply(sample.id, function(i)
        dGaussian.copula(
          n = data[, 1][[i]],
          Y = unlist(data[,-1][[i]]),
          x = x[i,],
          rho = rho
        ))))
    },
    grad = NULL, ui = ui, ci = ci, outer.eps = 1e-2)$par
    
    # Depart <- mle_rho[j,]
    
    print(paste0(j / nsim * 100,"%", " : ", list(mle_rho[j, ])))
  }
  mean_rho <- apply(mle_rho, 2, mean, na.rm = T)
  var_rho <- apply(mle_rho, 2, var)
  IC_rho <- matrix(c(mean_rho - qnorm(0.975) * sqrt(var_rho),
                     mean_rho + qnorm(0.975)* sqrt(var_rho)),
                   ncol=2)
  return(list("mean"=mean_rho, "variance"=var_rho, "IC"=IC_rho))
}


# Création d'une base de données avec les données aggrégées pour entraîner la copule.

x.mat_train <- model.matrix(losspaid ~ class4 + tgroup, data=SEV[train.id,])
LOSS <- aggregate(SEV[train.id,]$losspaid, by=list(SEV[train.id,]$pol_id), list)
N <- aggregate(SEV[train.id,]$losspaid, by=list(SEV[train.id,]$pol_id), length)
x.mat_train <- aggregate(x.mat_train, list(SEV[train.id,]$pol_id), data.table::first)[,-1]
# nrow(LOSS) ; nrow(N) ; nrow(x.mat_train)


(Depart <- spearmans.Rho(10, cbind(N$x, LOSS$x)))
# unique(copula::rho(normalCopula(P2p(Sig_rho1.2.z1(Depart$mean, k - 1)), k, dispstr = "un")))
# # On voit que le rho de Spearman empirique est un excellent point de départ pour l'optimisation
# # numérique puisque pour la copule gaussienne, celui-ci est pratiquement identique à ses
# # paramètres de dépendance.
rho2 <- 2 * sin(pi * Depart$mean[2] / 6)
k <- max(N[,2])
mle_rho <- optimize(function(rho1)
  - sum(log(sapply(1:length(N[,1]), function(i)
      dGaussian.copula(
        n = N$x[[i]],
        Y = unlist(LOSS$x[[i]]),
        x = x.mat_train[i,],
        rho = c(rho1, rho2)
      )))),
  lower = -0.2, upper = sqrt(((k - 1) * rho2 + 1) / k))$minimum

rho <- c(mle_rho, rho2)

# mle_rho <- optim_rho(cbind(N$x, LOSS$x), x.mat_train,
#                      Depart$mean, c(0.02, 0.07, 0.09, 0.2), nsim = 1, length.ech = 10000)
#
# # Il a été démontré dans les simulations de la section 7 qu'il était beaucoup plus rapide
# # sans qu'il n'y ait toutefois de perte significative de précision en utilisant la méthode
# # d'inversion du rho de spearman pour rho2 et de trouver rho par maximum de vraisemblance univariée.


# Reproduction du tableau 6 - Résultats
rGaussian.copula <- function(n_sim, rho, k = max(unlist(N$x)), z = 1) {
    # Fonction qui permet de simuler des réalisations d'une copule gaussienne.
    Sig_rho <- ifelse(z == 1, Sig_rho1.2.z1, Sig_rho1.2.z2)
    rCopula(n_sim, normalCopula(P2p(Sig_rho(rho, k)), k+1, dispstr = "un"))
}

simul.S <- function(n_sim, x, rho, k = max(unlist(N$x)), z = 1) {
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
AGG$losspaid[AGG$losspaid < 0] <- 0
summary(AGG)
head(AGG)


S.empirique <- aggregate(AGG$losspaid, by = list(AGG$pol_id, AGG$tgroup, AGG$class4), sum)
mean.S.empirique <- aggregate(S.empirique$x,
                              by = list(S.empirique$Group.2, S.empirique$Group.3),
                              mean)
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
risk_groups <- unique(x.mat_train)
for (i in 1:30) {
  proposed[i,] <- unlist(simul.S(5e+3, unlist(risk_groups[i,]), rho))
}
Class <- as.matrix(risk_groups[1:5]) %*% 0:4 + 1
Territory <- as.matrix(risk_groups[c(1,6:10)]) %*% 0:5 + 1

tbl_6.theorique <- cbind(Territory, Class, proposed)
colnames(tbl_6.theorique) <- c("Territory", "Class", "E[S]", "VaR_0.995")
tbl_6.theorique <- data.table::data.table(tbl_6.theorique)
tbl_6.theorique <- data.table::setorderv(tbl_6.theorique, c('Territory', 'Class'), c(1, 1))
tbl_6 <- cbind(tbl_6.empirique, tbl_6.theorique[,3:4])
# Le modèle surestime la somme des réclamations...
mean((tbl_6$`mean(S)` - tbl_6$`E[S]`)^2)
mean((tbl_6$VaR_n - tbl_6$VaR_0.995)^2)

#==== Modèle de Tweedie =============================
AGG_2 <- S.empirique
colnames(AGG_2) <- c("pol_id", "class4", "tgroup", "losspaid")
head(AGG_2)
summary(AGG_2)


train.id <- FREQ[1:1e+6, 1]
AGG_train <- AGG_2[pol_id %in% train.id, ]
AGG_test <- AGG_2[!(pol_id %in% train.id), ]


prof <- tweedie.profile(losspaid ~ class4 + tgroup, data = AGG_train, 
                p.vec=seq(1.2,1.8,0.1), method = "saddlepoint", control=list(maxit=400))
prof$p.max=1.2
Tweedie.mod <- glm(losspaid ~ class4 + tgroup, data = AGG_train,
                  family = tweedie(prof$p.max, 0))
summary(Tweedie.mod)


x.agg.train <- model.matrix(losspaid ~ class4 + tgroup,
                          data = AGG_train)
x.agg.test <- model.matrix(losspaid ~ class4 + tgroup,
                            data = AGG_test)


agg_loss <- aggregate(AGG_train$losspaid, list(AGG_train$pol_id, AGG_train$class4, AGG_train$tgroup), sum)
agg_loss <- aggregate(agg_loss$x, list(agg_loss$Group.2, agg_loss$Group.3), mean)
expected_loss <- aggregate(fitted(Tweedie.mod), list(AGG_train$class4, AGG_train$tgroup), mean)

MSE <-  mean((agg_loss$x - expected_loss$x)^2)

tbl_result <- cbind(agg_loss, expected_loss$x)
colnames(tbl_result) <- c("tgroup", "class4", "Mean(S)", "E[S]")




