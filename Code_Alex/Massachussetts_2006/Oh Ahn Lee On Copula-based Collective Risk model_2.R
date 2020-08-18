# Document R servant à essayer les exemples donnés dans l'article 
# "On Copula-based Collective Risk Models"
# écrit par Rosy Oh, Jae Youn Ahn et Woojoo Lee.

library(pscl)
library(statmod)
library(tweedie)
library(copula)
library(utils)
library(tidyverse)
library(tidymodels)


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

setwd("~/Github/Recherche_ete_2019/Code_Alex/Massachussetts_2006")
# Réarrangement des données en facteurs
data_train <- read.csv("training_dat_NY.csv")
data_test <- read.csv("test_dat_NY.csv")

summary(data_train)

data_train <-
  data_train %>%
  mutate(tgroup = factor(tgroup),
         cgroup = factor(cgroup)) %>%
  mutate(id=X, .before=1) %>%
  mutate(X=NULL)

data_test <-
  data_test %>%
  mutate(tgroup = factor(tgroup),
         cgroup = factor(cgroup)) %>%
  mutate(id=X, .before=1) %>%
  mutate(X=NULL)

head(data_train)
head(data_test)
# tgroup : 1 = least risky, 6 = most risky
# class4 : A = Adult, B = Business, I = <3 yrs experience, M = 3-6 yrs experience, S = Senior

#---- Reproduction du tableau 4 ----
data_train$freq %>% table() %>% barplot()
data_train %>% filter(freq>0) %>% select(freq) %>% table() %>% barplot()

n_obs <- nrow(data_train)

# Matrice d'aggrégation de la moyenne de la fréquence.
tbl_mean_freq <-
  data_train %>% group_by(tgroup, cgroup) %>%
  summarise(mean_freq = mean(freq)) 

tbl_mean_freq. <- 
  tbl_mean_freq[,3] %>% as.matrix() %>% matrix(ncol = 5, byrow = T)

dimnames_tbl4 <- list(unique(tbl_mean_freq$tgroup),
                      unique(tbl_mean_freq$cgroup))

dimnames(tbl_mean_freq.) <- dimnames_tbl4
round(tbl_mean_freq., 3)
# J'arrive à reproduire le pattern attendu par rapport à ce qui est donné
# dans le tableau 4 du papier de Oh, Ahn et Lee 2019 même si les valeurs
# ne sont pas les mêmes.

# Matrice d'aggrégation de la distribution des features
tbl_proportions_parClasse <-
  data_train %>% group_by(tgroup, cgroup) %>%
  summarise(proportions. = length(id)/n_obs * 100)
tbl_proportions_parClasse. <-
  tbl_proportions_parClasse[, 3] %>% as.matrix() %>% matrix(ncol = 5, byrow = T)
dimnames(tbl_proportions_parClasse.) <- list(unique(tbl_mean_freq$tgroup),
                                             unique(tbl_mean_freq$cgroup))
round(tbl_proportions_parClasse., 2) # En pourcentages.
# J'ai une répartion plutôt similaire à ce qui est présenté dans le tableau
# 4 du papier de Oh, Ahn et Lee.

all_loss <- data_train %>% filter(freq > 0) %>%
  pivot_longer(cols=contains("loss"), values_drop_na=T, values_to="loss")

# summary(all_loss)

data_train %>% filter(freq>0) %>% select(loss1) %>% unlist() %>% hist()
all_loss %>% select(loss) %>% unlist() %>% hist()

# Matrice d'aggrégation de la moyenne de la sévérité
tbl_mean_sev <-
  all_loss %>% group_by(tgroup, cgroup) %>%
  summarise(mean_sev = mean(loss)) 

tbl_mean_sev. <- 
  tbl_mean_sev[,3] %>% as.matrix() %>% matrix(ncol = 5, byrow = T)

dimnames(tbl_mean_sev.) <- dimnames_tbl4
round(tbl_mean_sev., 0)
# La moyenne des sévérités ne correspond pas au tableau 4 de l'article. 
# Toutefois certaines valeurs sont similaires.


#---- Analyse de la fréquence ===================
Poisson.model <- glm(freq ~ tgroup + cgroup, data=data_train,
                     family = poisson("log"))
(Summary_poisson <- summary(Poisson.model))
# J'ai les mêmes poids et les mêmes écarts types que dans l'article.

# hurdle.model <- hurdle(freq ~ tgroup + cgroup, data=data_train)
# (Summary_hurdle <- summary(hurdle.model))
# 
# AIC(hurdle.model, Poisson.model)
# 1-pchisq(Summary_poisson$aic-2*10 - 2*hurdle.model$loglik, 10)
# # Oh, Ahn et Lee avait pris le modèle Poisson simple. On reprend donc 
# # ce modèle aussi, même s'il s'ajuste moins bien aux données que le modèle
# # d'Hurdle.

F1 <- function(n, x) {
  mu <- predict(Poisson.model, newdata = x, type="response")
  ppois(n, mu)
}
F1.inv <- function(p, x) {
  mu <- predict(Poisson.model, newdata = x, type="response")
  qpois(p, mu)
}

#---- Analyse de la sévérité =================
gamma.model <- glm(loss ~ tgroup + cgroup, data = all_loss,
                   family = Gamma(link = "log"))


summary_gamma <- summary(gamma.model)
round(summary_gamma$coefficients, 3)
(Y.disp <- 1/MASS::gamma.dispersion(gamma.model))
(Y.coef <- summary_gamma$coefficients)

F2 <- function(Y, x) {
  mu <- predict(gamma.model, newdata = x, type="response")
  pgamma(Y, Y.disp, Y.disp/mu)
}
f2 <- function(Y, x) {
  mu <- predict(gamma.model, newdata = x, type="response")
  dgamma(Y, Y.disp, Y.disp/mu)
}
F2.inv <- function(p, x) {
  mu <- predict(gamma.model, newdata = x, type="response")
  qgamma(p, Y.disp, Y.disp/mu)
}

#--- Mesures de l'adéquation du modèle gamma 
predictions <- predict.glm(gamma.model, newdata = all_loss, "response")
tbl_mean_sev.pred <- cbind(all_loss, predictions) %>%
  group_by(tgroup, cgroup) %>% summarise(predictions = mean(predictions))

tbl_mean_sev.pred. <- 
  tbl_mean_sev.pred[,3] %>% as.matrix() %>% matrix(ncol = 5, byrow = T)

dimnames(tbl_mean_sev.pred.) <- dimnames_tbl4

round(tbl_mean_sev., 0)
round(tbl_mean_sev.pred., 0)
mean((tbl_mean_sev. - tbl_mean_sev.pred.)^2)


# Analyse graphique de l'adéquation
U <- ecdf(all_loss$loss)(all_loss$loss)

qqplot(all_loss$loss, F2.inv(U, all_loss))
abline(a=0, b=1, col="red")

# Kolmogorov-Smirnov test
ks.test(all_loss$loss,
        F2.inv(U, all_loss))
# Le modèle Gamma ne s'ajuste pas bien. Malgré tout, c'est celui que l'on va
# prendre pour comparer avec l'article.


#---- Mesure de la corrélation =================

# Calcul de la corrélation entre N et Y (Rho1)
(all_loss %>% select(freq, loss) %>% cor(method = "spearman"))[2]
(data_train %>% select(freq, loss1) %>% cor(use="complete.obs", method = "spearman"))[2]

cor_boostrap <- function(split) {
  dat <- analysis(split)
  tibble(
    term = "rho1",
    estimate = cor(dat$freq, dat$loss, method = "spearman"),
    std.err = NA_real_
  )
}
set.seed(2020)
bootstrap_rho1 <- all_loss %>%  select(freq, loss) %>%
  bootstraps(n = 1E+4, apparent = T)  %>%
  mutate(rho1 = map(splits, cor_boostrap)) %>%
  int_pctl(rho1)

bootstrap_rho1$.estimate
1-pchisq((length(all_loss$loss)+1) * bootstrap_rho1$.estimate^2, 1)
# Selon le test de Mantel-Haenszel, la corrélation n'est pas significativement
# non nulle. L'intervalle de confiance vient confirmer le résultat du test.


# Calcul de la corrélation entre les Y (Rho2)
data_train  %>% filter(freq>1) %>%
  select(contains("loss")) %>% as.matrix() %>%
  cor(use = "pairwise.complete.obs")

cor_boostrap <- function(split) {
  dat <- analysis(split)
  tibble(
    term="rho2",
    estimate = cor(dat$loss1, dat$loss2, use="complete.obs"),
    std.err = NA_real_
  )
}
set.seed(2020)
(bootstrap_rho2 <- data_train %>% filter(freq>1) %>%
    bootstraps(n=1E+4, apparent = T)  %>%
    mutate(rho2 = map(splits, cor_boostrap)) %>%
    int_pctl(rho2))

bootstrap_rho2$.estimate
1-pchisq((length(data_train$loss2)+1) * bootstrap_rho2$.estimate^2, 1)
# Selon le test de Mantel-Haenszel, la corrélation est significativement
# non nulle. Cependant, l'intervalle de confiance de la méthode Bootstrap
# inclut 0 au seuil de 5%

rho <- c(bootstrap_rho1$.estimate, bootstrap_rho2$.estimate)
# rho <- c(-0.018, 0.027)
# # Pour la suite de l'analyse, on va reprendre les résultats obtenus dans
# # l'article pour tenter de reproduire les résultats.


#==== Reproduction du tableau 6 - Résultats ----

VaR <- function(x) sort(x)[ceiling(length(x) * 0.995)]

rGaussian.copula <- function(n_sim, rho, k, z = 1) {
  # Fonction qui permet de simuler des réalisations d'une copule gaussienne.
  Sig_rho <- ifelse(z == 1, Sig_rho1.2.z1, Sig_rho1.2.z2)
  rCopula(n_sim, normalCopula(P2p(Sig_rho(rho, k)), k + 1, dispstr = "un"))
}

simul.S <- function(n_sim, x, rho, k = F1.inv(1-1e-6, x)+1, z = 1) {
    # Fonction qui permet de faire de la simulation Monte-Carlo
    # pour estimer les mesures de risque théorique du modèle.
    U <- rGaussian.copula(n_sim, rho, k, z)

    Couples <- cbind("N" = F1.inv(U[, 1], x),
                     "Y" = matrix(F2.inv(U[, -1], x),
                                  nrow = n_sim,
                                  byrow = T))
    
    S <- numeric(n_sim)
    for (i in 1:n_sim) {
      N <- Couples[i,1]
        if (N == 0)
            S[i] <- 0
        else
            S[i] <- sum(Couples[i, 2:(N + 1)])
    }
    
    mean_S <- mean(S)
    VaR_0.995 <- VaR(S)
    return(c(mean_S, VaR_0.995))
}


# Mesures empiriques
total_loss <- data_train %>%
  pivot_longer(cols=contains("loss"), values_drop_na=F, values_to="loss") %>%
  group_by(id, tgroup, cgroup) %>% summarise(S = sum(loss, na.rm = T))

tbl_6.empirique <- total_loss %>% group_by(tgroup, cgroup) %>%
  summarise(mean_S.emp=mean(S), VaR_S.emp=VaR(S))


# Mesures théoriques
proposed <- matrix(ncol=2, nrow=30)
features <- data_train %>% group_by(tgroup, cgroup) %>% summarise()
pb <- txtProgressBar(min=1, max=30, style=3)
set.seed(2020)
for (i in 1:30) {
  proposed[i,] <- simul.S(1E+4, features[i,], rho)
  setTxtProgressBar(pb, i)
}

tbl_6.theorique <- tibble(proposed[,1], proposed[,2],
                          .name_repair = ~c("mean_S.th", "VaR_S.th"))
tbl_6 <- cbind(tbl_6.empirique, tbl_6.theorique)


#---- Modèle de Tweedie ----

# prof <- tweedie.profile(S ~ tgroup + cgroup, data = total_loss,
#                         verbose=T, method="saddlepoint")
# xi = prof$xi.max
xi = 1.665306

Tweedie.mod <- glm(S ~ tgroup + cgroup, data = total_loss,
                  family = tweedie(xi, 0))
(Summary_Tweedie <- summary(Tweedie.mod))
phi = Summary_Tweedie$dispersion

mu <- predict.glm(Tweedie.mod, newdata = total_loss, type = "response")
U <- ecdf(total_loss$S)(total_loss$S)
ks.test(U, 
        ptweedie(total_loss$S, xi=xi, mu=mu, phi = phi))
# Le modèle ne s'ajuste pas bien, mais bon...on va faire comme dans l'article...


tbl_6.tweedie <- cbind(total_loss, fit=fitted(Tweedie.mod)) %>% group_by(tgroup, cgroup) %>%
  summarise(mean_S.tweed=mean(fit))

VaR_S.tweed <- sapply(1:30, function(i)
  qtweedie(
    0.995,
    xi = xi,
    mu = predict.glm(Tweedie.mod, newdata = tbl_6.tweedie[i, ], type = "response"),
    phi = phi
  ))

tbl_6.tweedie <- bind_cols(tbl_6.tweedie, VaR_S.tweed=VaR_S.tweed)
tbl_6 <- bind_cols(tbl_6, tbl_6.tweedie[,3:4])

#---- Calcul du MSE ====
# Avec les données d'entraînement
MSE <- matrix(ncol=4, nrow=2)
dimnames(MSE) <- list(c("train", "test"), 
                      c("mean.th", "mean.tweedie", "VaR.th", "VaR.tweedie"))
MSE[1,1] <- mean((unlist(tbl_6[,3] - tbl_6[,5]))^2)
MSE[1,2] <- mean((unlist(tbl_6[,3] - tbl_6[,7]))^2)
MSE[1,3] <- mean((unlist(tbl_6[,4] - tbl_6[,6]))^2)
MSE[1,4] <- mean((unlist(tbl_6[,4] - tbl_6[,8]))^2)

# Avec les données de test
total_loss <- data_test %>%
  pivot_longer(cols=contains("loss"), values_drop_na=F, values_to="loss") %>%
  group_by(id, tgroup, cgroup) %>% summarise(S = sum(loss, na.rm = T))

tbl_6.empirique <- total_loss %>% group_by(tgroup, cgroup) %>%
  summarise(mean_S.emp=mean(S), VaR_S.emp=VaR(S))

tbl_6_test <- bind_cols(tbl_6.empirique, tbl_6.theorique, tbl_6.tweedie[,3:4])

MSE[2,1] <- mean((unlist(tbl_6_test[,3] - tbl_6_test[,5]))^2)
MSE[2,2] <- mean((unlist(tbl_6_test[,3] - tbl_6_test[,7]))^2)
MSE[2,3] <- mean((unlist(tbl_6_test[,4] - tbl_6_test[,6]))^2)
MSE[2,4] <- mean((unlist(tbl_6_test[,4] - tbl_6_test[,8]))^2)


#---- Résultats ----
tbl_6
# tbl_6_test
MSE

#==== Appréciation graphique des résultats (Figure 3) ====
# Avec les données d'entraînement
plot.mean <- 
  bind_cols(tbl_6, Risk_group = 1:30) %>%
  pivot_longer(cols = contains("mean"), values_to = "mean") %>%
  ggplot(mapping = ) +
  geom_path(aes(
    x = Risk_group,
    y = mean,
    group = name,
    colour = name
  )) +
  scale_color_discrete(name = "Modèle",
                       labels = c("Train data", "Proposed", "Tweedie"))


plot.VaR <- 
  bind_cols(tbl_6, Risk_group = 1:30) %>%
    pivot_longer(cols = contains("VaR"), values_to = "VaR") %>%
    ggplot() + geom_path(aes(
      x = Risk_group,
      y = VaR,
      group = name,
      colour = name
    )) +
    scale_color_discrete(name = "Modèle",
                         labels = c("Train data", "Proposed", "Tweedie"))

gridExtra::grid.arrange(plot.mean, plot.VaR, nrow=2)

# Avec les données de test
plot.mean <- 
  bind_cols(tbl_6_test, Risk_group = 1:30) %>%
    pivot_longer(cols = contains("mean"), values_to = "mean") %>%
    ggplot(mapping = ) +
    geom_path(aes(
      x = Risk_group,
      y = mean,
      group = name,
      colour = name
    )) +
    scale_color_discrete(name = "Modèle",
                         labels = c("Test data", "Proposed", "Tweedie"))


plot.VaR <-
  bind_cols(tbl_6_test, Risk_group = 1:30) %>%
    pivot_longer(cols = contains("VaR"), values_to = "VaR") %>%
    ggplot() + geom_path(aes(
      x = Risk_group,
      y = VaR,
      group = name,
      colour = name
    )) +
    scale_color_discrete(name = "Modèle",
                         labels = c("Test data", "Proposed", "Tweedie"))

gridExtra::grid.arrange(plot.mean, plot.VaR, nrow=2)
