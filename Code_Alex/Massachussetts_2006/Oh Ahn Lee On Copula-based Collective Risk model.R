library(pscl)
library(rlist)

setwd("C:/Users/Alex/Desktop/Recherche_ete_2019/Code_Alex/Massachussetts_2006")
FREQ <- read.csv("FREQ.txt")
SEV <- read.csv("SEVERITY.txt")


# tgroup : 1 = least risky, 6 = most risky
# class4 : A = Adult, B = Business, I = <3 yrs experience, M = 3-6 yrs experience, S = Senior

#========================== Analyse de la fréquence ===================
# Réarrangement des données en facteurs
FREQ$class4 <- FREQ$class4 %% 10
FREQ$class4[FREQ$class4 == 1] <- "A"
FREQ$class4[FREQ$class4 == 2] <- "S"
FREQ$class4[FREQ$class4 %in% c(3, 4)] <- "M"
FREQ$class4[FREQ$class4 == 5] <- "B"
FREQ$class4[FREQ$class4 %in% (6:9)] <- "I"
FREQ$class4 <- factor(FREQ$class4)

FREQ$tgroup <- factor(FREQ$tgroup, 1:6)
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

# hurdle.model <-
#     hurdle(CompteDeclm_id ~ class4 + tgroup,
#            data = FREQ[train.id, ],
#            offset = log(earnexpo))
# (Sumry_1 <- summary(hurdle.model))


hurdle.model_2 <-
    hurdle(CompteDeclm_id ~ class4 + tgroup | 1,
           data = FREQ[train.id, ] ,
           offset = log(earnexpo))
(Sumry_2 <- summary(hurdle.model_2))


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


freq.pred <- predict(hurdle.model_2, newdata = FREQ[test.id, 3:5],
                     type = "response")
ecarts <- freq.pred - FREQ[test.id, 2]
head(data.frame(freq.pred, FREQ[test.id, 2:5], ecarts))
tail(data.frame(freq.pred, FREQ[test.id, 2:5], ecarts))

agg_pred_freq <- aggregate(freq.pred, by = list(class4[test.id], tgroup[test.id]), mean)
tbl_agg_freq <- cbind(agg_means_freq, agg_pred_freq$x)
colnames(tbl_agg_freq) <- c("Class", "Territory", "Empirical mean", "theorical mean")
tbl_agg_freq
mean((tbl_agg_freq$`Empirical mean` - tbl_agg_freq$`theorical mean`)^2)

# Test d'adéquation du Chi2 de pearson.
Q <- 0
for (n in unique(FREQ[test.id, 2])){
    Q <- Q + (sum(FREQ[test.id, 2] == n) - sum(freq.pred == n)) ^ 2 /
        sum(FREQ[test.id, 2] == n)
}
chi2 <- rbind(c(
    Q,
    qchisq(0.95, length(unique(FREQ[test.id, 2])) - 10, lower.tail = T),
    pchisq(Q, length(unique(FREQ[test.id, 2])) - 10, lower.tail = F)
))
colnames(chi2) <- c("Statistique", "valeur critique" , "P-value")
chi2
# Le modèle de fréquence n'est pas terrible...

p <- exp(hurdle.model_2$coefficients$zero)
p <- p / (1 + p)

#========================= Analyse de la sévérité =================
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
x.mat_train <- model.matrix(losspaid~ class4 + tgroup, data=SEV[train.id,])
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
                   family = Gamma)
(sumry_gam <- summary(gamma.model))
# Je suis loin des coefficients trouvés par Oh, Ahn et Lee... 


# Analyse de l'adéquation graphiquement
x.mat_test <- model.matrix(losspaid~ class4 + tgroup, data=SEV[-train.id,])
U <- pgamma(SEV[-train.id, ]$losspaid,
            1/sumry_gam$dispersion,
            (x.mat_test %*% gamma.model$coefficients)*sumry_gam$dispersion)
q <- qgamma(U, 
            1/sumry_gam$dispersion,
            (x.mat_test %*% gamma.model$coefficients)*sumry_gam$dispersion)
plot(SEV[-train.id, ]$losspaid, q, type = "l",
     xlab = "quantiles empiriques",
     ylab = "quantiles théoriques")
abline(a=0, b=1, col="red")
axis(2, tck = 1, lty = 2, col = "grey")
axis(1, tck=1, lty = 2, col = "grey")
# Graphiquement, le modèle semble adéquat.


# Mesure de l'adéquation du modèle gamma avec le Chi2 de pearson
sev.pred <- predict(gamma.model, newdata = SEV[-train.id, ], type = "response")
ecarts <- sev.pred - SEV[-train.id, ]$losspaid
head(data.frame(sev.pred, SEV[-train.id, 2:4], ecarts))
(Chi2_pearson <- sum(ecarts^2 / var(sev.pred)))
qchisq(0.99, nrow(SEV[-train.id, ]) - 10)
# Selon le chi2 de Pearson, le modèle n'est pas adéquat.

# Mesure de l'adéquation du modèle gamma avec le test de Kolmogorov-Smirnov
ks.test(sev.pred, SEV[-train.id, ]$losspaid)
# Selon le test de Kolmogorov-Smirnov, le modèle est adéquat.


# # Modèle de sévérité gamma construit sans glm
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

#======================= Modèle d'aggrégation =================

# Correlation Matrix for Frequency and Individual Severities
# avec z = 1 : Equicorrelation matrix 
Sig_rho2.z1 <- function(rho, k) {
     
    matrix(c(rep(c(1, rep(rho[2], k)), k-1),1), nrow = k, byrow = T)
}
Sig_rho1.2.z1 <- function(rho, k) {
    mat <- rbind(rep(rho[1], k), Sig_rho2.z1(rho, k))
    mat <- cbind(c(1, rep(rho[1], k)), mat)
    return(mat)
}

# Avec z = 2 : Autoregressive correlation matrix
Sig_rho2.z2 <- function(rho, k) {
    outer(0:(k-1), 0:(k-1), function(i,j) rho[2]^abs(i-j))
}
Sig_rho1.2.z2 <- function(rho, k) {
    mat <- rbind(rep(rho[1], k), Sig_rho2.z2(rho, k))
    mat <- cbind(c(1, rep(rho[1], k)), mat)
    return(mat)
}

mu.norm <- function(Y, x, F2, rho, z=1) {
    k <- length(Y)
    Sig_rho <- ifelse(z == 1, Sig_rho2.z1, Sig_rho2.z2)
    rho[1] * t(rep(1, k)) %*% solve(Sig_rho(rho, k)) %*% (qnorm(F2(Y, x)))
}
sig.norm <- function(Y, rho, z=1) {
    k <- length(Y)
    Sig_rho <- ifelse(z == 1, Sig_rho2.z1, Sig_rho2.z2)
    sqrt(1 - rho[1]^2 * t(rep(1, k)) %*% solve(Sig_rho(rho, k)) %*% rep(1, k))
}

F2 <- function(y, x) {
    x <- matrix(x, ncol = length(x))
    pgamma(y, 1/sumry_gam$dispersion, (x %*% (gamma.model$coefficients)) * sumry_gam$dispersion)
}
F1 <- function(n, x) {
    x <- matrix(x, ncol = length(x))
    lambda <- exp(x %*% hurdle.model_2$coefficients$count) / p
    (1 - p) + p * (ppois(n, lambda) - dpois(0, lambda)) / (1-dpois(0, lambda))
}


dGaussian.copula <- function(n, Y, x, F1, F2, rho, z=1) {
    if (n == 0)
        return(1-p)
    k <- length(Y)
    Sig <- sig.norm(Y, rho, z)
    mu <- mu.norm(Y, x, F2, rho, z)
    Sig_rho <- ifelse(z == 1, Sig_rho2.z1, Sig_rho2.z2)
    
    c_Y <- det(Sig_rho(rho, k)) ^ (-0.5) * exp(- t(qnorm(F2(Y, x))) %*%
                                                   (Sig_rho(rho, k) - diag(1, k)) %*%
                                                   qnorm(F2(Y, x)) / 2)
    c_Y <- c_Y * prod(F2(Y, x))
    
    d <- p * c_Y * (pnorm((qnorm(F1(n, x)) - mu) / Sig) - 
               pnorm((qnorm(F1(n-1, x)) - mu) / Sig))
    return(d)
}



# Création d'une base de données avec les données aggrégées pour entraîner la copule.
x.mat_train <- model.matrix(losspaid~ class4 + tgroup, data=SEV[train.id,])

LOSS <- aggregate(SEV[train.id,]$losspaid, by=list(SEV[train.id,]$pol_id), list)
N <- aggregate(SEV[train.id,]$losspaid, by=list(SEV[train.id,]$pol_id), length)


# Bornes d'optimisation
bounds <- matrix(c(-0.5, -1e-6,
                   1e-06, 0.5)
                   , ncol = 2, byrow = T)
# Convertir les constraintes en matrices ui et ci.
n_par <- nrow(bounds)
ui <- rbind( diag(n_par), -diag(n_par) )
ci <- c( bounds[,1], - bounds[,2] )
# Retirer les valeurs infinies
i <- as.vector(is.finite(bounds))
ui <- ui[i,]
ci <- ci[i]

mle_rho <- matrix(ncol=2, nrow=31)
mle_rho[1,] <- c(-0.018, 0.027)
for (j in 2:31) {
    sample.id <- sample(1:nrow(N), 100)
    
    mle_rho[j,] <- constrOptim(mle_rho[j-1,], function(rho) {
        -sum(log(sapply(sample.id, function(i)
            dGaussian.copula(N[i,]$x, unlist(LOSS[i,]$x), x.mat_train[i,], F1, F2, rho))))},
        grad = NULL, ui = ui, ci = ci, outer.eps = 1e-2)$par
    
    print(j)
}
mean_rho <- apply(mle_rho, 2, mean)
var_rho <- apply(mle_rho, 2, var)
IC_rho <- matrix(c(mean_rho - qnorm(0.975) * var_rho,
                   mean_rho + qnorm(0.975)* var_rho),
                 ncol=2)

# Reproduction du tableau 6 - Résultats


