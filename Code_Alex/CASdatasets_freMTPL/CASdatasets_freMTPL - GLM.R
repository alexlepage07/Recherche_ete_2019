# Dans ce projet R, on cherche à appliquer le modèle collectif
# du risque avec structure de dépendance. Cela implique de
# trouver la copule qui représente le mieux cette structure,
# À identifier les lois marginales et à estimer les paramètres
# du modèle.

library(CASdatasets)
library(ggplot2)
library(glmnet)
library(MASS)
library(pscl)


data("freMTPLfreq")
data("freMTPLsev")

summary(freMTPLfreq[, -1])

n_obs <- nrow(freMTPLfreq)
test <- sample(1:n_obs, 0.2 * n_obs)
train <- (1:n_obs)[-test]


# Test de la régression de Lasso. ----
x <- model.matrix(ClaimNb ~ . - Density + I(log(Density)) - Exposure +
                     offset(log(Exposure)), data = freMTPLfreq[train,-1])
y <- freMTPLfreq$ClaimNb[train]

x.test <- model.matrix(ClaimNb ~ . - Density + I(log(Density)) - Exposure +
                     offset(log(Exposure)), data = freMTPLfreq[test,-1])
y.test <- freMTPLfreq$ClaimNb[test]

grid <- 10 ^ seq(3, -2, leng = 100)
lasso.mod <-
    glmnet(
        x,
        y,
        alpha = 1,
        lambda = grid,
        family = "poisson"
    )
cv.lasso <-
    cv.glmnet(x,
              y,
              alpha = 1,
              family = "poisson")
plot(cv.lasso)

coef(cv.lasso, s="lambda.min", family="poisson")
coef(cv.lasso, s=exp(-7), family="poisson")
coef(cv.lasso, s="lambda.1se", family="poisson")


# lasso.pred <-
#     predict(lasso.mod, s = cv.lasso$lambda.min, newx = x.test)
# sum(((lasso.pred - y.test) ^ 2) / lasso.pred)
qchisq(0.99, length(y.test) - 31)

#----


# Comparaison de modèles ----
mod.poisson <-
    glm(
        ClaimNb ~ . - Density + I(log(Density)) - Exposure + offset(log(Exposure)),
        data = freMTPLfreq[train,-1],
        family = poisson
    )
sum_0 <- summary(mod.poisson)
mod.poisson1 <- glm(ClaimNb ~ 1 + offset(log(Exposure)),
                    data = freMTPLfreq[train,-1],
                    family = poisson)
sum_1 <- summary(mod.poisson1)

AIC(mod.poisson, mod.poisson1)
df <- sum_0$df[1] - sum_1$df[1]
# H0: Le modèle le plus simple est adéquat. H1: Le modèle le plus simple est moins bon.
sum_1$deviance - sum_0$deviance
qchisq(0.99, df)
# Comme la statistique est supérieure à la valeur critique, on rejette H0.


mod.nbinom <-
    glm.nb(ClaimNb ~ . - Density + I(log(Density)) - Exposure + offset(log(Exposure)),
           data = freMTPLfreq[train,-1])

mod.Pois.zeroinf <-
    zeroinfl(
        ClaimNb ~ DriverAge + CarAge + Gas + Power + Brand + I(log(Density)) +
            Region + offset(log(Exposure)),
        data = freMTPLfreq[train,-1], dist = "poisson"
    )

mod.nbinom.zeroinf <-
    zeroinfl(
        ClaimNb ~ DriverAge + CarAge + Gas + Power + Brand + I(log(Density)) +
            Region + offset(log(Exposure)),
        data = freMTPLfreq[train,-1], dist = "negbin"
    )

AIC(mod.poisson, mod.nbinom, mod.Pois.zeroinf, mod.nbinom.zeroinf)
qchisq(0.99, 1)


Y_predict <- predict(mod.nbinom, type="response", newdata=x.test)
mean((Y_predict - y.test)^2)


freq <- freMTPLfreq
freq$Power[freq$Power %in% c("k", "l", "m", "n", "o")] <- "k"
summary(freq[, -1])
