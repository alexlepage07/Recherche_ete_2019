# Normale
epsilon <- rnorm(1e+4)
mat <- matrix(1,nrow = 1e+4, ncol=1)

mu <- mean(epsilon)
sig <- sqrt(var(epsilon))
par <- c(mu,sig)
mle <- constrOptim(par,
                   function(par)
                       - sum(log(dnorm(epsilon, par[1], par[2]))),
                   grad = NULL, ui = c(0,1), ci = 0)
mle$par

# Poisson
lambda <- exp(2 + X * 1 + epsilon)
NN <- rpois(1e+4, lambda)
mat <- matrix(1,nrow = 1e+4, ncol=1)
lambda <- mean(NN)
log(lambda)

eta <- function(beta0, beta1) beta0 + beta1 * X + epsilon
mle <- optim(c(1,1), function(beta) -sum(log(dpois(NN, exp(eta(beta[1], beta[2]))))))
mle$par
mle$par[1] + mle$par[2] * X
X <- rnorm(1e+4)

modele <- glm(NN~X,family = poisson(), offset = 1/2)
summary(modele)


data("freMTPLfreq")
modele <- glmnet(model.matrix(ClaimNb~.-PolicyID,
              data = freMTPLfreq), freMTPLfreq$ClaimNb,
              family = "poisson")
plot(modele)
plot(modele$dev.ratio)
modele$lambda
cv <- cv.glmnet(model.matrix(ClaimNb~.-PolicyID,
                       data = freMTPLfreq), freMTPLfreq$ClaimNb,
          family = "poisson")
log(cv$lambda.min)

plot(cv)
summary(modele)
drop1(modele)
model.matrix(ClaimNb~.-PolicyID, freMTPLfreq)

library(glmnet)
library(mgcv)
# zip : zero inflated poisson
# gam : 
