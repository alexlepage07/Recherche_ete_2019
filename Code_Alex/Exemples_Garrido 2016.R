# Garrido et al., 2016 - Generalized linear models for dependent
# frequency and severity of insurance claims.
#
# Simulations

# Scénario 1:
x <- abs(rnorm(nsim <- 1e+4))
v <- exp(1 + x / 2)
N <- rpois(nsim, v)
mu.func <- function(theta) exp(5 + 5 * x + theta * N)
mu <- sapply(seq(-0.5, 0.5, length.out=101), mu.func)

S.func <- function(mu) {
    sapply(1:nsim, function(i)
        ifelse(N[i] == 0,
               0,
               sum(rgamma(N[i], 1 / 2, 1 / mu / 2))))
}
glm(S.func(mu[1])~
