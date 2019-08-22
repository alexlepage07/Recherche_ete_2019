x <- runif(nsim <- 1e+5)
w <- rbinom(nsim, 1, 0.5)
x <- sapply(1:nsim, function(i) qexp(x[i], exp(- t(c(0.5, 0.2)) %*% c(1, w[i]))))
summary(x)
u <- 5
x[x >= u] <- u
hist(x)


loglik <- function(beta, y, x, u){
    l <- 0
    for (i in 1:(n <- length(y))){
        l <- l + ifelse(y[i] < u, 
                        log(dexp(y[i], exp(-beta %*% c(1, x[i])))),
                        log(1-pexp(y[i], exp(-beta %*% c(1, x[i])))))
    }
    l
}
loglik(c(0.5, 0.2), x, w, u)

bounds <- matrix(c(-Inf, Inf,
                   -Inf, Inf), ncol = 2, byrow = T)
n_par <- nrow(bounds)
ui <- rbind( diag(n_par), -diag(n_par) )
ci <- c( bounds[,1], - bounds[,2] )
# Retirer les valeurs infinies
i <- as.vector(is.finite(bounds))
ui <- ui[i,]
ci <- ci[i]

constrOptim(c(0.5,0.2), function(beta) - loglik(beta, x, w, u),
            grad = NULL, ui = ui, ci = ci)$par

glm(x~w, family = Gamma("log"))
