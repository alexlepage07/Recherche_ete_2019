# Dans ce projet R, on cherche à appliquer le modèle collectif 
# du risque avec structure de dépendance. Cela implique de
# trouver la copule qui représente le mieux cette structure,
# À identifier les lois marginales et à estimer les paramètres
# du modèle.

library(CASdatasets)
library(ggplot2)
library(ReIns)
library(actuar)


data("freMTPLfreq")
data("freMTPLsev")
freq  <- freMTPLfreq$ClaimNb
loss <- freMTPLsev$ClaimAmount


#=============== Trouver la loi marginale de la v.a. de fréquence (N). =============
# Analyse des données
table(freq)
n_obs <- length(freq)
summary(freq)
max_freq <- max(freq)
plot(ecdf(freq), xlim=c(0,5), main="")
barplot(table(freq))


# Cas 1 : Loi de poisson
par <- mean(mean(freq), var(freq))
par_mle <- optimize(function(par) -sum(log(dpois(freq, par))), lower = 0.001, upper = 1)$minimum
rbind(par, par_mle)

datas <- data.frame(c(freq, qpois((0:100)/100, par)),
                    Source <- c(rep("Empirique", length(freq)),rep("Théorique", 101)))

ggplot(datas,  aes(datas[, 1], fill = Source)) + 
    geom_histogram(alpha = 0.3, aes(y = ..density..), position = 'identity', binwidth = 1)+
    xlab("n") + ylab("Probabilité") +
    theme(legend.title = element_blank())

Q <- 0
for (n in unique(freq)){
    Q <- Q + (sum(freq == n) - n_obs * dpois(n, par)) ^ 2 / sum(freq == n)
}
chi2 <- rbind(
              "Poisson" = c(
                  Q,
                  qchisq(0.95, length(unique(freq)) - 2, lower.tail = T),
                  pchisq(Q, length(unique(freq)) - 2, lower.tail = F)
              ))

# Loi de poisson inflationnée à zéro
par <- numeric(2)
par[1] <- mean(mean(freq), var(freq))
par[2] <- (1 - mean(freq==0)) / (1 - exp(-par[1]))


dpois.zero.inf <- function(n, lambda, prob){
    res <- numeric(length(n))
    for (i in 1:length(n)){
        if (n[i] == 0)
            res[i] <- (1 - prob) + prob * exp(-lambda)
        else
            res[i] <- prob * dpois(n[i], lambda)
        }
    return(res)
}

ppois.zero.inf <- function(x, lambda, prob){
    sum(dpois.zero.inf(0:x, lambda, prob))
}

qpois.zero.inf <- function(p, lambda, prob) {
    q <- numeric(length(p))
    for (i in 1:length(p)) {
        if (p[i] <= dpois.zero.inf(0, lambda, prob))
            q[i] <- 0
        else
            q[i] <- uniroot(function(x)
                ppois.zero.inf(x, lambda, prob) - p[i],
                c(0, qpois(p[i], lambda)))$root
    }
    return(ceiling(q))
}


mle <- optim(par, function(par) -sum(log(dpois.zero.inf(freq, par[1], par[2]))))
tbl_parametres <- rbind("moments"=par, "mle"=mle$par)
colnames(tbl_parametres) <- c("lambda", "pi")
tbl_parametres

datas <- data.frame(c(freq, qpois.zero.inf((0:99.99)/100, mle$par[1], mle$par[2])),
                    Source <- c(rep("Empirique", length(freq)),rep("Théorique", 100)))

ggplot(datas,  aes(datas[, 1], fill = Source)) + 
    geom_histogram(alpha = 0.3, aes(y = ..density..), position = 'identity', binwidth = 1)+
    xlab("n") + ylab("Probabilité") +
    theme(legend.title = element_blank())


Q <- 0
for (n in unique(freq)){
    Q <- Q + (sum(freq == n) - n_obs * dpois.zero.inf(n, mle$par[1], mle$par[2])) ^ 2 / sum(freq == n)
}
chi2 <- rbind(chi2,
    "Poisson inflationné à zéro" = c(
        Q,
        qchisq(0.95, length(unique(freq)) - 3, lower.tail = T),
        pchisq(Q, length(unique(freq)) - 3, lower.tail = F)
    ))

# Cas 2 : Loi binomiale négative
par <- numeric(2)
par[1] <- mean(freq) / var(freq)
par[2] <- mean(freq)

mle <- constrOptim(par,
                   function(par) -sum(log(dnbinom(freq, par[1], mu=par[2]))),
                   grad = NULL, ui=diag(2), ci = c(0,0))
par <- mle$par

datas <- data.frame(c(freq, qnbinom((0:100)/100, par[1], mu=par[2])),
                    Source <- c(rep("Empirique", length(freq)),rep("Théorique", 101)))

ggplot(datas,  aes(datas[, 1], fill = Source)) + 
    geom_histogram(alpha = 0.3, aes(y = ..density..), position = 'identity', binwidth = 1)+
    xlab("n") + ylab("Probabilité") +
    theme(legend.title = element_blank())


Q <- 0
for (n in unique(freq)){
    Q <- Q + (sum(freq == n) - n_obs * dnbinom(n, par[1], mu=par[2])) ^ 2 /
        sum(freq == n)
}

chi2 <- rbind(chi2,
              "Binomiale négative" = c(Q,
                                       qchisq(0.95, length(unique(freq)) - 3, lower.tail = T),
                                       pchisq(Q, length(unique(freq)) - 3, lower.tail = F))
              )

# Cas 3 : Loi binomiale
par <- numeric(2)
par[1] <- max_freq
par[2] <- mean(freq) / par[1]

datas <- data.frame(c(freq, qbinom((0:100)/100, par[1], par[2])),
                    Source <- c(rep("Empirique", length(freq)),rep("Théorique", 101)))

ggplot(datas,  aes(datas[, 1], fill = Source)) + 
    geom_histogram(alpha = 0.3, aes(y = ..density..), position = 'identity', binwidth = 1)+
    xlab("n") + ylab("Probabilité") +
    theme(legend.title = element_blank())

Q <- 0
for (n in unique(freq)){
    Q <- Q + (sum(freq == n) - n_obs * dbinom(n, par[1], par[2])) ^ 2 /
        sum(freq == n)
}


chi2 <- rbind(chi2,
              "Binomiale" = c(
                  Q,
                  qchisq(0.95, length(unique(freq)) - 3, lower.tail = T),
                  pchisq(Q, length(unique(freq)) - 3, lower.tail = F)
              ))
colnames(chi2) <- c("Statistique", "valeur critique" , "P-value")
chi2

# Malgré un test du chi-2 qui avantagerais la loi binomialle, le modèle de Poisson 
# est un choix plus judicieux puisque la loi binomiale possède un domaine finit
# qui ne s'apprète pas bien à notre contexte de nombre de sinistres automobiles.

#-----------------------------------------------------------------------------------
#=============== Trouver la loi marginale de la v.a. de sévérité (X). ==============
summary(loss)
loss <- sort(loss)
n_obs <- length(loss)

plot(ecdf(loss))
plot(ecdf(loss[1:(0.99 * n_obs)]))
plot(ecdf(loss[1:(0.9 * n_obs)]))
plot(ecdf(loss[1:(0.8 * n_obs)]))
plot(ecdf(loss[1:(0.3 * n_obs)]))

hist(loss[1:(0.999 * n_obs)], density = F, probability = T)
hist(loss[1:(0.99 * n_obs)], density = F, probability = T)
hist(loss[1:(0.9 * n_obs)], density = F, probability = T)
hist(loss[1:(0.8 * n_obs)], density = F, probability = T)
hist(loss[1:(0.3 * n_obs)], density = F, probability = T)

MeanExcess(loss)
MeanExcess(loss[1:(0.999*n_obs)])
MeanExcess(loss[1:(0.99*n_obs)])
MeanExcess(loss[1:(0.9*n_obs)])
MeanExcess(loss[1:(0.8*n_obs)])
MeanExcess(loss[1:(0.3*n_obs)])

pt_rupture <- 0.97*n_obs
loss[pt_rupture]
summary(loss[1:pt_rupture])

# Comparaison de la distribution des 97% premières données avec celle d'une v.a. suivant
# une loi Gamma. ----
par <- numeric(2)
par[1] <- mean(loss[1:pt_rupture]) / var(loss[1:pt_rupture])
par[2] <- par[1] * mean(loss[1:pt_rupture])
pgamma(2, par[2], par[1], lower.tail = F)

fct_vrais <- function(x, par) {
    sapply(x, function(x_i)
        dgamma(x_i, par[2], par[1]) / pgamma(2, par[2], par[1], lower.tail = F))
    }

mle_gamma <- constrOptim(par, function(par_) - sum(log(fct_vrais(loss[1:pt_rupture], par_))),
                         grad=NULL, ui=diag(2), ci=c(0,0))
mle_gamma$par
S_X.2 <- pgamma(2, mle_gamma$par[2], mle_gamma$par[1], lower.tail = F)
# MeanExcess(sort(rgamma(1e+5, mle_gamma$par[2], mle_gamma$par[1]))[1:(0.9*1e+5)])

p <- sapply(loss[1:pt_rupture],function(x){pgamma(x, mle_gamma$par[2], mle_gamma$par[1]) / S_X.2})
plot(loss[1:pt_rupture], qgamma(p, mle_gamma$par[2], mle_gamma$par[1]),
     xlab="Quantiles observes", ylab="Quantiles theoriques (gamma)", type="l")
lines(loss[1:pt_rupture], loss[1:pt_rupture], col='red')
axis(1, tck = 1, lty = 2, col = "grey")
axis(2, tck = 1, lty = 2, col = "grey")


# une loi Exponentielle ----
par <- 1 / mean(loss[1:pt_rupture])
pexp(2, par, lower.tail = F)

fct_vrais <- function(x, par) {
    sapply(x, function(x_i)
        dexp(x_i, par) / pexp(2, par, lower.tail = F))
}

mle_exp <- optimize(function(par_) - sum(log(fct_vrais(loss[1:pt_rupture], par_))),
                    interval = c(1e-5, 1e-3))
mle_exp$par <- mle_exp$minimum

S_X.2 <- pexp(2, mle_exp$par, lower.tail = F)
# MeanExcess(sort(rexp(1e+5, mle_exp$par))[1:(0.9*1e+5)])

p <- sapply(loss[1:pt_rupture],function(x){pexp(x, mle_exp$par) / S_X.2})
plot(loss[1:pt_rupture], qexp(p, mle_exp$par),
     xlab="Quantiles observes", ylab="Quantiles theoriques (Exponentielle)")
lines(loss[1:pt_rupture], loss[1:pt_rupture], col='red')
axis(1, tck = 1, lty = 2, col = "grey")
axis(2, tck = 1, lty = 2, col = "grey")

# H0: La loi expontielle est une simplification adéquate de la loi gamma.
# H1: La loi gamma représente mieux les données.
(TLR = 2 * (mle_exp$objective - mle_gamma$value))
qchisq(0.95, 1) ; pchisq(TLR, 1, lower.tail = F)
# Comme l'écart entre les vraisemblances est trop fort, on rejette 


# une loi lognormale ----
pt_rupture <- 0.97*n_obs

par <- numeric(2)
par[1] <- mean(log(loss[1:pt_rupture]))
par[2] <- var(log(loss[1:pt_rupture]))

plnorm(2, par[1], par[2], lower.tail = F)

fct_vrais <- function(x, par) {
    sapply(x, function(x_i)
        dlnorm(x_i, par[1], par[2]) / plnorm(2, par[1], par[2], lower.tail = F))
}

mle_lnorm <- constrOptim(par, function(par_) - sum(log(fct_vrais(loss[1:pt_rupture], par_))),
                         grad=NULL, ui=diag(2), ci=c(0,0))
mle_lnorm$par
S_X.2 <- plnorm(2, mle_lnorm$par[1], mle_lnorm$par[2], lower.tail = F)
# MeanExcess(sort(rlnorm(1e+5, mle_lnorm$par[1], mle_lnorm$par[2]))[1:(0.85*1e+5)])

p <- sapply(loss[1:pt_rupture],function(x){plnorm(x, mle_lnorm$par[1], mle_lnorm$par[2]) / S_X.2})
plot(loss[1:pt_rupture], qlnorm(p, mle_lnorm$par[2], mle_lnorm$par[1]),
     xlab="Quantiles observes", ylab="Quantiles theoriques (log-normale)")
lines(loss[1:pt_rupture], loss[1:pt_rupture], col='red')
axis(1, tck = 1, lty = 2, col = "grey")
axis(2, tck = 1, lty = 2, col = "grey")


# une loi Pareto (Lomax) ----
unloadNamespace("ReIns") # La pareto de Reins n'a pas la même forme que celle de actuar.
library(actuar)

pt_rupture <- 0.97*n_obs

par <- numeric(2)
par[1] <- mean(loss[1:pt_rupture])
par[2] <- (par[1] - 1) * mean(loss[1:pt_rupture])
ppareto(2, par[1], par[2], lower.tail = F)

fct_vrais <- function(x, par) {
    sapply(x, function(x_i)
        dpareto(x_i, par[1], par[2]) / ppareto(2, par[1], par[2], lower.tail = F))
}

mle_pareto_gen <- constrOptim(par, function(par_) - sum(log(fct_vrais(loss[1:pt_rupture], par_))),
                         grad=NULL, ui=diag(2), ci=c(0,0))
mle_pareto_gen$par
S_X.2 <- ppareto(2, mle_pareto_gen$par[1], mle_pareto_gen$par[2], lower.tail = F)
# MeanExcess(sort(rpareto(1e+5, mle_pareto_gen$par[1], mle_pareto_gen$par[2]))[1:(0.9*1e+5)])

p <- sapply(loss[1:pt_rupture], function(x)
    ppareto(x, mle_pareto_gen$par[1], mle_pareto_gen$par[2]) / S_X.2)
plot(loss[1:pt_rupture], qpareto(p, mle_pareto_gen$par[1], mle_pareto_gen$par[2]),
     xlab="Quantiles observes", ylab="Quantiles theoriques (Pareto - Lomax)", type="l")
lines(loss[1:pt_rupture], loss[1:pt_rupture], col='red')
axis(1, tck = 1, lty = 2, col = "grey")
axis(2, tck = 1, lty = 2, col = "grey")


# une loi Pareto simple ----
library(ReIns)

pt_rupture <- 0.97 * n_obs

par <- pt_rupture / sum(log(loss[1:pt_rupture]))
ppareto1(2, par, min=2, lower.tail = F)

fct_vrais <- function(x, par) {
    sapply(x, function(x_i)
        ppareto1(x_i, par, min=2) / ppareto1(2, par, min=2, lower.tail = F))
}

mle_pareto_1 <- optimize(function(par_) - sum(log(fct_vrais(loss[1:pt_rupture], par_))),
                    interval = c(1e-5, 1e-3))
mle_pareto_1$par <- mle_pareto_1$minimum

S_X.2 <- ppareto1(2, mle_pareto_1$par, 2, lower.tail = F)
# MeanExcess(sort(rexp(1e+5, mle_exp$par))[1:(0.9*1e+5)])

p <- sapply(loss[1:pt_rupture],function(x){ppareto1(x, mle_pareto_1$par, 2) / S_X.2})
plot(loss[1:pt_rupture], qpareto1(p, mle_pareto_1$par, 2),
     xlab="Quantiles observes", ylab="Quantiles theoriques (Pareto type 1)",
     type="l")
lines(loss[1:pt_rupture], loss[1:pt_rupture], col='red')
axis(1, tck = 1, lty = 2, col = "grey")
axis(2, tck = 1, lty = 2, col = "grey")

#-------------------------------------- Splicing -----------------------------------
# Gamma - Pareto simple ----
pt_rupture <- 0.97 * n_obs
w <- pt_rupture / n_obs

par <- length(loss[pt_rupture:n_obs]) / sum(log(loss[pt_rupture:n_obs]))


fct_vrais <- function(xx, par) {
    k <- 0
    d <- numeric(length(xx))
    for (x in xx) {
        if (x <= pt_rupture)
            d[k <- k + 1] <-
                w * (dgamma(x, mle_gamma$par[2], mle_gamma$par[1]) /
                     (
                         pgamma(pt_rupture, mle_gamma$par[2], mle_gamma$par[1]) -
                             pgamma(2, mle_gamma$par[2], mle_gamma$par[1])
                     ))
        else
            d[k <- k + 1] <-
                (1- w) * (dpareto1(x, par, pt_rupture) / 
                              ppareto1(pt_rupture, par, pt_rupture, lower.tail = F))
    }
    return(d)
}


mle_gam.pareto1 <- optimize(function(par_) - sum(log(fct_vrais(loss, par_))),
                                interval = c(0.001, 2))
(mle_gam.pareto1$par <- mle_gam.pareto1$minimum)


F_splice <- function(xx, par) {
    k <- 0
    P <- numeric(length(xx))
    for (x in xx) {
        if (x <= pt_rupture)
            P[k <- k + 1] <-
                w * (pgamma(x, mle_gamma$par[2], mle_gamma$par[1]) /
                         (
                             pgamma(pt_rupture, mle_gamma$par[2], mle_gamma$par[1]) -
                                 pgamma(2, mle_gamma$par[2], mle_gamma$par[1])
                         ))
        else
            P[k <- k + 1] <- w + (1 - w) *
                (
                    ppareto1(x, par, pt_rupture) /
                        ppareto1(pt_rupture, par, pt_rupture, lower.tail = F)
                )
    }
    return(P)
}

Q_splice <- function(uu, par) {
    k <- 0
    Q <- numeric(length(uu))
    for (u in uu){
        if (u <= w)
            Q[k <- k + 1] <- qgamma(u * pgamma(pt_rupture, mle_gamma$par[2], mle_gamma$par[1]) / w,
                                    mle_gamma$par[2], mle_gamma$par[1])
        else 
            Q[k <- k + 1] <- qpareto1((u - w) / (1 - w),
                                      par, pt_rupture)
    }
    return(Q)
}


p <- sapply(loss,function(x)F_splice(x, mle_gam.pareto1$par))
plot(loss[1:10], Q_splice(p[1:10], mle_gam.pareto1$par),
     xlab="Quantiles observes", ylab="Quantiles theoriques",
     type="l")
lines(loss, loss, col='red')
axis(1, tck = 1, lty = 2, col = "grey")
axis(2, tck = 1, lty = 2, col = "grey")




# =============== Mise en commun
# max_freq <- max(freq)
# ID <- sort(freMTPLfreq[freMTPLfreq$ClaimNb > 0, 1])
# DATA <- matrix(ncol = max_freq + 2, nrow = length(ID))
# for (i in 1:length(ID)){
#     DATA[i,] <- c(ID[i],
#                   N <- freMTPLfreq[freMTPLfreq$PolicyID == ID[i], 2],
#                   freMTPLsev[freMTPLsev$PolicyID == ID[i], 2],
#                   rep(NA, max_freq - N))
# }
# summary(DATA[,-1])

