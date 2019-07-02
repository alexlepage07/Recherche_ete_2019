# Dans ce projet R, on cherche à appliquer le modèle collectif 
# du risque avec structure de dépendance. Cela implique de
# trouver la copule qui représente le mieux cette structure,
# À identifier les lois marginales et à estimer les paramètres
# du modèle.

library(CASdatasets)
library(ggplot2)
library(ReIns)
library(actuar)

source("Mesure_dependance va mixtes.R")


data("freMTPLfreq")
data("freMTPLsev")
freq  <- freMTPLfreq$ClaimNb
loss <- freMTPLsev$ClaimAmount

n_obs <- length(freq)
max_freq <- max(freq)


#=============== Trouver la loi marginale de la v.a. de fréquence (N). =============
# Analyse des données --------------------------------------------------------------
table(freq)
summary(freq)
plot(ecdf(freq), xlim=c(0,5), main="")
barplot(table(freq))

# ----------------------------------------------------------------------------------
# Cas 1 : Loi de poisson -----------------------------------------------------------
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

# Cas 2 : Loi de poisson inflationnée à zéro -----------------------------------------------
par <- numeric(2)
par[1] <- mean(mean(freq), var(freq))
par[2] <- (1 - mean(freq==0)) / (1 - exp(-par[1]))


dpois.zero.inf <- function(x, lambda, prob){
    res <- numeric(length(x))
    for (i in 1:length(x)){
        if (x[i] == 0)
            res[i] <- (1 - prob) + prob * exp(-lambda)
        else
            res[i] <- prob * dpois(x[i], lambda)
        }
    return(res)
}

ppois.zero.inf <- function(q, lambda, prob) {
    sapply(q, function(n)
        sum(dpois.zero.inf(0:n, lambda, prob)))
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


mle_pois.zero.inf <-
    optim(par, function(par)
        - sum(log(dpois.zero.inf(freq, par[1], par[2]))))

tbl_parametres <- rbind("moments" = par, "mle" = mle_pois.zero.inf$par)
colnames(tbl_parametres) <- c("lambda", "pi")
tbl_parametres

datas <-
    data.frame(c(
        freq,
        qpois.zero.inf((0:99.99) / 100,
                       mle_pois.zero.inf$par[1],
                       mle_pois.zero.inf$par[2]
        )
    ),
    Source <-
        c(rep("Empirique", length(freq)), rep("Théorique", 100)))

ggplot(datas,  aes(datas[, 1], fill = Source)) +
    geom_histogram(
        alpha = 0.3,
        aes(y = ..density..),
        position = 'identity',
        binwidth = 1
    ) +
    xlab("n") + ylab("Probabilité") +
    theme(legend.title = element_blank())


Q <- 0
for (n in unique(freq)) {
    Q <-
        Q + (
            sum(freq == n) - n_obs * dpois.zero.inf(n, mle_pois.zero.inf$par[1],
                                                    mle_pois.zero.inf$par[2])
        ) ^ 2 /
        sum(freq == n)
}
chi2 <- rbind(chi2,
    "Poisson inflationné à zéro" = c(
        Q,
        qchisq(0.95, length(unique(freq)) - 3, lower.tail = T),
        pchisq(Q, length(unique(freq)) - 3, lower.tail = F)
    ))

# Cas 3 : Loi binomiale négative ---------------------------------------------------
par <- numeric(2)
par[1] <- mean(freq) / var(freq)
par[2] <- mean(freq)

mle_nbinom <- constrOptim(par,
                   function(par) -sum(log(dnbinom(freq, par[1], mu=par[2]))),
                   grad = NULL, ui=diag(2), ci = c(0,0))
par <- mle_nbinom$par

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

# Cas 4 : Loi binomiale ------------------------------------------------------------
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


#-----------------------------------------------------------------------------------
#=============== Trouver la loi marginale de la v.a. de sévérité (X). ==============
summary(loss)
loss <- sort(loss)
n_obs <- length(loss)
loss[200:300]
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
# Loi Gamma. ----
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


# Loi Exponentielle ----
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


# Loi lognormale ----
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


# Loi Pareto (Lomax) ----
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


# Loi Pareto simple ----
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

par <- numeric(4)
par[1] <- pt_rupture
par[2:3] <- mle_gamma$par
par[4] <- length(loss[pt_rupture:n_obs]) / sum(log(loss[pt_rupture:n_obs]))

w <- function(theta)theta / n_obs


fct_vrais <- function(xx, par) {
    k <- 0
    d <- numeric(length(xx))
    for (x in xx) {
        if (x <= par[1])
            d[k <- k + 1] <-
                w(par[1]) * (dgamma(x, par[3], par[2]) /
                     (
                         pgamma(par[1], par[3], par[2]) -
                             pgamma(2, par[3], par[2])
                     ))
        else
            d[k <- k + 1] <-
                (1 - w(par[1])) * (dpareto1(x, par[4], par[1]))
    }
    return(d)
}


mle_gam.pareto1 <- constrOptim(par, function(par_) - sum(log(fct_vrais(loss, par_))),
                                grad=NULL, ui=diag(4), ci=rep(0,4))
mle_gam.pareto1$par


F_splice <- function(xx, par) {
    k <- 0
    d <- numeric(length(xx))
    for (x in xx) {
        if (x <= par[1])
            d[k <- k + 1] <-
                w(par[1]) * (pgamma(x, par[3], par[2]) /
                                 (
                                     pgamma(par[1], par[3], par[2]) -
                                         pgamma(2, par[3], par[2])
                                 ))
        else
            d[k <- k + 1] <-
                w(par[1]) + (1 - w(par[1])) * (ppareto1(x, par[4], par[1]))
    }
    return(d)
}

Q_splice <- function(uu, par) {
    k <- 0
    Q <- numeric(length(uu))
    for (u in uu){
        if (u <= w(par[1]))
            Q[k <- k + 1] <- qgamma(u * pgamma(par[1], par[3], par[2]) / w(par[1]),
                                    par[3], par[2])
        else 
            Q[k <- k + 1] <- qpareto1((u - w(par[1])) / (1 - w(par[1])),
                                      par[4], par[1])
    }
    return(Q)
}


p <- sapply(loss,function(x)F_splice(x, mle_gam.pareto1$par))
plot(loss[], Q_splice(p[], mle_gam.pareto1$par),
     xlab="Quantiles observes", ylab="Quantiles theoriques",
     type="l")
lines(loss, loss, col='red')
axis(1, tck = 1, lty = 2, col = "grey")
axis(2, tck = 1, lty = 2, col = "grey")




# =========================== Calcul de la dépendance =====================
max_freq <- max(freq)
ID <- sort(freMTPLfreq[freMTPLfreq$ClaimNb > 0, 1])
DATA <- matrix(ncol = max_freq + 2, nrow = length(ID))
for (i in 1:length(ID)){
    DATA[i,] <- c(ID[i],
                  N <- freMTPLfreq[freMTPLfreq$PolicyID == ID[i], 2],
                  freMTPLsev[freMTPLsev$PolicyID == ID[i], 2],
                  rep(NA, max_freq - N))
}
summary(DATA[,-1])


# Test numéro 1 : Logarithmique-Gamma -----
{
    Copule0 <- frankCopula
    struct_copule <- function(alpha0, alpha1){
        LOG(1-exp(-alpha0), NULL, list(GAMMA(1/alpha1, 1:2, NULL)))
    }
    
    F_N <- ppois.zero.inf
    para <- list(N=list(lambda=mle_pois.zero.inf$par[1],
                        prob=mle_pois.zero.inf$par[2]))
    n_max <- qpois.zero.inf(1-1e-6, para$N$lambda, para$N$prob)
    
}# Les paramètres


set.seed(20190618)
tau_0 <- tau_NX(DATA[,-1], nsim=10, silent=F)

var_tau <- var(tau_0)
(mean_tau <- mean(tau_0))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))


{
    Tau_graph <- sapply(domaine <- c(-5:-1, 1:10), function(a)
        tau_kendall_theorique(F_N, para, Copule0, a, n_max))
    
    plot(domaine, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=domaine, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
}# Analyse graphique

{
    bornes <- c(6, 7)
    
    temps_0 <- system.time(
        alpha0_n <- inversion_tau_kendall(F_N, F_X,
                                          para,
                                          Copule0,
                                          bornes,
                                          mean_tau,
                                          n_max, x_max)
    )
    
    (tbl_0 <- rbind(
        "Vrai paramètre" = alpha0,
        "Paramètre trouvé" = alpha0_n,
        "Temps d'optimisation" = temps_0[[3]]
    ))
} # Alpha0

set.seed(20190618)
tau_1 <- tau_XX(DATA_train, nsim=30, silent=T)

(var_tau <- var(tau_1, na.rm = T))
(mean_tau <- mean(tau_1, na.rm = T))
(IC_tau <- c(mean_tau - sqrt(var_tau) * qnorm(0.975),
             mean_tau + sqrt(var_tau) * qnorm(0.975)))
tau_kendall_theorique_continues(struct_copule,
                                alpha0_n, alpha1)

{
    Tau_graph <- sapply(c(0.001,1:5), function(a)
        tau_kendall_theorique_continues(struct_copule,
                                        alpha0_n, alpha1=a))
    
    plot(0:5, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=0:5, tck=1, lty = 2, col = "grey", labels = F) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
    
    plot(0.001, Tau_graph[1],
         ylab="tau de kendall",
         xlab="alpha",
         type="p"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=(0:2)/10, tck=1, lty = 2, col = "grey", labels = F) # L'axe des abscisses
    abline(a=mean_tau, b=0, col="green")
}# Analyse graphique

bornes <- c(1e-12, 0.001)

temps_1 <- system.time(
    alpha1_n <- inversion_tau_kendall_XX(struct_copule, alpha0_n, mean_tau, bornes)
)

(tbl_1 <- rbind(
    "Vrai paramètre" = alpha1,
    "Paramètre trouvé" = alpha1_n,
    "Temps d'optimisation" = temps_1[[3]]
))

tbl_tot <- cbind(tbl_0, tbl_1)
colnames(tbl_tot) <- c("alpha0", "alpha1")
tbl_tot

