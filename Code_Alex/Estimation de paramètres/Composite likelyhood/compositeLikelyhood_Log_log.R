# Estimation de paramètres du modèle collectif du risque avec structure de dépendance
# sur des données simulées.
#
# Le modèle est construit avec une copule archimédienne hiérachique dont la v.a. mère
# suit une loi logarithmique de même que la v.a. enfant.
library(copula)
library(nCopula)
library(Deriv)
library(stringr)
library(xtable)
library(ggplot2)
library(psych)

source("Mesure_dependance va mixtes.R")

set.seed(20190702)
# ================================== Simulations des données d'entraînement ================
# Paramètres de simulation
nb_xi <- n <- 5
q <- 2/5
beta <- 1/100
alpha0 <- 5
alpha1 <- 7
nsim <- 1e+4

DATA_train <- rCompCop(nsim, LOG(1-exp(-alpha0), 1, list(LOG(1-exp(-alpha1), 1:5, NULL))))

for (i in 1:nsim){
    DATA_train[i, 1] <- N <- qbinom(DATA_train[i,1], n, q)
    if (N==0){ 
        DATA_train[i,-1] <- rep(NaN, nb_xi)
    }else{
        DATA_train[i,-1] <- c(qexp(DATA_train[i, 2:(N + 1)], beta), rep(NaN, nb_xi - N))}
}

colnames(DATA_train) <- c("N", sapply(1:(ncol(DATA_train)-1),function(i) paste0("X", i)))
head(DATA_train)

# Sommaire des résultats de simulation
(sommaire_train <- summary(DATA_train))
xtable(sommaire_train)  # Permet de  convertir un tableau de R à LaTeX.


#Graphiques de goodness of fit pour les données simulées.
datas <- data.frame(c(DATA_train[,1], qbinom((0:100)/100, n, q)),
                    Source <- c(rep("Empirique", length(DATA_train[,1])),rep("Théorique", 101)))

ggplot(datas,  aes(datas[, 1], fill = Source)) + 
    geom_histogram(alpha = 0.3, aes(y = ..density..), position = 'identity', binwidth = 1)+
    xlab("n") + ylab("Probabilité") +
    theme(legend.title = element_blank())


lst_XX <- array(DATA_train[,-1])
lst_XX <- lst_XX[!is.na(lst_XX)]
datas <- data.frame(c(lst_XX, qexp((0:100)/100, beta)),
                    Source <- c(rep("Simul", length(lst_XX)),
                                rep("théorique", 101)))

ggplot() + 
    geom_histogram(alpha = 0.3, aes(x= lst_XX, y = ..density.., fill = "Empirique"), position = 'identity')+
    geom_density(alpha = 0.3, aes(x= qexp((0:100)/100, beta), y = ..density.., fill = "Théorique")) + 
    xlab("x") + ylab("Densité") +
    theme(legend.title = element_blank())

# ====================== Estimation des paramètres des v.a. N et M ============================
para <- list("N"=list("size"=nb_xi, "prob"=q), "X"=list("beta"=beta), "alpha0"=alpha0, "alpha1"=alpha1)

Copule0 <- frankCopula
Copule1 <- function(alpha0, alpha1) {
    LOG(1-exp(-alpha0), NULL, list(LOG(1-exp(-alpha1), 1:2, NULL)))
}
F_N <- pbinom
F_X <- pexp
n_max <- nb_xi


F_N. <- function(n, prob) {
    # Paramétriser la fonction de répartition de N.
    do.call(F_N, list(q = n, size=nb_xi, prob))
}
F_X. <- function(x, beta) {
    # Paramétriser la fonction de répartition de X.
    do.call(F_X, list(q = x, beta))
}

pCopule0_NX <- function(n, x, para) {
    # Fonction de répartition conjointe de N et X.
    pCopula(c(F_N.(n, para[1]), F_X.(x, para[2])), Copule0(para[3]))
}
dCopule0_NX <- function(n, x, para) {
    # Fonction de densité conjointe de N et X.
    eps <- .Machine$double.eps^0.25
    
    f <- function(n,x,para){
        (pCopule0_NX(0, x + eps, para) - pCopule0_NX(0, x, para)) / eps}
    
    if (n==0)
        return(F_N.(0, para[1]))
    return(f(n, x, para) - f(n - 1, x, para))
}

fct_Score_NX <- function(para, Data = DATA_train){
    # Fonction de score: On cherchera à minimiser la log-vraisemblance négative
    neg_log_vrais <- 0
    for (i in 1:nrow(Data)) {
        N <- Data[i, 1]
        for (j in 1:N) {
            X <- Data[i, j + 1]
            neg_log_vrais <- neg_log_vrais - log(dCopule0_NX(N, X, para))
        }
    }
    return(neg_log_vrais)
}


{
    tau_0 <- mean(tau_NX(DATA_train, nsim=1, silent=F))


    Tau_graph <- sapply(domaine <- c(-5:-1, 1:10), function(a)
        tau_kendall_theorique(pbinom,
                              list(N=list(size=nb_xi, prob=q)),
                              Copule0, a, nb_xi))
    
    plot(domaine, Tau_graph,
         ylab="tau de kendall",
         xlab="alpha",
         type="l"
    )
    axis(2, tck = 1, lty = 2, col = "grey") # L'axe des ordonnées
    axis(1, at=domaine, tck=1, lty = 2, col = "grey",) # L'axe des abscisses
    abline(a=tau_0, b=0, col="green")


    bornes <- c(5.5, 6.5)
    
    temps_0 <- system.time(
        alpha0_n <- inversion_tau_kendall(pbinom,
                                          list(N=list(size=nb_xi, prob=q)),
                                          Copule0,
                                          bornes,
                                          tau_0,
                                          nb_xi)
    )
    
    (tbl_0 <- rbind(
        "Vrai paramètre" = alpha0,
        "Paramètre trouvé" = alpha0_n,
        "Temps d'optimisation" = temps_0[[3]]
    ))
} # Estimation du paramètre de dépendance entre N et X_i par la méthode des moments.


# Valeurs de départ de l'optimisation
val_depart <- c("q" = mean(DATA_train[,1]) / nb_xi,
                "beta" = 1/mean(DATA_train[,-1], na.rm = T),
                "alpha0" = alpha0_n)
# Bornes de l'optimisation
bounds <- matrix(c(0, 1,
                   0, Inf,
                   0, Inf), ncol = 2, byrow = T)
# Convert the constraints to the ui and ci matrices
n <- nrow(bounds)
ui <- rbind( diag(n), -diag(n) )
ci <- c( bounds[,1], - bounds[,2] )
# Remove the infinite values
i <- as.vector(is.finite(bounds))
ui <- ui[i,]
ci <- ci[i]

temps_solv <- system.time(
    # Estimation des paramètres
    mle_M <- constrOptim(val_depart, fct_Score, grad=NULL, 
                         ui=ui, ci=ci)
)

# Tableaux permettant de présenter les résultats
(resultats <- rbind("Valeurs de départ"=round(val_depart,4),
                    "Estimateurs" = round(mle$par, 4),
                    "Vrais paramètres" = round(para,4)))
(temps_tot <- rbind("Temps de dérivation"=temps_deriv,
                    "Temps d'estimation"=temps_solv[[3]]))
# Conversion en LaTeX
xtable(resultats, digits = 4)
xtable(temps_tot)
