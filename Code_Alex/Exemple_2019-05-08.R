# Exemple pour Diamilatou et Alexandre
# 
# Date : 2018-05-08
#
# ========================= Calculs demandés par Étiennes ====================================
#-------------------------- Modèle de dépendance ---------------------------------------------

Modele_dependance <- function(table){
    col_prob <- ncol(table)
    N <- unique(table[, 1])
    
    
    p_N <- function(n){
        return(sapply(n, function(k) sum(table[which(table[, 1] == k), col_prob])))
    }
    
    
    F_N <- function(n){
        return(sum(p_N(N[N <= n])))
    }
    
    
    X <- list()
    for (col in 1:(ncol(table)-2)){
        X[[col]] <- unique(table[, col+1])
    }
    
    
    p_x <- function(X_i, x) {
        col <- X_i + 1
        return(sapply(x, function(k) sum(table[which(table[, col] == k), col_prob])))
    }
    
    
    F_x <- function(X_i, x){
        col <- X_i
        return(sum(p_x(X_i, X[[col]][which(X[[col]] <= x)])))
    }
    
    
    S <- c()
    for(i in 1:nrow(table)){
        if (table[i, 1] == 0)
            S <- append(S, 0)
        else
            S <- append(S, sum(sapply(1:table[i, 1], function(j) table[i, j+1] )))
    }
    table_S <- cbind(table, S)
    S <- sort(unique(S))
    
    
    p_S <- function(s) {
        return(sapply(s, function(k) sum(table[which(table_S[, ncol(table_S)] == k), col_prob])))
    }
    
    
    F_S <- function(s){
        return(sapply(s, function(k) sum(p_S(S[S <= k]))))
    }
    
    
    stop_loss <- function(x){
        if (x == max(S)) return(0)
        return(sum(S[S > x] * p_S(S[S > x])) - x * (1 - F_S(x)))
    }
    
    
    VaR_S <- function(k) return(min(S[F_S(S) >= k]))
    
    
    TVaR_S <- function(k) return(stop_loss(VaR_S(k)) / (1 - k) + VaR_S(k))
    
    
    kappa <- F_S(S[-length(S)])
    
    rep <- list( 
        "S" = S,
        "p_S" = p_S(S),
        "F_S" = F_S(S),
        "E_S" = E_S <- sum(S * p_S(S)),
        "V_S" = sum((S - E_S)^2 * p_S(S)),
        "stop_loss" = sapply(S, stop_loss),
        "TVaR_S" = sapply(kappa, TVaR_S)
    )
    return(rep)
}

# ========================= Test ====================================
{
nlignes <- 27
ncolonnes <- 4
pNX1X2 <- matrix(0, nlignes, ncolonnes)
# pNX1X2
p1 <- 0.8
p2 <- 1 - p1
a1 <- 0.2
a2 <- 0.5
b1 <- 0.4
b2 <- 0.9

j <- 0
for (n in 0:2){
    for (k1 in 1:3){
        for (k2 in 1:3){
            j <- j+1
            pNX1X2[j,] <- c(n, k1, k2, 
                          p1 * dbinom(n,2,a1) * dbinom(k1-1,2,b1) * dbinom(k2-1,2,b1) +
                              p2 * dbinom(n,2,a2) * dbinom(k1-1,2,b2) * dbinom(k2-1,2,b2)
                          )
        }
    }
}

# pNX1X2
# sum(pNX1X2[,4])

} # Données fournies initialement par Étienne
Modele_dependance(pNX1X2)

#-------------------------- Modèle d'indépendance ---------------------------------------------

# Briser la structure de dépendance dans les probabilités conjointes
{table <- pNX1X2
col_prob <- ncol(table)
N <- unique(table[, 1])


p_N <- function(n){
    return(sapply(n, function(k) sum(table[which(table[, 1] == k), col_prob])))
}


F_N <- function(n){
    return(sum(p_N(N[N <= n])))
}


X <- list()
for (col in 1:(ncol(table)-2)){
    X[[col]] <- unique(table[, col+1])
}


p_x <- function(X_i, x) {
    col <- X_i + 1
    return(sapply(x, function(k) sum(table[which(table[, col] == k), col_prob])))
}


F_x <- function(X_i, x){
    col <- X_i
    return(sum(p_x(X_i, X[[col]][which(X[[col]] <= x)])))
}


S <- c()
for(i in 1:nrow(table)){
    if (table[i, 1] == 0)
        S <- append(S, 0)
    else
        S <- append(S, sum(sapply(1:table[i, 1], function(j) table[i, j+1] )))
}
table_S <- cbind(table, S)
}# Partir du modèle créé initialement
{
table_indep <- table

for (i in 1:nrow(table_indep)){
    prob_ind <- p_N(table_indep[i, 1]) * prod(sapply(1:(length(X)),
                                                     function(j) p_x(j, table_indep[i, j+1])))
    table_indep[i, col_prob] <- prob_ind
}
sum(table_indep[,col_prob])
}# Briser la structure de dépendance

# Calcul des mesures avec le modèle indépendant.
Modele_dependance(table_indep)
