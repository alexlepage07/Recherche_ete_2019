# Estimation des paramètres

library("erhcv")
#


lap_geo <- function(t, q){
  ans <- q * exp(-t) / (1 - (1 - q) * exp(-t))
  return(ans)
}

lap_geo_deriv <- function(t, q){
  (-q * exp(t)) / ((q + exp(t) - 1)^2)
}

lap_geo_deriv_deriv <- function(t, q){
  (q * exp(t) * (-q + exp(t) + 1)) / ((q + exp(t) - 1)^3)
}

lap_geo_inv <- function(y, q){
  ans <- -log(y / (q + y * (1 - q)))
  return(ans)
}

lap_geo_inv_deriv <- function(y, q){
  q / (y * (q * (y - 1) - y))
}

P_geo <- function(t, q){
  q * t / (1 - (1 - q) * t)
}

P_deriv <- function(t, q){
  q / ((q - 1) * t + 1)^2
}

P_deriv_deriv <- function(t, q){
  -2 * (q - 1) * q / (((q - 1) * t + 1)^3)
}

lap_gam <- function(t, al, be){
  ans <- (be / (be + t))^al
  return(ans)
}

density_cop <- function(u, pars){
  x <- sapply(u, function(i) lap_geo_inv(i, pars))
  lap_geo_deriv_deriv(sum(x), pars) * lap_geo_inv_deriv(u[1], pars) * lap_geo_inv_deriv(u[2], pars)
}


ll=function(u,x){
  s <- 0
  u1 <- u[1:2]
  u2 <- u[3:6]
  W <- expand.grid(u1,u2)#,u3)
  l <- dim(W)[1]
  for(i in 1:l) s <- s + log(density_cop(unlist(W[i,]),x))
  return(s)
}


structure <- GEO(0.5, NULL, list(GEO(0.1, c(1,2), NULL),
                                 GEO(0.1, NULL, list(GAMMA(0.1, c(3,4), NULL),
                                                     GAMMA(0.3, c(5,6), NULL))))) # Build structure
set.seed(2018)
UU <- rCompCop(10000, structure) # Random samples of u's
vect_spear <- cor(UU, method = "spearman")
faa <- dist(vect_spear, method = "maximum")
tree2plot(hclust2tree(hclust(faa)))
m <- 10000
logvrais <- function(x,UU) -sum(sapply(1:m,function(j) ll(UU[j,],x)))
gam_M <- optimize(logvrais, c(0.3, 0.7), UU = UU)$minimum


LapM01 <- stringr::str_replace_all(structure@PGF, "gamma", as.character(gam_M))
LapM01 <- stringr::str_replace_all(LapM01, "z", structure@structure[[1]]@Laplace) #L_N

LapInvM01 <- structure@structure[[1]]@LaplaceInv
LapInvM01 <- stringr::str_replace_all(LapInvM01, "z", stringr::str_replace_all(structure@PGFInv, "gamma", as.character(gam_M)))

# C = L_M01 ( L_M01^(-1)(u1) + L_M01^(-1)(u2))
ini <- "z1 + z2"
ini <- stringr::str_replace_all(ini, "z1", stringr::str_replace_all(LapInvM01, "z", "u1"))
ini <- stringr::str_replace_all(ini, "z2", stringr::str_replace_all(LapInvM01, "z", "u2"))

LapM01 <- stringr::str_replace_all(LapM01, "z", ini)

dM01 <- Deriv::Deriv(LapM01, "u1", cache.exp = F)
dM01 <- Deriv::Deriv(dM01, "u2", cache.exp = F)
dM01 <- parse(text = dM01)

logv_M01 <- function(alpha, data){
  u1 <- data[,1]
  u2 <- data[,2]
  -sum(log(eval(dM01)))
}

gam_M01 <- optimize(logv_M01, c(0, 1), data = UU)$minimum



LapM02 <- stringr::str_replace_all(structure@PGF, "gamma", as.character(gam_M))
LapM02 <- stringr::str_replace_all(LapM02, "z", structure@structure[[2]]@Laplace)

LapInvM02 <- structure@structure[[2]]@LaplaceInv
LapInvM02 <- stringr::str_replace_all(LapInvM02, "z", stringr::str_replace_all(structure@PGFInv, "gamma", as.character(gam_M)))

# C = L_M01 ( L_M01^(-1)(u1) + L_M01^(-1)(u2))
ini <- "z1 + z2"
ini <- stringr::str_replace_all(ini, "z1", stringr::str_replace_all(LapInvM02, "z", "u1"))
ini <- stringr::str_replace_all(ini, "z2", stringr::str_replace_all(LapInvM02, "z", "u2"))

LapM02 <- stringr::str_replace_all(LapM02, "z", ini)

dM02 <- Deriv::Deriv(LapM02, "u1", cache.exp = F)
dM02 <- Deriv::Deriv(dM02, "u2", cache.exp = F)
dM02 <- parse(text = dM02)

MAT <- rbind(c(3, 5),
             c(3, 6),
             c(4, 5),
             c(4, 6))

logv_M02 <- function(gamma, data){
  res <- rep(1, length(data[,1]))
  #ini <- paste("data[[i]][,", 1:length(n), "]", sep = "", collapse = ", ")
  for (i in 1:dim(MAT)[1]){
    for (b in 1:dim(MAT)[2])
      eval(parse(text = paste("u", b, " <- data[,", MAT[i,b], "]", sep = "")))
    res <- res * eval(dM02)
  }
  -sum(log(res))
}

gam_M02 <- optimize(logv_M02, c(0, 1), data = UU)$minimum



# THETA (0, 2, 1)

gamma0 <- gam_M
gamma1 <- gam_M02

PGF <- "(gamma)*(z) / (1 - (1-(gamma))*(z))"
PGFInv <- "1 / (((gamma)/(z)) + (1 - (gamma)))"
Lap <- "(gamma)*exp(-(z)) / (1 - (1 - (gamma)) * exp(-(z)))"
LapInv <- "-log(1 / (((gamma)/(z)) + (1 - (gamma))))"

L_B_gamma <- "(1 / (1 + (z)))^(alpha)"
L_B_gamma_inv <- "((z)^(-1/(alpha)) - 1)"

LapG1 <- stringr::str_replace_all(PGF, "gamma", as.character(gam_M)) #A
LapG1 <- stringr::str_replace_all(LapG1, "z", stringr::str_replace_all(PGF, "gamma", as.character(gam_M02)))

PM_inv02 <- stringr::str_replace_all(PGFInv, "gamma", as.character(gam_M02)) #B
PM_inv02 <- stringr::str_replace_all(PM_inv02, "z", stringr::str_replace_all(PGFInv, "gamma", as.character(gam_M)))

PM_inv_u1 <- stringr::str_replace_all(PM_inv02, "z", "u1") #C
PM_inv_u2 <- stringr::str_replace_all(PM_inv02, "z", "u2")

LB_inv_1 <- stringr::str_replace_all(L_B_gamma_inv, "z", PM_inv_u1)
LB_inv_2 <- stringr::str_replace_all(L_B_gamma_inv, "z", PM_inv_u2)

ini <- "z1 + z2"
ini <- stringr::str_replace_all(ini, "z1", LB_inv_1)
ini <- stringr::str_replace_all(ini, "z2", LB_inv_2)

L_BB <- stringr::str_replace_all(L_B_gamma, "z", ini)

LapG1 <- stringr::str_replace_all(LapG1, "z", L_BB)

dG01 <- Deriv::Deriv(LapG1, "u1", cache.exp = F)
dG01 <- Deriv::Deriv(dG01, "u2", cache.exp = F)
dG01 <- parse(text = dG01)

logv_G01 <- function(alpha, data){
  u1 <- data[,3]
  u2 <- data[,4]
  -sum(log(eval(dG01)))
}

(alpha021_test <- optimize(logv_G01, c(0, 1), data = UU)$minimum)


# Theta (0, 2, 2)

LapG2 <- stringr::str_replace_all(PGF, "gamma", as.character(gam_M)) #A
LapG2 <- stringr::str_replace_all(LapG1, "z", stringr::str_replace_all(PGF, "gamma", as.character(gam_M02)))

PM_inv02 <- stringr::str_replace_all(PGFInv, "gamma", as.character(gam_M02)) #B
PM_inv02 <- stringr::str_replace_all(PM_inv02, "z", stringr::str_replace_all(PGFInv, "gamma", as.character(gam_M)))

PM_inv_u1 <- stringr::str_replace_all(PM_inv02, "z", "u1") #C
PM_inv_u2 <- stringr::str_replace_all(PM_inv02, "z", "u2")

LB_inv_1 <- stringr::str_replace_all(L_B_gamma_inv, "z", PM_inv_u1)
LB_inv_2 <- stringr::str_replace_all(L_B_gamma_inv, "z", PM_inv_u2)

ini <- "z1 + z2"
ini <- stringr::str_replace_all(ini, "z1", LB_inv_1)
ini <- stringr::str_replace_all(ini, "z2", LB_inv_2)

L_BB <- stringr::str_replace_all(L_B_gamma, "z", ini)

LapG2 <- stringr::str_replace_all(LapG1, "z", L_BB)

dG02 <- Deriv::Deriv(LapG2, "u1", cache.exp = F)
dG02 <- Deriv::Deriv(dG02, "u2", cache.exp = F)
dG02 <- parse(text = dG02)

logv_G02 <- function(alpha, data){
  u1 <- data[,5]
  u2 <- data[,6]
  -sum(log(eval(dG02)))
}

(alpha022_test <- optimize(logv_G02, c(0, 1), data = UU)$minimum)

print(c(gam_M, gam_M01, gam_M02, alpha021_test, alpha022_test))

