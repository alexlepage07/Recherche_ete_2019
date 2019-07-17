library(nCopula)
library(stringr)
library(Deriv)


PGF_N <- "-(log(1 - (z) * (1 - exp(-a0))) / a0)"
PGFInv_N <- "(1 - exp(-a0)^(z)) / (1 - exp(-a0))"

Laplace_N <- "(1 / (1 + (z)))^(alpha)"
LaplaceInv_N <- "((z)^(-1/alpha) - 1)"

LapT <- stringr::str_replace_all(PGF_N, "z", Laplace_N)


#------------------------------------------------------------
func.TLS_ <- function(z, a0, alpha) eval(parse(text = LapT))
func.TLS_(2, -log(1 - 0.5), 0.7)
#------------------------------------------------------------


for (i in 1:5){
  LapT <- Deriv::Deriv(LapT, "z", cache.exp = F)
  print(i)
}
LapT_inv <- stringr::str_replace_all(LaplaceInv_N, "z", PGFInv_N)


#------------------------------------------------------------
func.TLS_ <- function(z, a0, alpha) eval(parse(text = LapT))
func.TLS_(2, -log(1 - 0.5), 0.7)
#------------------------------------------------------------
func.TLS_ <- function(z, a0, alpha) eval(parse(text = LapT_inv))
func.TLS_(2, -log(1 - 0.5), 0.7)
#------------------------------------------------------------


ini <- "z1 + z2 + z3 + z4 + z5"
ini <- stringr::str_replace_all(ini, "z1", stringr::str_replace_all(LapT_inv, "z", "u1"))
ini <- stringr::str_replace_all(ini, "z2", stringr::str_replace_all(LapT_inv, "z", "u2"))
ini <- stringr::str_replace_all(ini, "z3", stringr::str_replace_all(LapT_inv, "z", "u3"))
ini <- stringr::str_replace_all(ini, "z4", stringr::str_replace_all(LapT_inv, "z", "u4"))
ini <- stringr::str_replace_all(ini, "z5", stringr::str_replace_all(LapT_inv, "z", "u5"))

dG <- stringr::str_replace_all(LapT, "z", ini)


#-----------------------------------------------------------
func.TLS_ <- function(u1, u2, u3, u4, u5, a0, alpha) eval(parse(text=dG))
func.TLS_(0.4,0.5,0.3,0.6,0.5, -log(1 - 0.5), 0.7)
#-----------------------------------------------------------


ini2 <- "z1 * z2 * z3 * z4 * z5"
ini2 <- stringr::str_replace_all(ini2, "z1", Deriv::Deriv(stringr::str_replace_all(LapT_inv, "z", "u1"), "u1", cache.exp = F))
ini2 <- stringr::str_replace_all(ini2, "z2", Deriv::Deriv(stringr::str_replace_all(LapT_inv, "z", "u2"), "u2", cache.exp = F))
ini2 <- stringr::str_replace_all(ini2, "z3", Deriv::Deriv(stringr::str_replace_all(LapT_inv, "z", "u3"), "u3", cache.exp = F))
ini2 <- stringr::str_replace_all(ini2, "z4", Deriv::Deriv(stringr::str_replace_all(LapT_inv, "z", "u4"), "u4", cache.exp = F))
ini2 <- stringr::str_replace_all(ini2, "z5", Deriv::Deriv(stringr::str_replace_all(LapT_inv, "z", "u5"), "u5", cache.exp = F))

dG_final <- "z1 * z2"
dG_final <- stringr::str_replace_all(dG_final, "z1", dG)
dG_final <- stringr::str_replace_all(dG_final, "z2", ini2)
dG_final <- parse(text = dG_final)


#-----------------------------------------------------------
func.TLS_ <- function(u1, u2, u3, u4, u5, a0, alpha) eval(dG_final)
func.TLS_(0.4,0.5,0.3,0.6,0.5, -log(1 - 0.5), 0.7)
#-----------------------------------------------------------


eval(parse(text = paste("a0 <- -(log(1 - 0.5))", sep = "")))
structure <- LOG(0.5, NULL, list(GAMMA(0.7, 1:5, NULL)))
UU <- rCompCop(10000, structure = structure)

logv_M <- function(alpha){
  res <- rep(1, length(UU[,1]))
  eval(parse(text = paste("u1 <- UU[,1]", sep = "")))
  eval(parse(text = paste("u2 <- UU[,2]", sep = "")))
  eval(parse(text = paste("u3 <- UU[,3]", sep = "")))
  eval(parse(text = paste("u4 <- UU[,4]", sep = "")))
  eval(parse(text = paste("u5 <- UU[,5]", sep = "")))
  res <- res * eval(dG_final)
  -sum(log(res))
}


plot(domaine <- (6:9)/10, 
     sapply(domaine, logv_M), type="l")
axis(2, tck = 1, lty = 2, col = "grey")
axis(1, tck=1, lty = 2, col = "grey")

gam_M <- optimize(logv_M, c(0.6, 0.8))$minimum
