library(actuar)
library(copula)
library(nCopula)
library(xtable)
library(ggplot2)
library(Deriv)
library(stringr)
library(pracma)

########################################################## Alpha0
# log == frank et geo == amh
#calcul du tau de kendall pour variables mixtes
## Ici binom et exp a 4 fct

###### Calcul tau de kendall empirique
####### calcule le tau de kendall entre deux va mixte
## prends en entrée une matrice a deux colones avec 1 discrete et l'autre continu
## sort le tau de kendall
# non definit si N vaut 0

ChoixN<- function(obs){
  obs<- obs[!(obs[,1]==0),]
  sortie<- matrix(0,nrow(obs),2)
  sortie[,1]<- obs[,1]
  pos<-sapply(obs[,1],function(x) sample.int(x,1,replace = FALSE))
  for (i in 1:nrow(obs)) {
    sortie[i,2]<- obs[i,pos[i]+1]
  }
  sortie
}


TauKendallMixte<- function(obs1){
  knd<- numeric(100)
  for (k in 1:100){ 
    obs<- ChoixN(obs1)
    n<- nrow(obs)
    nbrconc<- 0 
    nbregal<-0
    for (i in 1:(n-2)) {
      d<-obs[-1,]
      d[,1]<- d[,1]-obs[1,1]
      d[,2]<- d[,2]-obs[1,2]
      nbregal<- nbregal+sum(d[,1]==0)
      nbrconc<-nbrconc+sum(sign(d[,1])==sign(d[,2]))
      obs<- obs[-1,]
    }
    knd[k]<-( 4*nbrconc + 2*nbregal )/(n*(n-1)) -1
    print(k)
    print(knd[k])
  }
  kendall<- mean(knd)
  print(kendall)
  return(kendall)
}
#################################### Fonctions

####### Clayton 
#Fr bivarié d'une copule de clayton,amh et frank
FrClayton <- function(u1,u2,alpha){
  (u1^(-alpha) + u2^(-alpha) - 1)^(-1/alpha)
}

FrFrank <- function(u1,u2,alpha){
  -log( 1+ (exp(-alpha*u1)-1)*(exp(-alpha*u2)-1)/(exp(-alpha)-1))/alpha
}

FrAMH <- function(u1,u2,alpha){
  u1*u2/(1- alpha *(1-u1)*(1-u2))
}

####### Joe alpha = 1/gamma 
######## alpha sur 1,Inf
FrJoe <- function(u1,u2,alpha){
  1-(1-exp(log(1-(1-u1)^(alpha))+log(1-(1-u2)^(alpha))))^(1/alpha)
}

## densité bivarié mixte d'une copule de clayton,frank,amh avec marginales 
###########################################################################exp 1/100
fClayton <- function(u1,x,alpha){
  #expr<- "(u1^(-alpha) + (1-exp(-x/100))^(-alpha) - 1)^(-1/alpha)"
  #expr2<-parse(text = Deriv(expr, "x"))
  expr2<-expression({.e2 <- exp(-(x/100)); .e3 <- 1 - .e2; 0.01 * (.e2/(.e3^(1 + alpha) * (1/.e3^alpha + 1/u1^alpha - 1)^(1 +      1/alpha)))})
  eval(expr2)
}

fFrank <- function(u1,x,alpha){
  #expr<- "-log( 1+ ((exp(-alpha*u1)-1)*(exp(-alpha*(1-exp(-x/100)))-1)/(exp(-alpha)-1)))/alpha"
  #expr2<-parse(text = Deriv(expr, "x"))
  expr2<-expression({.e2 <- exp(-(x/100)); .e6 <- exp(-(alpha * (1 - .e2))); .e8 <- exp(-(alpha * u1)) - 1; .e10 <- exp(-alpha) - 1; 0.01 * (.e6 * .e8 * .e2/(((.e6 - 1) * .e8/.e10 + 1) * .e10))})
  eval(expr2)
}

fAMH <- function(u1,x,alpha){
  #expr<- " u1*(1-exp(-x/100))/(1- alpha *(1-u1)*(exp(-x/100)))"
  #expr2<-parse(text = Deriv(expr, "x"))
  expr2<-expression({.e2 <- exp(-(x/100)); .e3 <- 1 - u1; .e4 <- 1 - alpha * .e3 * .e2; u1 * (0.01 - 0.01 * (alpha * (1 - .e2) * .e3/.e4)) * .e2/.e4})
  eval(expr2)
}

fJoe <- function(u1,x,alpha){
  #expr<- " 1-(1-(1-(1-u1)^(alpha))*(1-(1-(1-exp(-x/100)))^(alpha)))^(1/alpha)"
  #expr2<-parse(text = Deriv(expr, "x"))
  expr2<-expression({.e3 <- 1 - (1 - u1)^alpha; .e5 <- exp(-(x/100))^alpha; 0.01 * ((1 - .e3 * (1 - .e5))^(1/alpha - 1) * .e3 * .e5)})
  eval(expr2)
}

#### fonction a interger pour une copule de clayton,frank,amh avec 
#####################################################################marginales binomiales et exponentielle 
fonctionc <- function(n,x,m,q,alpha){
  FrClayton(pbinom(n-1,m,q),pexp(x,1/100),alpha) * (fClayton(pbinom(n,m,q),x,alpha)-fClayton(pbinom(n-1,m,q),x,alpha))
} 

fonctionk <- function(n,x,m,q,alpha){
  FrFrank(pbinom(n-1,m,q),pexp(x,1/100),alpha) * (fFrank(pbinom(n,m,q),x,alpha)-fFrank(pbinom(n-1,m,q),x,alpha))
}

fonctiona <- function(n,x,m,q,alpha){
  FrAMH(pbinom(n-1,m,q),pexp(x,1/100),alpha) * (fAMH(pbinom(n,m,q),x,alpha)-fAMH(pbinom(n-1,m,q),x,alpha))
}

fonctionj <- function(n,x,m,q,alpha){
  FrJoe(pbinom(n-1,m,q),pexp(x,1/100),alpha) * (fJoe(pbinom(n,m,q),x,alpha)-fJoe(pbinom(n-1,m,q),x,alpha))
}

######### Fonction qui calcule le tau de Kendall theorique 
#selon la loi de la copule mere et 
###########################################################les parametres de la binomiale
taukndmixte <- function(m,q,alpha,loi=c("frank","AMH","Clayton","Joe")){
  loi <- match.arg(loi)
  somme<-0
  if(loi == "Clayton") {
    for (i in 1:m) {
      aa<-try(integrate(function(x) fonctionc(i,x,m,q,alpha),0,1000)$value,silent = TRUE)
      if (!is.numeric(aa))    print(i)
      else somme <- somme + aa
    } }
  else if(loi == "frank") {
    for (i in 1:m) {
      aa<-try(integrate(function(x) fonctionk(i,x,m,q,alpha),0,1000)$value,silent = TRUE)
      if (!is.numeric(aa))    print(i)
      else somme <- somme + aa
    } }
  else if(loi == "AMH") {
    for (i in 1:m) {
      aa<-try(integrate(function(x) fonctiona(i,x,m,q,alpha),0,1000)$value,silent = TRUE)
      if (!is.numeric(aa))    print(i)
      else somme <- somme + aa
    } }
  else {
    for (i in 1:m) {
      aa<-try(integrate(function(x) fonctionj(i,x,m,q,alpha),0,1000)$value,silent = TRUE)
      if (!is.numeric(aa))    print(i)
      else somme <- somme + aa
    } }
  
  
  4*somme + sum(dbinom(1:m,m,q)^2) -1 
}


#####################################################################
## Fonction qui estime le parametre alpha0 selon la loi mere de la copule
## prend en entree la matrice des obs et en sort la valeur estimé de alpha0

estalpha<- function(mat,loi=c("log","geo","gamma","sibuya")){
  loi <- match.arg(loi)
  #tau<-TauKendallMixte2(mat,3)
  tau<- TauKendallMixte(mat)
  if(loi == "log") Val<- optimize(function(x) abs(taukndmixte(m,q,x,"frank")-tau),c(0,15))$minimum
  if(loi == "geo") Val<- optimize(function(x) abs(taukndmixte(m,q,x,"AMH")-tau),c(0,1))$minimum
  if(loi == "gamma") Val<- optimize(function(x) abs(taukndmixte(m,q,x,"clayton")-tau),c(0,15))$minimum
  if(loi == "sibuya") Val<- optimize(function(x) abs(taukndmixte(m,q,x,"Joe")-tau),c(1,15))$minimum
  sortie<- c(tau,Val)
  sortie
}


############################################### Alpha1

### Fonction qui estime alpha1 pour une copule hierarchique
# Prend en entrée une matrice des observation uniformes, loi mere, loi enfant, et alpha0 estimé
# retourne le tau calculé et alpha1

ChoixX<- function(obs){
  obs<- obs[!(obs[,1]<2),]
  sortie<- matrix(0,nrow(obs),2)
  pos<-sapply(obs[,1],function(x) sample.int(x,2,replace = FALSE))
  for (i in 1:nrow(obs)) {
    sortie[i,]<- obs[i,(pos[,i]+1)]
  }
  sortie
}

TauxRep<- function(obs1){
  knd<- numeric(100)
  for (k in 1:100) {
    obs<-ChoixX(obs1)
    knd[k]<-cor(obs,method = "kendall")[1,2]
    print(k)
    print(knd[k])
  }
  tau<- mean(knd)
  tau
}

ConstrCopule<-function(mere=c("geo","log"),enfant=c("geo","log","gamma")){
  mere<- match.arg(mere)
  enfant<- match.arg(enfant)
  
  if (mere=="geo") strm <- c("(1-alpha0)/(exp(u)-alpha0)","-log(1/(((1-alpha0)/(u))+alpha0))")
  else strm <- c(" -log(1-(1-exp(-alpha0))*exp(-(u)))/alpha0","-log((1-exp(-alpha0*(u)))/(1-exp(-alpha0)))")
  
  if (enfant=="geo") stre <- c(" (1-alpha1)/(exp(u)-alpha1)","-log(1/(((1-alpha1)/(u))+alpha1))")
  else if (enfant == "log") stre<- c(" -log(1-(1-exp(-alpha1))*exp(-(u)))/alpha1","-log((1-exp(-alpha1*(u)))/(1-exp(-alpha1)))")
  else stre <- c(" 1/((1+(u))^(1/alpha1))","(u)^(-alpha1) -1 ")
  
  Laptheta<- str_replace_all(strm[1],"u",paste0("-log(",as.character(stre[1]),")"))
  Lapthetainv<- str_replace_all(stre[2],"u",paste0("exp(-",as.character(strm[2]),")"))
  FrCop<- "chr"
  FrCop<-str_replace_all(Laptheta,"u","u1+u2")
  Lapthetainv<- str_replace_all(Lapthetainv,"u","u1")
  FrCop<-str_replace_all(FrCop,"u1",as.character(Lapthetainv))
  Lapthetainv<- str_replace_all(Lapthetainv,"u1","u2")
  FrCop<-str_replace_all(FrCop,"u2",as.character(Lapthetainv))
  fCop<- Deriv(Deriv(as.character(FrCop),"u1"),"u2")
  expres<- paste0("(",as.character(FrCop),")*(",as.character(fCop),")")
  expr<-parse(text = expres)
  return(expr)
}

fct <- function(obs,mere=c("geo","log"),enfant=c("geo","log","gamma"),alpha0,expr){
  mere<- match.arg(mere)
  enfant<- match.arg(enfant)
  
 #expr<- ConstrCopule(mere,enfant)
  
  PEval<- function(u1,u2,alpha0,alpha1,expr){
    eval(expr)
  }
  
  tau<-TauxRep(obs)
  
  taukndcont<- function(alpha0,alpha1,expr){
    probb<-integral2(function(u1,u2) PEval(u1,u2,alpha0,alpha1,expr),0,1,0,1)$`Q`
    4*probb-1
  }
  
  if( enfant == "geo") Val<- optimize(function(x) abs(taukndcont(alpha0,x,expr)-tau),c(0,1))$minimum
  else Val<- optimize(function(x) abs(taukndcont(alpha0,x,expr)-tau),c(0,15))$minimum
  sortie<- c(tau,Val)
  sortie
}


# ### met les couples (N,X) ensemble
# ## marche moins bien
# CoupleNX<- function(obs){
#   nbr<- sum(obs[,1])
#   k<-0
#   mat<- matrix(0,nrow=nbr,2)
#   for (i in 1:n) {
#     nn<- obs[i,1]
#     if(nn!=0){ 
#       for (j in 1:nn) {
#         k<-k+1
#         mat[k,]<- obs[i,c(1,j+1)]
#       }}
#   }
#   mat
# }
# 
# ######################### 
# 
# n<-10000  # nbr de simulation
# m<-7      # param de la binmiale et nbr de Xi simulés 
# q<-0.4    # parm de la binomiale
# 
# 
# ### geo-geo
# alpha0gg<-0.6
# alpha1gg<-0.4
# UGEOGEO<- rCompCop(n,GEO(1-alpha0gg, 1, list(GEO(1-alpha1gg, 1:m, NULL))))
# 
# obs1<- matrix(NA,n,m+1)
# obs1[,1]<- qbinom(UGEOGEO[,1],m,q)
# for (i in 1:n) {
#   count <- obs1[i,1]
#   if (count !=0)  obs1[i,2:(count+1)]<- qexp(UGEOGEO[i,2:(count+1)],1/100)
# }

#mat<- CoupleNX(obs1)

# obss1<- obs1
# for (i in 1:n) {
#   count <- obs1[i,1]
#   obss1[i,1] <- count
#   if (count !=0)  obss1[i,2:(count+1)]<- UGEOGEO[i,2:(count+1)]
# }
# 
# tempsalpha01<-system.time(alpha0geogeo<-estalpha(obs1,"geo"))
# exprgg<- ConstrCopule("geo","geo")
# tempsalpha11<- system.time(alpha1geogeo<-fct(obss1,"geo","geo",alpha0geogeo[2],exprgg))
# 
# 
# ##geo-log
# alpha0gl<-0.6
# alpha1gl<-7
# UGEOLOG<- rCompCop(n,GEO(1-alpha0gl, 1, list(LOG(1-exp(-alpha1gl), 1:m, NULL))))
# 
# obs2<- matrix(NA,n,m+1)
# obs2[,1]<- qbinom(UGEOLOG[,1],m,q)
# for (i in 1:n) {
#   count <- obs2[i,1]
#   if (count !=0)  obs2[i,2:(count+1)]<- qexp(UGEOLOG[i,2:(count+1)],1/100)
# }
# 
# obss2<- obs2
# for (i in 1:n) {
#   count <- obs2[i,1]
#   obss2[i,1] <- count
#   if (count !=0)  obss2[i,2:(count+1)]<- UGEOLOG[i,2:(count+1)]
# }
# 
# tempsalpha02<-system.time(alpha0geolog<-estalpha(obs2,"geo"))
# tempsalpha12<- system.time(alpha1geolog<-fct(obss2,"geo","log",alpha0geolog[2]))
# 
# 
# ## geo-gamma
# alpha0ga<-0.6
# alpha1ga<-4
# 
# samplet<- matrix(0,100,2)
# simplet<- matrix(0,100,2)
# for (l in 1:100) {
#   cat(l,"iteration")
#   UGEOGAMMA<- rCompCop(n,GEO(1-alpha0ga, 1, list(GAMMA(1/(alpha1ga), 1:m, NULL))))
#   
#   obs3<- matrix(NA,n,m+1)
#   obs3[,1]<- qbinom(UGEOGAMMA[,1],m,q)
#   for (i in 1:n) {
#     count <- obs3[i,1]
#     if (count !=0) obs3[i,2:(count+1)]<- qexp(UGEOGAMMA[i,2:(count+1)],1/100)
#   }
#   
#   obss3<- obs3
#   for (i in 1:n) {
#     count <- obs3[i,1]
#     obss3[i,1] <- count
#     if (count !=0) obss3[i,2:(count+1)]<- UGEOGAMMA[i,2:(count+1)]
#   }
#   
#   samplet[l,1]<-estalpha(obs3,"geo")[2]
#   expr<- ConstrCopule("geo","gamma")
#   samplet[l,2]<- fct(obss3,"geo","gamma",samplet[l,1],expr)[2]
#   simplet[l,1]<-estalpha2(obs3,"geo",0,1)[2]
#   expr2<- ConstrCopule2("geo","gamma")
#   simplet[l,2]<- fct2(obss3,"geo","gamma","opt",simplet[l,1],0,15,expr2)[2]
#   
#   
# }
# 
# 
# UGEOGAMMA<- rCompCop(n,GEO(1-alpha0ga, 1, list(GAMMA(1/(alpha1ga), 1:m, NULL))))
# 
# obs3<- matrix(NA,n,m+1)
# obs3[,1]<- qbinom(UGEOGAMMA[,1],m,q)
# for (i in 1:n) {
#   count <- obs3[i,1]
#   if (count !=0) obs3[i,2:(count+1)]<- qexp(UGEOGAMMA[i,2:(count+1)],1/100)
# }
# 
# obss3<- obs3
# for (i in 1:n) {
#   count <- obs3[i,1]
#   obss3[i,1] <- count
#   if (count !=0) obss3[i,2:(count+1)]<- UGEOGAMMA[i,2:(count+1)]
# }
# 
# tempsalpha03<-system.time(alpha0geogamma<-estalpha(obs3,"geo"))
# tempsalpha13<- system.time(alpha1geogamma<-fct(obss3,"geo","gamma",alpha0geogamma[2]))
# 
# 
# ### log-log
# alpha0ll<-4
# alpha1ll<-7
# ULOGLOG<- rCompCop(n,LOG(1-exp(-alpha0ll), 1, list(LOG(1-exp(-alpha1ll), 1:m, NULL))))
# 
# obs4<- matrix(NA,n,m+1)
# obs4[,1]<- qbinom(ULOGLOG[,1],m,q)
# for (i in 1:n) {
#   count <- obs4[i,1]
#   if (count !=0) obs4[i,2:(count+1)]<- qexp(ULOGLOG[i,2:(count+1)],1/100)
# }
# 
# obss4<- matrix(NA,n,m+1)
# for (i in 1:n) {
#   count <- obs4[i,1]
#   obss4[i,1] <- count
#   if (count !=0) obss4[i,2:(count+1)]<- ULOGLOG[i,2:(count+1)]
# }
# 
# tempsalpha04<-system.time(alpha0loglog<-estalpha(obs4,"log"))
# tempsalpha14<- system.time(alpha1loglog<-fct(obss4,"log","log",alpha0loglog[2]))
# 
# ##log-geo
# alpha0lg<-5
# alpha1lg<-0.6
# ULOGGEO<- rCompCop(n,LOG(1-exp(-alpha0lg), 1, list(GEO(1-alpha1lg, 1:m, NULL))))
# 
# obs5<- matrix(NA,n,m+1)
# obs5[,1]<- qbinom(ULOGGEO[,1],m,q)
# for (i in 1:n) {
#   count <- obs5[i,1]
#   if (count !=0)  obs5[i,2:(count+1)]<- qexp(ULOGGEO[i,2:(count+1)],1/100)
# }
# 
# obss5<- matrix(NA,n,m+1)
# for (i in 1:n) {
#   count <- obs5[i,1]
#   obss5[i,1] <- count
#   if (count !=0) obss5[i,2:(count+1)]<- ULOGGEO[i,2:(count+1)]
# }
# 
# tempsalpha05<-system.time(alpha0loggeo<-estalpha(obs5,"log"))
# tempsalpha15<- system.time(alpha1loggeo<-fct(obss5,"log","geo",alpha0loggeo[2]))
# 
# 
# ## log-gamma
# alpha0la<-6
# alpha1la<-4
# ULOGGAMMA<- rCompCop(n,LOG(1-exp(-alpha0la), 1, list(GAMMA(1/(alpha1la), 1:m, NULL))))
# 
# obs6<- matrix(NA,n,m+1)
# obs6[,1]<- qbinom(ULOGGAMMA[,1],m,q)
# for (i in 1:n) {
#   count <- obs6[i,1]
#   if (count !=0) obs6[i,2:(count+1)]<- qexp(ULOGGAMMA[i,2:(count+1)],1/100)
# }
# 
# obss6<- obs6
# for (i in 1:n) {
#   count <- obs6[i,1]
#   obss6[i,1] <- count
#   if (count !=0) obss6[i,2:(count+1)]<- ULOGGAMMA[i,2:(count+1)]
# }
# 
# tempsalpha06<-system.time(alpha0loggamma<-estalpha(obs6,"log"))
# tempsalpha16<- system.time(alpha1loggamma<-fct(obss6,"log","gamma",alpha0loggamma))
# 
# ################# Resultats
# 
# 
# res<-matrix(0,6,5)
# as.data.frame(res)
# rownames(res)<- c("Géométrique-Géométrique","Géométrique-Logarithmique","Géométrique-Gamma","Logarithmique-Logarithmique","Logarithmique-Géométrique","Logarithmique-Gamma")
# colnames(res)<- c("Paramètre alpha0","Estimateur alpha0","Paramètre alpha1","Estimateur alpha1","temps")
# 
# res[1,]<- c(alpha0gg,alpha0geogeo[2],alpha1gg,alpha1geogeo[2],tempsalpha01[1]+tempsalpha11[1])
# res[2,]<- c(alpha0gl,alpha0geolog[2],alpha1gl,alpha1geolog[2],tempsalpha02[1]+tempsalpha12[1])
# res[3,]<- c(alpha0ga,alpha0geogamma[2],alpha1ga,alpha1geogamma[2],tempsalpha03[1]+tempsalpha13[1])
# 
# res[4,]<- c(alpha0ll,alpha0loglog[2],alpha1ll,alpha1loglog[2],tempsalpha04[1]+tempsalpha14[1])
# res[5,]<- c(alpha0lg,alpha0loggeo[2],alpha1lg,alpha1loggeo[2],tempsalpha05[1]+tempsalpha15[1])
# res[6,]<- c(alpha0la,alpha0loggamma[2],alpha1la,alpha1loggamma[2],tempsalpha06[1]+tempsalpha16[1])
# res
# #ressample<- res
# xtable(res,digits = 4)
