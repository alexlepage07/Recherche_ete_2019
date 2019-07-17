library(pscl)

setwd("C:/Users/Alex/Desktop/Recherche_ete_2019/Code_Alex/Massachussetts_2006")
FREQ <- read.csv("FREQ.txt")


# tgroup : 1 = least risky, 6 = most risky
# class4 : A = Adult, B = Business, I = <3 yrs experience, M = 3-6 yrs experience, S = Senior

FREQ$class4 <- FREQ$class4 %% 10
FREQ$class4[FREQ$class4 == 1] <- "A"
FREQ$class4[FREQ$class4 == 2] <- "S"
FREQ$class4[FREQ$class4 %in% c(3, 4)] <- "M"
FREQ$class4[FREQ$class4 == 5] <- "B"
FREQ$class4[FREQ$class4 %in% (6:9)] <- "I"
FREQ$class4 <- factor(FREQ$class4)

FREQ$tgroup <- factor(FREQ$tgroup, 1:6)
FREQ$clm_id <- factor(FREQ$clm_id)
FREQ$pol_id <- factor(FREQ$pol_id)
summary(FREQ)
head(FREQ)
attach(FREQ)

agg_means <- aggregate(CompteDeclm_id, by=list(class4, tgroup), mean)
mat_means <- matrix(agg_means$x, ncol = 5, byrow = T)
colnames(mat_means) <- agg_means$Group.1[1:5]
rownames(mat_means) <- 1:6
mat_means

n_obs <- nrow(FREQ)
train.id <- 1:1e+6
test.id <- (max(train.id)+1):n_obs

hurdle.model <- hurdle(CompteDeclm_id ~ class4 + tgroup, data = FREQ, offset = earnexpo)
(Sumry_1 <- summary(hurdle.model))

hurdle.model_2 <-
    hurdle(CompteDeclm_id ~ class4 + tgroup | 1,
           data = FREQ ,
           offset = earnexpo)
(Sumry_2 <- summary(hurdle.model_2))

Sumry_1$loglik - Sumry_2$loglik
qchisq(0.95, 11)
