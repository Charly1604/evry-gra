
# Packages
library(tidyverse)
library(readxl)
library(quadprog)

# Data importation GMV and TP
(workshop1 <- read_xlsx("data.xlsx", sheet = "workshop1", skip = 4))
fin_return1 <- workshop1 %>% select(-"TRADE DATE")
(workshop2 <- read_xlsx("data.xlsx", sheet = "workshop2", skip = 3))
fin_return2 <- workshop2 %>% select(-"TRADE DATE")

# Modelling setup GMV
n <- ncol(fin_return1)
T <- nrow(fin_return1)
e <- rep(1, n)
perio <- 12

# Empirical GMV
Sigma <- cov(fin_return1) * (T - 1) / (T - n - 2) * perio
C <- t(e) %*% solve(Sigma) %*% e
sigmag <- sqrt(1 / C)
omega <- 1 / as.numeric(C) * solve(Sigma) %*% e
barplot(as.numeric(omega), col = 'black')

# Empirical constrained GMV with no short selling
# min(-d^T b + 1/2 b^T D b) with the constraints A^T b >= b_0
dvec <-  numeric(n)
# Constraints of sum of the weights equal to 1 and positive weights
Amat <-  cbind(e, diag(1, n, n))
bvec <-  cbind(1, t(numeric(n)))
sigmag_constraint <- sqrt(t(omega_constraint) %*% Sigma %*% omega_constraint)
omega_constraint <- solve.QP(Sigma, dvec, Amat, bvec, meq = 1)$solution
barplot(as.numeric(omega_constraint), col = 'black')

# Modelling setup TP
n <- ncol(fin_return2)
T <- nrow(fin_return2)
e <- rep(1, n)
rf <- 0

# Strategic TP
mu <- colMeans(fin_return2) * perio - rf
Sigma <- cov(fin_return2) * (T - 1) / (T - n - 2) * perio
A <- t(e) %*% solve(Sigma) %*% mu
omega <- 1 / as.numeric(A) * solve(Sigma) %*% mu
sigma <- sqrt(t(omega) %*% Sigma %*% omega)
barnames <- c('France Equity', 'BRIC Equity', 'US Corporate Bond')
barplot(as.numeric(omega), col = 'black', names.arg = barnames, ylim = c(0, 1))

# Strategic constrained portfolio with TP target return and no short selling
# min(-d^T b + 1/2 b^T D b) with the constraints A^T b >= b_0
dvec <-  numeric(n)
# Constraints of TP target return, sum of the weights equal to 1 and positive weights
Amat <-  cbind(mu, e, diag(1, n, n))
target <- as.numeric(t(mu) %*% solve(Sigma) %*% mu/A)
bvec <-  cbind(target, 1, t(numeric(n)))
omega_constraint <-  solve.QP(Sigma, dvec, Amat, bvec, meq = 2)$solution
sigma_constraint <- sqrt(t(omega_constraint) %*% Sigma %*% omega_constraint)
barplot(as.numeric(omega_constraint), col = 'black', names.arg = barnames,ylim = c(0, 1))

