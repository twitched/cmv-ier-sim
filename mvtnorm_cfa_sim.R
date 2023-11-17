library(mvtnorm)
library(MASS)
library(semPlot)
library(lavaan)

J <- 1000
I <- 6
K <- 2
psi <- matrix(c(1, 0.5,
                0.5, 0.8), nrow = K)  
beta <- seq(1, 2, by = .2)

# loading matrix
Lambda <- cbind(c(1, 1.5, 2, 0, 0, 0), c(0, 0, 0, 1, 1.5, 2))

# error covariance
Theta <- diag(0.3, nrow = I)

# factor scores
eta <- mvrnorm(J, mu = c(0, 0), Sigma = psi)

# error term
epsilon <- mvrnorm(J, mu = rep(0, ncol(Theta)),Sigma = Theta)

dat <- tcrossprod(eta, Lambda) + epsilon
dat_cfa  <-  dat %>% as.data.frame() %>% setNames(c("Y1", "Y2", "Y3", "Y4", "Y5", "Y6"))

lavaan_cfa <- 'eta1 =~ Y1 + Y2 + Y3
               eta2 =~ Y4 + Y5 + Y6'

semPaths(semPlotModel_lavaanModel(lavaan_cfa))
lav_cfa_fit <- cfa(lavaan_cfa, data = dat_cfa, meanstructure = TRUE, std.lv=TRUE)
summary(lav_cfa_fit, fit.measures = TRUE, standardized=TRUE)
inspect(lav_cfa_fit,"cov.ov")
standardizedsolution(lav_cfa_fit)
