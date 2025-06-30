library(Matrix)
library(MASS)
library(mgcv)
library(psych)
library(splines2)
library(fda)

# Generate the Data Sample

N <- 100

u_grid <- seq(0,1,length.out = 500)

cov_X1 <- rnorm(N, 0, 2)
cov_X2 <- rbinom(N, 100, 0.5)
cov_X <- cbind(cov_X1, cov_X2)
cov_Z <- 2 + matrix(rep(u_grid, each = N), nrow = N) + rnorm(N,0,1) * u_grid^2

# Scenario I

alpha <- c(1,0.5)
beta_0_init <- 3 + sin(2 * pi * u_grid)+2*cos(2 * pi * u_grid)
beta_1_init <- 1 + 2 * u_grid^3
beta_2_init <- 2 + 5 * u_grid
beta_0 <- (beta_0_init - mean(beta_0_init)) / sd(beta_0_init)
beta_1 <- (beta_1_init - mean(beta_1_init)) / sd(beta_1_init)
beta_2 <- (beta_2_init - mean(beta_2_init)) / sd(beta_2_init)

g_0_fun <- function(s) { 3 * s + 4 * s^2 }
g_1_fun <- function(s) { -2 * s - 3 * s^3 }
g_2_fun <- function(s) { 2 * s^2 + 5 * s^3 }

basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = 10)
cov_Z_fd <- Data2fd(argvals = u_grid, y = cov_Z, basisobj = basis)

integral_beta_0 <- apply(cov_Z, 1, function(Z_i) {
  cov_Z_fd_i <- Data2fd(argvals = u_grid, y = Z_i, basisobj = basis)
  integral_value <- eval.fd(u_grid, cov_Z_fd) %*% beta_0  
  return(integral_value)
})

integral_beta_1 <- apply(cov_Z, 1, function(Z_i) {
  cov_Z_fd_i <- Data2fd(argvals = u_grid, y = Z_i, basisobj = basis)
  integral_value <- eval.fd(u_grid, cov_Z_fd) %*% beta_1  
  return(integral_value)
})

integral_beta_2 <- apply(cov_Z, 1, function(Z_i) {
  cov_Z_fd_i <- Data2fd(argvals = u_grid, y = Z_i, basisobj = basis)
  integral_value <- eval.fd(u_grid, cov_Z_fd_i) %*% beta_2  
  return(integral_value)
})

g_0_values <- g_0_fun(integral_beta_0)
g_1_values <- g_1_fun(integral_beta_1)
g_2_values <- g_2_fun(integral_beta_2)

respon_Y <- cov_X %*% alpha + g_0_values + cov_X[,1] * g_1_values + cov_X[,2] * g_2_values  + rnorm(N,0,1)

# Scenario II

alpha <- 8 * tau * c(1,0.5)
beta_0_init <- 3 + 3 * tau * sin(2 * pi * u_grid)+ 6 * tau * cos(2 * pi * u_grid)
beta_1_init <- 1 + 6 * tau * u_grid^3
beta_2_init <- 2 + 10  * tau * u_grid
beta_0 <- (beta_0_init - mean(beta_0_init)) / sd(beta_0_init)
beta_1 <- (beta_1_init - mean(beta_1_init)) / sd(beta_1_init)
beta_2 <- (beta_2_init - mean(beta_2_init)) / sd(beta_2_init)

g_0_fun <- function(s) { 3 * s + 4 * s^2 }
g_1_fun <- function(s) { -2 * s - 3 * s^3 }
g_2_fun <- function(s) { 2 * s^2 + 5 * s^3 }

basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = 10)
cov_Z_fd <- Data2fd(argvals = u_grid, y = cov_Z, basisobj = basis)

integral_beta_0 <- apply(cov_Z, 1, function(Z_i) {
  cov_Z_fd_i <- Data2fd(argvals = u_grid, y = Z_i, basisobj = basis)
  integral_value <- eval.fd(u_grid, cov_Z_fd) %*% beta_0  
  return(integral_value)
})

integral_beta_1 <- apply(cov_Z, 1, function(Z_i) {
  cov_Z_fd_i <- Data2fd(argvals = u_grid, y = Z_i, basisobj = basis)
  integral_value <- eval.fd(u_grid, cov_Z_fd) %*% beta_1  
  return(integral_value)
})

integral_beta_2 <- apply(cov_Z, 1, function(Z_i) {
  cov_Z_fd_i <- Data2fd(argvals = u_grid, y = Z_i, basisobj = basis)
  integral_value <- eval.fd(u_grid, cov_Z_fd_i) %*% beta_2  
  return(integral_value)
})

g_0_values <- g_0_fun(integral_beta_0)
g_1_values <- g_1_fun(integral_beta_1)
g_2_values <- g_2_fun(integral_beta_2)

respon_Y <- cov_X %*% alpha + g_0_values + cov_X[,1] * g_1_values + cov_X[,2] * g_2_values  + rnorm(N,0,1)
