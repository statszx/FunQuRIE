library(Matrix)
library(MASS)
library(mgcv)
library(psych)
library(splines2)
library(fda)

# Hypothesis Test for the Interaction Effects

compute_rho_prime <- function(tau) {
  return(function(x) { dnorm(x, mean = tau, sd = 1) })  
}

compute_kernel <- function(h) {
  return(function(x) { exp(-x^2 / (2 * h^2)) }) 
}

convolution_derivative <- function(rho_tau_prime, K_h, u_grid) {
  return(sapply(u_grid, function(u) rho_tau_prime(u) * K_h(u)))
}

compute_PZ <- function(X, Psi_0) {
  Z <- cbind(X, Psi_0) 
  PZ <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
  return(PZ)
}

compute_R <- function(PZ, W) {
  I <- diag(nrow(PZ)) 
  R <- (I - PZ) %*% W  
  return(R)
}

objective_function <- function(alpha, xi_list, gamma_list, xi_mat, X, Z, Y, h, tau, u_grid) {
  N <- nrow(X)
  loss <- 0
  
  rho_tau_prime <- compute_rho_prime(tau)  
  K_h <- compute_kernel(h)  
  conv_deriv <- convolution_derivative(rho_tau_prime, K_h, u_grid)
  
  for (i in 1:N) {
    residual <- Y[i] - (X[i, ] %*% alpha + sum(xi_list[[1]] * gamma_list[[1]]))
    loss <- loss + residual^2 
  }
  
  loss <- loss + sum(conv_deriv^2)  
  return(loss)
}


gradient_descent <- function(alpha, xi_list, gamma_list, xi_mat, X, Z, Y, tau, h, lr = 0.01, max_iter = 50) {
  loss_trace <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    loss <- objective_function(alpha, xi_list, gamma_list, xi_mat, X, Z, Y, h, tau, u_grid)
    loss_trace[iter] <- loss
    cat("Iter:", iter, "Loss:", loss, "\n")
    
    grad_alpha <- rep(0, length(alpha))
    for (i in 1:nrow(X)) {
      Ai <- X[i, ]
      residual <- Y[i] - sum(Ai * alpha) - sum(xi_list[[1]] * eval.basis(sum(gamma_list[[1]] * Ai), beta_basis, Lfdobj = 1))
      grad_alpha <- grad_alpha - 2 * residual * Ai 
    }
    
    grad_gamma <- rep(0, length(gamma_list[[1]]))
    for (i in 1:nrow(X)) {
      Ai <- X[i, ]
      residual <- Y[i] - sum(Ai * alpha) - sum(xi_list[[1]] * eval.basis(sum(gamma_list[[1]] * Ai), beta_basis, Lfdobj = 1))
      grad_gamma <- grad_gamma - 2 * residual * eval.basis(sum(gamma_list[[1]] * Ai), beta_basis, deriv = 1) * xi_list[[1]]
    }
    
    grad_xi <- rep(0, length(xi_list[[1]]))
    for (i in 1:nrow(X)) {
      Ai <- X[i, ]
      residual <- Y[i] - sum(Ai * alpha) - sum(xi_list[[1]] * eval.basis(sum(gamma_list[[1]] * Ai), beta_basis, Lfdobj = 1))
      grad_xi <- grad_xi - 2 * residual * eval.basis(sum(gamma_list[[1]] * Ai), beta_basis)
    }
    
    alpha <- alpha - lr * grad_alpha
    gamma_list[[1]] <- gamma_list[[1]] - lr * grad_gamma
    xi_list[[1]] <- xi_list[[1]] - lr * grad_xi
  }
  
  return(list(alpha = alpha, xi_list = xi_list, gamma_list = gamma_list, loss_trace = loss_trace))
}


compute_VN <- function(R, Y, X, Z, alpha, gamma_list, xi_list, u_grid, h, tau) {
  V_N <- numeric(length(R))
  
  rho_tau_prime <- compute_rho_prime(tau)  
  K_h <- compute_kernel(h)
  conv_deriv <- convolution_derivative(rho_tau_prime, K_h, u_grid)
  
  for (i in 1:length(R)) {
    Ai <- X[i, ]
    residual <- Y[i] - sum(Ai * alpha) - sum(xi_list[[1]] * eval.basis(sum(gamma_list[[1]] * Ai), beta_basis, Lfdobj = 1))
    V_N[i] <- sqrt(length(R)) * sum(R[i, ]) * conv_deriv[i] * residual
  }
  
  return(V_N)
}

compute_SigmaN <- function(R, Y, X, Z, alpha, gamma_list, xi_list, u_grid, h, tau) {
  Sigma_N <- matrix(0, nrow = length(R), ncol = length(R))
  
  rho_tau_prime <- compute_rho_prime(tau)
  K_h <- compute_kernel(h)
  conv_deriv <- convolution_derivative(rho_tau_prime, K_h, u_grid)
  
  for (i in 1:length(R)) {
    Ai <- X[i, ]
    residual <- Y[i] - sum(Ai * alpha) - sum(xi_list[[1]] * eval.basis(sum(gamma_list[[1]] * Ai), beta_basis, Lfdobj = 1))
    Sigma_N[i, i] <- sqrt(length(R)) * sum(R[i, ]) * conv_deriv[i] * residual^2
  }
  
  return(Sigma_N)
}


compute_test_statistic <- function(V_N, Sigma_N, p, m) {
  Lambda_N <- t(V_N) %*% solve(Sigma_N) %*% V_N
  S_N <- (Lambda_N - p * m) / sqrt(2 * p * m)
  return(S_N)
}

run_test <- function(X, Z, Y, tau, h, u_grid, p, m, max_iter = 50) {

  alpha_init <- rep(0, ncol(X)) 
  xi_list_init <- list(rep(0, ncol(X)))  
  gamma_list_init <- list(rep(0, ncol(X)))  

  Psi_0 <- cbind(rep(1, nrow(X)))  
  PZ <- compute_PZ(X, Psi_0) 

  W <- matrix(rnorm(nrow(X) * p), nrow = nrow(X), ncol = p)  
  R <- compute_R(PZ, W)  

  fit_result <- gradient_descent(alpha_init, xi_list_init, gamma_list_init, NULL, X, Y, Z, phi, tau, h, u_grid, lr = 0.01, max_iter = max_iter)
  
  alpha_hat <- fit_result$alpha
  xi_list_hat <- fit_result$xi_list
  gamma_list_hat <- fit_result$gamma_list
  
  V_N <- compute_VN(R, Y, X, Z, alpha_hat, gamma_list_hat, xi_list_hat, u_grid, h, tau)
  Sigma_N <- compute_SigmaN(R, Y, X, Z, alpha_hat, gamma_list_hat, xi_list_hat, u_grid, h, tau)

  S_N <- compute_test_statistic(V_N, Sigma_N, p, m)

  return(S_N)
}

S_N_result <- run_test(X, Z, Y, tau, h, u_grid, p, m)


