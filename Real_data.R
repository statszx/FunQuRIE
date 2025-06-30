library(Matrix)
library(MASS)
library(mgcv)
library(psych)
library(splines2)
library(fda)
library(fds)
library(fda.usc)


# Biscuit Data Analysis

respon_Y <- cbind(labp[2,],labc[2,])
cov_X <- c(1,cbind(labp[1,],labc[1,]))
cov_Z <- rbind(nirp$y,nirc$y)
u_grid <- (nirp$x - min(nirp$x)) / (max(nirp$x) - min(nirp$x))

knots <- seq(range(u_grid)[1],range(u_grid)[2], length.out=M+1)
nknots <- length(knots)
nbasis <- nknots+norder-2
beta_basis <- create.bspline.basis(knots,nbasis,norder)
cov_Z_fd <- Data2fd(argvals = u_grid, y = cov_Z, basisobj = basis)

convolution_derivative <- function(rho_tau, K_h, u_grid) {
  
  convolved_derivative <- sapply(u_grid, function(u) {
    integral_value <- sum(rho_tau * sapply(u_grid, function(v) K_h(v - u)))
    return(integral_value)
  })
  return(convolved_derivative)
}


objective_function <- function(alpha, xi_list, gamma_list, xi_mat, X, Z, Y, h, tau) {
  N <- nrow(X)
  loss <- 0
  
  rho_tau_prime <- compute_rho_prime(tau)  
  K_h <- compute_kernel(h)  
  conv_deriv <- convolution_derivative(rho_tau_prime, K_h, u_grid)
  
  for (i in 1:N) {
    par_mat <- X[i, ]
    Zi <- Z[[i]] 
    
    integral_term <- sum(Zi * phi) 
    g0_i <- sum(xi_list[[1]] * eval.basis(integral_term * gamma_list[[1]], beta_basis, Lfdobj = 1)) 
    
    residual <- Y[i] - sum(par_mat * alpha) - g0_i
    loss <- loss + residual^2  
  }
  
  loss <- loss + sum(conv_deriv^2)  
  return(loss)
}


gradient_descent <- function(alpha, xi_list, gamma_list, xi_mat, X, Z, Y, tau = 0.5, h = 0.05, lr = 0.01, max_iter = 50) {
  loss_trace <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    loss <- objective_function(alpha, xi_list, gamma_list, xi_mat, X, Y, h, tau)
    loss_trace[iter] <- loss
    cat("Iter:", iter, "Loss:", loss, "\n")
    
    for (j in 1:length(alpha)) {
      alpha_eps <- alpha
      alpha_eps[j] <- alpha_eps[j] + 1e-5
      grad <- (objective_function(alpha_eps, xi_list, gamma_list, xi_mat, X, Z, Y, h, tau) - loss) / 1e-5
      alpha[j] <- alpha[j] - lr * grad
    }
    
    for (j in 1:length(gamma_list)) {
      gamma_eps <- gamma_list
      gamma_eps[[j]] <- gamma_eps[[j]] + 1e-5
      grad <- (objective_function(alpha, xi_list, gamma_eps, xi_mat, X, Z, Y, h, tau) - loss) / 1e-5
      gamma_list[[j]] <- gamma_list[[j]] - lr * grad
    }
    
    for (j in 1:length(xi_list)) {
      for (k in 1:length(xi_list[[j]])) {
        xi_eps <- xi_list
        xi_eps[[j]][k] <- xi_eps[[j]][k] + 1e-5
        grad <- (objective_function(alpha, xi_eps, gamma_list, xi_mat, X, Z, Y, h, tau) - loss) / 1e-5
        xi_list[[j]][k] <- xi_list[[j]][k] - lr * grad
      }
    }
  }
  
  return(list(alpha = alpha, xi_list = xi_list, gamma_list = gamma_list, loss_trace = loss_trace))
}


iterated_estimation <- function(X, Z, Y, u_grid, alpha_init, xi_list_init, gamma_list_init, xi_mat_init, tau = 0.5, h = 0.05, max_iter = 50, lr = 0.01, convergence_threshold = 1e-5) {
  
  alpha <- alpha_init
  xi_list <- xi_list_init
  gamma_list <- gamma_list_init
  xi_mat <- xi_mat_init
  
  loss_trace <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    cat("Iteration", iter, "\n")
    
    grad_results <- compute_gradients(alpha, xi_list, gamma_list, xi_mat, X, Z, Y, tau, h)
    
    grad_gamma <- grad_results$grad_alpha  
    alpha <- alpha - lr * grad_gamma 
    
    grad_xi <- grad_results$grad_xi_mat  
    xi_mat <- xi_mat - lr * grad_xi  
    
    grad_xi_list <- grad_results$grad_xi_list
    grad_gamma_list <- grad_results$grad_gamma_list
    
    
    for (j in 1:length(xi_list)) {
      xi_list[[j]] <- xi_list[[j]] - lr * grad_xi_list[[j]]
    }
    for (j in 1:length(gamma_list)) {
      gamma_list[[j]] <- gamma_list[[j]] - lr * grad_gamma_list[[j]]
    }
    
    loss <- objective_function(alpha, xi_list, gamma_list, xi_mat, X, Z, Y, h, tau)
    loss_trace[iter] <- loss
    cat("Loss:", loss, "\n")
    
    
    if (iter > 1) {
      alpha_diff <- max(abs(alpha - prev_alpha))
      xi_diff <- max(abs(xi_mat - prev_xi_mat))
      xi_list_diff <- max(sapply(1:length(xi_list), function(j) max(abs(xi_list[[j]] - prev_xi_list[[j]]))))
      gamma_list_diff <- max(sapply(1:length(gamma_list), function(j) max(abs(gamma_list[[j]] - prev_gamma_list[[j]]))))
      
      if (alpha_diff < convergence_threshold && xi_diff < convergence_threshold && xi_list_diff < convergence_threshold && gamma_list_diff < convergence_threshold) {
        cat("Convergence achieved.\n")
        break
      }
    }
    
    prev_alpha <- alpha
    prev_xi_mat <- xi_mat
    prev_xi_list <- xi_list
    prev_gamma_list <- gamma_list
  }
  
  return(list(alpha = alpha, xi_list = xi_list, gamma_list = gamma_list, xi_mat = xi_mat, loss_trace = loss_trace))
}

get_beta_g0_gj <- function(alpha, gamma_list, xi_list) {
  beta <- alpha 
  g0 <- gamma_list[[1]] 
  gj <- xi_list[[1]]
  
  return(list(beta = beta, g0 = g0, gj = gj))
}

alpha_init <- rep(0, p)  
xi_list_init <- list(rep(0, K))  
gamma_list_init <- list(rep(0, K))  
xi_mat_init <- matrix(0, N, p)  

result <- iterated_estimation(cov_X, cov_Z, response_Y, u_grid, alpha_init, xi_list_init, gamma_list_init, xi_mat_init, tau = 0.5, h = 0.05)
beta_g0_gj <- get_beta_g0_gj(result$alpha, result$gamma_list, result$xi_list)


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
  
  fit_result <- gradient_descent(alpha_init, xi_list_init, gamma_list_init, NULL, X, Y, Z, tau, h, u_grid, lr = 0.01, max_iter = max_iter)
  
  alpha_hat <- fit_result$alpha
  xi_list_hat <- fit_result$xi_list
  gamma_list_hat <- fit_result$gamma_list
  
  V_N <- compute_VN(R, Y, X, Z, alpha_hat, gamma_list_hat, xi_list_hat, u_grid, h, tau)
  Sigma_N <- compute_SigmaN(R, Y, X, Z, alpha_hat, gamma_list_hat, xi_list_hat, u_grid, h, tau)
  
  S_N <- compute_test_statistic(V_N, Sigma_N, p, m)
  
  return(S_N)
}

S_N_result <- run_test(cov_X, cov_Z, response_Y, tau=0.05, h, u_grid, p, m)


# Tecator Data Analysis

respon_Y <- tecator$y$Fat
cov_X <- c(1,tecator$y$Protein[,1])
cov_Z <- tecator$absorp.fdata$data
u_point <- tecator$absorp.fdata$argvals
u_grid <- (u_point - min(u_point)) / (max(u_point) - min(u_point))

knots <- seq(range(u_grid)[1],range(u_grid)[2], length.out=M+1)
nknots <- length(knots)
nbasis <- nknots+norder-2
beta_basis <- create.bspline.basis(knots,nbasis,norder)
cov_Z_fd <- Data2fd(argvals = u_grid, y = cov_Z, basisobj = basis)

convolution_derivative <- function(rho_tau, K_h, u_grid) {
  
  convolved_derivative <- sapply(u_grid, function(u) {
    integral_value <- sum(rho_tau * sapply(u_grid, function(v) K_h(v - u)))
    return(integral_value)
  })
  return(convolved_derivative)
}


objective_function <- function(alpha, xi_list, gamma_list, xi_mat, X, Z, Y, h, tau) {
  N <- nrow(X)
  loss <- 0
  
  rho_tau_prime <- compute_rho_prime(tau)  
  K_h <- compute_kernel(h)  
  conv_deriv <- convolution_derivative(rho_tau_prime, K_h, u_grid)
  
  for (i in 1:N) {
    par_mat <- X[i, ]
    Zi <- Z[[i]] 
    
    integral_term <- sum(Zi * phi) 
    g0_i <- sum(xi_list[[1]] * eval.basis(integral_term * gamma_list[[1]], beta_basis, Lfdobj = 1)) 
    
    residual <- Y[i] - sum(par_mat * alpha) - g0_i
    loss <- loss + residual^2  
  }
  
  loss <- loss + sum(conv_deriv^2)  
  return(loss)
}


gradient_descent <- function(alpha, xi_list, gamma_list, xi_mat, X, Z, Y, tau = 0.5, h = 0.05, lr = 0.01, max_iter = 50) {
  loss_trace <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    loss <- objective_function(alpha, xi_list, gamma_list, xi_mat, X, Y, h, tau)
    loss_trace[iter] <- loss
    cat("Iter:", iter, "Loss:", loss, "\n")
    
    for (j in 1:length(alpha)) {
      alpha_eps <- alpha
      alpha_eps[j] <- alpha_eps[j] + 1e-5
      grad <- (objective_function(alpha_eps, xi_list, gamma_list, xi_mat, X, Z, Y, h, tau) - loss) / 1e-5
      alpha[j] <- alpha[j] - lr * grad
    }
    
    for (j in 1:length(gamma_list)) {
      gamma_eps <- gamma_list
      gamma_eps[[j]] <- gamma_eps[[j]] + 1e-5
      grad <- (objective_function(alpha, xi_list, gamma_eps, xi_mat, X, Z, Y, h, tau) - loss) / 1e-5
      gamma_list[[j]] <- gamma_list[[j]] - lr * grad
    }
    
    for (j in 1:length(xi_list)) {
      for (k in 1:length(xi_list[[j]])) {
        xi_eps <- xi_list
        xi_eps[[j]][k] <- xi_eps[[j]][k] + 1e-5
        grad <- (objective_function(alpha, xi_eps, gamma_list, xi_mat, X, Z, Y, h, tau) - loss) / 1e-5
        xi_list[[j]][k] <- xi_list[[j]][k] - lr * grad
      }
    }
  }
  
  return(list(alpha = alpha, xi_list = xi_list, gamma_list = gamma_list, loss_trace = loss_trace))
}


iterated_estimation <- function(X, Z, Y, u_grid, alpha_init, xi_list_init, gamma_list_init, xi_mat_init, tau = 0.5, h = 0.05, max_iter = 50, lr = 0.01, convergence_threshold = 1e-5) {
  
  alpha <- alpha_init
  xi_list <- xi_list_init
  gamma_list <- gamma_list_init
  xi_mat <- xi_mat_init
  
  loss_trace <- numeric(max_iter)
  
  for (iter in 1:max_iter) {
    cat("Iteration", iter, "\n")
    
    grad_results <- compute_gradients(alpha, xi_list, gamma_list, xi_mat, X, Z, Y, tau, h)
    
    grad_gamma <- grad_results$grad_alpha  
    alpha <- alpha - lr * grad_gamma 
    
    grad_xi <- grad_results$grad_xi_mat  
    xi_mat <- xi_mat - lr * grad_xi  
    
    grad_xi_list <- grad_results$grad_xi_list
    grad_gamma_list <- grad_results$grad_gamma_list
    
    
    for (j in 1:length(xi_list)) {
      xi_list[[j]] <- xi_list[[j]] - lr * grad_xi_list[[j]]
    }
    for (j in 1:length(gamma_list)) {
      gamma_list[[j]] <- gamma_list[[j]] - lr * grad_gamma_list[[j]]
    }
    
    loss <- objective_function(alpha, xi_list, gamma_list, xi_mat, X, Z, Y, h, tau)
    loss_trace[iter] <- loss
    cat("Loss:", loss, "\n")
    
    
    if (iter > 1) {
      alpha_diff <- max(abs(alpha - prev_alpha))
      xi_diff <- max(abs(xi_mat - prev_xi_mat))
      xi_list_diff <- max(sapply(1:length(xi_list), function(j) max(abs(xi_list[[j]] - prev_xi_list[[j]]))))
      gamma_list_diff <- max(sapply(1:length(gamma_list), function(j) max(abs(gamma_list[[j]] - prev_gamma_list[[j]]))))
      
      if (alpha_diff < convergence_threshold && xi_diff < convergence_threshold && xi_list_diff < convergence_threshold && gamma_list_diff < convergence_threshold) {
        cat("Convergence achieved.\n")
        break
      }
    }
    
    prev_alpha <- alpha
    prev_xi_mat <- xi_mat
    prev_xi_list <- xi_list
    prev_gamma_list <- gamma_list
  }
  
  return(list(alpha = alpha, xi_list = xi_list, gamma_list = gamma_list, xi_mat = xi_mat, loss_trace = loss_trace))
}

get_beta_g0_gj <- function(alpha, gamma_list, xi_list) {
  beta <- alpha 
  g0 <- gamma_list[[1]] 
  gj <- xi_list[[1]]
  
  return(list(beta = beta, g0 = g0, gj = gj))
}

alpha_init <- rep(0, p)  
xi_list_init <- list(rep(0, K))  
gamma_list_init <- list(rep(0, K))  
xi_mat_init <- matrix(0, N, p)  

result <- iterated_estimation(cov_X, cov_Z, response_Y, u_grid, alpha_init, xi_list_init, gamma_list_init, xi_mat_init, tau = 0.5, h = 0.05)
beta_g0_gj <- get_beta_g0_gj(result$alpha, result$gamma_list, result$xi_list)


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
  
  fit_result <- gradient_descent(alpha_init, xi_list_init, gamma_list_init, NULL, X, Y, Z, tau, h, u_grid, lr = 0.01, max_iter = max_iter)
  
  alpha_hat <- fit_result$alpha
  xi_list_hat <- fit_result$xi_list
  gamma_list_hat <- fit_result$gamma_list
  
  V_N <- compute_VN(R, Y, X, Z, alpha_hat, gamma_list_hat, xi_list_hat, u_grid, h, tau)
  Sigma_N <- compute_SigmaN(R, Y, X, Z, alpha_hat, gamma_list_hat, xi_list_hat, u_grid, h, tau)
  
  S_N <- compute_test_statistic(V_N, Sigma_N, p, m)
  
  return(S_N)
}

S_N_result <- run_test(cov_X, cov_Z, response_Y, tau=0.05, h, u_grid, p, m)
