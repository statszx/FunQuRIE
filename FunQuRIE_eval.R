library(Matrix)
library(MASS)
library(mgcv)
library(psych)
library(fdapace)
library(splines2)
library(fda)

# Evaluate Estimations from the FunQuRIE and FunQuR Models

rho_tau <- function(u, tau) {
  u * (tau - as.numeric(u < 0))
}

K_h <- function(t, h) {
  dnorm(t / h) / h  
}

gaussian_kernel <- function(u, h) {
  return((1 / (h * sqrt(2 * pi))) * exp(-0.5 * (u / h)^2))  
}

smooth_rho_tau_convolution <- function(u, h, tau) {
  delta <- mean(diff(u))  
  kernel_values <- sapply(u, function(v) {
    sum(sapply(u, function(v_prime) {
      rho_tau(v_prime, tau) * gaussian_kernel(v_prime - v, h)  
    })) * delta  
  })
  return(kernel_values)
}


knots <- seq(range(u_grid)[1],range(u_grid)[2], length.out=M+1)
nknots <- length(knots)
nbasis <- nknots +norder-2
beta_basis <- create.bspline.basis(knots,nbasis,norder)

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
    loss <- objective_function(alpha, xi_list, gamma_list, xi_mat, X, Z, Y, h, tau)
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

single_simulation <- function() {

  alpha_true <- rep(2, p)
  gamma_true <- rnorm(K, mean = 0, sd = 1)
  xi_true <- rnorm(K, mean = 0, sd = 1)
  
  Phi <- eval.basis(u_grid, beta_basis)
  Psi <- Phi  
  
  beta0_true_vals <- as.vector(Phi %*% gamma_true)
  g0_true_vals <- as.vector(Psi %*% xi_true)
  
  alpha_init <- rep(0, p)
  xi_list_init <- list(rep(0, K))  
  gamma_list_init <- list(rep(0, K))  
  xi_mat_init <- matrix(0, N, p)
  
  result <- iterated_estimation(
    X, Z, Y, u_grid, alpha_init, xi_list_init, gamma_list_init, xi_mat_init,
    tau = 0.5, h = 0.05, max_iter = 50
  )
  
  beta0_est_vals <- as.vector(Phi %*% result$gamma_list[[1]])
  g0_est_vals <- as.vector(Psi %*% result$xi_list[[1]])
  
  rmse_alpha <- sqrt(mean((result$alpha - alpha_true)^2))
  rmse_beta0_pointwise <- sqrt((beta0_est_vals - beta0_true_vals)^2) 
  rmse_g0_pointwise <- sqrt((g0_est_vals - g0_true_vals)^2)           
  
  return(list(
    rmse_alpha = rmse_alpha,
    rmse_beta0_pointwise = rmse_beta0_pointwise,
    rmse_g0_pointwise = rmse_g0_pointwise
  ))
}


T_len <- length(u_grid)
rmse_alpha_all <- numeric(sim.num)
rmse_beta0_all <- matrix(NA, nrow = sim.num, ncol = T_len)
rmse_g0_all <- matrix(NA, nrow = sim.num, ncol = T_len)

for (s in 1:sim.num) {
  sim_result <- single_simulation()
  rmse_alpha_all[s] <- sim_result$rmse_alpha
  rmse_beta0_all[s, ] <- sim_result$rmse_beta0_pointwise
  rmse_g0_all[s, ] <- sim_result$rmse_g0_pointwise
}

mean_rmse_alpha <- mean(rmse_alpha_all)
mean_rmse_beta0_u <- colMeans(rmse_beta0_all)
mean_rmse_g0_u <- colMeans(rmse_g0_all)
