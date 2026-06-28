SubbaggingCore <- function(data, k_N, m_N, y_name = "y", model = "Linear") {
  
  ############################################################
  ## data:  data frame containing response y and covariates x
  ## k_N:   subsample size
  ## m_N:   number of subsamples
  ## model: "Linear" or "Logistic"
  ############################################################
  
  y <- data[[y_name]]
  x <- as.matrix(data[, setdiff(names(data), y_name)])
  
  N <- nrow(x)
  p <- ncol(x)
  
  ############################################################
  ## 1. Solve original estimating equation
  ############################################################
  
  fit_original <- function(x_sub, y_sub) {
    
    if (model == "Linear") {
      beta_hat <- solve(t(x_sub) %*% x_sub) %*% t(x_sub) %*% y_sub
      beta_hat <- as.numeric(beta_hat)
    }
    
    if (model == "Logistic") {
      fit <- glm(y_sub ~ x_sub - 1, family = binomial(link = "logit"))
      beta_hat <- as.numeric(coef(fit))
    }
    
    return(beta_hat)
  }
  
  ############################################################
  ## 2. Compute psi, V, and estimated bias B_hat
  ############################################################
  
  get_quantities <- function(x_sub, y_sub, beta) {
    
    k <- nrow(x_sub)
    
    if (model == "Linear") {
      
      ## psi_i(beta) = x_i (y_i - x_i^T beta)
      residual <- as.numeric(y_sub - x_sub %*% beta)
      psi <- x_sub * residual
      
      ## V = average derivative of psi
      V <- -t(x_sub) %*% x_sub / k
      
      ## For linear regression, second derivative is zero
      second_part <- rep(0, p)
      
      ## a_i = V^{-1} psi_i
      A <- t(solve(V, t(psi)))   # k by p
      
      ## D_i a_i = - x_i x_i^T a_i
      xa <- rowSums(x_sub * A)
      D_a <- -x_sub * xa
      
      first_part <- -colMeans(D_a - psi)
    }
    
    if (model == "Logistic") {
      
      eta <- as.numeric(x_sub %*% beta)
      prob <- 1 / (1 + exp(-eta))
      w <- prob * (1 - prob)
      
      ## psi_i(beta) = x_i (y_i - p_i)
      psi <- x_sub * as.numeric(y_sub - prob)
      
      ## V = average derivative of psi
      V <- -t(x_sub) %*% (x_sub * w) / k
      
      ## a_i = V^{-1} psi_i
      A <- t(solve(V, t(psi)))   # k by p
      
      ## D_i a_i = - w_i x_i x_i^T a_i
      xa <- rowSums(x_sub * A)
      D_a <- -x_sub * as.numeric(w * xa)
      
      first_part <- -colMeans(D_a - psi)
      
      ## second derivative part
      second_part <- -0.5 * colMeans(
        x_sub * as.numeric(w * (1 - 2 * prob) * xa^2)
      )
    }
    
    ## Estimated bias term B_hat
    B_hat <- -solve(V, first_part + second_part)
    B_hat <- as.numeric(B_hat)
    
    return(list(
      psi = psi,
      V = V,
      B_hat = B_hat
    ))
  }
  
  ############################################################
  ## 3. Bias correction by adding correction term
  ############################################################
  
  bias_correct_add <- function(x_sub, y_sub, beta_hat) {
    
    q <- get_quantities(x_sub, y_sub, beta_hat)
    
    ## Paper form:
    ## beta_bc = beta_hat - B_hat / k_N
    beta_bc <- beta_hat - q$B_hat / k_N
    
    return(beta_bc)
  }
  
  ############################################################
  ## 4. Bias correction by solving adjusted equation
  ############################################################
  
  bias_correct_equation <- function(x_sub, y_sub, beta_start, max_iter = 20) {
    
    beta <- beta_start
    
    for (iter in 1:max_iter) {
      
      q <- get_quantities(x_sub, y_sub, beta)
      
      mean_psi <- colMeans(q$psi)
      
      ## Adjusted equation:
      ## sum psi + V B = 0
      ##
      ## Divide by k_N:
      ## mean psi + V B / k_N = 0
      adjusted_score <- mean_psi + as.numeric(q$V %*% q$B_hat) / k_N
      
      step <- solve(q$V, adjusted_score)
      
      beta_new <- beta - as.numeric(step)
      
      if (max(abs(beta_new - beta)) < 1e-8) {
        beta <- beta_new
        break
      }
      
      beta <- beta_new
    }
    
    return(beta)
  }
  
  ############################################################
  ## Main subbagging loop
  ############################################################
  
  beta_simple_list <- list()
  beta_bc_add_list <- list()
  beta_bc_equation_list <- list()
  
  for (s in 1:m_N) {
    
    subsample_id <- sample(1:N, size = k_N, replace = FALSE)
    
    x_sub <- x[subsample_id, , drop = FALSE]
    y_sub <- y[subsample_id]
    
    ## 1. Original subsample estimator
    beta_hat <- fit_original(x_sub, y_sub)
    
    ## 2. Bias correction by adding correction term
    beta_bc_add <- bias_correct_add(x_sub, y_sub, beta_hat)
    
    ## 3. Bias correction by solving adjusted equation
    beta_bc_equation <- bias_correct_equation(x_sub, y_sub, beta_bc_add)
    
    beta_simple_list[[s]] <- beta_hat
    beta_bc_add_list[[s]] <- beta_bc_add
    beta_bc_equation_list[[s]] <- beta_bc_equation
  }
  
  ############################################################
  ## Aggregate over m_N subsamples
  ############################################################
  
  beta_simple_average <- Reduce("+", beta_simple_list) / m_N
  beta_bc_add_average <- Reduce("+", beta_bc_add_list) / m_N
  beta_bc_equation_average <- Reduce("+", beta_bc_equation_list) / m_N
  
  return(list(
    simple_average = beta_simple_average,
    bias_correction_add = beta_bc_add_average,
    bias_correction_equation = beta_bc_equation_average,
    
    subsample_simple = beta_simple_list,
    subsample_bc_add = beta_bc_add_list,
    subsample_bc_equation = beta_bc_equation_list
  ))
}
