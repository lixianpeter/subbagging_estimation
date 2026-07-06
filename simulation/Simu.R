############################################################
## Subbagging simulation based on the paper setting
## Uses:
##   1. simple subbagging average
##   2. bias correction by subtracting estimated bias
##   3. bias correction by solving adjusted equation
############################################################

############################################################
## 1. Core subbagging function
############################################################

SubbaggingCore <- function(data, k_N, m_N, y_name = "y",
                           model = c("Linear", "Logistic"),
                           return_subsamples = FALSE) {
  
  model <- match.arg(model)
  
  y <- data[[y_name]]
  x <- as.matrix(data[, setdiff(names(data), y_name), drop = FALSE])
  
  N <- nrow(x)
  p <- ncol(x)
  
  if (k_N > N) {
    stop("k_N cannot be larger than N.")
  }
  
  ############################################################
  ## 1. Solve original estimating equation
  ############################################################
  
  fit_original <- function(x_sub, y_sub) {
    
    if (model == "Linear") {
      beta_hat <- qr.solve(crossprod(x_sub), crossprod(x_sub, y_sub))
      beta_hat <- as.numeric(beta_hat)
    }
    
    if (model == "Logistic") {
      fit <- suppressWarnings(
        glm.fit(
          x = x_sub,
          y = y_sub,
          family = binomial(link = "logit")
        )
      )
      beta_hat <- as.numeric(coef(fit))
      
      if (any(!is.finite(beta_hat))) {
        stop("Logistic glm failed. Possible separation or singular design.")
      }
    }
    
    return(beta_hat)
  }
  
  ############################################################
  ## 2. Compute psi, V, and estimated bias B_hat
  ############################################################
  
  get_quantities <- function(x_sub, y_sub, beta) {
    
    k <- nrow(x_sub)
    
    if (model == "Linear") {
      
      residual <- as.numeric(y_sub - x_sub %*% beta)
      psi <- x_sub * residual
      
      V <- -crossprod(x_sub) / k
      
      second_part <- rep(0, p)
      
      A <- t(qr.solve(V, t(psi)))   # k by p
      
      xa <- rowSums(x_sub * A)
      D_a <- -x_sub * xa
      
      first_part <- -colMeans(D_a - psi)
    }
    
    if (model == "Logistic") {
      
      eta <- as.numeric(x_sub %*% beta)
      prob <- plogis(eta)
      w <- prob * (1 - prob)
      
      psi <- x_sub * as.numeric(y_sub - prob)
      
      V <- -crossprod(x_sub, x_sub * w) / k
      
      A <- t(qr.solve(V, t(psi)))   # k by p
      
      xa <- rowSums(x_sub * A)
      D_a <- -x_sub * as.numeric(w * xa)
      
      first_part <- -colMeans(D_a - psi)
      
      second_part <- -0.5 * colMeans(
        x_sub * as.numeric(w * (1 - 2 * prob) * xa^2)
      )
    }
    
    B_hat <- -qr.solve(V, first_part + second_part)
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
    
    k <- nrow(x_sub)
    q <- get_quantities(x_sub, y_sub, beta_hat)
    
    beta_bc <- beta_hat - q$B_hat / k
    
    return(beta_bc)
  }
  
  ############################################################
  ## 4. Bias correction by solving adjusted equation
  ############################################################
  
  bias_correct_equation <- function(x_sub, y_sub, beta_start,
                                    max_iter = 20, tol = 1e-8) {
    
    beta <- beta_start
    k <- nrow(x_sub)
    
    for (iter in 1:max_iter) {
      
      q <- get_quantities(x_sub, y_sub, beta)
      
      mean_psi <- colMeans(q$psi)
      
      adjusted_score <- mean_psi + as.numeric(q$V %*% q$B_hat) / k
      
      step <- qr.solve(q$V, adjusted_score)
      
      beta_new <- beta - as.numeric(step)
      
      if (max(abs(beta_new - beta)) < tol) {
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
  
  sum_simple <- rep(0, p)
  sum_bc_add <- rep(0, p)
  sum_bc_equation <- rep(0, p)
  
  if (return_subsamples) {
    beta_simple_mat <- matrix(NA_real_, nrow = m_N, ncol = p)
    beta_bc_add_mat <- matrix(NA_real_, nrow = m_N, ncol = p)
    beta_bc_equation_mat <- matrix(NA_real_, nrow = m_N, ncol = p)
  }
  
  for (s in 1:m_N) {
    
    subsample_id <- sample.int(N, size = k_N, replace = FALSE)
    
    x_sub <- x[subsample_id, , drop = FALSE]
    y_sub <- y[subsample_id]
    
    beta_hat <- fit_original(x_sub, y_sub)
    
    beta_bc_add <- bias_correct_add(x_sub, y_sub, beta_hat)
    
    beta_bc_equation <- bias_correct_equation(x_sub, y_sub, beta_bc_add)
    
    sum_simple <- sum_simple + beta_hat
    sum_bc_add <- sum_bc_add + beta_bc_add
    sum_bc_equation <- sum_bc_equation + beta_bc_equation
    
    if (return_subsamples) {
      beta_simple_mat[s, ] <- beta_hat
      beta_bc_add_mat[s, ] <- beta_bc_add
      beta_bc_equation_mat[s, ] <- beta_bc_equation
    }
  }
  
  out <- list(
    simple_average = sum_simple / m_N,
    bias_correction_add = sum_bc_add / m_N,
    bias_correction_equation = sum_bc_equation / m_N
  )
  
  if (return_subsamples) {
    out$subsample_simple <- beta_simple_mat
    out$subsample_bc_add <- beta_bc_add_mat
    out$subsample_bc_equation <- beta_bc_equation_mat
  }
  
  return(out)
}


############################################################
## 2. Full-sample estimator
############################################################

FitFullEstimator <- function(data, y_name = "y",
                             model = c("Linear", "Logistic")) {
  
  model <- match.arg(model)
  
  y <- data[[y_name]]
  x <- as.matrix(data[, setdiff(names(data), y_name), drop = FALSE])
  
  if (model == "Linear") {
    beta_hat <- qr.solve(crossprod(x), crossprod(x, y))
  }
  
  if (model == "Logistic") {
    fit <- suppressWarnings(
      glm.fit(
        x = x,
        y = y,
        family = binomial(link = "logit")
      )
    )
    beta_hat <- coef(fit)
  }
  
  beta_hat <- as.numeric(beta_hat)
  
  if (any(!is.finite(beta_hat))) {
    stop("Full-sample estimator failed.")
  }
  
  return(beta_hat)
}


############################################################
## 3. Paper data-generating process
############################################################

GeneratePaperData <- function(N,
                              model = c("Logistic", "Linear"),
                              theta0 = c(0, 1),
                              sigma = 1) {
  
  model <- match.arg(model)
  
  x1 <- rnorm(N)
  X <- cbind(intercept = 1, x1 = x1)
  
  eta <- as.numeric(X %*% theta0)
  
  if (model == "Logistic") {
    y <- rbinom(N, size = 1, prob = plogis(eta))
  }
  
  if (model == "Linear") {
    y <- eta + rnorm(N, mean = 0, sd = sigma)
  }
  
  data <- data.frame(
    y = y,
    intercept = X[, 1],
    x1 = X[, 2]
  )
  
  return(data)
}


############################################################
## 4. One simulation replication
############################################################

OneReplication <- function(rep_id,
                           N,
                           k_exponent,
                           alpha,
                           m_rule = c("alphaN", "alphaN43"),
                           model = c("Logistic", "Linear"),
                           theta0 = c(0, 1),
                           seed = 12345,
                           max_m_N = Inf) {
  
  model <- match.arg(model)
  m_rule <- match.arg(m_rule)
  
  set.seed(seed + rep_id)
  
  data <- GeneratePaperData(
    N = N,
    model = model,
    theta0 = theta0
  )
  
  k_N <- floor(N^k_exponent)
  
  if (m_rule == "alphaN") {
    m_N_exact <- floor(alpha * N / k_N)
  }
  
  if (m_rule == "alphaN43") {
    m_N_exact <- floor(alpha * N^(4 / 3) / k_N)
  }
  
  m_N <- m_N_exact
  
  if (is.finite(max_m_N)) {
    m_N <- min(m_N, max_m_N)
  }
  
  m_N <- as.integer(max(1, m_N))
  
  time_start <- Sys.time()
  
  beta_full <- FitFullEstimator(
    data = data,
    model = model
  )
  
  subbagging_result <- SubbaggingCore(
    data = data,
    k_N = k_N,
    m_N = m_N,
    y_name = "y",
    model = model,
    return_subsamples = FALSE
  )
  
  time_used <- as.numeric(Sys.time() - time_start, units = "secs")
  
  est_list <- list(
    full_sample = beta_full,
    simple_average = subbagging_result$simple_average,
    bias_correction_add = subbagging_result$bias_correction_add,
    bias_correction_equation = subbagging_result$bias_correction_equation
  )
  
  out <- do.call(
    rbind,
    lapply(names(est_list), function(method_name) {
      est <- as.numeric(est_list[[method_name]])
      data.frame(
        rep = rep_id,
        model = model,
        N = N,
        alpha = alpha,
        k_exponent = k_exponent,
        k_N = k_N,
        m_rule = m_rule,
        m_N_exact = m_N_exact,
        m_N_used = m_N,
        method = method_name,
        parameter = paste0("beta", seq_along(theta0) - 1),
        estimate = est,
        truth = theta0,
        time = time_used,
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(out) <- NULL
  return(out)
}


############################################################
## 5. Summarise simulation results
############################################################

SummariseSimulation <- function(raw_result) {
  
  split_id <- interaction(
    raw_result$method,
    raw_result$parameter,
    drop = TRUE
  )
  
  out <- do.call(
    rbind,
    lapply(split(raw_result, split_id), function(d) {
      
      bias <- mean(d$estimate - d$truth)
      
      ## Paper-style SD: denominator R, not R - 1
      sd_paper <- sqrt(mean((d$estimate - mean(d$estimate))^2))
      
      rmse <- sqrt(bias^2 + sd_paper^2)
      
      data.frame(
        model = d$model[1],
        N = d$N[1],
        alpha = d$alpha[1],
        k_exponent = d$k_exponent[1],
        k_N = d$k_N[1],
        m_rule = d$m_rule[1],
        m_N_exact = d$m_N_exact[1],
        m_N_used = d$m_N_used[1],
        method = d$method[1],
        parameter = d$parameter[1],
        truth = d$truth[1],
        mean_estimate = mean(d$estimate),
        BIAS = bias,
        SD = sd_paper,
        RMSE = rmse,
        avg_time = mean(d$time),
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(out) <- NULL
  return(out)
}


############################################################
## 6. Run one setting
############################################################

RunSubbaggingSimulation <- function(R = 1000,
                                    N = 20000,
                                    alpha = 1,
                                    k_exponent = 5 / 12,
                                    m_rule = c("alphaN", "alphaN43"),
                                    model = c("Logistic", "Linear"),
                                    theta0 = c(0, 1),
                                    seed = 12345,
                                    max_m_N = Inf,
                                    verbose = TRUE) {
  
  model <- match.arg(model)
  m_rule <- match.arg(m_rule)
  
  raw_list <- vector("list", R)
  
  for (r in 1:R) {
    
    if (verbose && (r == 1 || r %% max(1, floor(R / 10)) == 0)) {
      message("Running replication ", r, " / ", R)
    }
    
    raw_list[[r]] <- OneReplication(
      rep_id = r,
      N = N,
      k_exponent = k_exponent,
      alpha = alpha,
      m_rule = m_rule,
      model = model,
      theta0 = theta0,
      seed = seed,
      max_m_N = max_m_N
    )
  }
  
  raw_result <- do.call(rbind, raw_list)
  summary_result <- SummariseSimulation(raw_result)
  
  return(list(
    raw = raw_result,
    summary = summary_result
  ))
}


############################################################
## 7. Paper simulation grid
############################################################

MakePaperGrid <- function() {
  
  grid <- expand.grid(
    N = c(20000, 100000, 500000),
    alpha = c(1, 1 / 3),
    k_exponent = c(5 / 12, 6 / 12, 7 / 12, 8 / 12),
    m_rule = c("alphaN", "alphaN43"),
    stringsAsFactors = FALSE
  )
  
  return(grid)
}


RunPaperGrid <- function(R = 1000,
                         grid = MakePaperGrid(),
                         model = "Logistic",
                         theta0 = c(0, 1),
                         seed = 12345,
                         max_m_N = Inf,
                         output_dir = "subbagging_results") {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  all_summary <- list()
  
  for (g in 1:nrow(grid)) {
    
    N_g <- grid$N[g]
    alpha_g <- grid$alpha[g]
    k_exp_g <- grid$k_exponent[g]
    m_rule_g <- grid$m_rule[g]
    
    message("============================================")
    message("Grid ", g, " / ", nrow(grid))
    message("N = ", N_g,
            ", alpha = ", alpha_g,
            ", k_exponent = ", k_exp_g,
            ", m_rule = ", m_rule_g)
    message("============================================")
    
    res <- RunSubbaggingSimulation(
      R = R,
      N = N_g,
      alpha = alpha_g,
      k_exponent = k_exp_g,
      m_rule = m_rule_g,
      model = model,
      theta0 = theta0,
      seed = seed + 100000 * g,
      max_m_N = max_m_N,
      verbose = TRUE
    )
    
    tag <- paste0(
      model,
      "_N", N_g,
      "_alpha", signif(alpha_g, 4),
      "_kexp", gsub("\\.", "_", as.character(round(k_exp_g, 6))),
      "_", m_rule_g
    )
    
    write.csv(
      res$raw,
      file = file.path(output_dir, paste0("raw_", tag, ".csv")),
      row.names = FALSE
    )
    
    write.csv(
      res$summary,
      file = file.path(output_dir, paste0("summary_", tag, ".csv")),
      row.names = FALSE
    )
    
    all_summary[[g]] <- res$summary
  }
  
  all_summary <- do.call(rbind, all_summary)
  
  write.csv(
    all_summary,
    file = file.path(output_dir, "all_summary.csv"),
    row.names = FALSE
  )
  
  return(all_summary)
}


############################################################
## 8. Example runs
############################################################

## Quick test first
## This is NOT the full paper setting. It is just to check the code runs.

test_res <- RunSubbaggingSimulation(
  R = 20,
  N = 20000,
  alpha = 1,
  k_exponent = 5 / 12,
  m_rule = "alphaN",
  model = "Logistic",
  theta0 = c(0, 1),
  max_m_N = Inf
)

print(test_res$summary)


## Full paper-style grid
## Warning: R = 1000 and m_rule = "alphaN43" can be very slow.
## Run this only when you are ready.

# paper_summary <- RunPaperGrid(
#   R = 1000,
#   model = "Logistic",
#   theta0 = c(0, 1),
#   max_m_N = Inf,
#   output_dir = "subbagging_results"
# )
