############################################################
## Safe Subbagging Simulation for Logistic Regression
## Based on the paper setting
##
## Methods:
##   1. full_sample
##   2. simple_average
##   3. bias_correction_add      = bc2 in paper style
##   4. bias_correction_equation = bc3 in paper style
############################################################

rm(list = ls())

############################################################
## 0. Safe matrix solver
############################################################

safe_solve <- function(A, B, tol = 1e-10) {
  
  if (any(!is.finite(A)) || any(!is.finite(B))) {
    stop("Non-finite value in matrix solve.")
  }
  
  if (nrow(A) != ncol(A)) {
    stop("Matrix A is not square.")
  }
  
  qrA <- qr(A, tol = tol)
  
  if (qrA$rank < ncol(A)) {
    stop("Singular or rank-deficient matrix.")
  }
  
  out <- tryCatch(
    solve(A, B),
    error = function(e) {
      stop("Matrix solve failed.")
    }
  )
  
  return(out)
}


############################################################
## 1. Core subbagging function
############################################################

SubbaggingCore <- function(data,
                           k_N,
                           m_N,
                           y_name = "y",
                           model = c("Linear", "Logistic"),
                           max_bad_draws = 10000,
                           verbose = TRUE) {
  
  model <- match.arg(model)
  
  y <- data[[y_name]]
  x <- as.matrix(data[, setdiff(names(data), y_name), drop = FALSE])
  
  N <- nrow(x)
  p <- ncol(x)
  
  if (k_N > N) stop("k_N cannot be larger than N.")
  if (k_N <= p) stop("k_N must be larger than p.")
  
  ############################################################
  ## 1. Original estimator on one subsample
  ############################################################
  
  fit_original <- function(x_sub, y_sub) {
    
    if (qr(x_sub)$rank < ncol(x_sub)) {
      stop("x_sub is rank deficient.")
    }
    
    if (model == "Linear") {
      beta_hat <- safe_solve(crossprod(x_sub), crossprod(x_sub, y_sub))
      beta_hat <- as.numeric(beta_hat)
    }
    
    if (model == "Logistic") {
      fit <- suppressWarnings(
        glm.fit(
          x = x_sub,
          y = y_sub,
          family = binomial(link = "logit"),
          control = glm.control(maxit = 50)
        )
      )
      
      beta_hat <- as.numeric(coef(fit))
      
      if (any(!is.finite(beta_hat))) {
        stop("Logistic glm produced non-finite coefficients.")
      }
    }
    
    return(beta_hat)
  }
  
  ############################################################
  ## 2. Compute psi, V, and B_hat
  ############################################################
  
  get_quantities <- function(x_sub, y_sub, beta) {
    
    k <- nrow(x_sub)
    
    if (model == "Linear") {
      
      residual <- as.numeric(y_sub - x_sub %*% beta)
      
      ## psi_i(beta) = x_i (y_i - x_i^T beta)
      psi <- x_sub * residual
      
      ## V = average derivative of psi
      V <- -crossprod(x_sub) / k
      
      ## A_i = V^{-1} psi_i
      A <- t(safe_solve(V, t(psi)))
      
      ## D_i A_i = -x_i x_i^T A_i
      xa <- rowSums(x_sub * A)
      D_a <- -x_sub * xa
      
      first_part <- -colMeans(D_a - psi)
      
      ## Linear regression second derivative is zero
      second_part <- rep(0, ncol(x_sub))
    }
    
    if (model == "Logistic") {
      
      eta <- as.numeric(x_sub %*% beta)
      prob <- plogis(eta)
      w <- prob * (1 - prob)
      
      ## If all weights are basically zero, V will be singular
      if (min(w) < 1e-14 && mean(w) < 1e-8) {
        stop("Logistic weights are too small. Possible separation.")
      }
      
      ## psi_i(beta) = x_i (y_i - p_i)
      psi <- x_sub * as.numeric(y_sub - prob)
      
      ## V = average derivative of psi
      V <- -crossprod(x_sub, x_sub * w) / k
      
      ## A_i = V^{-1} psi_i
      A <- t(safe_solve(V, t(psi)))
      
      ## D_i A_i = -w_i x_i x_i^T A_i
      xa <- rowSums(x_sub * A)
      D_a <- -x_sub * as.numeric(w * xa)
      
      first_part <- -colMeans(D_a - psi)
      
      ## second derivative part
      second_part <- -0.5 * colMeans(
        x_sub * as.numeric(w * (1 - 2 * prob) * xa^2)
      )
    }
    
    ## B_hat
    B_hat <- -safe_solve(V, first_part + second_part)
    B_hat <- as.numeric(B_hat)
    
    return(list(
      psi = psi,
      V = V,
      B_hat = B_hat
    ))
  }
  
  ############################################################
  ## 3. Bias correction by subtracting estimated bias
  ############################################################
  
  bias_correct_add <- function(x_sub, y_sub, beta_hat) {
    
    k <- nrow(x_sub)
    q <- get_quantities(x_sub, y_sub, beta_hat)
    
    beta_bc <- beta_hat - q$B_hat / k
    
    if (any(!is.finite(beta_bc))) {
      stop("Non-finite beta_bc.")
    }
    
    return(beta_bc)
  }
  
  ############################################################
  ## 4. Bias correction by solving adjusted equation
  ############################################################
  
  bias_correct_equation <- function(x_sub,
                                    y_sub,
                                    beta_start,
                                    max_iter = 20,
                                    tol = 1e-8) {
    
    beta <- beta_start
    k <- nrow(x_sub)
    
    for (iter in 1:max_iter) {
      
      q <- get_quantities(x_sub, y_sub, beta)
      
      mean_psi <- colMeans(q$psi)
      
      ## Adjusted equation:
      ## sum psi + V B = 0
      ## divide by k:
      ## mean psi + V B / k = 0
      adjusted_score <- mean_psi + as.numeric(q$V %*% q$B_hat) / k
      
      step <- safe_solve(q$V, adjusted_score)
      step <- as.numeric(step)
      
      beta_new <- beta - step
      
      if (any(!is.finite(beta_new))) {
        stop("Non-finite beta in adjusted equation.")
      }
      
      if (max(abs(beta_new - beta)) < tol) {
        beta <- beta_new
        break
      }
      
      beta <- beta_new
    }
    
    return(beta)
  }
  
  ############################################################
  ## 5. Main subbagging loop
  ## Bad subsamples are redrawn
  ############################################################
  
  sum_simple <- rep(0, p)
  sum_bc_add <- rep(0, p)
  sum_bc_equation <- rep(0, p)
  
  accepted <- 0
  bad_draws <- 0
  total_draws <- 0
  
  while (accepted < m_N) {
    
    total_draws <- total_draws + 1
    
    one_draw <- tryCatch({
      
      subsample_id <- sample.int(N, size = k_N, replace = FALSE)
      
      x_sub <- x[subsample_id, , drop = FALSE]
      y_sub <- y[subsample_id]
      
      beta_hat <- fit_original(x_sub, y_sub)
      
      beta_bc_add <- bias_correct_add(x_sub, y_sub, beta_hat)
      
      beta_bc_equation <- bias_correct_equation(
        x_sub = x_sub,
        y_sub = y_sub,
        beta_start = beta_bc_add
      )
      
      list(
        beta_hat = beta_hat,
        beta_bc_add = beta_bc_add,
        beta_bc_equation = beta_bc_equation
      )
      
    }, error = function(e) {
      NULL
    })
    
    if (is.null(one_draw)) {
      bad_draws <- bad_draws + 1
      
      if (bad_draws > max_bad_draws) {
        stop(
          paste0(
            "Too many bad subsamples. ",
            "Check k_N, p, collinearity, or logistic separation. ",
            "bad_draws = ", bad_draws
          )
        )
      }
      
      next
    }
    
    accepted <- accepted + 1
    
    sum_simple <- sum_simple + one_draw$beta_hat
    sum_bc_add <- sum_bc_add + one_draw$beta_bc_add
    sum_bc_equation <- sum_bc_equation + one_draw$beta_bc_equation
    
    if (verbose && accepted %% max(1, floor(m_N / 5)) == 0) {
      message("Accepted subsamples: ", accepted, " / ", m_N)
    }
  }
  
  return(list(
    simple_average = sum_simple / m_N,
    bias_correction_add = sum_bc_add / m_N,
    bias_correction_equation = sum_bc_equation / m_N,
    accepted_subsamples = accepted,
    bad_draws = bad_draws,
    total_draws = total_draws
  ))
}


############################################################
## 2. Full sample estimator
############################################################

FitFullEstimator <- function(data,
                             y_name = "y",
                             model = c("Linear", "Logistic")) {
  
  model <- match.arg(model)
  
  y <- data[[y_name]]
  x <- as.matrix(data[, setdiff(names(data), y_name), drop = FALSE])
  
  if (model == "Linear") {
    beta_hat <- safe_solve(crossprod(x), crossprod(x, y))
    beta_hat <- as.numeric(beta_hat)
  }
  
  if (model == "Logistic") {
    fit <- suppressWarnings(
      glm.fit(
        x = x,
        y = y,
        family = binomial(link = "logit"),
        control = glm.control(maxit = 50)
      )
    )
    
    beta_hat <- as.numeric(coef(fit))
  }
  
  if (any(!is.finite(beta_hat))) {
    stop("Full sample estimator failed.")
  }
  
  return(beta_hat)
}


############################################################
## 3. Paper logistic DGP
############################################################

GeneratePaperLogisticData <- function(N, theta0 = c(0, 1)) {
  
  x1 <- rnorm(N)
  
  X <- cbind(
    intercept = 1,
    x1 = x1
  )
  
  eta <- as.numeric(X %*% theta0)
  prob <- plogis(eta)
  
  y <- rbinom(N, size = 1, prob = prob)
  
  data <- data.frame(
    y = y,
    intercept = X[, 1],
    x1 = X[, 2]
  )
  
  return(data)
}


############################################################
## 4. One replication
############################################################

OneReplication <- function(rep_id,
                           N,
                           alpha,
                           k_exponent,
                           m_rule = c("alphaN", "alphaN43"),
                           theta0 = c(0, 1),
                           seed = 12345,
                           max_m_N = Inf,
                           verbose = FALSE) {
  
  m_rule <- match.arg(m_rule)
  
  set.seed(seed + rep_id)
  
  data <- GeneratePaperLogisticData(
    N = N,
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
  
  if (verbose) {
    message("N = ", N)
    message("k_N = ", k_N)
    message("m_N_exact = ", m_N_exact)
    message("m_N_used = ", m_N)
  }
  
  time_start <- Sys.time()
  
  beta_full <- FitFullEstimator(
    data = data,
    model = "Logistic"
  )
  
  subbagging_result <- SubbaggingCore(
    data = data,
    k_N = k_N,
    m_N = m_N,
    y_name = "y",
    model = "Logistic",
    verbose = FALSE
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
        N = N,
        alpha = alpha,
        k_exponent = k_exponent,
        k_N = k_N,
        m_rule = m_rule,
        m_N_exact = m_N_exact,
        m_N_used = m_N,
        method = method_name,
        parameter = paste0("theta", seq_along(theta0)),
        estimate = est,
        truth = theta0,
        bad_draws = subbagging_result$bad_draws,
        total_draws = subbagging_result$total_draws,
        time_seconds = time_used,
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(out) <- NULL
  return(out)
}


############################################################
## 5. Run simulation for one setting
############################################################

RunOneSetting <- function(R = 100,
                          N = 20000,
                          alpha = 1,
                          k_exponent = 5 / 12,
                          m_rule = c("alphaN", "alphaN43"),
                          theta0 = c(0, 1),
                          seed = 12345,
                          max_m_N = Inf,
                          verbose = TRUE) {
  
  m_rule <- match.arg(m_rule)
  
  raw_list <- vector("list", R)
  
  for (r in 1:R) {
    
    if (verbose) {
      message("Replication ", r, " / ", R)
    }
    
    raw_list[[r]] <- OneReplication(
      rep_id = r,
      N = N,
      alpha = alpha,
      k_exponent = k_exponent,
      m_rule = m_rule,
      theta0 = theta0,
      seed = seed,
      max_m_N = max_m_N,
      verbose = FALSE
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
## 6. Summarise results
############################################################

SummariseSimulation <- function(raw_result) {
  
  group_id <- interaction(
    raw_result$N,
    raw_result$alpha,
    raw_result$k_exponent,
    raw_result$k_N,
    raw_result$m_rule,
    raw_result$m_N_used,
    raw_result$method,
    raw_result$parameter,
    drop = TRUE
  )
  
  out <- do.call(
    rbind,
    lapply(split(raw_result, group_id), function(d) {
      
      bias <- mean(d$estimate - d$truth)
      
      ## Paper-style SD uses denominator R, not R - 1
      sd_paper <- sqrt(mean((d$estimate - mean(d$estimate))^2))
      
      rmse <- sqrt(bias^2 + sd_paper^2)
      
      data.frame(
        N = d$N[1],
        alpha = d$alpha[1],
        k_exponent = d$k_exponent[1],
        k_N = d$k_N[1],
        m_rule = d$m_rule[1],
        m_N_used = d$m_N_used[1],
        method = d$method[1],
        parameter = d$parameter[1],
        truth = d$truth[1],
        mean_estimate = mean(d$estimate),
        BIAS = bias,
        SD = sd_paper,
        RMSE = rmse,
        BIAS_x100 = 100 * bias,
        SD_x100 = 100 * sd_paper,
        RMSE_x100 = 100 * rmse,
        avg_bad_draws = mean(d$bad_draws),
        avg_total_draws = mean(d$total_draws),
        avg_time_seconds = mean(d$time_seconds),
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(out) <- NULL
  return(out)
}


############################################################
## 7. Paper grid
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
    message("N = ", N_g)
    message("alpha = ", alpha_g)
    message("k_exponent = ", k_exp_g)
    message("m_rule = ", m_rule_g)
    message("============================================")
    
    res <- RunOneSetting(
      R = R,
      N = N_g,
      alpha = alpha_g,
      k_exponent = k_exp_g,
      m_rule = m_rule_g,
      theta0 = theta0,
      seed = seed + 100000 * g,
      max_m_N = max_m_N,
      verbose = TRUE
    )
    
    tag <- paste0(
      "N", N_g,
      "_alpha", gsub("\\.", "_", as.character(alpha_g)),
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
## 8. Copy-paste runnable test
############################################################

## This quick test should run first.
## It uses the paper DGP but only R = 20 replications.

test_res <- RunOneSetting(
  R = 20,
  N = 20000,
  alpha = 1,
  k_exponent = 5 / 12,
  m_rule = "alphaN",
  theta0 = c(0, 1),
  seed = 12345,
  max_m_N = Inf,
  verbose = TRUE
)

print(test_res$summary)


############################################################
## 9. Full paper-style run
############################################################

## This is expensive.
## Uncomment only after the quick test works.

# paper_summary <- RunPaperGrid(
#   R = 1000,
#   theta0 = c(0, 1),
#   seed = 12345,
#   max_m_N = Inf,
#   output_dir = "subbagging_results"
# )
#
# print(paper_summary)
