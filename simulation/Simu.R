############################################################
## Required package
############################################################

if (!requireNamespace("nleqslv", quietly = TRUE)) {
  install.packages("nleqslv")
}


############################################################
## 1. Subbagging core
############################################################

SubbaggingCore <- function(data, k_N, alpha,
                           model = c("linear", "logistic"),
                           seed = 1) {
  
  model <- match.arg(model)
  set.seed(seed)
  
  y <- as.numeric(data[, 1])
  X <- as.matrix(data[, -1])
  X <- cbind(Intercept = 1, X)
  
  N <- nrow(X)
  d <- ncol(X)
  
  m_N <- floor(alpha * N / k_N)
  alpha_N <- k_N * m_N / N
  
  ############################################################
  ## Z-estimator on one subsample
  ############################################################
  
  fit_theta <- function(Xs, ys) {
    
    if (model == "linear") {
      beta <- lm.fit(Xs, ys)$coefficients
    } else {
      fit <- suppressWarnings(
        glm.fit(Xs, ys, family = binomial())
      )
      
      if (!isTRUE(fit$converged)) {
        stop("logistic fit did not converge")
      }
      
      beta <- fit$coefficients
    }
    
    beta <- as.numeric(beta)
    
    if (any(!is.finite(beta))) {
      stop("bad fit")
    }
    
    beta
  }
  
  ############################################################
  ## Compute psi, V, and B_hat
  ############################################################
  
  moments <- function(theta, Xs, ys) {
    
    n <- nrow(Xs)
    eta <- as.numeric(Xs %*% theta)
    
    if (model == "linear") {
      mu <- eta
      w <- rep(1, n)
      h <- rep(0, n)
    } else {
      mu <- plogis(eta)
      w <- mu * (1 - mu)
      h <- w * (1 - 2 * mu)
    }
    
    psi <- Xs * as.numeric(ys - mu)
    V <- -crossprod(Xs, Xs * w) / n
    
    list(
      psi = psi,
      V = V,
      w = w,
      h = h
    )
  }
  
  Bhat <- function(theta, Xs, ys) {
    
    n <- nrow(Xs)
    d <- ncol(Xs)
    
    M <- moments(theta, Xs, ys)
    V <- M$V
    
    if (!is.finite(rcond(V)) || rcond(V) < 1e-12) {
      stop("V is singular or numerically ill-conditioned")
    }
    
    Vinv <- solve(V)
    
    psi <- M$psi
    A <- psi %*% t(Vinv)
    
    XA <- rowSums(Xs * A)
    
    D_part <- colMeans(Xs * (-M$w * XA))
    V_part <- as.numeric(V %*% colMeans(A))
    
    term1 <- -(D_part - V_part)
    
    Q <- crossprod(A) / n
    term2 <- rep(0, d)
    
    if (model == "logistic") {
      for (j in seq_len(d)) {
        Hj <- -crossprod(
          Xs,
          Xs * (M$h * Xs[, j])
        ) / n
        
        term2[j] <- 0.5 * sum(Hj * Q)
      }
    }
    
    B <- -as.numeric(Vinv %*% (term1 + term2))
    
    if (any(!is.finite(B))) {
      stop("bad Bhat")
    }
    
    B
  }
  
  ############################################################
  ## Adjusted estimating equation
  ############################################################
  
  adjusted_score <- function(theta, Xs, ys) {
    
    M <- moments(theta, Xs, ys)
    score <- colSums(M$psi)
    B <- Bhat(theta, Xs, ys)
    
    adjusted_value <- as.numeric(score + M$V %*% B)
    
    if (any(!is.finite(adjusted_value))) {
      stop("non-finite adjusted score")
    }
    
    adjusted_value
  }
  
  ############################################################
  ## Solve adjusted_score(theta) = 0 directly using nleqslv
  ############################################################
  
  solve_bc_equation <- function(theta_start, Xs, ys) {
    
    result <- nleqslv::nleqslv(
      x = theta_start,
      fn = function(theta) {
        adjusted_score(theta, Xs, ys)
      },
      method = "Broyden",
      global = "dbldog",
      control = list(
        ftol = 1e-8,
        xtol = 1e-8,
        maxit = 200,
        allowSingular = FALSE
      )
    )
    
    theta_bc <- as.numeric(result$x)
    
    if (any(!is.finite(theta_bc))) {
      stop("non-finite adjusted-equation estimate")
    }
    
    final_score <- adjusted_score(theta_bc, Xs, ys)
    
    if (!(result$termcd %in% c(1, 2))) {
      stop(
        paste0(
          "adjusted equation solver did not converge; termcd = ",
          result$termcd
        )
      )
    }
    
    if (max(abs(final_score)) > 1e-6) {
      stop(
        paste0(
          "adjusted equation residual is too large: ",
          max(abs(final_score))
        )
      )
    }
    
    theta_bc
  }
  
  ############################################################
  ## Storage
  ############################################################
  
  theta_simple <- matrix(0, m_N, d)
  theta_bc_bias <- matrix(0, m_N, d)
  theta_bc_equation <- matrix(0, m_N, d)
  
  bad_draws <- 0
  total_draws <- 0
  
  ############################################################
  ## Timing accumulators, in seconds
  ############################################################
  
  time_subsampling_seconds <- 0
  time_simple_seconds <- 0
  time_bc_bias_seconds <- 0
  time_bc_equation_seconds <- 0
  
  ############################################################
  ## Draw and process subsamples
  ############################################################
  
  for (b in seq_len(m_N)) {
    
    valid_draw <- FALSE
    
    while (!valid_draw) {
      
      total_draws <- total_draws + 1
      
      ##########################################################
      ## Time subsample drawing and extraction
      ##########################################################
      
      time_start <- proc.time()[["elapsed"]]
      
      id <- sample(seq_len(N), k_N, replace = FALSE)
      Xs <- X[id, , drop = FALSE]
      ys <- y[id]
      
      time_subsampling_seconds <- time_subsampling_seconds +
        (proc.time()[["elapsed"]] - time_start)
      
      ##########################################################
      ## Time simple subsample estimator
      ##########################################################
      
      time_start <- proc.time()[["elapsed"]]
      
      theta_result <- try(
        fit_theta(Xs, ys),
        silent = TRUE
      )
      
      time_simple_seconds <- time_simple_seconds +
        (proc.time()[["elapsed"]] - time_start)
      
      if (inherits(theta_result, "try-error")) {
        bad_draws <- bad_draws + 1
        next
      }
      
      theta_hat <- theta_result
      
      ##########################################################
      ## Time bias-correction calculation
      ##########################################################
      
      time_start <- proc.time()[["elapsed"]]
      
      B_result <- try(
        Bhat(theta_hat, Xs, ys),
        silent = TRUE
      )
      
      time_bc_bias_seconds <- time_bc_bias_seconds +
        (proc.time()[["elapsed"]] - time_start)
      
      if (inherits(B_result, "try-error")) {
        bad_draws <- bad_draws + 1
        next
      }
      
      theta_bias <- theta_hat - B_result / k_N
      
      ##########################################################
      ## Time adjusted estimating-equation solution
      ##########################################################
      
      time_start <- proc.time()[["elapsed"]]
      
      theta_equation_result <- try(
        solve_bc_equation(
          theta_start = theta_hat,
          Xs = Xs,
          ys = ys
        ),
        silent = TRUE
      )
      
      time_bc_equation_seconds <- time_bc_equation_seconds +
        (proc.time()[["elapsed"]] - time_start)
      
      if (inherits(theta_equation_result, "try-error")) {
        bad_draws <- bad_draws + 1
        next
      }
      
      theta_simple[b, ] <- theta_hat
      theta_bc_bias[b, ] <- theta_bias
      theta_bc_equation[b, ] <- theta_equation_result
      
      valid_draw <- TRUE
    }
  }
  
  colnames(theta_simple) <- colnames(X)
  colnames(theta_bc_bias) <- colnames(X)
  colnames(theta_bc_equation) <- colnames(X)
  
  ############################################################
  ## Include final averaging time for each estimator
  ############################################################
  
  time_start <- proc.time()[["elapsed"]]
  estimate_simple <- colMeans(theta_simple)
  time_simple_seconds <- time_simple_seconds +
    (proc.time()[["elapsed"]] - time_start)
  
  time_start <- proc.time()[["elapsed"]]
  estimate_bc_bias <- colMeans(theta_bc_bias)
  time_bc_bias_seconds <- time_bc_bias_seconds +
    (proc.time()[["elapsed"]] - time_start)
  
  time_start <- proc.time()[["elapsed"]]
  estimate_bc_equation <- colMeans(theta_bc_equation)
  time_bc_equation_seconds <- time_bc_equation_seconds +
    (proc.time()[["elapsed"]] - time_start)
  
  estimate <- rbind(
    simple = estimate_simple,
    bc_bias = estimate_bc_bias,
    bc_equation = estimate_bc_equation
  )
  
  ############################################################
  ## Standard error estimates
  ############################################################
  
  se_fun <- function(A) {
    Omega <- cov(A) * (m_N - 1) / m_N
    sqrt((1 + 1 / alpha_N) * k_N * diag(Omega) / N)
  }
  
  se <- rbind(
    simple = se_fun(theta_simple),
    bc_bias = se_fun(theta_bc_bias),
    bc_equation = se_fun(theta_bc_equation)
  )
  
  colnames(se) <- colnames(X)
  
  list(
    N = N,
    k_N = k_N,
    m_N = m_N,
    alpha = alpha,
    estimate = estimate,
    se = se,
    bad_draws = bad_draws,
    total_draws = total_draws,
    timing = list(
      subsampling_seconds = time_subsampling_seconds,
      simple_seconds = time_simple_seconds,
      bc_bias_seconds = time_bc_bias_seconds,
      bc_equation_seconds = time_bc_equation_seconds
    ),
    subsample_estimates = list(
      simple = theta_simple,
      bc_bias = theta_bc_bias,
      bc_equation = theta_bc_equation
    )
  )
}


############################################################
## 2. Generate simulated data
############################################################

GenerateSimulatedData <- function(N = 1000,
                                  theta = c(0, 1),
                                  model = c("linear", "logistic"),
                                  sigma = 1,
                                  seed = 1) {
  
  model <- match.arg(model)
  set.seed(seed)
  
  p <- length(theta) - 1
  
  X <- matrix(rnorm(N * p), nrow = N, ncol = p)
  colnames(X) <- paste0("x", seq_len(p))
  
  eta <- as.numeric(cbind(1, X) %*% theta)
  
  if (model == "linear") {
    y <- eta + rnorm(N, sd = sigma)
  }
  
  if (model == "logistic") {
    prob <- plogis(eta)
    y <- rbinom(N, size = 1, prob = prob)
  }
  
  data <- data.frame(y = y, X)
  attr(data, "theta_true") <- theta
  
  data
}


############################################################
## 3. Run simulation
############################################################

RunSubbaggingSimulation <- function(R = 100,
                                    seed_start = 1,
                                    N = 1000,
                                    theta = c(0, 1),
                                    k_N = 100,
                                    alpha = 1,
                                    model = c("linear", "logistic"),
                                    sigma = 1,
                                    output_dir = ".") {
  
  model <- match.arg(model)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  m_N <- floor(alpha * N / k_N)
  methods <- c("simple", "bc_bias", "bc_equation")
  
  results <- list()
  
  for (method in methods) {
    results[[method]] <- NULL
  }
  
  ############################################################
  ## Helper: avoid scientific notation in file names
  ############################################################
  
  plain_number <- function(x) {
    x <- format(x, scientific = FALSE, trim = TRUE)
    x <- gsub("\\.", "p", x)
    x
  }
  
  ############################################################
  ## Run replications
  ############################################################
  
  for (r in seq_len(R)) {
    
    seed_r <- seed_start + r - 1
    
    data_r <- GenerateSimulatedData(
      N = N,
      theta = theta,
      model = model,
      sigma = sigma,
      seed = seed_r
    )
    
    fit_r <- SubbaggingCore(
      data = data_r,
      k_N = k_N,
      alpha = alpha,
      model = model,
      seed = seed_r
    )
    
    coef_names <- colnames(fit_r$estimate)
    
    for (method in methods) {
      
      estimate_r <- fit_r$estimate[method, ]
      se_r <- fit_r$se[method, ]
      
      row_r <- data.frame(
        replication = r,
        seed = seed_r,
        N = N,
        k_N = k_N,
        m_N = fit_r$m_N,
        alpha = alpha,
        model = model,
        method = method,
        bad_draws = fit_r$bad_draws,
        total_draws = fit_r$total_draws,
        time_subsampling_seconds = fit_r$timing$subsampling_seconds,
        time_simple_seconds = fit_r$timing$simple_seconds,
        time_bc_bias_seconds = fit_r$timing$bc_bias_seconds,
        time_bc_equation_seconds = fit_r$timing$bc_equation_seconds
      )
      
      for (j in seq_along(coef_names)) {
        row_r[[paste0("estimate_", coef_names[j])]] <- estimate_r[j]
      }
      
      for (j in seq_along(coef_names)) {
        row_r[[paste0("se_", coef_names[j])]] <- se_r[j]
      }
      
      results[[method]] <- rbind(results[[method]], row_r)
    }
    
    cat(
      "Finished replication", r,
      "with seed", seed_r,
      "| bad draws:", fit_r$bad_draws,
      "| total draws:", fit_r$total_draws,
      "| subsampling:", round(fit_r$timing$subsampling_seconds, 4),
      "| simple:", round(fit_r$timing$simple_seconds, 4),
      "| bc bias:", round(fit_r$timing$bc_bias_seconds, 4),
      "| bc equation:", round(fit_r$timing$bc_equation_seconds, 4),
      "seconds\n"
    )
  }
  
  ############################################################
  ## Save CSV files
  ############################################################
  
  N_name <- plain_number(N)
  k_N_name <- plain_number(k_N)
  m_N_name <- plain_number(m_N)
  alpha_name <- plain_number(alpha)
  
  file_paths <- c()
  
  for (method in methods) {
    
    file_name <- paste0(
      "simulation_",
      "model=", model,
      "_N=", N_name,
      "_k_N=", k_N_name,
      "_m_N=", m_N_name,
      "_alpha=", alpha_name,
      "_method=", method,
      ".csv"
    )
    
    file_path <- file.path(output_dir, file_name)
    
    write.csv(
      results[[method]],
      file_path,
      row.names = FALSE
    )
    
    file_paths[method] <- file_path
  }
  
  list(
    settings = list(
      R = R,
      seed_start = seed_start,
      N = N,
      theta = theta,
      k_N = k_N,
      m_N = m_N,
      alpha = alpha,
      model = model,
      sigma = sigma
    ),
    results = results,
    file_paths = file_paths
  )
}


############################################################
## 4. Linear regression implementation
############################################################

R <- 100
seed_start <- 1
theta <- c(1, 2)
model <- "linear"
sigma <- 1
output_dir <- "./Subbagging new"
N <- 20000
k_N <- floor(N^(1 / 2 + 1 / 4))
alpha <- 1

sim_linear <- RunSubbaggingSimulation(
  R = R,
  seed_start = seed_start,
  N = N,
  theta = theta,
  k_N = k_N,
  alpha = alpha,
  model = model,
  sigma = sigma,
  output_dir = output_dir
)


############################################################
## 5. Logistic regression implementation
############################################################

R <- 100
seed_start <- 1
theta <- c(1, 2)
model <- "logistic"
sigma <- 1
output_dir <- "./Subbagging new"
N <- 20000
k_N <- floor(N^(1 / 2 + 1 / 4))
alpha <- 1

sim_logistic <- RunSubbaggingSimulation(
  R = R,
  seed_start = seed_start,
  N = N,
  theta = theta,
  k_N = k_N,
  alpha = alpha,
  model = model,
  sigma = sigma,
  output_dir = output_dir
)
