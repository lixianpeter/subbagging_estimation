############################################################
## 1. Subbagging core
############################################################
SubbaggingCore <- function(data, k_N, alpha,
                           model = c("linear", "logistic"),
                           seed = 1) {
  # model settings
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
  ## M-estimator on one subsample
  ############################################################
  fit_theta <- function(Xs, ys) {
    if (model == "linear") {
      beta <- lm.fit(Xs, ys)$coefficients
    } else {
      beta <- glm.fit(Xs, ys, family = binomial())$coefficients
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
      w  <- rep(1, n)
      h  <- rep(0, n)
    } else {
      mu <- 1 / (1 + exp(-eta))
      w  <- mu * (1 - mu)
      h  <- w * (1 - 2 * mu)
    }
    
    psi <- Xs * as.numeric(ys - mu)
    V   <- -crossprod(Xs, Xs * w) / n
    
    list(psi = psi, V = V, w = w, h = h)
  }
  
  Bhat <- function(theta, Xs, ys) {
    n <- nrow(Xs)
    d <- ncol(Xs)
    
    M <- moments(theta, Xs, ys)
    V <- M$V
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
      for (j in 1:d) {
        Hj <- -crossprod(Xs, Xs * (M$h * Xs[, j])) / n
        term2[j] <- 0.5 * sum(Hj * Q)
      }
    }
    
    B <- -as.numeric(Vinv %*% (term1 + term2))
    
    if (any(!is.finite(B))) {
      stop("bad Bhat")
    }
    
    B
  }
  
  adjusted_score <- function(theta, Xs, ys) {
    M <- moments(theta, Xs, ys)
    score <- colSums(M$psi)
    B <- Bhat(theta, Xs, ys)
    
    as.numeric(score + M$V %*% B)
  }
  
  solve_bc_equation <- function(theta_start, Xs, ys) {
    obj <- function(theta) {
      score <- adjusted_score(theta, Xs, ys)
      sum(score^2)
    }
    
    theta_bc <- optim(theta_start, obj, method = "BFGS")$par
    
    if (any(!is.finite(theta_bc))) {
      stop("bad bc equation")
    }
    
    theta_bc
  }
  
  theta_simple      <- matrix(0, m_N, d)
  theta_bc_bias     <- matrix(0, m_N, d)
  theta_bc_equation <- matrix(0, m_N, d)
  
  bad_draws <- 0
  total_draws <- 0
  
  for (b in 1:m_N) {
    
    valid_draw <- FALSE
    
    while (!valid_draw) {
      
      total_draws <- total_draws + 1
      
      id <- sample(1:N, k_N, replace = FALSE)
      
      Xs <- X[id, , drop = FALSE]
      ys <- y[id]
      
      tmp <- try({
        
        theta_hat <- fit_theta(Xs, ys)
        
        theta_simple[b, ] <- theta_hat
        
        ############################################################
        ## bc2: theta_hat - B_hat / k_N
        ############################################################
        B <- Bhat(theta_hat, Xs, ys)
        theta_bc_bias[b, ] <- theta_hat - B / k_N
        
        ############################################################
        ## bc3: solve adjusted estimating equation
        ############################################################
        theta_bc_equation[b, ] <- solve_bc_equation(theta_hat, Xs, ys)
        
      }, silent = TRUE)
      
      if (!inherits(tmp, "try-error")) {
        valid_draw <- TRUE
      } else {
        bad_draws <- bad_draws + 1
      }
    }
  }
  
  colnames(theta_simple)      <- colnames(X)
  colnames(theta_bc_bias)     <- colnames(X)
  colnames(theta_bc_equation) <- colnames(X)
  
  estimate <- rbind(
    simple      = colMeans(theta_simple),
    bc_bias     = colMeans(theta_bc_bias),
    bc_equation = colMeans(theta_bc_equation)
  )
  
  se_fun <- function(A) {
    Omega <- cov(A) * (m_N - 1) / m_N
    sqrt((1 + 1 / alpha_N) * k_N * diag(Omega) / N)
  }
  
  se <- rbind(
    simple      = se_fun(theta_simple),
    bc_bias     = se_fun(theta_bc_bias),
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
  colnames(X) <- paste0("x", 1:p)
  
  eta <- as.numeric(cbind(1, X) %*% theta)
  
  if (model == "linear") {
    y <- eta + rnorm(N, sd = sigma)
  }
  
  if (model == "logistic") {
    prob <- 1 / (1 + exp(-eta))
    y <- rbinom(N, size = 1, prob = prob)
  }
  
  data <- data.frame(y = y, X)
  
  attr(data, "theta_true") <- theta
  
  data
}


############################################################
## 3. Simulation 
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
  
  for (r in 1:R) {
    
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
      se_r       <- fit_r$se[method, ]
      
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
        total_draws = fit_r$total_draws
      )
      
      for (j in seq_along(coef_names)) {
        row_r[[paste0("estimate_", coef_names[j])]] <- estimate_r[j]
      }
      
      for (j in seq_along(coef_names)) {
        row_r[[paste0("se_", coef_names[j])]] <- se_r[j]
      }
      
      results[[method]] <- rbind(results[[method]], row_r)
    }
    
    cat("Finished replication", r, "with seed", seed_r,
        "| bad draws:", fit_r$bad_draws,
        "| total draws:", fit_r$total_draws, "\n")
  }
  
  ############################################################
  ## Save CSV files
  ## No scientific notation will appear in file names
  ############################################################
  N_name     <- plain_number(N)
  k_N_name   <- plain_number(k_N)
  m_N_name   <- plain_number(m_N)
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
    
    write.csv(results[[method]], file_path, row.names = FALSE)
    
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
## 4. Implementation 
############################################################

R = 100
seed_start = 1
theta = c(1, 2)
model = "linear"
sigma = 1
output_dir = "./Subbagging new"
N = 500000
k_N = floor(N^(1/2 + 1/4))
alpha = 1

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


R = 100
seed_start = 1
theta = c(1, 2)
model = "logistic"
sigma = 1
output_dir = "./Subbagging new"
N = 500000
k_N = floor(N^(1/2 + 1/4))
alpha = 1

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
