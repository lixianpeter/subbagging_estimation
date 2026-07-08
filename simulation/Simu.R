############################################################
## Clean Subbagging Simulation Script
##
## One manually defined setting only.
## Supports:
##   1. model_type = "logistic"
##   2. model_type = "linear"
##
## For each setting, save one raw CSV file.
## The CSV file name includes model_type, N, alpha, k_N, and m_N.
############################################################

rm(list = ls())

############################################################
## 0. Safe solver
############################################################

safe_solve <- function(A, B, tol = 1e-10) {
  
  if (any(!is.finite(A)) || any(!is.finite(B))) {
    stop("Non-finite value in solve.")
  }
  
  if (nrow(A) != ncol(A)) {
    stop("A is not square.")
  }
  
  qrA <- qr(A, tol = tol)
  
  if (qrA$rank < ncol(A)) {
    stop("Singular matrix.")
  }
  
  solve(A, B)
}


############################################################
## 0b. Model type checker
############################################################

CheckModelType <- function(model_type) {
  
  model_type <- match.arg(
    arg = model_type,
    choices = c("logistic", "linear")
  )
  
  model_type
}


############################################################
## 1. Subbagging core
############################################################

SubbaggingCore <- function(data,
                           k_N,
                           m_N,
                           model_type = "logistic",
                           y_name = "y",
                           max_bad_draws = NULL) {
  
  model_type <- CheckModelType(model_type)
  
  y <- data[[y_name]]
  x <- as.matrix(data[, setdiff(names(data), y_name), drop = FALSE])
  
  N <- nrow(x)
  p <- ncol(x)
  
  if (k_N > N) stop("k_N cannot be larger than N.")
  if (k_N <= p) stop("k_N must be larger than p.")
  
  if (is.null(max_bad_draws)) {
    max_bad_draws <- max(10000, 10 * m_N)
  }
  
  ############################################################
  ## M-estimator on one subsample
  ############################################################
  
  fit_original <- function(x_sub, y_sub) {
    
    if (qr(x_sub)$rank < ncol(x_sub)) {
      stop("Rank deficient x_sub.")
    }
    
    if (model_type == "logistic") {
      
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
        stop("Non-finite logistic coefficient.")
      }
      
      return(beta_hat)
    }
    
    if (model_type == "linear") {
      
      fit <- lm.fit(x = x_sub, y = y_sub)
      beta_hat <- as.numeric(coef(fit))
      
      if (any(!is.finite(beta_hat))) {
        stop("Non-finite linear coefficient.")
      }
      
      return(beta_hat)
    }
    
    stop("Unknown model_type.")
  }
  
  ############################################################
  ## Compute psi, V, and B_hat
  ############################################################
  
  get_quantities <- function(x_sub, y_sub, beta) {
    
    k <- nrow(x_sub)
    p <- ncol(x_sub)
    
    eta <- as.numeric(x_sub %*% beta)
    
    if (model_type == "logistic") {
      
      prob <- plogis(eta)
      w <- prob * (1 - prob)
      
      if (min(w) < 1e-14 && mean(w) < 1e-8) {
        stop("Possible logistic separation.")
      }
      
      psi <- x_sub * as.numeric(y_sub - prob)
      
      V <- -crossprod(x_sub, x_sub * w) / k
      
      A <- t(safe_solve(V, t(psi)))
      
      xa <- rowSums(x_sub * A)
      D_a <- -x_sub * as.numeric(w * xa)
      
      first_part <- -colMeans(D_a - psi)
      
      second_part <- -0.5 * colMeans(
        x_sub * as.numeric(w * (1 - 2 * prob) * xa^2)
      )
      
      B_hat <- -safe_solve(V, first_part + second_part)
      B_hat <- as.numeric(B_hat)
      
      return(list(
        psi = psi,
        V = V,
        B_hat = B_hat
      ))
    }
    
    if (model_type == "linear") {
      
      resid <- as.numeric(y_sub - eta)
      
      psi <- x_sub * resid
      
      V <- -crossprod(x_sub) / k
      
      A <- t(safe_solve(V, t(psi)))
      
      xa <- rowSums(x_sub * A)
      D_a <- -x_sub * as.numeric(xa)
      
      first_part <- -colMeans(D_a - psi)
      
      ## Linear regression has zero second derivative in beta.
      second_part <- rep(0, p)
      
      B_hat <- -safe_solve(V, first_part + second_part)
      B_hat <- as.numeric(B_hat)
      
      return(list(
        psi = psi,
        V = V,
        B_hat = B_hat
      ))
    }
    
    stop("Unknown model_type.")
  }
  
  ############################################################
  ## bc2: theta_hat - B_hat / k_N
  ############################################################
  
  bias_correct_add <- function(x_sub, y_sub, beta_hat) {
    
    k <- nrow(x_sub)
    q <- get_quantities(x_sub, y_sub, beta_hat)
    
    beta_bc <- beta_hat - q$B_hat / k
    
    if (any(!is.finite(beta_bc))) {
      stop("Non-finite bc2.")
    }
    
    beta_bc
  }
  
  ############################################################
  ## bc3: solve adjusted estimating equation
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
      
      adjusted_score <- mean_psi + as.numeric(q$V %*% q$B_hat) / k
      
      step <- safe_solve(q$V, adjusted_score)
      step <- as.numeric(step)
      
      beta_new <- beta - step
      
      if (any(!is.finite(beta_new))) {
        stop("Non-finite bc3.")
      }
      
      if (max(abs(beta_new - beta)) < tol) {
        beta <- beta_new
        break
      }
      
      beta <- beta_new
    }
    
    beta
  }
  
  ############################################################
  ## Main loop
  ############################################################
  
  sum_simple <- rep(0, p)
  sum_bc2 <- rep(0, p)
  sum_bc3 <- rep(0, p)
  
  sumsq_simple <- rep(0, p)
  sumsq_bc2 <- rep(0, p)
  sumsq_bc3 <- rep(0, p)
  
  accepted <- 0
  bad_draws <- 0
  total_draws <- 0
  
  while (accepted < m_N) {
    
    total_draws <- total_draws + 1
    
    one_draw <- tryCatch({
      
      id <- sample.int(N, size = k_N, replace = FALSE)
      
      x_sub <- x[id, , drop = FALSE]
      y_sub <- y[id]
      
      beta_hat <- fit_original(x_sub, y_sub)
      beta_bc2 <- bias_correct_add(x_sub, y_sub, beta_hat)
      beta_bc3 <- bias_correct_equation(x_sub, y_sub, beta_bc2)
      
      list(
        simple = beta_hat,
        bc2 = beta_bc2,
        bc3 = beta_bc3
      )
      
    }, error = function(e) {
      NULL
    })
    
    if (is.null(one_draw)) {
      
      bad_draws <- bad_draws + 1
      
      if (bad_draws > max_bad_draws) {
        stop("Too many bad subsamples.")
      }
      
      next
    }
    
    accepted <- accepted + 1
    
    sum_simple <- sum_simple + one_draw$simple
    sum_bc2 <- sum_bc2 + one_draw$bc2
    sum_bc3 <- sum_bc3 + one_draw$bc3
    
    sumsq_simple <- sumsq_simple + one_draw$simple^2
    sumsq_bc2 <- sumsq_bc2 + one_draw$bc2^2
    sumsq_bc3 <- sumsq_bc3 + one_draw$bc3^2
  }
  
  beta_simple <- sum_simple / m_N
  beta_bc2 <- sum_bc2 / m_N
  beta_bc3 <- sum_bc3 / m_N
  
  sd_simple <- sqrt(pmax(sumsq_simple / m_N - beta_simple^2, 0))
  sd_bc2 <- sqrt(pmax(sumsq_bc2 / m_N - beta_bc2^2, 0))
  sd_bc3 <- sqrt(pmax(sumsq_bc3 / m_N - beta_bc3^2, 0))
  
  list(
    simple = beta_simple,
    bc2 = beta_bc2,
    bc3 = beta_bc3,
    
    simple_sd = sd_simple,
    bc2_sd = sd_bc2,
    bc3_sd = sd_bc3,
    
    bad_draws = bad_draws,
    total_draws = total_draws
  )
}


############################################################
## 2. Full-sample estimators
############################################################

FitFullLogistic <- function(data, y_name = "y") {
  
  y <- data[[y_name]]
  x <- as.matrix(data[, setdiff(names(data), y_name), drop = FALSE])
  
  fit <- suppressWarnings(
    glm.fit(
      x = x,
      y = y,
      family = binomial(link = "logit"),
      control = glm.control(maxit = 50)
    )
  )
  
  beta <- as.numeric(coef(fit))
  
  if (any(!is.finite(beta))) {
    stop("Full sample logistic failed.")
  }
  
  prob <- fitted(fit)
  w <- prob * (1 - prob)
  
  fisher <- crossprod(x, x * w)
  cov_hat <- safe_solve(fisher, diag(ncol(x)))
  se <- sqrt(diag(cov_hat))
  
  list(
    beta = beta,
    se = se
  )
}


FitFullLinear <- function(data, y_name = "y") {
  
  y <- data[[y_name]]
  x <- as.matrix(data[, setdiff(names(data), y_name), drop = FALSE])
  
  if (qr(x)$rank < ncol(x)) {
    stop("Full sample linear design is rank deficient.")
  }
  
  fit <- lm.fit(x = x, y = y)
  beta <- as.numeric(coef(fit))
  
  if (any(!is.finite(beta))) {
    stop("Full sample linear failed.")
  }
  
  resid <- as.numeric(y - x %*% beta)
  n <- nrow(x)
  p <- ncol(x)
  sigma2_hat <- sum(resid^2) / max(n - p, 1)
  
  xtx_inv <- safe_solve(crossprod(x), diag(p))
  se <- sqrt(diag(sigma2_hat * xtx_inv))
  
  list(
    beta = beta,
    se = se
  )
}


FitFullEstimator <- function(data,
                             model_type = "logistic",
                             y_name = "y") {
  
  model_type <- CheckModelType(model_type)
  
  if (model_type == "logistic") {
    return(FitFullLogistic(data = data, y_name = y_name))
  }
  
  if (model_type == "linear") {
    return(FitFullLinear(data = data, y_name = y_name))
  }
  
  stop("Unknown model_type.")
}


############################################################
## 3. Data generation
############################################################

GenerateData <- function(N,
                         theta0 = c(0, 1),
                         model_type = "logistic",
                         noise_sd = 1) {
  
  model_type <- CheckModelType(model_type)
  
  if (length(theta0) != 2) {
    stop("This simple DGP expects theta0 to have length 2: intercept and x1 coefficient.")
  }
  
  x1 <- rnorm(N)
  
  X <- cbind(
    intercept = 1,
    x1 = x1
  )
  
  eta <- as.numeric(X %*% theta0)
  
  if (model_type == "logistic") {
    
    prob <- plogis(eta)
    y <- rbinom(N, size = 1, prob = prob)
    
    return(data.frame(
      y = y,
      intercept = X[, 1],
      x1 = X[, 2]
    ))
  }
  
  if (model_type == "linear") {
    
    if (!is.finite(noise_sd) || noise_sd <= 0) {
      stop("noise_sd must be positive for linear regression.")
    }
    
    y <- eta + rnorm(N, mean = 0, sd = noise_sd)
    
    return(data.frame(
      y = y,
      intercept = X[, 1],
      x1 = X[, 2]
    ))
  }
  
  stop("Unknown model_type.")
}


############################################################
## 4. One replication, one row output
############################################################

OneReplication <- function(rep_id,
                           seed_base,
                           N,
                           alpha,
                           k_N,
                           m_N,
                           theta0 = c(0, 1),
                           model_type = "logistic",
                           noise_sd = 1) {
  
  model_type <- CheckModelType(model_type)
  
  seed_index <- seed_base + rep_id
  set.seed(seed_index)
  
  data <- GenerateData(
    N = N,
    theta0 = theta0,
    model_type = model_type,
    noise_sd = noise_sd
  )
  
  k_N <- as.integer(k_N)
  if (k_N <= 0) stop("k_N must be positive.")
  
  m_N <- as.integer(m_N)
  if (m_N <= 0) stop("m_N must be positive.")
  
  full <- FitFullEstimator(
    data = data,
    model_type = model_type
  )
  
  sub <- SubbaggingCore(
    data = data,
    k_N = k_N,
    m_N = m_N,
    model_type = model_type
  )
  
  data.frame(
    seed_index = seed_index,
    rep = rep_id,
    model_type = model_type,
    noise_sd = ifelse(model_type == "linear", noise_sd, NA_real_),
    
    N = N,
    alpha = alpha,
    k_N = k_N,
    m_N = m_N,
    
    full_theta1 = full$beta[1],
    full_theta2 = full$beta[2],
    full_sd_theta1 = full$se[1],
    full_sd_theta2 = full$se[2],
    
    simple_theta1 = sub$simple[1],
    simple_theta2 = sub$simple[2],
    simple_sd_theta1 = sub$simple_sd[1],
    simple_sd_theta2 = sub$simple_sd[2],
    
    bc2_theta1 = sub$bc2[1],
    bc2_theta2 = sub$bc2[2],
    bc2_sd_theta1 = sub$bc2_sd[1],
    bc2_sd_theta2 = sub$bc2_sd[2],
    
    bc3_theta1 = sub$bc3[1],
    bc3_theta2 = sub$bc3[2],
    bc3_sd_theta1 = sub$bc3_sd[1],
    bc3_sd_theta2 = sub$bc3_sd[2],
    
    bad_draws = sub$bad_draws,
    total_draws = sub$total_draws,
    
    stringsAsFactors = FALSE
  )
}


############################################################
## 6. Folder name
############################################################

SettingName <- function(model_type, N, alpha, k_N, m_N, noise_sd = 1) {
  
  model_type <- CheckModelType(model_type)
  
  alpha_name <- ifelse(
    abs(alpha - 1) < 1e-12,
    "alpha_1",
    ifelse(
      abs(alpha - 1 / 3) < 1e-12,
      "alpha_1over3",
      paste0("alpha_", alpha)
    )
  )
  
  noise_name <- ifelse(
    model_type == "linear",
    paste0("__noise_sd_", noise_sd),
    ""
  )
  
  paste0(
    "model_", model_type,
    "__N_", N,
    "__", alpha_name,
    "__k_N_", k_N,
    "__m_N_", m_N,
    noise_name
  )
}


############################################################
## 7. Run one setting and write raw file only
############################################################

RunOneSetting <- function(R,
                          N,
                          alpha,
                          k_N,
                          m_N,
                          model_type = "logistic",
                          theta0 = c(0, 1),
                          noise_sd = 1,
                          seed_base = seed,
                          output_root = "subbagging_results") {
  
  model_type <- CheckModelType(model_type)
  
  setting_name <- SettingName(
    model_type = model_type,
    N = N,
    alpha = alpha,
    k_N = k_N,
    m_N = m_N,
    noise_sd = noise_sd
  )
  setting_dir <- file.path(output_root, setting_name)
  
  dir.create(setting_dir, recursive = TRUE, showWarnings = FALSE)
  
  raw_list <- vector("list", R)
  
  for (r in 1:R) {
    
    message("[", setting_name, "] replication ", r, " / ", R)
    
    raw_list[[r]] <- OneReplication(
      rep_id = r,
      seed_base = seed_base,
      N = N,
      alpha = alpha,
      k_N = k_N,
      m_N = m_N,
      theta0 = theta0,
      model_type = model_type,
      noise_sd = noise_sd
    )
  }
  
  raw <- do.call(rbind, raw_list)
  
  noise_file_part <- ifelse(
    model_type == "linear",
    paste0("_noise_sd=", noise_sd),
    ""
  )
  
  raw_file_name <- paste0(
    "model=", model_type,
    "_N=", N,
    "_alpha=", alpha,
    "_k_N=", k_N,
    "_m_N=", m_N,
    noise_file_part,
    ".csv"
  )
  
  raw_file_path <- file.path(setting_dir, raw_file_name)
  
  write.csv(
    raw,
    file = raw_file_path,
    row.names = FALSE
  )
  
  message("Saved raw file to: ", raw_file_path)
  
  invisible(list(
    raw = raw,
    raw_file_path = raw_file_path,
    setting_dir = setting_dir
  ))
}


############################################################
## Global run values: one manually defined setting only
############################################################

seed <- 12345
R <- 20
N <- 20000
alpha <- 1

## Choose one:
model_type <- "logistic"
## model_type <- "linear"

theta0 <- c(0, 1)

## Used only when model_type = "linear".
## Ignored when model_type = "logistic".
noise_sd <- 1

## Choose the subsample size k_N directly through the formula you want.
k_N <- floor(N^(5 / 12))

## Choose the number of subsamples m_N directly through the formula you want.
## Paper option 1:
m_N <- floor(alpha * N / k_N)

## Paper option 2, if needed later:
## m_N <- floor(alpha * N^(4 / 3) / k_N)

output_root <- "subbagging_one_setting_results"


############################################################
## Run one setting
############################################################

RunOneSetting(
  R = R,
  N = N,
  alpha = alpha,
  k_N = k_N,
  m_N = m_N,
  model_type = model_type,
  theta0 = theta0,
  noise_sd = noise_sd,
  seed_base = seed,
  output_root = output_root
)
