#########################################################---
# Required package ----
#########################################################---
if (!requireNamespace("nleqslv", quietly = TRUE)) {
  install.packages("nleqslv")
}

#########################################################---
# 1. Simulation settings ----
#########################################################---
## Choose one model: "linear" or "logistic"
model <- "linear"

## Choose one estimator:
## "simple", "bc_bias", or "bc_equation"
method <- "bc_equation"

## Number of simulation replications
R <- 100

## Seed for the first replication
seed_start <- 1

## True parameter, including the intercept
true_theta <- c(1, 2)

## Error standard deviation for linear regression
sigma <- 1

## Full sample size
N <- 20000

## Subsample size
k_N <- floor(N^(1 / 2 + 1 / 4))

## Controls the number of subsamples:
## m_N = floor(alpha * N / k_N)
alpha <- 1

## Output folder
output_dir <- "./Subbagging new"


## model and estimation method.
model <- match.arg(model, c("linear", "logistic"))
method <- match.arg(method, c("simple", "bc_bias", "bc_equation"))

m_N <- floor(alpha * N / k_N)

## Effective value after applying floor() to m_N.
alpha_N <- k_N * m_N / N

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#########################################################---
# 2. Functions used in the simulation ----
## 2.1 Fit the ordinary estimator on one data set ----
#########################################################---


FitTheta <- function(Xs, ys, model) {
  
  if (model == "linear") {
    theta_hat <- lm.fit(Xs, ys)$coefficients
  } else {
    fit <- suppressWarnings(
      glm.fit(Xs, ys, family = binomial())
    )
    
    if (!isTRUE(fit$converged)) {
      stop("Logistic regression did not converge.")
    }
    
    theta_hat <- fit$coefficients
  }
  
  theta_hat <- as.numeric(theta_hat)
  
  ## A singular regression fit can return non-finite
  ## coefficient estimates.
  if (any(!is.finite(theta_hat))) {
    stop("The estimator is not finite.")
  }
  
  theta_hat
}

#########################################################---
## 2.2 Compute the estimating functions and derivative matrix ----
#########################################################---
Moments <- function(theta, Xs, ys, model) {
  
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
  
  ## Individual estimating functions psi_i(theta).
  psi <- Xs * as.numeric(ys - mu)
  
  ## Sample derivative matrix V_hat(theta).
  V <- -crossprod(Xs, Xs * w) / n
  
  list(
    psi = psi,
    V = V,
    w = w,
    h = h
  )
}

#########################################################---
## 2.3 Estimate the bias term on one subsample ----
#########################################################---
Bhat <- function(theta, Xs, ys, model) {
  
  n <- nrow(Xs)
  d <- ncol(Xs)
  
  M <- Moments(theta, Xs, ys, model)
  V <- M$V
  
  ## The bias formula requires V to be invertible.
  if (!is.finite(rcond(V)) || rcond(V) < 1e-12) {
    stop("V is singular or numerically ill-conditioned.")
  }
  
  V_inverse <- solve(V)
  
  psi <- M$psi
  A <- psi %*% t(V_inverse)
  XA <- rowSums(Xs * A)
  
  D_part <- colMeans(Xs * (-M$w * XA))
  V_part <- as.numeric(V %*% colMeans(A))
  term1 <- -(D_part - V_part)
  
  Q <- crossprod(A) / n
  term2 <- rep(0, d)
  
  if (model == "logistic") {
    for (j in seq_len(d)) {
      H_j <- -crossprod(
        Xs,
        Xs * (M$h * Xs[, j])
      ) / n
      
      term2[j] <- 0.5 * sum(H_j * Q)
    }
  }
  
  B <- -as.numeric(
    V_inverse %*% (term1 + term2)
  )
  
  if (any(!is.finite(B))) {
    stop("The estimated bias term is not finite.")
  }
  
  B
}

#########################################################---
## 2.4 Estimating equation ----
#########################################################---
AdjustedScore <- function(theta, Xs, ys, model) {
  
  M <- Moments(theta, Xs, ys, model)
  score <- colSums(M$psi)
  B <- Bhat(theta, Xs, ys, model)
  
  adjusted_score <- as.numeric(
    score + M$V %*% B
  )
  
  if (any(!is.finite(adjusted_score))) {
    stop("The adjusted score is not finite.")
  }
  
  adjusted_score
}

#########################################################---
## 2.5 Solve the estimating equation ----
#########################################################---
SolveBCEquation <- function(theta_start, Xs, ys, model) {
  
  solution <- nleqslv(
    x = theta_start,
    fn = function(theta) {
      AdjustedScore(theta, Xs, ys, model)
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
  
  theta_bc <- as.numeric(solution$x)
  
  ## termcd 1 or 2 indicates convergence in nleqslv.
  if (!(solution$termcd %in% c(1, 2))) {
    stop("The adjusted estimating equation did not converge.")
  }
  
  if (any(!is.finite(theta_bc))) {
    stop("The adjusted estimating-equation estimate is not finite.")
  }
  
  final_score <- AdjustedScore(theta_bc, Xs, ys, model)
  
  if (max(abs(final_score)) > 1e-6) {
    stop("The adjusted estimating-equation residual is too large.")
  }
  
  theta_bc
}

#########################################################---
## 2.6 Calculate the full-data asymptotic standard deviation ----
#########################################################---
CalculateASD <- function(X, y, true_theta, model) {
  
  ## This is the simulation benchmark based on the full-data
  ## sandwich variance evaluated at the known true parameter.
  n <- nrow(X)
  M <- Moments(true_theta, X, y, model)
  
  V_hat <- M$V
  Sigma_hat <- crossprod(M$psi) / n
  
  if (!is.finite(rcond(V_hat)) || rcond(V_hat) < 1e-12) {
    stop("The full-data derivative matrix is singular.")
  }
  
  V_inverse <- solve(V_hat)
  Xi_hat <- V_inverse %*%
    Sigma_hat %*%
    t(V_inverse)
  
  sqrt(diag(Xi_hat) / n)
}

#########################################################---
## 2.7 Calculate the subbagging standard error ----
#########################################################---
CalculateSSE <- function(theta_subsample, N, k_N) {
  
  ## Omega_hat uses the denominator m_N, matching the paper's
  ## subbagging variance estimator.
  theta_average <- colMeans(theta_subsample)
  centered_theta <- sweep(
    theta_subsample,
    MARGIN = 2,
    STATS = theta_average,
    FUN = "-"
  )
  
  m_N <- nrow(theta_subsample)
  Omega_hat <- crossprod(centered_theta) / m_N
  
  sqrt(
    k_N * diag(Omega_hat) / N
  )
}

#########################################################---
## 2.8 Format numbers used in the output filename ----
#########################################################---
PlainNumber <- function(x) {
  x <- format(x, scientific = FALSE, trim = TRUE)
  gsub("\\.", "p", x)
}

#########################################################---
# 3. Formatting on the output filename ----
#########################################################---
file_name <- paste0(
  "simulation_",
  "model=", model,
  "_N=", PlainNumber(N),
  "_k_N=", PlainNumber(k_N),
  "_m_N=", PlainNumber(m_N),
  "_alpha=", PlainNumber(alpha),
  "_method=", method,
  ".csv"
)

file_path <- file.path(output_dir, file_name)

#########################################################---
# 4. Main simulation loop ----
#########################################################---
for (r in seq_len(R)) {
  
  seed_r <- seed_start + r - 1
  set.seed(seed_r)
  
  #########################################################---
  ## 4.1 Generate one full dataset ----
  #########################################################---
  p <- length(true_theta) - 1
  
  X_covariates <- matrix(
    rnorm(N * p),
    nrow = N,
    ncol = p
  )
  
  colnames(X_covariates) <- paste0("x", seq_len(p))
  
  X <- cbind(
    Intercept = 1,
    X_covariates
  )
  
  eta <- as.numeric(X %*% true_theta)
  
  if (model == "linear") {
    y <- eta + rnorm(N, mean = 0, sd = sigma)
  } else {
    probability <- plogis(eta)
    y <- rbinom(N, size = 1, prob = probability)
  }
  
  ########################################################---
  # 4.2 Full-data ASD evaluated at the true parameter ----
  ########################################################---
  asd <- CalculateASD(
    X = X,
    y = y,
    true_theta = true_theta,
    model = model
  )
  
  ########################################################---
  # 4.3 Storage for the selected subbagging estimator ----
  ########################################################---
  coefficient_names <- colnames(X)
  d <- ncol(X)
  
  theta_subsample <- matrix(
    NA_real_,
    nrow = m_N,
    ncol = d,
    dimnames = list(NULL, coefficient_names)
  )
  
  bad_draws <- 0
  total_draws <- 0
  
  ########################################################---
  # 4.4 Start timing the selected subbagging method ----
  ########################################################---
  ## Resetting the seed ensures that separate runs of the
  ## three methods use the same sequence of subsamples.
  set.seed(seed_r)
  
  time_start <- proc.time()[["elapsed"]]
  
  ########################################################---
  # 4.5 Draw and process m_N valid subsamples ----
  ########################################################---
  for (b in seq_len(m_N)) {
    
    valid_draw <- FALSE
    
    while (!valid_draw) {
      
      total_draws <- total_draws + 1
      
      subsample_id <- sample(
        seq_len(N),
        size = k_N,
        replace = FALSE
      )
      
      Xs <- X[subsample_id, , drop = FALSE]
      ys <- y[subsample_id]
      
    #########################################################---
    ## 4.5.1 Apply the selected estimator to this subsample ----
    #########################################################---
      theta_result <- try({
        
        theta_hat <- FitTheta(Xs, ys, model)
        
        if (method == "simple") {
          
          theta_hat
          
        } else if (method == "bc_bias") {
          
          B <- Bhat(theta_hat, Xs, ys, model)
          theta_hat - B / k_N
          
        } else {
          
          SolveBCEquation(
            theta_start = theta_hat,
            Xs = Xs,
            ys = ys,
            model = model
          )
        }
        
      }, silent = TRUE)
      
      #########################################################---
      ## 4.5.2 Redraw only when the selected estimator fails ----
      #########################################################---
      if (inherits(theta_result, "try-error")) {
        bad_draws <- bad_draws + 1
        next
      }
      
      theta_subsample[b, ] <- theta_result
      valid_draw <- TRUE
    }
  }
  
  #########################################################---
  ## 4.6 Average the m_N subsample estimates ----
  #########################################################---
  estimate <- colMeans(theta_subsample)
  
  #########################################################---
  ## 4.7 End timing ----
  #########################################################---
  time_end <- proc.time()[["elapsed"]]
  time_total_seconds <- time_end - time_start
  
  #########################################################---
  ## 4.8 Calculate SSE and variance-inflation adjustments ----
  #########################################################---
  sse <- CalculateSSE(
    theta_subsample = theta_subsample,
    N = N,
    k_N = k_N
  )
  
  adjustment_factor <- sqrt(1 + 1 / alpha_N)
  adjusted_asd <- adjustment_factor * asd
  adjusted_sse <- adjustment_factor * sse
  
  #########################################################---
  ## 4.9 Save this replication ----
  #########################################################---
  result <- data.frame(
    replication = r,
    seed = seed_r,
    N = N,
    k_N = k_N,
    m_N = m_N,
    alpha = alpha,
    alpha_N = alpha_N,
    model = model,
    method = method,
    bad_draws = bad_draws,
    total_draws = total_draws,
    time_total_seconds = time_total_seconds,
    check.names = FALSE
  )
  
  for (j in seq_along(coefficient_names)) {
    result[[paste0("estimate_", coefficient_names[j])]] <- estimate[j]
  }
  
  for (j in seq_along(coefficient_names)) {
    result[[paste0("asd_", coefficient_names[j])]] <- asd[j]
  }
  
  for (j in seq_along(coefficient_names)) {
    result[[paste0("adjusted_asd_", coefficient_names[j])]] <- adjusted_asd[j]
  }
  
  for (j in seq_along(coefficient_names)) {
    result[[paste0("sse_", coefficient_names[j])]] <- sse[j]
  }
  
  for (j in seq_along(coefficient_names)) {
    result[[paste0("adjusted_sse_", coefficient_names[j])]] <- adjusted_sse[j]
  }
  
  write.table(
    result,
    file = file_path,
    sep = ",",
    row.names = FALSE,
    col.names = (r == 1),
    append = (r > 1),
    quote = FALSE
  )
  
  cat(
    "Finished replication", r,
    "| seed:", seed_r,
    "| model:", model,
    "| method:", method,
    "| bad draws:", bad_draws,
    "| total draws:", total_draws,
    "| time:", round(time_total_seconds, 4),
    "seconds\n"
  )
}

#########################################################---
# 5. Completion message ----
#########################################################---
cat("\nSimulation completed.\n")
cat("Output file:", file_path, "\n")

