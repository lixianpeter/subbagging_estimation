Simu <- function(a) {
  
  library(MASS)
  library(glmnet)
  ## Input Parameters
  id <- a
  list <- read.csv("Paralist.csv", stringsAsFactors = FALSE)
  
  para <- list[id, ]
  
  model <- as.character(para$model)
  p <- para$p
  N <- para$N
  delta <- para$delta
  alpha <- para$alpha
  SNR <- para$SNR
  p0 <- 4
  
  LoopStart <- para$LoopStart
  LoopEnd <- para$LoopEnd
  
  rm(list, para)
  
  if(model == 'Linear'){
    filename1 <- file.path("result", paste0(model, "_SNR", SNR, "_p", p, "_N", sprintf("%.0f", N), "_delta", sprintf("%.2f", delta), "_alpha", alpha, ".csv"))
    filename2 <- file.path("result", paste0(model, "_loglogp_SNR", SNR, "_p", p, "_N", sprintf("%.0f", N), "_delta", sprintf("%.2f", delta), "_alpha", alpha, ".csv"))
    filename3 <- file.path("result", paste0(model, "_dflambda_SNR", SNR, "_p", p, "_N", sprintf("%.0f", N), "_delta", sprintf("%.2f", delta), "_alpha", alpha, ".csv"))
    filename4 <- file.path("result", paste0(model, "_loglogp_dflambda_SNR", SNR, "_p", p, "_N", sprintf("%.0f", N), "_delta", sprintf("%.2f", delta), "_alpha", alpha, ".csv"))
  }else{
    filename1 <- file.path("result", paste0(model, "_p", p, "_N", sprintf("%.0f", N), "_delta", sprintf("%.2f", delta), "_alpha", alpha, ".csv"))
    filename2 <- file.path("result", paste0(model, "_loglogp_p", p, "_N", sprintf("%.0f", N), "_delta", sprintf("%.2f", delta), "_alpha", alpha, ".csv"))
    filename3 <- file.path("result", paste0(model, "_dflambda_p", p, "_N", sprintf("%.0f", N), "_delta", sprintf("%.2f", delta), "_alpha", alpha, ".csv"))
    filename4 <- file.path("result", paste0(model, "_loglogp_dflambda_p", p, "_N", sprintf("%.0f", N), "_delta", sprintf("%.2f", delta), "_alpha", alpha, ".csv"))
  }
  beta_true <- c(seq(-1, -0.5, length.out = p0/2), seq(0.5, 1, length.out = p0/2), rep(0, p - p0))
  
  for (loop in LoopStart:LoopEnd) {
    set.seed(loop)
    
    rho <- 0.5
    idx <- 0:(p - 1)
    
    Sigma <- rho ^ abs(outer(idx, idx, "-"))
    
    ## Generate correlated covariates
    x <- mvrnorm(
      n = N,
      mu = rep(0, p),
      Sigma = Sigma
    )   # N x p
    
    e_sigma <- t(beta_true) %*% Sigma %*% beta_true / SNR
      
    if (model == "Linear") {
      
      # code for Linear model
      error <- rnorm(n = N, mean = 0, sd = sqrt(e_sigma))
      y <- x %*% beta_true + error
      
    } else {
      
      # code for non-Linear model
      y <- rbinom(n = N, 1, plogis(x %*% beta_true))
    }
    
    beta_subsample_list <- list()
    second_derivative_subsample_list <- list()
    Sigma_hat_variance_subsample_list <- list()
    V_subsample_list <- list()
    
    kn <- floor(N^(1/2+delta))
    mn <- floor(alpha*N/kn)
    
    t0 <- Sys.time()  
    for (i in 1:mn) {
      
      subsample_idx <- sample(N, size = kn, replace = FALSE)
      
      y_subsample <- y[subsample_idx]
      x_subsample <- x[subsample_idx,]
      
      
      if (model == "Linear") {
        
        # code for Linear model
        beta_subsample <- solve(t(x_subsample) %*% x_subsample) %*% t(x_subsample) %*% y_subsample
        second_derivative_subsample <- 2 * t(x_subsample) %*% x_subsample / kn
        grad <- -2 * x_subsample * matrix(y_subsample - x_subsample %*% beta_subsample, nrow=kn, ncol=p, byrow = FALSE)
        Sigma_hat_variance_subsample <- t(grad) %*% grad / kn
        
      } else {
        
        fit <- glm(y_subsample ~ x_subsample -1, family = binomial(link = "logit"))
        beta_subsample <- as.numeric(coef(fit)) 
        
        prob <- plogis(x_subsample %*% beta_subsample) 
        second_derivative_subsample <- t(x_subsample) %*% diag(as.numeric(prob * (1 - prob))) %*% x_subsample / kn
        grad <- x_subsample *  matrix(prob - y_subsample, nrow=kn, ncol=p, byrow = FALSE)
        Sigma_hat_variance_subsample <- t(grad) %*% grad / kn
        
      }
      
      beta_subsample_list[[i]] <- beta_subsample
      second_derivative_subsample_list[[i]] <- second_derivative_subsample
      Sigma_hat_variance_subsample_list[[i]] <- Sigma_hat_variance_subsample
    }
    
    beta_average <- Reduce("+", beta_subsample_list) / mn
    # LSA <- Reduce(
    #   `+`,
    #   Map(
    #     function(beta_i, H_i) {
    #       diff <- beta_average - beta_i
    #       as.numeric(t(diff) %*% H_i %*% diff)
    #     },
    #     beta_subsample_list,
    #     second_derivative_subsample_list
    #   )
    # ) / m_N
    
    # BIC_min = k_N * LSA + df * log(N)
    
    alpha0 = (kn * mn)/N 
    
    glm_X=NULL
    glm_Y=NULL
    for (k in seq_len(mn)) {
      #Vk_half <- chol(second_derivative_subsample_list[[k]])
      eig <- eigen(second_derivative_subsample_list[[k]])
      
      Vk_half <- eig$vectors %*%
        diag(sqrt(eig$values)) %*%
        t(eig$vectors)
      
      glm_Y <- rbind(glm_Y, Vk_half%*%beta_subsample_list[[k]])
      glm_X <-rbind(glm_X, Vk_half)
    }
    glm_X=glm_X/sqrt(mn)
    glm_Y=glm_Y/sqrt(mn)
    cy = sd(glm_Y)*sqrt(length(glm_Y)-1)/sqrt(length(glm_Y))
    
    gridLambda=10^seq(0,log10(log(N)/N*1e-1),length=100)
    glmnet_fit <- glmnet(glm_X/cy, glm_Y/cy,
                         family = "gaussian", 
                         alpha=1, 
                         standardize= FALSE, 
                         intercept = FALSE, 
                         penalty.factor= 1/(abs(beta_average)/cy^2/2/length(glm_Y)),
                         lambda=gridLambda)
    
    beta_hat = predict(glmnet_fit,type="coefficients")[-1, ]
    
    BIC_vec <- apply(beta_hat, 2, FUN = function(x) {
      
      df <- sum(x != 0)
      LSA <- 0
      
      for (k in seq_len(mn)) {
        diff <- x - beta_subsample_list[[k]]
        LSA <- LSA +
          as.numeric(
            t(diff) %*%
              second_derivative_subsample_list[[k]] %*%
              diff
          )
      }
      LSA <- LSA/mn
      kn * LSA + df * log(N)
    })
    
    
    beta_hat_optimal  = beta_hat[,which.min(BIC_vec)]    
    lambda_min  = glmnet_fit$lambda[which.min(BIC_vec)]
    BIC_min  = min(BIC_vec)   
    
    Sigma_hat <- Reduce("+", Sigma_hat_variance_subsample_list) / mn
    V_hat <- Reduce("+", second_derivative_subsample_list) / mn
    
    SE_hat <- sqrt(diag((1 + 1 / alpha0) / N * solve(V_hat[1:p0, 1:p0]) %*% Sigma_hat[1:p0, 1:p0] %*% solve(V_hat[1:p0, 1:p0])))
    
    CI_indicator <- as.numeric((beta_hat_optimal[1:p0] + qnorm(0.975) * SE_hat > beta_true[1:p0]) & (beta_hat_optimal[1:p0] - qnorm(0.975) * SE_hat < beta_true[1:p0]))
    
    time <- as.numeric(Sys.time() - t0, units = "secs")
    
    result <- c(loop, beta_hat_optimal, SE_hat, CI_indicator, BIC_min, lambda_min, time)
    
    write.table(t(result), filename1, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  
    
    BIC_vec_loglogp <- apply(beta_hat, 2, FUN = function(x) {
      
      df <- sum(x != 0)
      LSA <- 0
      
      for (k in seq_len(mn)) {
        diff <- x - beta_subsample_list[[k]]
        LSA <- LSA +
          as.numeric(
            t(diff) %*%
              second_derivative_subsample_list[[k]] %*%
              diff
          )
      }
      
      ## My changes
      LSA <- LSA/mn
      
      ## My changes
      kn * LSA + df * log(N) * log(log(p))
    })
    
    beta_hat_optimal_loglogp  = beta_hat[,which.min(BIC_vec_loglogp)]    
    
    lambda_min_loglogp  = glmnet_fit$lambda[which.min(BIC_vec_loglogp)]
    
    BIC_min_loglogp  = min(BIC_vec_loglogp)  
    
    CI_indicator_loglogp <- as.numeric((beta_hat_optimal_loglogp[1:p0] + qnorm(0.975) * SE_hat > beta_true[1:p0]) & (beta_hat_optimal_loglogp[1:p0] - qnorm(0.975) * SE_hat < beta_true[1:p0]))
    
    result <- c(loop, beta_hat_optimal_loglogp, SE_hat, CI_indicator_loglogp, BIC_min_loglogp, lambda_min_loglogp, time)
    
    write.table(t(result), filename2, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    
    
    glmnet_fit <- glmnet(glm_X/cy, glm_Y/cy,
                         family = "gaussian", 
                         alpha=1, 
                         standardize= FALSE, 
                         intercept = FALSE, 
                         penalty.factor= 1/(abs(beta_average)/cy^2/2/length(glm_Y)))
    
    beta_hat = predict(glmnet_fit,type="coefficients")[-1, ]
    
    BIC_vec <- apply(beta_hat, 2, FUN = function(x) {
      
      df <- sum(x != 0)
      LSA <- 0
      
      for (k in seq_len(mn)) {
        diff <- x - beta_subsample_list[[k]]
        LSA <- LSA +
          as.numeric(
            t(diff) %*%
              second_derivative_subsample_list[[k]] %*%
              diff
          )
      }
      LSA <- LSA/mn
      kn * LSA + df * log(N)
    })
    
    
    beta_hat_optimal_dflambda  = beta_hat[,which.min(BIC_vec)]    
    lambda_min_dflambda  = glmnet_fit$lambda[which.min(BIC_vec)]
    BIC_min_dflambda  = min(BIC_vec)   
    
    CI_indicator_dflambda <- as.numeric((beta_hat_optimal_dflambda[1:p0] + qnorm(0.975) * SE_hat > beta_true[1:p0]) & (beta_hat_optimal_dflambda[1:p0] - qnorm(0.975) * SE_hat < beta_true[1:p0]))
    
    result <- c(loop, beta_hat_optimal_dflambda, SE_hat, CI_indicator_dflambda, BIC_min_dflambda, lambda_min_dflambda, time)
    
    write.table(t(result), filename3, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
    
    
    BIC_vec_loglogp_dflambda <- apply(beta_hat, 2, FUN = function(x) {
      
      df <- sum(x != 0)
      LSA <- 0
      
      for (k in seq_len(mn)) {
        diff <- x - beta_subsample_list[[k]]
        LSA <- LSA +
          as.numeric(
            t(diff) %*%
              second_derivative_subsample_list[[k]] %*%
              diff
          )
      }
      
      ## My changes
      LSA <- LSA/mn
      
      ## My changes
      kn * LSA + df * log(N) * log(log(p))
    })
    
    beta_hat_optimal_loglogp_dflambda  = beta_hat[,which.min(BIC_vec_loglogp_dflambda)]    
    
    lambda_min_loglogp_dflambda  = glmnet_fit$lambda[which.min(BIC_vec_loglogp_dflambda)]
    
    BIC_min_loglogp_dflambda  = min(BIC_vec_loglogp_dflambda)  
    
    CI_indicator_loglogp_dflambda <- as.numeric((beta_hat_optimal_loglogp_dflambda[1:p0] + qnorm(0.975) * SE_hat > beta_true[1:p0]) & (beta_hat_optimal_loglogp_dflambda[1:p0] - qnorm(0.975) * SE_hat < beta_true[1:p0]))
    
    result <- c(loop, beta_hat_optimal_loglogp_dflambda, SE_hat, CI_indicator_loglogp_dflambda, BIC_min_loglogp_dflambda, lambda_min_loglogp_dflambda, time)
    
    write.table(t(result), filename4, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}
