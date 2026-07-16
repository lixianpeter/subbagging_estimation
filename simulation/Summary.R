############################################################
## Summarise simulation CSV files in one folder
############################################################

SummariseSimulationFolder <- function(
    folder = ".",
    theta = c(1, 2),
    output_file = "all_simulation_results_summary.csv",
    z_value = 1.96
) {
  
  ############################################################
  ## 1. Find simulation CSV files ----
  ############################################################
  
  files <- list.files(
    path = folder,
    pattern = "^simulation_.*\\.csv$",
    full.names = TRUE
  )
  
  files <- files[
    !grepl(
      pattern = "summary",
      x = basename(files),
      ignore.case = TRUE
    )
  ]
  
  if (length(files) == 0) {
    stop("No simulation CSV files were found in the selected folder.")
  }
  
  
  ############################################################
  ## 2. Summarise one simulation file ----
  ############################################################
  
  SummariseOneFile <- function(file) {
    
    df <- read.csv(
      file = file,
      check.names = FALSE
    )
    
    estimate_cols <- grep(
      pattern = "^estimate_",
      x = names(df),
      value = TRUE
    )
    
    if (length(estimate_cols) == 0) {
      stop(
        paste0(
          "No estimate columns were found in file: ",
          basename(file)
        )
      )
    }
    
    if (length(estimate_cols) != length(theta)) {
      stop(
        paste0(
          "The number of estimate columns does not match ",
          "the length of theta in file: ",
          basename(file),
          "\nNumber of estimate columns = ",
          length(estimate_cols),
          "\nLength of theta = ",
          length(theta)
        )
      )
    }
    
    
    ############################################################
    ## 2.1 Helper: extract one constant setting ----
    ############################################################
    
    GetSetting <- function(column_name) {
      
      if (!column_name %in% names(df)) {
        return(NA)
      }
      
      values <- unique(df[[column_name]])
      values <- values[!is.na(values)]
      
      if (length(values) == 0) {
        return(NA)
      }
      
      if (length(values) > 1) {
        warning(
          paste0(
            "More than one value of ",
            column_name,
            " was found in file: ",
            basename(file),
            ". The first value will be used."
          )
        )
      }
      
      values[1]
    }
    
    
    ############################################################
    ## 2.2 Helper: safely calculate a column mean ----
    ############################################################
    
    MeanColumn <- function(column_name) {
      
      if (!column_name %in% names(df)) {
        return(NA_real_)
      }
      
      values <- suppressWarnings(
        as.numeric(df[[column_name]])
      )
      
      values <- values[is.finite(values)]
      
      if (length(values) == 0) {
        return(NA_real_)
      }
      
      mean(values)
    }
    
    
    ############################################################
    ## 2.3 Helper: safely calculate seed minimum and maximum ----
    ############################################################
    
    SeedRange <- function() {
      
      if (!"seed" %in% names(df)) {
        return(
          list(
            minimum = NA_real_,
            maximum = NA_real_
          )
        )
      }
      
      seed_values <- suppressWarnings(
        as.numeric(df$seed)
      )
      
      seed_values <- seed_values[is.finite(seed_values)]
      
      if (length(seed_values) == 0) {
        return(
          list(
            minimum = NA_real_,
            maximum = NA_real_
          )
        )
      }
      
      list(
        minimum = min(seed_values),
        maximum = max(seed_values)
      )
    }
    
    
    ############################################################
    ## 2.4 Helper: calculate ASE and coverage probability ----
    ############################################################
    
    SummariseStandardError <- function(
    estimate_column,
    standard_error_column,
    true_value
    ) {
      
      if (!standard_error_column %in% names(df)) {
        return(
          list(
            ASE = NA_real_,
            CP = NA_real_,
            n_valid = 0
          )
        )
      }
      
      estimates <- suppressWarnings(
        as.numeric(df[[estimate_column]])
      )
      
      standard_errors <- suppressWarnings(
        as.numeric(df[[standard_error_column]])
      )
      
      valid <- is.finite(estimates) &
        is.finite(standard_errors) &
        standard_errors >= 0
      
      if (!any(valid)) {
        return(
          list(
            ASE = NA_real_,
            CP = NA_real_,
            n_valid = 0
          )
        )
      }
      
      estimates <- estimates[valid]
      standard_errors <- standard_errors[valid]
      
      lower_bound <- estimates - z_value * standard_errors
      upper_bound <- estimates + z_value * standard_errors
      
      coverage <- (
        lower_bound <= true_value &
          upper_bound >= true_value
      )
      
      list(
        ASE = mean(standard_errors),
        CP = mean(coverage),
        n_valid = length(estimates)
      )
    }
    
    
    ############################################################
    ## 2.5 Simulation settings and timing summaries ----
    ############################################################
    
    model_name <- GetSetting("model")
    method_name <- GetSetting("method")
    
    N_value <- GetSetting("N")
    k_N_value <- GetSetting("k_N")
    m_N_value <- GetSetting("m_N")
    alpha_value <- GetSetting("alpha")
    alpha_N_value <- GetSetting("alpha_N")
    
    mean_bad_draws <- MeanColumn("bad_draws")
    mean_time_subsampling_seconds <- MeanColumn("time_subsampling_seconds")
    mean_time_estimation_seconds <- MeanColumn("time_estimation_seconds")
    mean_time_sse_seconds <- MeanColumn("time_sse_seconds")
    mean_time_asd_seconds <- MeanColumn("time_asd_seconds")
    
    seed_range <- SeedRange()
    
    
    ############################################################
    ## 2.6 Summarise each parameter ----
    ############################################################
    
    rows <- vector(
      mode = "list",
      length = length(estimate_cols)
    )
    
    for (j in seq_along(estimate_cols)) {
      
      estimate_col <- estimate_cols[j]
      
      parameter_name <- sub(
        pattern = "^estimate_",
        replacement = "",
        x = estimate_col
      )
      
      true_value <- theta[j]
      
      estimates <- suppressWarnings(
        as.numeric(df[[estimate_col]])
      )
      
      estimates <- estimates[is.finite(estimates)]
      
      if (length(estimates) == 0) {
        next
      }
      
      
      ############################################################
      ## Point-estimation performance ----
      ############################################################
      
      mean_estimate <- mean(estimates)
      bias <- mean_estimate - true_value
      
      if (length(estimates) >= 2) {
        sd_value <- sd(estimates)
      } else {
        sd_value <- NA_real_
      }
      
      rmse <- sqrt(
        mean(
          (estimates - true_value)^2
        )
      )
      
      
      ############################################################
      ## ASD-based inference ----
      ############################################################
      
      asd_summary <- SummariseStandardError(
        estimate_column = estimate_col,
        standard_error_column = paste0(
          "asd_",
          parameter_name
        ),
        true_value = true_value
      )
      
      
      ############################################################
      ## Adjusted-ASD-based inference ----
      ############################################################
      
      adjusted_asd_summary <- SummariseStandardError(
        estimate_column = estimate_col,
        standard_error_column = paste0(
          "adjusted_asd_",
          parameter_name
        ),
        true_value = true_value
      )
      
      
      ############################################################
      ## SSE-based inference ----
      ############################################################
      
      sse_summary <- SummariseStandardError(
        estimate_column = estimate_col,
        standard_error_column = paste0(
          "sse_",
          parameter_name
        ),
        true_value = true_value
      )
      
      
      ############################################################
      ## Adjusted-SSE-based inference ----
      ############################################################
      
      adjusted_sse_summary <- SummariseStandardError(
        estimate_column = estimate_col,
        standard_error_column = paste0(
          "adjusted_sse_",
          parameter_name
        ),
        true_value = true_value
      )
      
      
      ############################################################
      ## Create one summary row ----
      ############################################################
      
      rows[[j]] <- data.frame(
        model = model_name,
        N = N_value,
        k_N = k_N_value,
        m_N = m_N_value,
        alpha = alpha_value,
        alpha_N = alpha_N_value,
        method = method_name,
        parameter_index = j,
        parameter = parameter_name,
        true_value = true_value,
        n_rep = nrow(df),
        n_valid_estimates = length(estimates),
        seed_min = seed_range$minimum,
        seed_max = seed_range$maximum,
        mean_estimate = mean_estimate,
        BIAS = bias,
        SD = sd_value,
        RMSE = rmse,
        ASD = asd_summary$ASE,
        CP_asd = asd_summary$CP,
        adjusted_ASD = adjusted_asd_summary$ASE,
        CP_adjusted_asd = adjusted_asd_summary$CP,
        SSE = sse_summary$ASE,
        CP_sse = sse_summary$CP,
        adjusted_SSE = adjusted_sse_summary$ASE,
        CP_adjusted_sse = adjusted_sse_summary$CP,
        mean_bad_draws = mean_bad_draws,
        mean_time_subsampling_seconds = mean_time_subsampling_seconds,
        mean_time_estimation_seconds = mean_time_estimation_seconds,
        mean_time_sse_seconds = mean_time_sse_seconds,
        mean_time_asd_seconds = mean_time_asd_seconds,
        source_file = basename(file),
        check.names = FALSE
      )
    }
    
    rows <- rows[
      !vapply(
        rows,
        is.null,
        logical(1)
      )
    ]
    
    if (length(rows) == 0) {
      return(NULL)
    }
    
    do.call(rbind, rows)
  }
  
  
  ############################################################
  ## 3. Apply the summary function to every file ----
  ############################################################
  
  all_rows <- vector(
    mode = "list",
    length = length(files)
  )
  
  for (i in seq_along(files)) {
    
    cat(
      "Processing",
      i,
      "out of",
      length(files),
      ":",
      basename(files[i]),
      "\n"
    )
    
    all_rows[[i]] <- SummariseOneFile(files[i])
  }
  
  all_rows <- all_rows[
    !vapply(
      all_rows,
      is.null,
      logical(1)
    )
  ]
  
  if (length(all_rows) == 0) {
    stop("No valid simulation results could be summarised.")
  }
  
  summary_df <- do.call(
    rbind,
    all_rows
  )
  
  
  ############################################################
  ## 4. Sort the summary table ----
  ############################################################
  
  summary_df <- summary_df[
    order(
      summary_df$model,
      summary_df$N,
      summary_df$alpha,
      summary_df$k_N,
      summary_df$m_N,
      summary_df$method,
      summary_df$parameter_index
    ),
  ]
  
  rownames(summary_df) <- NULL
  
  
  ############################################################
  ## 5. Save the summary CSV file ----
  ############################################################
  
  output_path <- file.path(
    folder,
    output_file
  )
  
  write.csv(
    summary_df,
    file = output_path,
    row.names = FALSE
  )
  
  
  ############################################################
  ## 6. Completion message ----
  ############################################################
  
  cat("\nSummary completed.\n")
  cat("Files processed:", length(files), "\n")
  cat("Summary rows:", nrow(summary_df), "\n")
  cat("Output saved to:", output_path, "\n")
  
  summary_df
}


############################################################
## Implementation ----
############################################################

summary_df <- SummariseSimulationFolder(
  folder = "./Subbagging new",
  theta = c(1, 2),
  output_file = "all_simulation_results_summary.csv",
  z_value = 1.96
)

summary_df
