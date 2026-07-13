############################################################
## Summary simulation CSV files in one folder
############################################################

SummariseSimulationFolder <- function(folder = ".",
                                      theta = c(0, 1),
                                      output_file = "all_simulation_results_summary.csv",
                                      z_value = 1.96) {

  ############################################################
  ## Find simulation CSV files
  ############################################################
  files <- list.files(
    path = folder,
    pattern = "^simulation_.*\\.csv$",
    full.names = TRUE
  )

  files <- files[
    !grepl("summary", basename(files), ignore.case = TRUE)
  ]

  if (length(files) == 0) {
    stop("No simulation CSV files found in this folder.")
  }

  ############################################################
  ## Summarise one file
  ############################################################
  summarise_one_file <- function(file) {

    df <- read.csv(
      file,
      check.names = FALSE
    )

    estimate_cols <- grep(
      "^estimate_",
      names(df),
      value = TRUE
    )

    if (length(estimate_cols) != length(theta)) {
      stop(
        paste0(
          "Number of estimate columns does not match length(theta) in file: ",
          basename(file),
          "\nNumber of estimate columns = ",
          length(estimate_cols),
          "\nlength(theta) = ",
          length(theta)
        )
      )
    }

    ############################################################
    ## Helper: safely calculate the mean of a timing column
    ############################################################
    mean_time <- function(column_name) {

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
    ## Timing summaries across replications
    ##
    ## Each timing column is summarised separately.
    ## No timing columns are added together.
    ############################################################
    mean_time_subsampling_seconds <- mean_time(
      "time_subsampling_seconds"
    )

    mean_time_simple_seconds <- mean_time(
      "time_simple_seconds"
    )

    mean_time_bc_bias_seconds <- mean_time(
      "time_bc_bias_seconds"
    )

    mean_time_bc_equation_seconds <- mean_time(
      "time_bc_equation_seconds"
    )

    ############################################################
    ## Method represented by the current file
    ############################################################
    if (!"method" %in% names(df)) {
      stop(
        paste0(
          "Column 'method' is missing from file: ",
          basename(file)
        )
      )
    }

    method_values <- unique(
      as.character(df$method[!is.na(df$method)])
    )

    if (length(method_values) == 0) {
      stop(
        paste0(
          "No valid method value found in file: ",
          basename(file)
        )
      )
    }

    if (length(method_values) > 1) {
      stop(
        paste0(
          "More than one method found in file: ",
          basename(file),
          "\nMethods found: ",
          paste(method_values, collapse = ", ")
        )
      )
    }

    method_name <- method_values[1]

    if (!method_name %in% c(
      "simple",
      "bc_bias",
      "bc_equation"
    )) {
      stop(
        paste0(
          "Unknown method in file: ",
          basename(file),
          "\nMethod found: ",
          method_name
        )
      )
    }

    rows <- list()

    ############################################################
    ## Summarise each parameter
    ############################################################
    for (j in seq_along(estimate_cols)) {

      estimate_col <- estimate_cols[j]

      se_col <- paste0(
        "se_",
        sub("^estimate_", "", estimate_col)
      )

      true_value <- theta[j]

      estimates_all <- suppressWarnings(
        as.numeric(df[[estimate_col]])
      )

      estimates <- estimates_all[
        is.finite(estimates_all)
      ]

      if (length(estimates) == 0) {
        next
      }

      mean_estimate <- mean(estimates)
      bias <- mean_estimate - true_value

      if (length(estimates) >= 2) {
        sd_value <- sd(estimates)
      } else {
        sd_value <- NA_real_
      }

      rmse <- sqrt(
        mean((estimates - true_value)^2)
      )

      ############################################################
      ## ASE and coverage probability
      ############################################################
      if (se_col %in% names(df)) {

        se_values <- suppressWarnings(
          as.numeric(df[[se_col]])
        )

        valid_se <- is.finite(se_values)

        if (sum(valid_se) > 0) {
          ase <- mean(se_values[valid_se])
        } else {
          ase <- NA_real_
        }

        valid <- is.finite(estimates_all) &
          is.finite(se_values)

        if (sum(valid) > 0) {

          lower <- estimates_all[valid] -
            z_value * se_values[valid]

          upper <- estimates_all[valid] +
            z_value * se_values[valid]

          cp <- mean(
            lower <= true_value &
              upper >= true_value
          )

        } else {
          cp <- NA_real_
        }

      } else {
        ase <- NA_real_
        cp <- NA_real_
      }

      ############################################################
      ## Bad draw summaries
      ############################################################
      if ("bad_draws" %in% names(df)) {

        bad_draw_values <- suppressWarnings(
          as.numeric(df$bad_draws)
        )

        bad_draw_values <- bad_draw_values[
          is.finite(bad_draw_values)
        ]

        if (length(bad_draw_values) > 0) {
          mean_bad_draws <- mean(bad_draw_values)
        } else {
          mean_bad_draws <- NA_real_
        }

      } else {
        mean_bad_draws <- NA_real_
      }

      if ("total_draws" %in% names(df)) {

        total_draw_values <- suppressWarnings(
          as.numeric(df$total_draws)
        )

        total_draw_values <- total_draw_values[
          is.finite(total_draw_values)
        ]

        if (length(total_draw_values) > 0) {
          mean_total_draws <- mean(total_draw_values)
        } else {
          mean_total_draws <- NA_real_
        }

      } else {
        mean_total_draws <- NA_real_
      }

      ############################################################
      ## Seed range
      ############################################################
      if ("seed" %in% names(df)) {

        seed_values <- suppressWarnings(
          as.numeric(df$seed)
        )

        seed_values <- seed_values[
          is.finite(seed_values)
        ]

        if (length(seed_values) > 0) {
          seed_min <- min(seed_values)
          seed_max <- max(seed_values)
        } else {
          seed_min <- NA_real_
          seed_max <- NA_real_
        }

      } else {
        seed_min <- NA_real_
        seed_max <- NA_real_
      }

      ############################################################
      ## Create one summary row
      ############################################################
      row <- data.frame(
        model = unique(df$model)[1],
        N = unique(df$N)[1],
        k_N = unique(df$k_N)[1],
        m_N = unique(df$m_N)[1],
        alpha = unique(df$alpha)[1],
        method = method_name,
        parameter_index = j,
        true_value = true_value,
        n_rep = nrow(df),
        seed_min = seed_min,
        seed_max = seed_max,
        mean_estimate = mean_estimate,
        BIAS = bias,
        SD = sd_value,
        ASE = ase,
        RMSE = rmse,
        CP = cp,
        mean_bad_draws = mean_bad_draws,
        mean_total_draws = mean_total_draws,

        mean_time_subsampling_seconds =
          mean_time_subsampling_seconds,

        mean_time_simple_seconds =
          mean_time_simple_seconds,

        mean_time_bc_bias_seconds =
          mean_time_bc_bias_seconds,

        mean_time_bc_equation_seconds =
          mean_time_bc_equation_seconds,

        source_file = basename(file),
        check.names = FALSE
      )

      rows[[length(rows) + 1]] <- row
    }

    if (length(rows) == 0) {
      return(NULL)
    }

    do.call(rbind, rows)
  }

  ############################################################
  ## Apply to all files
  ############################################################
  all_rows <- list()

  for (i in seq_along(files)) {

    cat(
      "Processing:",
      basename(files[i]),
      "\n"
    )

    all_rows[[i]] <- summarise_one_file(files[i])
  }

  all_rows <- all_rows[
    !vapply(all_rows, is.null, logical(1))
  ]

  if (length(all_rows) == 0) {
    stop("No valid simulation results could be summarised.")
  }

  summary_df <- do.call(
    rbind,
    all_rows
  )

  ############################################################
  ## Sort output
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
  ## Save output
  ############################################################
  output_path <- file.path(
    folder,
    output_file
  )

  write.csv(
    summary_df,
    output_path,
    row.names = FALSE
  )

  cat("\nDONE.\n")
  cat("Files processed:", length(files), "\n")
  cat("Summary rows:", nrow(summary_df), "\n")
  cat("Output saved to:", output_path, "\n")

  summary_df
}


############################################################
## Implementation
############################################################

summary_df <- SummariseSimulationFolder(
  folder = "./Subbagging new",
  theta = c(1, 2),
  output_file = "all_simulation_results_summary.csv"
)

summary_df
