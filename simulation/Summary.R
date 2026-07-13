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

  files <- files[!grepl("summary", basename(files), ignore.case = TRUE)]

  if (length(files) == 0) {
    stop("No simulation CSV files found in this folder.")
  }

  ############################################################
  ## Summarise one file
  ############################################################
  summarise_one_file <- function(file) {

    df <- read.csv(file, check.names = FALSE)

    estimate_cols <- grep("^estimate_", names(df), value = TRUE)

    if (length(estimate_cols) != length(theta)) {
      stop(
        paste0(
          "Number of estimate columns does not match length(theta) in file: ",
          basename(file),
          "\nNumber of estimate columns = ", length(estimate_cols),
          "\nlength(theta) = ", length(theta)
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

      values <- suppressWarnings(as.numeric(df[[column_name]]))
      values <- values[is.finite(values)]

      if (length(values) == 0) {
        return(NA_real_)
      }

      mean(values)
    }

    ############################################################
    ## Timing summaries across replications
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
    ## Mean total time for the method represented by this file
    ##
    ## simple:
    ##   subsampling + simple
    ##
    ## bc_bias:
    ##   subsampling + simple + bias correction
    ##
    ## bc_equation:
    ##   subsampling + simple + adjusted-equation solution
    ############################################################
    method_name <- as.character(unique(df$method)[1])

    if (method_name == "simple") {

      required_times <- c(
        mean_time_subsampling_seconds,
        mean_time_simple_seconds
      )

    } else if (method_name == "bc_bias") {

      required_times <- c(
        mean_time_subsampling_seconds,
        mean_time_simple_seconds,
        mean_time_bc_bias_seconds
      )

    } else if (method_name == "bc_equation") {

      required_times <- c(
        mean_time_subsampling_seconds,
        mean_time_simple_seconds,
        mean_time_bc_equation_seconds
      )

    } else {

      required_times <- numeric(0)
    }

    if (length(required_times) > 0 &&
        all(is.finite(required_times))) {
      mean_method_total_seconds <- sum(required_times)
    } else {
      mean_method_total_seconds <- NA_real_
    }

    rows <- list()

    for (j in seq_along(estimate_cols)) {

      estimate_col <- estimate_cols[j]
      se_col <- paste0("se_", sub("^estimate_", "", estimate_col))

      true_value <- theta[j]

      estimates_all <- df[[estimate_col]]
      estimates <- estimates_all[!is.na(estimates_all)]

      if (length(estimates) == 0) {
        next
      }

      mean_estimate <- mean(estimates)
      bias <- mean_estimate - true_value
      sd_value <- sd(estimates)
      rmse <- sqrt(mean((estimates - true_value)^2))

      if (se_col %in% names(df)) {

        se_values <- df[[se_col]]
        ase <- mean(se_values, na.rm = TRUE)

        lower <- df[[estimate_col]] - z_value * df[[se_col]]
        upper <- df[[estimate_col]] + z_value * df[[se_col]]

        valid <- !is.na(lower) & !is.na(upper)

        if (sum(valid) > 0) {
          cp <- mean(
            lower[valid] <= true_value &
              upper[valid] >= true_value
          )
        } else {
          cp <- NA
        }

      } else {
        ase <- NA
        cp <- NA
      }

      ############################################################
      ## Bad draw summaries
      ############################################################
      if ("bad_draws" %in% names(df)) {
        total_bad_draws <- sum(df$bad_draws, na.rm = TRUE)
        mean_bad_draws <- mean(df$bad_draws, na.rm = TRUE)
      } else {
        total_bad_draws <- NA
        mean_bad_draws <- NA
      }

      if ("total_draws" %in% names(df)) {
        total_draws_sum <- sum(df$total_draws, na.rm = TRUE)
        mean_total_draws <- mean(df$total_draws, na.rm = TRUE)
      } else {
        total_draws_sum <- NA
        mean_total_draws <- NA
      }

      if (!is.na(total_bad_draws) &&
          !is.na(total_draws_sum) &&
          total_draws_sum > 0) {
        bad_draw_rate <- total_bad_draws / total_draws_sum
      } else {
        bad_draw_rate <- NA
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
        seed_min = min(df$seed, na.rm = TRUE),
        seed_max = max(df$seed, na.rm = TRUE),
        mean_estimate = mean_estimate,
        BIAS = bias,
        SD = sd_value,
        ASE = ase,
        RMSE = rmse,
        CP = cp,
        total_bad_draws = total_bad_draws,
        mean_bad_draws = mean_bad_draws,
        total_draws_sum = total_draws_sum,
        mean_total_draws = mean_total_draws,
        bad_draw_rate = bad_draw_rate,
        mean_time_subsampling_seconds =
          mean_time_subsampling_seconds,
        mean_time_simple_seconds =
          mean_time_simple_seconds,
        mean_time_bc_bias_seconds =
          mean_time_bc_bias_seconds,
        mean_time_bc_equation_seconds =
          mean_time_bc_equation_seconds,
        mean_method_total_seconds =
          mean_method_total_seconds,
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
    cat("Processing:", basename(files[i]), "\n")
    all_rows[[i]] <- summarise_one_file(files[i])
  }

  all_rows <- all_rows[!vapply(all_rows, is.null, logical(1))]

  if (length(all_rows) == 0) {
    stop("No valid simulation results could be summarised.")
  }

  summary_df <- do.call(rbind, all_rows)

  ############################################################
  ## Sort output
  ############################################################
  summary_df <- summary_df[order(
    summary_df$model,
    summary_df$N,
    summary_df$alpha,
    summary_df$k_N,
    summary_df$m_N,
    summary_df$method,
    summary_df$parameter_index
  ), ]

  rownames(summary_df) <- NULL

  ############################################################
  ## Save output
  ############################################################
  output_path <- file.path(folder, output_file)
  write.csv(summary_df, output_path, row.names = FALSE)

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
