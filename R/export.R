#' Export spike data to CSV
#'
#' Exports spike train data to CSV format for use in other analysis tools.
#' Supports multiple channels and optional metadata inclusion.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}}
#' @param output_file Character. Output file path (must end in .csv)
#' @param include_metadata Logical. Whether to include metadata columns like X, Y coordinates
#'   (default: TRUE)
#' @return Invisibly returns the output file path
#' @import rlang
#' @importFrom readr write_csv
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Export to CSV
#' exportToCSV(spikes, "spike_data.csv")
#'
#' h5$close()
#' }
exportToCSV <- function(spike_data, output_file, include_metadata = TRUE) {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  if (!is.character(output_file) || length(output_file) != 1) {
    rlang::abort("output_file must be a single character string")
  }

  if (!grepl("\\.csv$", output_file, ignore.case = TRUE)) {
    rlang::abort("output_file must have .csv extension")
  }

  # Prepare data for export
  if (!include_metadata) {
    # Keep only essential columns
    essential_cols <- c("spike_times", "spike_chid", "spike_units", "peak_amplitude")
    export_data <- spike_data[, essential_cols[essential_cols %in% names(spike_data)], drop = FALSE]
  } else {
    export_data <- spike_data
  }

  # Write to CSV
  tryCatch({
    readr::write_csv(export_data, output_file)
    message("Successfully exported ", nrow(export_data), " spikes to ", output_file)
  }, error = function(e) {
    rlang::abort(paste("Failed to write CSV file:", e$message))
  })

  invisible(output_file)
}


#' Export data to Parquet format
#'
#' Exports MEA data to Apache Parquet format for efficient storage and
#' compatibility with big data tools. Parquet provides better compression
#' and faster read speeds than CSV.
#'
#' @param data data.frame. Data to export (spike data, burst data, etc.)
#' @param output_file Character. Output file path (must end in .parquet)
#' @param compression Character. Compression algorithm: "snappy", "gzip", "brotli"
#'   (default: "snappy")
#' @return Invisibly returns the output file path
#' @import rlang arrow
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Export to Parquet
#' exportToParquet(spikes, "spike_data.parquet")
#'
#' h5$close()
#' }
exportToParquet <- function(data, output_file, compression = c("snappy", "gzip", "brotli")) {
  # Validate inputs
  if (!is.data.frame(data)) {
    rlang::abort("data must be a data.frame")
  }

  if (!is.character(output_file) || length(output_file) != 1) {
    rlang::abort("output_file must be a single character string")
  }

  if (!grepl("\\.parquet$", output_file, ignore.case = TRUE)) {
    rlang::abort("output_file must have .parquet extension")
  }

  compression <- match.arg(compression)

  # Write to Parquet
  tryCatch({
    arrow::write_parquet(data, output_file, compression = compression)
    message("Successfully exported ", nrow(data), " rows to ", output_file)
  }, error = function(e) {
    rlang::abort(paste("Failed to write Parquet file:", e$message))
  })

  invisible(output_file)
}


#' Export analysis summary report
#'
#' Generates a comprehensive analysis summary report in text format including
#' spike statistics, burst metrics, and network properties.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}}
#' @param burst_data data.frame. Burst data from \code{\link{bxrSpikeBursts}} or
#'   \code{\link{detectBursts}} (optional)
#' @param output_file Character. Output file path (default: "analysis_report.txt")
#' @param include_plots Logical. Whether to generate and save plots (default: FALSE)
#' @param plot_dir Character. Directory for saving plots if include_plots=TRUE
#'   (default: "plots")
#' @return Invisibly returns the report content as a character vector
#' @import rlang
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#' bursts <- bxrSpikeBursts(h5)
#'
#' # Generate summary report
#' exportSummaryReport(spikes, bursts, output_file = "report.txt")
#'
#' h5$close()
#' }
exportSummaryReport <- function(spike_data, burst_data = NULL,
                                output_file = "analysis_report.txt",
                                include_plots = FALSE, plot_dir = "plots") {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  report <- c()

  # Header
  report <- c(report, "=" %>% rep(70) %>% paste(collapse = ""))
  report <- c(report, "MEA ANALYSIS SUMMARY REPORT")
  report <- c(report, paste("Generated:", Sys.time()))
  report <- c(report, "=" %>% rep(70) %>% paste(collapse = ""))
  report <- c(report, "")

  # Spike Statistics
  report <- c(report, "1. SPIKE TRAIN STATISTICS")
  report <- c(report, "-" %>% rep(70) %>% paste(collapse = ""))

  tryCatch({
    stats <- getSpikeStatistics(spike_data)

    report <- c(report, sprintf("Total spikes: %d", stats$total_spikes))
    report <- c(report, sprintf("Active channels: %d", stats$n_active_channels))
    report <- c(report, sprintf("Recording duration: %.2f seconds", stats$recording_duration))
    report <- c(report, sprintf("Mean firing rate: %.2f Hz", stats$mean_firing_rate))
    report <- c(report, sprintf("Median firing rate: %.2f Hz", stats$median_firing_rate))
    report <- c(report, sprintf("Global firing rate: %.2f Hz", stats$global_firing_rate))
    report <- c(report, sprintf("Most active channel: %d (%.2f Hz)",
                               stats$most_active_channel, stats$max_firing_rate))
  }, error = function(e) {
    report <<- c(report, "Error calculating spike statistics")
  })

  report <- c(report, "")

  # ISI Statistics
  report <- c(report, "2. INTER-SPIKE INTERVAL (ISI) STATISTICS")
  report <- c(report, "-" %>% rep(70) %>% paste(collapse = ""))

  tryCatch({
    isi_stats <- calculateISI(spike_data)

    report <- c(report, sprintf("Mean ISI: %.4f s (%.2f ms)",
                               mean(isi_stats$mean_isi, na.rm = TRUE),
                               mean(isi_stats$mean_isi, na.rm = TRUE) * 1000))
    report <- c(report, sprintf("Median ISI: %.4f s (%.2f ms)",
                               mean(isi_stats$median_isi, na.rm = TRUE),
                               mean(isi_stats$median_isi, na.rm = TRUE) * 1000))
    report <- c(report, sprintf("Mean CV: %.3f", mean(isi_stats$cv_isi, na.rm = TRUE)))
  }, error = function(e) {
    report <<- c(report, "Error calculating ISI statistics")
  })

  report <- c(report, "")

  # Burst Statistics (if provided)
  if (!is.null(burst_data)) {
    report <- c(report, "3. BURST STATISTICS")
    report <- c(report, "-" %>% rep(70) %>% paste(collapse = ""))

    tryCatch({
      burst_stats <- calculateBurstStatistics(burst_data)

      report <- c(report, sprintf("Total bursts: %d", nrow(burst_data)))
      report <- c(report, sprintf("Channels with bursts: %d", nrow(burst_stats)))
      report <- c(report, sprintf("Mean burst duration: %.4f s (%.2f ms)",
                                 mean(burst_stats$mean_burst_duration),
                                 mean(burst_stats$mean_burst_duration) * 1000))
      report <- c(report, sprintf("Mean inter-burst interval: %.4f s",
                                 mean(burst_stats$mean_ibi, na.rm = TRUE)))
      report <- c(report, sprintf("Mean burst frequency: %.2f bursts/min",
                                 mean(burst_stats$burst_frequency)))
    }, error = function(e) {
      report <<- c(report, "Error calculating burst statistics")
    })

    report <- c(report, "")
  }

  # Firing Rate Distribution
  report <- c(report, "4. FIRING RATE DISTRIBUTION")
  report <- c(report, "-" %>% rep(70) %>% paste(collapse = ""))

  tryCatch({
    rates <- calculateFiringRate(spike_data)

    quantiles <- quantile(rates$firing_rate, probs = c(0.25, 0.5, 0.75, 0.9, 0.95))
    report <- c(report, sprintf("25th percentile: %.2f Hz", quantiles[1]))
    report <- c(report, sprintf("50th percentile (median): %.2f Hz", quantiles[2]))
    report <- c(report, sprintf("75th percentile: %.2f Hz", quantiles[3]))
    report <- c(report, sprintf("90th percentile: %.2f Hz", quantiles[4]))
    report <- c(report, sprintf("95th percentile: %.2f Hz", quantiles[5]))
  }, error = function(e) {
    report <<- c(report, "Error calculating firing rate distribution")
  })

  report <- c(report, "")
  report <- c(report, "=" %>% rep(70) %>% paste(collapse = ""))
  report <- c(report, "END OF REPORT")
  report <- c(report, "=" %>% rep(70) %>% paste(collapse = ""))

  # Write report to file
  tryCatch({
    writeLines(report, output_file)
    message("Report saved to: ", output_file)
  }, error = function(e) {
    rlang::abort(paste("Failed to write report file:", e$message))
  })

  # Generate plots if requested
  if (include_plots) {
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir, recursive = TRUE)
    }

    message("Generating plots in: ", plot_dir)

    # Save raster plot
    tryCatch({
      p <- plotSpikeRaster(spike_data, time_range = c(0, min(60, max(spike_data$spike_times))))
      ggplot2::ggsave(file.path(plot_dir, "raster_plot.png"), p, width = 10, height = 6)
    }, error = function(e) {
      warning("Could not generate raster plot")
    })

    # Save firing rate heatmap if spatial data available
    if (all(c("X", "Y") %in% names(spike_data))) {
      tryCatch({
        p <- plotFiringRateHeatmap(spike_data)
        ggplot2::ggsave(file.path(plot_dir, "firing_rate_heatmap.png"), p, width = 8, height = 8)
      }, error = function(e) {
        warning("Could not generate heatmap")
      })
    }

    message("Plots saved to: ", plot_dir)
  }

  invisible(report)
}


#' Batch export multiple data types
#'
#' Exports multiple MEA data types (spikes, bursts, field potentials) to
#' specified format in a single operation.
#'
#' @param data_list Named list. List of data frames to export (e.g., list(spikes=..., bursts=...))
#' @param output_dir Character. Output directory for files
#' @param format Character. Output format: "csv" or "parquet" (default: "csv")
#' @param prefix Character. Prefix for output file names (default: "mea_data")
#' @return Invisibly returns vector of output file paths
#' @import rlang
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#' bursts <- bxrSpikeBursts(h5)
#'
#' # Batch export
#' data_list <- list(spikes = spikes, bursts = bursts)
#' batchExport(data_list, output_dir = "exported_data", format = "parquet")
#'
#' h5$close()
#' }
batchExport <- function(data_list, output_dir, format = c("csv", "parquet"),
                        prefix = "mea_data") {
  # Validate inputs
  if (!is.list(data_list) || is.null(names(data_list))) {
    rlang::abort("data_list must be a named list")
  }

  if (!is.character(output_dir) || length(output_dir) != 1) {
    rlang::abort("output_dir must be a single character string")
  }

  format <- match.arg(format)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  output_files <- character(length(data_list))
  names(output_files) <- names(data_list)

  # Export each data frame
  for (name in names(data_list)) {
    data <- data_list[[name]]

    if (!is.data.frame(data)) {
      warning("Skipping '", name, "': not a data frame")
      next
    }

    # Generate output filename
    if (format == "csv") {
      output_file <- file.path(output_dir, paste0(prefix, "_", name, ".csv"))
      exportToCSV(data, output_file)
    } else {
      output_file <- file.path(output_dir, paste0(prefix, "_", name, ".parquet"))
      exportToParquet(data, output_file)
    }

    output_files[name] <- output_file
  }

  message("Batch export complete. ", length(output_files), " files exported.")

  invisible(output_files)
}
