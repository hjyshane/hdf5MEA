#' Detect inactive/dead channels
#'
#' Identifies channels with insufficient activity based on spike count and
#' firing rate thresholds. Dead channels may indicate hardware issues or
#' poor electrode contact.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}}
#' @param min_spike_count Integer. Minimum number of spikes required (default: 10)
#' @param min_firing_rate Numeric. Minimum firing rate in Hz (default: 0.1)
#' @param all_channels Integer vector. All expected channel IDs. If NULL, uses
#'   only channels present in data (default: NULL)
#' @return List containing:
#'   \describe{
#'     \item{dead_channels}{Vector of channel IDs classified as dead}
#'     \item{low_activity_channels}{Vector of channel IDs with low activity}
#'     \item{channel_summary}{data.frame with per-channel statistics}
#'     \item{n_dead}{Number of dead channels}
#'     \item{n_low_activity}{Number of low activity channels}
#'   }
#' @import dplyr rlang
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Detect dead channels
#' qc <- detectDeadChannels(spikes, min_spike_count = 10, min_firing_rate = 0.1)
#' print(qc$n_dead)
#'
#' h5$close()
#' }
detectDeadChannels <- function(spike_data, min_spike_count = 10,
                               min_firing_rate = 0.1, all_channels = NULL) {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  if (!all(c("spike_times", "spike_chid") %in% names(spike_data))) {
    rlang::abort("spike_data must contain 'spike_times' and 'spike_chid' columns")
  }

  if (!is.numeric(min_spike_count) || min_spike_count < 0) {
    rlang::abort("min_spike_count must be non-negative")
  }

  if (!is.numeric(min_firing_rate) || min_firing_rate < 0) {
    rlang::abort("min_firing_rate must be non-negative")
  }

  # Calculate duration
  duration <- max(spike_data$spike_times) - min(spike_data$spike_times)

  # Calculate per-channel statistics
  channel_stats <- spike_data %>%
    dplyr::group_by(spike_chid) %>%
    dplyr::summarise(
      spike_count = dplyr::n(),
      firing_rate = dplyr::n() / duration,
      .groups = "drop"
    ) %>%
    dplyr::rename(channel_id = spike_chid)

  # Identify dead and low activity channels
  dead_channels <- channel_stats$channel_id[channel_stats$spike_count < min_spike_count]
  low_activity_channels <- channel_stats$channel_id[
    channel_stats$spike_count >= min_spike_count &
    channel_stats$firing_rate < min_firing_rate
  ]

  # Add completely silent channels if all_channels provided
  if (!is.null(all_channels)) {
    active_channels <- channel_stats$channel_id
    silent_channels <- setdiff(all_channels, active_channels)

    # Add silent channels to summary
    if (length(silent_channels) > 0) {
      silent_stats <- data.frame(
        channel_id = silent_channels,
        spike_count = 0,
        firing_rate = 0
      )
      channel_stats <- rbind(channel_stats, silent_stats)
      dead_channels <- c(dead_channels, silent_channels)
    }
  }

  # Add classification to summary
  channel_stats$status <- "active"
  channel_stats$status[channel_stats$channel_id %in% low_activity_channels] <- "low_activity"
  channel_stats$status[channel_stats$channel_id %in% dead_channels] <- "dead"

  result <- list(
    dead_channels = sort(dead_channels),
    low_activity_channels = sort(low_activity_channels),
    channel_summary = channel_stats[order(channel_stats$channel_id), ],
    n_dead = length(dead_channels),
    n_low_activity = length(low_activity_channels),
    n_active = nrow(channel_stats) - length(dead_channels) - length(low_activity_channels)
  )

  return(result)
}


#' Calculate signal-to-noise ratio (SNR)
#'
#' Estimates SNR for each channel based on spike amplitudes and background
#' noise levels. High SNR indicates good recording quality.
#'
#' @param spike_data data.frame. Spike data with amplitude information
#' @param voltage_data List. Raw voltage traces from \code{\link{brwtimeseriesConvert}}
#'   (optional, for noise estimation)
#' @param method Character. SNR calculation method: "amplitude" or "voltage" (default: "amplitude")
#' @return data.frame with columns:
#'   \describe{
#'     \item{channel_id}{Channel index}
#'     \item{signal_mean}{Mean signal amplitude}
#'     \item{noise_std}{Noise standard deviation}
#'     \item{snr_db}{Signal-to-noise ratio in dB}
#'     \item{quality}{Quality classification: "excellent", "good", "fair", "poor"}
#'   }
#' @import dplyr rlang
#' @importFrom stats sd mad median
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Calculate SNR from spike amplitudes
#' snr <- calculateSNR(spikes, method = "amplitude")
#'
#' h5$close()
#' }
calculateSNR <- function(spike_data, voltage_data = NULL,
                        method = c("amplitude", "voltage")) {
  method <- match.arg(method)

  if (method == "amplitude") {
    # Calculate SNR from spike amplitudes
    if (!is.data.frame(spike_data)) {
      rlang::abort("spike_data must be a data.frame")
    }

    if (!"peak_amplitude" %in% names(spike_data)) {
      rlang::abort("spike_data must contain 'peak_amplitude' column for amplitude method")
    }

    if (!"spike_chid" %in% names(spike_data)) {
      rlang::abort("spike_data must contain 'spike_chid' column")
    }

    # Calculate per-channel SNR
    snr_data <- spike_data %>%
      dplyr::group_by(spike_chid) %>%
      dplyr::summarise(
        signal_mean = mean(abs(peak_amplitude)),
        noise_std = stats::sd(abs(peak_amplitude)) / sqrt(2),  # Approximate noise from spike variability
        n_spikes = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::rename(channel_id = spike_chid)

  } else {  # voltage method
    if (is.null(voltage_data)) {
      rlang::abort("voltage_data required for voltage method")
    }

    if (!is.list(voltage_data)) {
      rlang::abort("voltage_data must be a list from brwtimeseriesConvert()")
    }

    snr_list <- list()

    for (ch_name in names(voltage_data)) {
      ch_data <- voltage_data[[ch_name]]

      if (!is.data.frame(ch_data) || !"voltage" %in% names(ch_data)) {
        next
      }

      voltage <- ch_data$voltage

      # Estimate noise using median absolute deviation (robust to spikes)
      noise_est <- stats::mad(voltage, center = stats::median(voltage))

      # Estimate signal from spike-like events (threshold crossings)
      threshold <- -4 * noise_est
      spikes <- voltage[voltage < threshold]

      if (length(spikes) > 0) {
        signal_mean <- mean(abs(spikes))
      } else {
        signal_mean <- 0
      }

      snr_list[[ch_name]] <- data.frame(
        channel_id = as.numeric(ch_name),
        signal_mean = signal_mean,
        noise_std = noise_est,
        n_samples = length(voltage)
      )
    }

    snr_data <- do.call(rbind, snr_list)
  }

  # Calculate SNR in dB
  snr_data$snr_db <- 20 * log10(snr_data$signal_mean / snr_data$noise_std)

  # Classify quality
  snr_data$quality <- cut(snr_data$snr_db,
                          breaks = c(-Inf, 6, 12, 20, Inf),
                          labels = c("poor", "fair", "good", "excellent"))

  return(as.data.frame(snr_data))
}


#' Check data completeness and integrity
#'
#' Performs comprehensive quality checks on MEA data including temporal
#' coverage, channel availability, and data consistency.
#'
#' @param h5 H5File object. Open file handle from \code{\link{openBRW}} or \code{\link{openBXR}}
#' @param spike_data data.frame. Spike data (optional)
#' @param expected_duration Numeric. Expected recording duration in seconds (optional)
#' @param expected_channels Integer. Expected number of channels (optional)
#' @return List containing:
#'   \describe{
#'     \item{file_info}{Basic file information}
#'     \item{temporal_coverage}{Recording duration and gaps}
#'     \item{channel_coverage}{Number and distribution of channels}
#'     \item{data_quality}{Overall quality assessment}
#'     \item{warnings}{Vector of warning messages}
#'     \item{passed}{Logical indicating if all checks passed}
#'   }
#' @import hdf5r rlang
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Check data completeness
#' qc <- checkDataCompleteness(h5, spike_data = spikes,
#'                               expected_duration = 600,
#'                               expected_channels = 2304)
#' print(qc$warnings)
#'
#' h5$close()
#' }
checkDataCompleteness <- function(h5, spike_data = NULL,
                                  expected_duration = NULL,
                                  expected_channels = NULL) {
  # Validate inputs
  if (!inherits(h5, "H5File")) {
    rlang::abort("h5 must be an H5File object")
  }

  warnings <- character(0)
  file_info <- list()

  # Get basic file attributes
  tryCatch({
    file_info$sampling_rate <- hdf5r::h5attr(h5, "SamplingRate")
    file_info$file_format <- hdf5r::h5attr(h5, "FileFormat")

    if ("MeaName" %in% names(hdf5r::h5attributes(h5))) {
      file_info$mea_name <- hdf5r::h5attr(h5, "MeaName")
    }
  }, error = function(e) {
    warnings <<- c(warnings, "Could not read file attributes")
  })

  # Check temporal coverage
  temporal_coverage <- list()

  if (!is.null(spike_data)) {
    if (is.data.frame(spike_data) && "spike_times" %in% names(spike_data)) {
      temporal_coverage$start_time <- min(spike_data$spike_times)
      temporal_coverage$end_time <- max(spike_data$spike_times)
      temporal_coverage$duration <- temporal_coverage$end_time - temporal_coverage$start_time

      # Check against expected duration
      if (!is.null(expected_duration)) {
        duration_diff <- abs(temporal_coverage$duration - expected_duration)
        if (duration_diff > expected_duration * 0.05) {  # More than 5% difference
          warnings <- c(warnings,
                       sprintf("Duration mismatch: expected %.1f s, got %.1f s",
                              expected_duration, temporal_coverage$duration))
        }
      }

      # Check for temporal gaps
      time_diffs <- diff(sort(spike_data$spike_times))
      max_gap <- max(time_diffs)
      temporal_coverage$max_gap <- max_gap

      if (max_gap > 1.0) {  # Gap larger than 1 second
        warnings <- c(warnings,
                     sprintf("Large temporal gap detected: %.2f seconds", max_gap))
      }
    }
  }

  # Check channel coverage
  channel_coverage <- list()

  if (!is.null(spike_data)) {
    if (is.data.frame(spike_data) && "spike_chid" %in% names(spike_data)) {
      active_channels <- unique(spike_data$spike_chid)
      channel_coverage$n_active_channels <- length(active_channels)
      channel_coverage$active_channels <- sort(active_channels)

      # Check against expected channels
      if (!is.null(expected_channels)) {
        if (length(active_channels) < expected_channels * 0.8) {  # Less than 80%
          warnings <- c(warnings,
                       sprintf("Low channel coverage: %d/%d active (%.1f%%)",
                              length(active_channels), expected_channels,
                              100 * length(active_channels) / expected_channels))
        }
      }

      # Check for channel clustering (potential hardware issues)
      channel_density <- hist(active_channels, breaks = 20, plot = FALSE)$counts
      if (max(channel_density) > 2 * mean(channel_density)) {
        warnings <- c(warnings, "Uneven channel distribution detected")
      }
    }
  }

  # Data quality assessment
  data_quality <- list()

  if (!is.null(spike_data)) {
    total_spikes <- nrow(spike_data)
    data_quality$total_spikes <- total_spikes

    if (total_spikes < 100) {
      warnings <- c(warnings, "Very low spike count (<100)")
      data_quality$overall <- "poor"
    } else if (total_spikes < 1000) {
      data_quality$overall <- "fair"
    } else if (total_spikes < 10000) {
      data_quality$overall <- "good"
    } else {
      data_quality$overall <- "excellent"
    }

    # Check spike rate consistency
    if ("spike_times" %in% names(spike_data)) {
      bin_size <- 10  # seconds
      time_range <- range(spike_data$spike_times)
      bins <- seq(time_range[1], time_range[2], by = bin_size)

      if (length(bins) > 2) {
        spike_counts <- hist(spike_data$spike_times, breaks = bins, plot = FALSE)$counts
        cv_spike_rate <- sd(spike_counts) / mean(spike_counts)

        data_quality$rate_cv <- cv_spike_rate

        if (cv_spike_rate > 2) {
          warnings <- c(warnings, "High variability in spike rate across recording")
        }
      }
    }
  }

  # Determine if all checks passed
  passed <- length(warnings) == 0

  result <- list(
    file_info = file_info,
    temporal_coverage = temporal_coverage,
    channel_coverage = channel_coverage,
    data_quality = data_quality,
    warnings = warnings,
    passed = passed
  )

  return(result)
}


#' Generate quality control summary
#'
#' Creates a comprehensive QC report combining multiple quality metrics
#' into a single summary with recommendations.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}}
#' @param h5 H5File object. Open file handle (optional)
#' @param voltage_data List. Raw voltage traces (optional)
#' @return List containing:
#'   \describe{
#'     \item{dead_channels}{Dead channel analysis}
#'     \item{snr}{Signal-to-noise analysis}
#'     \item{completeness}{Data completeness checks}
#'     \item{recommendations}{Vector of recommended actions}
#'     \item{overall_quality}{Overall quality score (0-100)}
#'   }
#' @import rlang
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Generate QC summary
#' qc_summary <- generateQCSummary(spikes, h5 = h5)
#' print(qc_summary$overall_quality)
#' print(qc_summary$recommendations)
#'
#' h5$close()
#' }
generateQCSummary <- function(spike_data, h5 = NULL, voltage_data = NULL) {
  recommendations <- character(0)
  quality_score <- 100

  # Check for dead channels
  dead_analysis <- detectDeadChannels(spike_data)

  if (dead_analysis$n_dead > 0) {
    recommendations <- c(recommendations,
                        sprintf("Remove or investigate %d dead channels", dead_analysis$n_dead))
    quality_score <- quality_score - min(20, dead_analysis$n_dead * 2)
  }

  if (dead_analysis$n_low_activity > 0) {
    recommendations <- c(recommendations,
                        sprintf("Review %d low-activity channels", dead_analysis$n_low_activity))
    quality_score <- quality_score - min(10, dead_analysis$n_low_activity)
  }

  # Calculate SNR
  snr_analysis <- NULL
  if ("peak_amplitude" %in% names(spike_data)) {
    snr_analysis <- calculateSNR(spike_data, voltage_data, method = "amplitude")

    poor_snr <- sum(snr_analysis$quality == "poor", na.rm = TRUE)
    if (poor_snr > 0) {
      recommendations <- c(recommendations,
                          sprintf("%d channels have poor SNR", poor_snr))
      quality_score <- quality_score - min(15, poor_snr * 2)
    }
  }

  # Check completeness
  completeness_analysis <- NULL
  if (!is.null(h5)) {
    completeness_analysis <- checkDataCompleteness(h5, spike_data = spike_data)

    if (!completeness_analysis$passed) {
      recommendations <- c(recommendations, completeness_analysis$warnings)
      quality_score <- quality_score - length(completeness_analysis$warnings) * 5
    }
  }

  # Ensure quality score is within bounds
  quality_score <- max(0, min(100, quality_score))

  # Add general recommendations based on quality
  if (quality_score < 50) {
    recommendations <- c(recommendations,
                        "Data quality is poor. Consider re-recording or reviewing experimental setup.")
  } else if (quality_score < 70) {
    recommendations <- c(recommendations,
                        "Data quality is fair. Some channels may need attention.")
  }

  result <- list(
    dead_channels = dead_analysis,
    snr = snr_analysis,
    completeness = completeness_analysis,
    recommendations = recommendations,
    overall_quality = quality_score
  )

  return(result)
}
