#' Calculate firing rate per channel
#'
#' Computes firing rate (spikes per second) for each channel with optional
#' time binning for temporal analysis. Supports both overall and binned rate calculations.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}} with columns:
#'   spike_times, spike_chid
#' @param bin_size Numeric. Time bin size in seconds for rate calculation.
#'   If NULL, calculates overall firing rate (default: NULL)
#' @param time_range Numeric vector of length 2. Time range [start, end] in seconds.
#'   If NULL, uses full spike data range (default: NULL)
#' @return data.frame with columns:
#'   \describe{
#'     \item{channel_id}{Channel index}
#'     \item{firing_rate}{Mean firing rate in Hz}
#'     \item{spike_count}{Total number of spikes}
#'     \item{duration}{Recording duration in seconds}
#'     \item{time_bin}{Time bin center (only if bin_size specified)}
#'   }
#' @import dplyr rlang
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Overall firing rate per channel
#' rates <- calculateFiringRate(spikes)
#'
#' # Firing rate with 1-second bins
#' rates_binned <- calculateFiringRate(spikes, bin_size = 1)
#'
#' h5$close()
#' }
calculateFiringRate <- function(spike_data, bin_size = NULL, time_range = NULL) {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  if (!all(c("spike_times", "spike_chid") %in% names(spike_data))) {
    rlang::abort("spike_data must contain 'spike_times' and 'spike_chid' columns")
  }

  if (nrow(spike_data) == 0) {
    rlang::abort("spike_data is empty")
  }

  # Determine time range
  if (is.null(time_range)) {
    time_range <- c(min(spike_data$spike_times), max(spike_data$spike_times))
  } else {
    if (!is.numeric(time_range) || length(time_range) != 2) {
      rlang::abort("time_range must be a numeric vector of length 2")
    }
    # Filter spike data to time range
    spike_data <- spike_data[spike_data$spike_times >= time_range[1] &
                               spike_data$spike_times <= time_range[2], ]
  }

  duration <- time_range[2] - time_range[1]

  if (duration <= 0) {
    rlang::abort("Invalid time range: duration must be positive")
  }

  # Calculate overall firing rate per channel
  if (is.null(bin_size)) {
    firing_rates <- spike_data %>%
      dplyr::group_by(spike_chid) %>%
      dplyr::summarise(
        spike_count = dplyr::n(),
        duration = duration,
        firing_rate = spike_count / duration,
        .groups = "drop"
      ) %>%
      dplyr::rename(channel_id = spike_chid)

    return(as.data.frame(firing_rates))

  } else {
    # Validate bin_size
    if (!is.numeric(bin_size) || length(bin_size) != 1 || bin_size <= 0) {
      rlang::abort("bin_size must be a positive numeric value")
    }

    # Create time bins
    bin_edges <- seq(time_range[1], time_range[2], by = bin_size)
    if (tail(bin_edges, 1) < time_range[2]) {
      bin_edges <- c(bin_edges, time_range[2])
    }

    # Assign spikes to bins
    spike_data$time_bin <- cut(spike_data$spike_times,
                                breaks = bin_edges,
                                labels = FALSE,
                                include.lowest = TRUE)

    # Calculate bin centers
    bin_centers <- (bin_edges[-1] + bin_edges[-length(bin_edges)]) / 2

    # Calculate firing rate per channel per bin
    firing_rates_binned <- spike_data %>%
      dplyr::group_by(spike_chid, time_bin) %>%
      dplyr::summarise(
        spike_count = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        time_bin_center = bin_centers[time_bin],
        firing_rate = spike_count / bin_size
      ) %>%
      dplyr::rename(channel_id = spike_chid) %>%
      dplyr::select(channel_id, time_bin_center, spike_count, firing_rate)

    return(as.data.frame(firing_rates_binned))
  }
}


#' Calculate inter-spike intervals (ISI)
#'
#' Computes inter-spike intervals for each channel with summary statistics.
#' ISI is a fundamental measure of spike train regularity and firing patterns.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}} with columns:
#'   spike_times, spike_chid
#' @param channel_id Integer. Specific channel to analyze. If NULL, analyzes all channels
#'   (default: NULL)
#' @return data.frame with columns:
#'   \describe{
#'     \item{channel_id}{Channel index}
#'     \item{isi}{Inter-spike interval in seconds (if single channel)}
#'     \item{mean_isi}{Mean ISI in seconds}
#'     \item{median_isi}{Median ISI in seconds}
#'     \item{sd_isi}{Standard deviation of ISI}
#'     \item{cv_isi}{Coefficient of variation (SD/mean)}
#'     \item{n_intervals}{Number of ISI intervals}
#'   }
#' @import dplyr rlang
#' @importFrom stats sd median
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # ISI for all channels (summary statistics)
#' isi_all <- calculateISI(spikes)
#'
#' # ISI for specific channel (full distribution)
#' isi_ch1 <- calculateISI(spikes, channel_id = 100)
#'
#' h5$close()
#' }
calculateISI <- function(spike_data, channel_id = NULL) {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  if (!all(c("spike_times", "spike_chid") %in% names(spike_data))) {
    rlang::abort("spike_data must contain 'spike_times' and 'spike_chid' columns")
  }

  if (nrow(spike_data) == 0) {
    rlang::abort("spike_data is empty")
  }

  # Filter to specific channel if requested
  if (!is.null(channel_id)) {
    if (!is.numeric(channel_id) || length(channel_id) != 1) {
      rlang::abort("channel_id must be a single numeric value")
    }

    spike_data <- spike_data[spike_data$spike_chid == channel_id, ]

    if (nrow(spike_data) < 2) {
      rlang::abort(paste("Channel", channel_id, "has fewer than 2 spikes"))
    }

    # Sort by time
    spike_data <- spike_data[order(spike_data$spike_times), ]

    # Calculate ISI
    isi_values <- diff(spike_data$spike_times)

    # Return full ISI distribution for single channel
    result <- data.frame(
      channel_id = channel_id,
      isi = isi_values,
      mean_isi = mean(isi_values),
      median_isi = stats::median(isi_values),
      sd_isi = stats::sd(isi_values),
      cv_isi = stats::sd(isi_values) / mean(isi_values),
      n_intervals = length(isi_values)
    )

    return(result)

  } else {
    # Calculate ISI statistics for all channels
    isi_stats <- spike_data %>%
      dplyr::arrange(spike_chid, spike_times) %>%
      dplyr::group_by(spike_chid) %>%
      dplyr::summarise(
        n_spikes = dplyr::n(),
        mean_isi = if (dplyr::n() > 1) mean(diff(spike_times)) else NA_real_,
        median_isi = if (dplyr::n() > 1) stats::median(diff(spike_times)) else NA_real_,
        sd_isi = if (dplyr::n() > 1) stats::sd(diff(spike_times)) else NA_real_,
        cv_isi = if (dplyr::n() > 1) stats::sd(diff(spike_times)) / mean(diff(spike_times)) else NA_real_,
        n_intervals = dplyr::n() - 1,
        .groups = "drop"
      ) %>%
      dplyr::rename(channel_id = spike_chid)

    return(as.data.frame(isi_stats))
  }
}


#' Get spike train statistics
#'
#' Computes comprehensive summary statistics for spike trains including
#' firing rates, spike counts, active channels, and temporal coverage.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}} with columns:
#'   spike_times, spike_chid
#' @return List containing:
#'   \describe{
#'     \item{total_spikes}{Total number of spikes across all channels}
#'     \item{n_active_channels}{Number of channels with at least one spike}
#'     \item{recording_duration}{Total recording duration in seconds}
#'     \item{mean_firing_rate}{Mean firing rate across all active channels (Hz)}
#'     \item{median_firing_rate}{Median firing rate across active channels (Hz)}
#'     \item{global_firing_rate}{Overall firing rate (total spikes / duration) (Hz)}
#'     \item{spike_density}{Spikes per channel per second}
#'     \item{most_active_channel}{Channel ID with highest firing rate}
#'     \item{max_firing_rate}{Maximum firing rate across all channels (Hz)}
#'   }
#' @import dplyr rlang
#' @importFrom stats median
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Get comprehensive statistics
#' stats <- getSpikeStatistics(spikes)
#' print(stats)
#'
#' h5$close()
#' }
getSpikeStatistics <- function(spike_data) {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  if (!all(c("spike_times", "spike_chid") %in% names(spike_data))) {
    rlang::abort("spike_data must contain 'spike_times' and 'spike_chid' columns")
  }

  if (nrow(spike_data) == 0) {
    rlang::abort("spike_data is empty")
  }

  # Calculate duration
  duration <- max(spike_data$spike_times) - min(spike_data$spike_times)

  if (duration <= 0) {
    rlang::abort("Invalid spike data: all spikes at same time")
  }

  # Calculate per-channel firing rates
  channel_rates <- spike_data %>%
    dplyr::group_by(spike_chid) %>%
    dplyr::summarise(
      spike_count = dplyr::n(),
      firing_rate = spike_count / duration,
      .groups = "drop"
    )

  # Find most active channel
  max_rate_idx <- which.max(channel_rates$firing_rate)

  # Compile statistics
  stats <- list(
    total_spikes = nrow(spike_data),
    n_active_channels = length(unique(spike_data$spike_chid)),
    recording_duration = duration,
    mean_firing_rate = mean(channel_rates$firing_rate),
    median_firing_rate = stats::median(channel_rates$firing_rate),
    global_firing_rate = nrow(spike_data) / duration,
    spike_density = nrow(spike_data) / (length(unique(spike_data$spike_chid)) * duration),
    most_active_channel = channel_rates$spike_chid[max_rate_idx],
    max_firing_rate = channel_rates$firing_rate[max_rate_idx]
  )

  return(stats)
}


#' Calculate spike time autocorrelation
#'
#' Computes autocorrelation of spike times for assessing temporal patterns
#' and rhythmicity in spike trains.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}}
#' @param channel_id Integer. Channel to analyze
#' @param max_lag Numeric. Maximum lag in seconds (default: 1)
#' @param bin_size Numeric. Bin size for autocorrelation in seconds (default: 0.001)
#' @return data.frame with columns:
#'   \describe{
#'     \item{lag}{Time lag in seconds}
#'     \item{autocorr}{Autocorrelation value}
#'   }
#' @import rlang
#' @importFrom stats acf
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Calculate autocorrelation for channel 100
#' autocorr <- calculateSpikeAutocorrelation(spikes, channel_id = 100)
#' plot(autocorr$lag, autocorr$autocorr, type = "l")
#'
#' h5$close()
#' }
calculateSpikeAutocorrelation <- function(spike_data, channel_id,
                                          max_lag = 1, bin_size = 0.001) {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  if (!all(c("spike_times", "spike_chid") %in% names(spike_data))) {
    rlang::abort("spike_data must contain 'spike_times' and 'spike_chid' columns")
  }

  if (!is.numeric(channel_id) || length(channel_id) != 1) {
    rlang::abort("channel_id must be a single numeric value")
  }

  if (!is.numeric(max_lag) || max_lag <= 0) {
    rlang::abort("max_lag must be positive")
  }

  if (!is.numeric(bin_size) || bin_size <= 0) {
    rlang::abort("bin_size must be positive")
  }

  # Filter to specific channel
  channel_spikes <- spike_data[spike_data$spike_chid == channel_id, ]

  if (nrow(channel_spikes) < 2) {
    rlang::abort(paste("Channel", channel_id, "has fewer than 2 spikes"))
  }

  # Create binned spike train
  time_range <- range(channel_spikes$spike_times)
  bins <- seq(time_range[1], time_range[2], by = bin_size)
  spike_train <- hist(channel_spikes$spike_times, breaks = bins, plot = FALSE)$counts

  # Calculate autocorrelation
  max_lag_bins <- ceiling(max_lag / bin_size)
  if (max_lag_bins >= length(spike_train)) {
    max_lag_bins <- length(spike_train) - 1
  }

  acf_result <- stats::acf(spike_train, lag.max = max_lag_bins, plot = FALSE)

  # Format output
  result <- data.frame(
    lag = acf_result$lag * bin_size,
    autocorr = as.vector(acf_result$acf)
  )

  return(result)
}
