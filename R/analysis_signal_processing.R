#' Detect spikes from raw voltage traces
#'
#' Identifies spike events in continuous voltage recordings using threshold-based
#' detection with optional adaptive thresholding and dead-time enforcement.
#'
#' @param voltage_data data.frame or list. Voltage traces from \code{\link{brwtimeseriesConvert}}
#' @param sampling_rate Numeric. Sampling rate in Hz
#' @param threshold Numeric. Detection threshold in voltage units. Can be:
#'   - Absolute voltage value (e.g., -100)
#'   - NULL for automatic threshold (default: NULL, uses -4.5 * MAD)
#' @param method Character. Threshold method: "absolute", "std", or "mad" (default: "mad")
#' @param dead_time Numeric. Minimum time between spikes in seconds (default: 0.001)
#' @param polarity Character. Spike polarity: "negative", "positive", or "both" (default: "negative")
#' @return data.frame with columns:
#'   \describe{
#'     \item{channel_id}{Channel index}
#'     \item{spike_time}{Spike time in seconds}
#'     \item{spike_amplitude}{Peak amplitude at detection}
#'     \item{threshold_used}{Threshold value used for this channel}
#'   }
#' @import rlang
#' @importFrom stats sd mad median
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBRW("data.brw")
#' raw_data <- get_brw_data(h5, start = 0, duration = 10)
#' parsed <- brwdataParse(raw_data$binary_chunk)
#' timeseries <- brwtimeseriesConvert(parsed, raw_data$sampling_rate,
#'                                    raw_data$start_frame)
#'
#' # Detect spikes with automatic threshold
#' spikes <- detectSpikes(timeseries, sampling_rate = raw_data$sampling_rate)
#'
#' h5$close()
#' }
detectSpikes <- function(voltage_data, sampling_rate, threshold = NULL,
                         method = c("mad", "std", "absolute"),
                         dead_time = 0.001, polarity = c("negative", "positive", "both")) {
  # Validate inputs
  if (!is.list(voltage_data)) {
    rlang::abort("voltage_data must be a list from brwtimeseriesConvert()")
  }

  if (!is.numeric(sampling_rate) || sampling_rate <= 0) {
    rlang::abort("sampling_rate must be positive")
  }

  method <- match.arg(method)
  polarity <- match.arg(polarity)

  if (!is.numeric(dead_time) || dead_time < 0) {
    rlang::abort("dead_time must be non-negative")
  }

  dead_time_samples <- ceiling(dead_time * sampling_rate)

  # Process each channel
  all_spikes <- list()

  for (ch_name in names(voltage_data)) {
    ch_data <- voltage_data[[ch_name]]

    if (!is.data.frame(ch_data) || !"voltage" %in% names(ch_data) || !"time" %in% names(ch_data)) {
      next
    }

    if (nrow(ch_data) < 10) {
      next  # Skip channels with insufficient data
    }

    voltage <- ch_data$voltage
    time <- ch_data$time

    # Calculate threshold if not provided
    if (is.null(threshold)) {
      if (method == "mad") {
        noise_est <- stats::mad(voltage, center = stats::median(voltage))
        ch_threshold <- -4.5 * noise_est
      } else if (method == "std") {
        noise_est <- stats::sd(voltage)
        ch_threshold <- -4.5 * noise_est
      } else {
        rlang::abort("threshold must be specified for 'absolute' method")
      }
    } else {
      ch_threshold <- threshold
    }

    # Detect threshold crossings based on polarity
    if (polarity == "negative") {
      crossings <- which(voltage < ch_threshold)
    } else if (polarity == "positive") {
      crossings <- which(voltage > abs(ch_threshold))
    } else {  # both
      crossings <- which(abs(voltage) > abs(ch_threshold))
    }

    if (length(crossings) == 0) {
      next
    }

    # Enforce dead time and find local peaks
    spike_idx <- c()
    last_spike <- -Inf

    for (i in crossings) {
      if (i - last_spike > dead_time_samples) {
        # Find local peak around crossing
        search_window <- max(1, i - 5):min(length(voltage), i + 5)
        if (polarity == "negative") {
          peak_idx <- search_window[which.min(voltage[search_window])]
        } else if (polarity == "positive") {
          peak_idx <- search_window[which.max(voltage[search_window])]
        } else {
          peak_idx <- search_window[which.max(abs(voltage[search_window]))]
        }

        spike_idx <- c(spike_idx, peak_idx)
        last_spike <- peak_idx
      }
    }

    if (length(spike_idx) > 0) {
      ch_spikes <- data.frame(
        channel_id = as.numeric(ch_name),
        spike_time = time[spike_idx],
        spike_amplitude = voltage[spike_idx],
        threshold_used = ch_threshold
      )

      all_spikes[[ch_name]] <- ch_spikes
    }
  }

  # Combine results
  if (length(all_spikes) == 0) {
    return(data.frame(channel_id = numeric(), spike_time = numeric(),
                      spike_amplitude = numeric(), threshold_used = numeric()))
  }

  result <- do.call(rbind, all_spikes)
  rownames(result) <- NULL

  return(result)
}


#' Apply bandpass filter to voltage traces
#'
#' Filters voltage data to remove low-frequency drift and high-frequency noise
#' using a butterworth bandpass filter.
#'
#' @param voltage_data Numeric vector or data.frame with 'voltage' column
#' @param sampling_rate Numeric. Sampling rate in Hz
#' @param low_freq Numeric. Low cutoff frequency in Hz (default: 300)
#' @param high_freq Numeric. High cutoff frequency in Hz (default: 3000)
#' @param order Integer. Filter order (default: 4)
#' @return Filtered voltage data in same format as input
#' @import rlang
#' @importFrom stats filter
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBRW("data.brw")
#' raw_data <- get_brw_data(h5, start = 0, duration = 10)
#' parsed <- brwdataParse(raw_data$binary_chunk)
#' timeseries <- brwtimeseriesConvert(parsed, raw_data$sampling_rate,
#'                                    raw_data$start_frame)
#'
#' # Apply bandpass filter
#' filtered <- applyBandpassFilter(timeseries[[1]]$voltage,
#'                                  sampling_rate = raw_data$sampling_rate,
#'                                  low_freq = 300, high_freq = 3000)
#'
#' h5$close()
#' }
applyBandpassFilter <- function(voltage_data, sampling_rate,
                                low_freq = 300, high_freq = 3000, order = 4) {
  # Validate inputs
  if (!is.numeric(sampling_rate) || sampling_rate <= 0) {
    rlang::abort("sampling_rate must be positive")
  }

  if (!is.numeric(low_freq) || low_freq <= 0) {
    rlang::abort("low_freq must be positive")
  }

  if (!is.numeric(high_freq) || high_freq <= 0) {
    rlang::abort("high_freq must be positive")
  }

  if (high_freq >= sampling_rate / 2) {
    rlang::abort("high_freq must be less than Nyquist frequency (sampling_rate/2)")
  }

  if (low_freq >= high_freq) {
    rlang::abort("low_freq must be less than high_freq")
  }

  # Extract voltage vector
  if (is.data.frame(voltage_data)) {
    if (!"voltage" %in% names(voltage_data)) {
      rlang::abort("voltage_data must contain 'voltage' column")
    }
    voltage <- voltage_data$voltage
    is_df <- TRUE
  } else if (is.numeric(voltage_data)) {
    voltage <- voltage_data
    is_df <- FALSE
  } else {
    rlang::abort("voltage_data must be numeric vector or data.frame with 'voltage' column")
  }

  if (length(voltage) < order * 3) {
    rlang::abort("Signal too short for specified filter order")
  }

  # Simple butterworth approximation using cascaded filters
  # Normalize frequencies
  nyquist <- sampling_rate / 2
  low_norm <- low_freq / nyquist
  high_norm <- high_freq / nyquist

  # High-pass filter (remove low frequencies)
  alpha_high <- 1 / (1 + 2 * pi * low_norm)
  high_passed <- voltage
  for (i in 2:length(voltage)) {
    high_passed[i] <- alpha_high * (high_passed[i-1] + voltage[i] - voltage[i-1])
  }

  # Low-pass filter (remove high frequencies)
  alpha_low <- 2 * pi * high_norm / (1 + 2 * pi * high_norm)
  filtered <- high_passed
  for (i in 2:length(high_passed)) {
    filtered[i] <- alpha_low * high_passed[i] + (1 - alpha_low) * filtered[i-1]
  }

  # Return in same format as input
  if (is_df) {
    voltage_data$voltage <- filtered
    return(voltage_data)
  } else {
    return(filtered)
  }
}


#' Remove baseline drift from voltage traces
#'
#' Removes slow baseline drift using median filtering or linear detrending.
#'
#' @param voltage_data Numeric vector or data.frame with 'voltage' column
#' @param method Character. Detrending method: "median" or "linear" (default: "median")
#' @param window_size Numeric. Window size for median filter in seconds
#'   (default: 1, only used for median method)
#' @param sampling_rate Numeric. Sampling rate in Hz (required for median method)
#' @return Detrended voltage data in same format as input
#' @import rlang
#' @importFrom stats median lm fitted
#' @export
#' @examples
#' \dontrun{
#' # Remove baseline drift
#' detrended <- removeBaselineDrift(voltage_trace, method = "median",
#'                                   window_size = 1, sampling_rate = 10000)
#' }
removeBaselineDrift <- function(voltage_data, method = c("median", "linear"),
                                window_size = 1, sampling_rate = NULL) {
  method <- match.arg(method)

  # Extract voltage vector
  if (is.data.frame(voltage_data)) {
    if (!"voltage" %in% names(voltage_data)) {
      rlang::abort("voltage_data must contain 'voltage' column")
    }
    voltage <- voltage_data$voltage
    is_df <- TRUE
  } else if (is.numeric(voltage_data)) {
    voltage <- voltage_data
    is_df <- FALSE
  } else {
    rlang::abort("voltage_data must be numeric vector or data.frame")
  }

  if (method == "median") {
    if (is.null(sampling_rate)) {
      rlang::abort("sampling_rate required for median method")
    }

    if (!is.numeric(window_size) || window_size <= 0) {
      rlang::abort("window_size must be positive")
    }

    # Convert window size to samples
    window_samples <- ceiling(window_size * sampling_rate)
    if (window_samples > length(voltage)) {
      window_samples <- length(voltage)
    }

    # Apply running median
    baseline <- rep(0, length(voltage))
    half_window <- floor(window_samples / 2)

    for (i in 1:length(voltage)) {
      start_idx <- max(1, i - half_window)
      end_idx <- min(length(voltage), i + half_window)
      baseline[i] <- stats::median(voltage[start_idx:end_idx])
    }

    corrected <- voltage - baseline

  } else {  # linear
    # Fit linear trend and remove
    time_index <- 1:length(voltage)
    trend <- stats::lm(voltage ~ time_index)
    baseline <- stats::fitted(trend)
    corrected <- voltage - baseline
  }

  # Return in same format as input
  if (is_df) {
    voltage_data$voltage <- corrected
    return(voltage_data)
  } else {
    return(corrected)
  }
}


#' Detect and remove artifacts
#'
#' Identifies and removes electrical artifacts based on amplitude thresholds
#' and duration criteria.
#'
#' @param voltage_data data.frame or list. Voltage traces
#' @param artifact_threshold Numeric. Amplitude threshold for artifact detection
#'   (default: 1000, in voltage units)
#' @param min_artifact_duration Numeric. Minimum duration in seconds to consider
#'   as artifact (default: 0.001)
#' @param max_artifact_duration Numeric. Maximum duration in seconds to consider
#'   as artifact (default: 0.1)
#' @param sampling_rate Numeric. Sampling rate in Hz
#' @param action Character. Action to take: "remove", "interpolate", or "flag"
#'   (default: "interpolate")
#' @return List containing:
#'   \describe{
#'     \item{cleaned_data}{Voltage data with artifacts handled}
#'     \item{artifact_times}{data.frame with artifact locations}
#'     \item{n_artifacts}{Number of artifacts detected}
#'   }
#' @import rlang
#' @importFrom stats approx
#' @export
#' @examples
#' \dontrun{
#' # Detect and remove artifacts
#' cleaned <- detectArtifacts(voltage_data, artifact_threshold = 1000,
#'                             sampling_rate = 10000, action = "interpolate")
#' }
detectArtifacts <- function(voltage_data, artifact_threshold = 1000,
                            min_artifact_duration = 0.001,
                            max_artifact_duration = 0.1,
                            sampling_rate, action = c("interpolate", "remove", "flag")) {
  action <- match.arg(action)

  if (!is.numeric(sampling_rate) || sampling_rate <= 0) {
    rlang::abort("sampling_rate must be positive")
  }

  if (!is.numeric(artifact_threshold) || artifact_threshold <= 0) {
    rlang::abort("artifact_threshold must be positive")
  }

  # Convert durations to samples
  min_samples <- ceiling(min_artifact_duration * sampling_rate)
  max_samples <- ceiling(max_artifact_duration * sampling_rate)

  # Process based on input type
  if (is.list(voltage_data) && !is.data.frame(voltage_data)) {
    # Process each channel
    cleaned_data <- list()
    all_artifacts <- list()

    for (ch_name in names(voltage_data)) {
      ch_result <- detectArtifacts(voltage_data[[ch_name]],
                                    artifact_threshold = artifact_threshold,
                                    min_artifact_duration = min_artifact_duration,
                                    max_artifact_duration = max_artifact_duration,
                                    sampling_rate = sampling_rate,
                                    action = action)
      cleaned_data[[ch_name]] <- ch_result$cleaned_data
      all_artifacts[[ch_name]] <- ch_result$artifact_times
    }

    return(list(
      cleaned_data = cleaned_data,
      artifact_times = do.call(rbind, all_artifacts),
      n_artifacts = sum(sapply(all_artifacts, nrow))
    ))
  }

  # Single channel processing
  if (!is.data.frame(voltage_data) || !all(c("voltage", "time") %in% names(voltage_data))) {
    rlang::abort("voltage_data must contain 'voltage' and 'time' columns")
  }

  voltage <- voltage_data$voltage
  time <- voltage_data$time

  # Detect artifacts (extreme amplitude deviations)
  artifact_mask <- abs(voltage) > artifact_threshold

  # Find continuous artifact segments
  artifact_starts <- which(diff(c(0, artifact_mask)) == 1)
  artifact_ends <- which(diff(c(artifact_mask, 0)) == -1)

  if (length(artifact_starts) == 0) {
    return(list(
      cleaned_data = voltage_data,
      artifact_times = data.frame(),
      n_artifacts = 0
    ))
  }

  # Filter by duration
  artifact_durations <- artifact_ends - artifact_starts + 1
  valid_artifacts <- artifact_durations >= min_samples & artifact_durations <= max_samples

  artifact_starts <- artifact_starts[valid_artifacts]
  artifact_ends <- artifact_ends[valid_artifacts]

  if (length(artifact_starts) == 0) {
    return(list(
      cleaned_data = voltage_data,
      artifact_times = data.frame(),
      n_artifacts = 0
    ))
  }

  # Handle artifacts based on action
  cleaned_voltage <- voltage

  if (action == "interpolate") {
    for (i in 1:length(artifact_starts)) {
      start_idx <- max(1, artifact_starts[i] - 2)
      end_idx <- min(length(voltage), artifact_ends[i] + 2)

      # Linear interpolation
      good_indices <- c(start_idx, end_idx)
      artifact_indices <- artifact_starts[i]:artifact_ends[i]

      cleaned_voltage[artifact_indices] <- stats::approx(
        x = time[good_indices],
        y = voltage[good_indices],
        xout = time[artifact_indices]
      )$y
    }
  } else if (action == "remove") {
    # Mark for removal
    remove_mask <- rep(FALSE, length(voltage))
    for (i in 1:length(artifact_starts)) {
      remove_mask[artifact_starts[i]:artifact_ends[i]] <- TRUE
    }
    voltage_data <- voltage_data[!remove_mask, ]
    return(list(
      cleaned_data = voltage_data,
      artifact_times = data.frame(
        start_time = time[artifact_starts],
        end_time = time[artifact_ends],
        duration = time[artifact_ends] - time[artifact_starts]
      ),
      n_artifacts = length(artifact_starts)
    ))
  }
  # If action == "flag", just return data with artifact info

  voltage_data$voltage <- cleaned_voltage

  artifact_info <- data.frame(
    start_time = time[artifact_starts],
    end_time = time[artifact_ends],
    duration = time[artifact_ends] - time[artifact_starts],
    max_amplitude = sapply(1:length(artifact_starts), function(i) {
      max(abs(voltage[artifact_starts[i]:artifact_ends[i]]))
    })
  )

  return(list(
    cleaned_data = voltage_data,
    artifact_times = artifact_info,
    n_artifacts = nrow(artifact_info)
  ))
}
