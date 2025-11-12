#' Calculate pairwise cross-correlation between channels
#'
#' Computes spike time cross-correlation between channel pairs to assess
#' functional connectivity. Identifies temporal relationships and potential
#' synaptic connections.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}} with columns:
#'   spike_times, spike_chid
#' @param channel_pairs data.frame. Pairs of channels to analyze with columns: ch1, ch2.
#'   If NULL, analyzes all unique pairs (default: NULL, limited to 100 random pairs)
#' @param max_lag Numeric. Maximum lag in seconds (default: 0.05)
#' @param bin_size Numeric. Bin size for correlation in seconds (default: 0.001)
#' @return data.frame with columns:
#'   \describe{
#'     \item{channel_1}{First channel ID}
#'     \item{channel_2}{Second channel ID}
#'     \item{peak_correlation}{Maximum correlation value}
#'     \item{peak_lag}{Lag at maximum correlation (seconds)}
#'     \item{correlation_values}{List column with full correlation values}
#'     \item{lags}{List column with lag values}
#'   }
#' @import dplyr rlang
#' @importFrom stats ccf
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Analyze specific pairs
#' pairs <- data.frame(ch1 = c(1, 1), ch2 = c(2, 3))
#' corr <- calculateCrossCorrelation(spikes, channel_pairs = pairs)
#'
#' h5$close()
#' }
calculateCrossCorrelation <- function(spike_data, channel_pairs = NULL,
                                      max_lag = 0.05, bin_size = 0.001) {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  if (!all(c("spike_times", "spike_chid") %in% names(spike_data))) {
    rlang::abort("spike_data must contain 'spike_times' and 'spike_chid' columns")
  }

  if (!is.numeric(max_lag) || max_lag <= 0) {
    rlang::abort("max_lag must be positive")
  }

  if (!is.numeric(bin_size) || bin_size <= 0) {
    rlang::abort("bin_size must be positive")
  }

  # Get all channels
  all_channels <- unique(spike_data$spike_chid)

  # Generate channel pairs if not provided
  if (is.null(channel_pairs)) {
    # Limit to 100 random pairs to avoid excessive computation
    if (length(all_channels) > 15) {
      warning("Too many channels. Randomly selecting 100 pairs. Specify channel_pairs for custom selection.")
      n_pairs <- min(100, choose(length(all_channels), 2))
      all_pairs <- t(combn(all_channels, 2))
      sample_idx <- sample(1:nrow(all_pairs), n_pairs)
      channel_pairs <- data.frame(ch1 = all_pairs[sample_idx, 1],
                                   ch2 = all_pairs[sample_idx, 2])
    } else {
      all_pairs <- t(combn(all_channels, 2))
      channel_pairs <- data.frame(ch1 = all_pairs[, 1], ch2 = all_pairs[, 2])
    }
  }

  if (!all(c("ch1", "ch2") %in% names(channel_pairs))) {
    rlang::abort("channel_pairs must contain 'ch1' and 'ch2' columns")
  }

  # Create binned spike trains
  time_range <- range(spike_data$spike_times)
  bins <- seq(time_range[1], time_range[2], by = bin_size)

  # Pre-compute binned spike trains for all channels
  spike_trains <- lapply(all_channels, function(ch) {
    ch_spikes <- spike_data$spike_times[spike_data$spike_chid == ch]
    hist(ch_spikes, breaks = bins, plot = FALSE)$counts
  })
  names(spike_trains) <- as.character(all_channels)

  # Calculate cross-correlation for each pair
  max_lag_bins <- ceiling(max_lag / bin_size)

  results <- lapply(1:nrow(channel_pairs), function(i) {
    ch1 <- channel_pairs$ch1[i]
    ch2 <- channel_pairs$ch2[i]

    train1 <- spike_trains[[as.character(ch1)]]
    train2 <- spike_trains[[as.character(ch2)]]

    if (is.null(train1) || is.null(train2)) {
      return(NULL)
    }

    # Compute cross-correlation
    ccf_result <- tryCatch({
      stats::ccf(train1, train2, lag.max = max_lag_bins, plot = FALSE)
    }, error = function(e) {
      return(NULL)
    })

    if (is.null(ccf_result)) {
      return(NULL)
    }

    # Extract results
    lags <- as.vector(ccf_result$lag) * bin_size
    correlations <- as.vector(ccf_result$acf)

    # Find peak
    peak_idx <- which.max(abs(correlations))

    data.frame(
      channel_1 = ch1,
      channel_2 = ch2,
      peak_correlation = correlations[peak_idx],
      peak_lag = lags[peak_idx],
      correlation_values = I(list(correlations)),
      lags = I(list(lags))
    )
  })

  # Combine results
  results <- results[!sapply(results, is.null)]

  if (length(results) == 0) {
    return(data.frame())
  }

  result <- do.call(rbind, results)
  return(result)
}


#' Calculate network-wide synchrony index
#'
#' Computes synchrony metrics for the entire network to assess coordinated
#' activity across channels. High synchrony indicates network-wide events.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}} with columns:
#'   spike_times, spike_chid
#' @param bin_size Numeric. Time bin size in seconds (default: 0.01)
#' @param method Character. Method for synchrony calculation: "correlation" or "coincidence"
#'   (default: "coincidence")
#' @return List containing:
#'   \describe{
#'     \item{synchrony_index}{Overall network synchrony (0-1)}
#'     \item{temporal_profile}{data.frame with time bins and synchrony over time}
#'     \item{active_channels_per_bin}{Mean number of active channels per bin}
#'   }
#' @import dplyr rlang
#' @importFrom stats cor
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Calculate network synchrony
#' synchrony <- calculateNetworkSynchrony(spikes, bin_size = 0.01)
#' print(synchrony$synchrony_index)
#'
#' h5$close()
#' }
calculateNetworkSynchrony <- function(spike_data, bin_size = 0.01,
                                      method = c("coincidence", "correlation")) {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  if (!all(c("spike_times", "spike_chid") %in% names(spike_data))) {
    rlang::abort("spike_data must contain 'spike_times' and 'spike_chid' columns")
  }

  if (!is.numeric(bin_size) || bin_size <= 0) {
    rlang::abort("bin_size must be positive")
  }

  method <- match.arg(method)

  # Create time bins
  time_range <- range(spike_data$spike_times)
  bins <- seq(time_range[1], time_range[2], by = bin_size)
  n_bins <- length(bins) - 1

  # Get all channels
  all_channels <- unique(spike_data$spike_chid)
  n_channels <- length(all_channels)

  if (n_channels < 2) {
    rlang::abort("Need at least 2 channels for synchrony calculation")
  }

  if (method == "coincidence") {
    # Count active channels per bin
    spike_data$time_bin <- cut(spike_data$spike_times,
                                breaks = bins,
                                labels = FALSE,
                                include.lowest = TRUE)

    bin_activity <- spike_data %>%
      dplyr::group_by(time_bin) %>%
      dplyr::summarise(
        n_active = dplyr::n_distinct(spike_chid),
        n_spikes = dplyr::n(),
        .groups = "drop"
      )

    # Calculate synchrony as proportion of channels active together
    bin_activity$synchrony <- bin_activity$n_active / n_channels

    # Overall synchrony index
    synchrony_index <- mean(bin_activity$synchrony)

    # Temporal profile
    bin_centers <- (bins[-1] + bins[-length(bins)]) / 2
    temporal_profile <- data.frame(
      time = bin_centers[bin_activity$time_bin],
      synchrony = bin_activity$synchrony,
      n_active = bin_activity$n_active
    )

    result <- list(
      synchrony_index = synchrony_index,
      temporal_profile = temporal_profile,
      active_channels_per_bin = mean(bin_activity$n_active)
    )

  } else {  # correlation method
    # Create spike train matrix
    spike_trains <- matrix(0, nrow = n_channels, ncol = n_bins)
    rownames(spike_trains) <- as.character(all_channels)

    for (i in 1:n_channels) {
      ch <- all_channels[i]
      ch_spikes <- spike_data$spike_times[spike_data$spike_chid == ch]
      spike_trains[i, ] <- hist(ch_spikes, breaks = bins, plot = FALSE)$counts
    }

    # Calculate pairwise correlations
    cor_matrix <- stats::cor(t(spike_trains))
    diag(cor_matrix) <- NA  # Exclude self-correlations

    # Synchrony index as mean correlation
    synchrony_index <- mean(cor_matrix, na.rm = TRUE)

    # Temporal profile based on population activity
    pop_activity <- colSums(spike_trains > 0)
    bin_centers <- (bins[-1] + bins[-length(bins)]) / 2

    temporal_profile <- data.frame(
      time = bin_centers,
      n_active = pop_activity,
      synchrony = pop_activity / n_channels
    )

    result <- list(
      synchrony_index = synchrony_index,
      temporal_profile = temporal_profile,
      active_channels_per_bin = mean(pop_activity)
    )
  }

  return(result)
}


#' Detect spatial patterns in spike activity
#'
#' Identifies spatial clusters and propagating waves in MEA recordings
#' using spatial statistics on grid-mapped spike data.
#'
#' @param spike_data data.frame. Spike data with columns: spike_times, spike_chid, X, Y
#' @param time_window Numeric. Time window in seconds for pattern detection (default: 0.1)
#' @param spatial_threshold Numeric. Maximum spatial distance to consider channels as
#'   neighbors (default: 3)
#' @return List containing:
#'   \describe{
#'     \item{spatial_events}{data.frame with detected spatial events}
#'     \item{hotspots}{data.frame with persistent high-activity regions}
#'     \item{wave_events}{data.frame with potential wave propagation events}
#'   }
#' @import dplyr rlang
#' @importFrom stats dist
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Detect spatial patterns
#' patterns <- detectSpatialPatterns(spikes, time_window = 0.1)
#'
#' h5$close()
#' }
detectSpatialPatterns <- function(spike_data, time_window = 0.1, spatial_threshold = 3) {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  required_cols <- c("spike_times", "spike_chid", "X", "Y")
  if (!all(required_cols %in% names(spike_data))) {
    rlang::abort("spike_data must contain 'spike_times', 'spike_chid', 'X', and 'Y' columns")
  }

  if (nrow(spike_data) == 0) {
    return(list(spatial_events = data.frame(),
                hotspots = data.frame(),
                wave_events = data.frame()))
  }

  # Remove rows with missing spatial coordinates
  spike_data <- spike_data[!is.na(spike_data$X) & !is.na(spike_data$Y), ]

  if (nrow(spike_data) == 0) {
    rlang::abort("No spikes with valid spatial coordinates")
  }

  # Sort by time
  spike_data <- spike_data[order(spike_data$spike_times), ]

  # Detect spatial events (clustered activity in space and time)
  time_range <- range(spike_data$spike_times)
  time_bins <- seq(time_range[1], time_range[2], by = time_window)

  spike_data$time_bin <- cut(spike_data$spike_times,
                              breaks = time_bins,
                              labels = FALSE,
                              include.lowest = TRUE)

  # Analyze each time bin
  spatial_events <- spike_data %>%
    dplyr::group_by(time_bin) %>%
    dplyr::summarise(
      time = mean(spike_times),
      n_spikes = dplyr::n(),
      n_channels = dplyr::n_distinct(spike_chid),
      center_x = mean(X),
      center_y = mean(Y),
      spread_x = stats::sd(X),
      spread_y = stats::sd(Y),
      .groups = "drop"
    )

  # Identify hotspots (persistent high-activity regions)
  channel_activity <- spike_data %>%
    dplyr::group_by(spike_chid, X, Y) %>%
    dplyr::summarise(
      n_spikes = dplyr::n(),
      firing_rate = dplyr::n() / (max(spike_times) - min(spike_times)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(firing_rate))

  # Define hotspots as top 10% most active channels
  n_hotspots <- max(1, ceiling(nrow(channel_activity) * 0.1))
  hotspots <- channel_activity[1:n_hotspots, ]

  # Detect wave-like propagation
  # Look for sequential activation across neighboring channels
  wave_events <- list()
  wave_id <- 1

  for (tb in unique(spike_data$time_bin)) {
    bin_spikes <- spike_data[spike_data$time_bin == tb, ]

    if (nrow(bin_spikes) < 3) next  # Need at least 3 channels

    # Sort by spike time within bin
    bin_spikes <- bin_spikes[order(bin_spikes$spike_times), ]

    # Check for spatial progression
    coords <- as.matrix(bin_spikes[, c("X", "Y")])

    # Calculate center of mass progression
    n <- nrow(bin_spikes)
    if (n >= 5) {
      # Split into early and late phases
      early_idx <- 1:floor(n/2)
      late_idx <- (floor(n/2)+1):n

      early_center <- colMeans(coords[early_idx, , drop = FALSE])
      late_center <- colMeans(coords[late_idx, , drop = FALSE])

      # Calculate displacement
      displacement <- sqrt(sum((late_center - early_center)^2))

      # If significant displacement, mark as potential wave
      if (displacement >= spatial_threshold) {
        wave_events[[wave_id]] <- data.frame(
          wave_id = wave_id,
          time = mean(bin_spikes$spike_times),
          n_channels = n,
          start_x = early_center[1],
          start_y = early_center[2],
          end_x = late_center[1],
          end_y = late_center[2],
          displacement = displacement,
          direction = atan2(late_center[2] - early_center[2],
                           late_center[1] - early_center[1]) * 180 / pi
        )
        wave_id <- wave_id + 1
      }
    }
  }

  # Combine wave events
  if (length(wave_events) > 0) {
    wave_events <- do.call(rbind, wave_events)
  } else {
    wave_events <- data.frame()
  }

  result <- list(
    spatial_events = as.data.frame(spatial_events),
    hotspots = as.data.frame(hotspots),
    wave_events = as.data.frame(wave_events)
  )

  return(result)
}


#' Build functional connectivity matrix
#'
#' Constructs a functional connectivity matrix based on spike time correlations
#' or mutual information between all channel pairs.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}}
#' @param method Character. Connectivity metric: "correlation" or "granger" (default: "correlation")
#' @param bin_size Numeric. Time bin size in seconds (default: 0.01)
#' @param threshold Numeric. Threshold for considering channels connected (default: 0.3)
#' @return List containing:
#'   \describe{
#'     \item{connectivity_matrix}{Matrix of connectivity strengths}
#'     \item{adjacency_matrix}{Binary adjacency matrix (thresholded)}
#'     \item{network_density}{Proportion of possible connections present}
#'     \item{hub_channels}{Channels with highest connectivity}
#'   }
#' @import dplyr rlang
#' @importFrom stats cor
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Build connectivity matrix
#' connectivity <- buildConnectivityMatrix(spikes, threshold = 0.3)
#' heatmap(connectivity$connectivity_matrix)
#'
#' h5$close()
#' }
buildConnectivityMatrix <- function(spike_data, method = c("correlation", "granger"),
                                    bin_size = 0.01, threshold = 0.3) {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  if (!all(c("spike_times", "spike_chid") %in% names(spike_data))) {
    rlang::abort("spike_data must contain 'spike_times' and 'spike_chid' columns")
  }

  method <- match.arg(method)

  if (!is.numeric(bin_size) || bin_size <= 0) {
    rlang::abort("bin_size must be positive")
  }

  if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
    rlang::abort("threshold must be between 0 and 1")
  }

  # Get all channels
  all_channels <- sort(unique(spike_data$spike_chid))
  n_channels <- length(all_channels)

  # Create time bins
  time_range <- range(spike_data$spike_times)
  bins <- seq(time_range[1], time_range[2], by = bin_size)
  n_bins <- length(bins) - 1

  # Create spike train matrix
  spike_trains <- matrix(0, nrow = n_channels, ncol = n_bins)
  rownames(spike_trains) <- as.character(all_channels)

  for (i in 1:n_channels) {
    ch <- all_channels[i]
    ch_spikes <- spike_data$spike_times[spike_data$spike_chid == ch]
    spike_trains[i, ] <- hist(ch_spikes, breaks = bins, plot = FALSE)$counts
  }

  # Calculate connectivity matrix
  if (method == "correlation") {
    connectivity_matrix <- stats::cor(t(spike_trains))
    diag(connectivity_matrix) <- 0  # Remove self-connections
  } else {  # granger causality (simplified)
    # Simplified Granger causality using lagged correlation
    connectivity_matrix <- matrix(0, n_channels, n_channels)

    for (i in 1:n_channels) {
      for (j in 1:n_channels) {
        if (i != j) {
          # Lag j by 1 bin and correlate with i
          lagged_j <- c(0, spike_trains[j, -n_bins])
          connectivity_matrix[i, j] <- stats::cor(spike_trains[i, ], lagged_j)
        }
      }
    }

    rownames(connectivity_matrix) <- as.character(all_channels)
    colnames(connectivity_matrix) <- as.character(all_channels)
  }

  # Create binary adjacency matrix
  adjacency_matrix <- abs(connectivity_matrix) >= threshold
  diag(adjacency_matrix) <- FALSE

  # Calculate network density
  n_possible <- n_channels * (n_channels - 1)
  n_actual <- sum(adjacency_matrix)
  network_density <- n_actual / n_possible

  # Identify hub channels (high degree)
  degree <- rowSums(adjacency_matrix)
  hub_threshold <- quantile(degree, 0.8)
  hub_channels <- all_channels[degree >= hub_threshold]

  result <- list(
    connectivity_matrix = connectivity_matrix,
    adjacency_matrix = adjacency_matrix,
    network_density = network_density,
    hub_channels = hub_channels,
    degree = data.frame(channel_id = all_channels, degree = degree)
  )

  return(result)
}
