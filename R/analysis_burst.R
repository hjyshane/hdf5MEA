#' Detect bursts in spike trains
#'
#' Identifies burst events in spike trains using threshold-based detection.
#' A burst is defined as a sequence of spikes with ISI below a threshold.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}} with columns:
#'   spike_times, spike_chid
#' @param isi_threshold Numeric. Maximum ISI in seconds to consider spikes as part of burst
#'   (default: 0.1)
#' @param min_spikes Integer. Minimum number of spikes required for a burst (default: 5)
#' @param channel_id Integer. Specific channel to analyze. If NULL, analyzes all channels
#'   (default: NULL)
#' @return data.frame with columns:
#'   \describe{
#'     \item{channel_id}{Channel index}
#'     \item{burst_id}{Sequential burst identifier within channel}
#'     \item{start_time}{Burst start time in seconds}
#'     \item{end_time}{Burst end time in seconds}
#'     \item{duration}{Burst duration in seconds}
#'     \item{n_spikes}{Number of spikes in burst}
#'     \item{mean_isi}{Mean ISI within burst}
#'     \item{intra_burst_rate}{Firing rate within burst (Hz)}
#'   }
#' @import dplyr rlang
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Detect bursts with default parameters
#' bursts <- detectBursts(spikes, isi_threshold = 0.1, min_spikes = 5)
#'
#' # Detect bursts for specific channel
#' bursts_ch <- detectBursts(spikes, channel_id = 100)
#'
#' h5$close()
#' }
detectBursts <- function(spike_data, isi_threshold = 0.1, min_spikes = 5, channel_id = NULL) {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  if (!all(c("spike_times", "spike_chid") %in% names(spike_data))) {
    rlang::abort("spike_data must contain 'spike_times' and 'spike_chid' columns")
  }

  if (!is.numeric(isi_threshold) || isi_threshold <= 0) {
    rlang::abort("isi_threshold must be positive")
  }

  if (!is.numeric(min_spikes) || min_spikes < 2) {
    rlang::abort("min_spikes must be at least 2")
  }

  # Filter to specific channel if requested
  if (!is.null(channel_id)) {
    spike_data <- spike_data[spike_data$spike_chid == channel_id, ]
  }

  if (nrow(spike_data) == 0) {
    return(data.frame())
  }

  # Sort by channel and time
  spike_data <- spike_data[order(spike_data$spike_chid, spike_data$spike_times), ]

  # Process each channel
  channels <- unique(spike_data$spike_chid)
  all_bursts <- list()

  for (ch in channels) {
    ch_spikes <- spike_data[spike_data$spike_chid == ch, ]

    if (nrow(ch_spikes) < min_spikes) {
      next
    }

    # Calculate ISI
    isi <- diff(ch_spikes$spike_times)

    # Identify burst segments
    in_burst <- isi < isi_threshold
    burst_starts <- which(c(TRUE, in_burst) & !c(in_burst, FALSE))
    burst_ends <- which(c(in_burst, FALSE) & !c(FALSE, in_burst))

    if (length(burst_starts) == 0) {
      next
    }

    # Extract burst information
    bursts <- data.frame(
      channel_id = ch,
      burst_id = seq_along(burst_starts),
      start_idx = burst_starts,
      end_idx = burst_ends + 1  # Include the last spike
    )

    # Calculate burst properties
    bursts$start_time <- ch_spikes$spike_times[bursts$start_idx]
    bursts$end_time <- ch_spikes$spike_times[bursts$end_idx]
    bursts$duration <- bursts$end_time - bursts$start_time
    bursts$n_spikes <- bursts$end_idx - bursts$start_idx + 1

    # Calculate mean ISI within each burst
    bursts$mean_isi <- sapply(1:nrow(bursts), function(i) {
      spike_times_in_burst <- ch_spikes$spike_times[bursts$start_idx[i]:bursts$end_idx[i]]
      if (length(spike_times_in_burst) > 1) {
        mean(diff(spike_times_in_burst))
      } else {
        NA_real_
      }
    })

    # Calculate intra-burst firing rate
    bursts$intra_burst_rate <- ifelse(bursts$duration > 0,
                                       bursts$n_spikes / bursts$duration,
                                       NA_real_)

    # Filter by minimum spike count
    bursts <- bursts[bursts$n_spikes >= min_spikes, ]

    # Remove helper columns
    bursts$start_idx <- NULL
    bursts$end_idx <- NULL

    all_bursts[[as.character(ch)]] <- bursts
  }

  # Combine results
  if (length(all_bursts) == 0) {
    return(data.frame())
  }

  result <- do.call(rbind, all_bursts)
  rownames(result) <- NULL

  return(result)
}


#' Calculate burst statistics
#'
#' Computes comprehensive statistics for burst data including inter-burst
#' intervals, burst frequency, and burst characteristics.
#'
#' @param burst_data data.frame. Burst data from \code{\link{detectBursts}} or
#'   \code{\link{bxrSpikeBursts}} with columns: start_time, end_time, duration,
#'   channel_id or burst_chidx
#' @return data.frame with columns:
#'   \describe{
#'     \item{channel_id}{Channel index}
#'     \item{n_bursts}{Total number of bursts}
#'     \item{mean_burst_duration}{Mean burst duration in seconds}
#'     \item{median_burst_duration}{Median burst duration in seconds}
#'     \item{mean_ibi}{Mean inter-burst interval in seconds}
#'     \item{median_ibi}{Median inter-burst interval in seconds}
#'     \item{burst_frequency}{Bursts per minute}
#'     \item{total_burst_time}{Total time spent bursting in seconds}
#'     \item{burst_percentage}{Percentage of time spent bursting}
#'   }
#' @import dplyr rlang
#' @importFrom stats median sd
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#' bursts <- detectBursts(spikes)
#'
#' # Calculate burst statistics
#' burst_stats <- calculateBurstStatistics(bursts)
#'
#' h5$close()
#' }
calculateBurstStatistics <- function(burst_data) {
  # Validate inputs
  if (!is.data.frame(burst_data)) {
    rlang::abort("burst_data must be a data.frame")
  }

  if (nrow(burst_data) == 0) {
    return(data.frame())
  }

  # Check for required columns (handle both naming conventions)
  channel_col <- if ("channel_id" %in% names(burst_data)) {
    "channel_id"
  } else if ("burst_chidx" %in% names(burst_data)) {
    "burst_chidx"
  } else {
    rlang::abort("burst_data must contain 'channel_id' or 'burst_chidx' column")
  }

  if (!all(c("start_time", "end_time", "duration") %in% names(burst_data))) {
    rlang::abort("burst_data must contain 'start_time', 'end_time', and 'duration' columns")
  }

  # Rename to standard column name
  if (channel_col == "burst_chidx") {
    burst_data$channel_id <- burst_data$burst_chidx
  }

  # Calculate recording duration
  recording_duration <- max(burst_data$end_time) - min(burst_data$start_time)

  # Calculate statistics per channel
  burst_stats <- burst_data %>%
    dplyr::arrange(channel_id, start_time) %>%
    dplyr::group_by(channel_id) %>%
    dplyr::summarise(
      n_bursts = dplyr::n(),
      mean_burst_duration = mean(duration),
      median_burst_duration = stats::median(duration),
      sd_burst_duration = stats::sd(duration),
      mean_ibi = if (dplyr::n() > 1) mean(diff(start_time)) else NA_real_,
      median_ibi = if (dplyr::n() > 1) stats::median(diff(start_time)) else NA_real_,
      burst_frequency = dplyr::n() / (max(end_time) - min(start_time)) * 60,  # per minute
      total_burst_time = sum(duration),
      .groups = "drop"
    )

  # Calculate burst percentage
  burst_stats$burst_percentage <- (burst_stats$total_burst_time / recording_duration) * 100

  return(as.data.frame(burst_stats))
}


#' Analyze burst propagation across channels
#'
#' Examines spatial and temporal propagation of bursts across the MEA grid.
#' Identifies coordinated burst activity and propagation patterns.
#'
#' @param burst_data data.frame. Burst data with columns: start_time, channel_id, X, Y
#' @param time_window Numeric. Time window in seconds to consider bursts as coordinated
#'   (default: 0.05)
#' @return data.frame with columns:
#'   \describe{
#'     \item{propagation_id}{Sequential ID for coordinated burst events}
#'     \item{start_time}{Start time of coordinated burst}
#'     \item{end_time}{End time of coordinated burst}
#'     \item{n_channels}{Number of channels participating}
#'     \item{duration}{Duration of propagation in seconds}
#'     \item{center_x}{Center X coordinate of burst activity}
#'     \item{center_y}{Center Y coordinate of burst activity}
#'     \item{spatial_spread}{Spatial spread (max distance from center)}
#'   }
#' @import dplyr rlang
#' @importFrom stats dist
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' bursts <- bxrSpikeBursts(h5)
#'
#' # Analyze burst propagation
#' propagation <- analyzeBurstPropagation(bursts, time_window = 0.05)
#'
#' h5$close()
#' }
analyzeBurstPropagation <- function(burst_data, time_window = 0.05) {
  # Validate inputs
  if (!is.data.frame(burst_data)) {
    rlang::abort("burst_data must be a data.frame")
  }

  if (nrow(burst_data) == 0) {
    return(data.frame())
  }

  required_cols <- c("start_time", "channel_id")
  if (!all(required_cols %in% names(burst_data))) {
    # Check alternative naming
    if ("burst_chidx" %in% names(burst_data)) {
      burst_data$channel_id <- burst_data$burst_chidx
    } else {
      rlang::abort("burst_data must contain 'start_time' and 'channel_id' columns")
    }
  }

  if (!is.numeric(time_window) || time_window <= 0) {
    rlang::abort("time_window must be positive")
  }

  # Check if spatial coordinates are available
  has_spatial <- all(c("X", "Y") %in% names(burst_data))

  # Sort by start time
  burst_data <- burst_data[order(burst_data$start_time), ]

  # Group bursts by time proximity
  propagation_events <- list()
  current_group <- list()
  group_id <- 1

  for (i in 1:nrow(burst_data)) {
    if (length(current_group) == 0) {
      current_group[[1]] <- i
    } else {
      # Check if current burst is within time window of group
      group_start <- burst_data$start_time[current_group[[1]]]
      current_time <- burst_data$start_time[i]

      if (current_time - group_start <= time_window) {
        current_group[[length(current_group) + 1]] <- i
      } else {
        # Save current group if it has multiple channels
        if (length(current_group) > 1) {
          propagation_events[[group_id]] <- current_group
          group_id <- group_id + 1
        }
        # Start new group
        current_group <- list(i)
      }
    }
  }

  # Save last group
  if (length(current_group) > 1) {
    propagation_events[[group_id]] <- current_group
  }

  if (length(propagation_events) == 0) {
    return(data.frame())
  }

  # Extract propagation properties
  propagation_summary <- lapply(1:length(propagation_events), function(i) {
    indices <- unlist(propagation_events[[i]])
    group_bursts <- burst_data[indices, ]

    result <- data.frame(
      propagation_id = i,
      start_time = min(group_bursts$start_time),
      end_time = max(group_bursts$end_time),
      n_channels = length(unique(group_bursts$channel_id)),
      duration = max(group_bursts$end_time) - min(group_bursts$start_time)
    )

    # Add spatial information if available
    if (has_spatial) {
      result$center_x <- mean(group_bursts$X, na.rm = TRUE)
      result$center_y <- mean(group_bursts$Y, na.rm = TRUE)

      # Calculate spatial spread
      if (nrow(group_bursts) > 1) {
        coords <- as.matrix(group_bursts[, c("X", "Y")])
        center <- c(result$center_x, result$center_y)
        distances <- sqrt(rowSums((sweep(coords, 2, center))^2))
        result$spatial_spread <- max(distances, na.rm = TRUE)
      } else {
        result$spatial_spread <- 0
      }
    }

    return(result)
  })

  result <- do.call(rbind, propagation_summary)
  return(result)
}


#' Calculate burst synchrony index
#'
#' Computes synchrony index for burst activity across channels, measuring
#' how coordinated burst events are across the network.
#'
#' @param burst_data data.frame. Burst data with columns: start_time, channel_id
#' @param bin_size Numeric. Time bin size in seconds for synchrony calculation
#'   (default: 0.1)
#' @return List containing:
#'   \describe{
#'     \item{synchrony_index}{Overall synchrony index (0-1)}
#'     \item{burst_coincidence}{data.frame with time bins and coincident burst counts}
#'     \item{mean_coincident_bursts}{Mean number of simultaneous bursts per bin}
#'   }
#' @import dplyr rlang
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' bursts <- bxrSpikeBursts(h5)
#'
#' # Calculate burst synchrony
#' synchrony <- calculateBurstSynchrony(bursts, bin_size = 0.1)
#' print(synchrony$synchrony_index)
#'
#' h5$close()
#' }
calculateBurstSynchrony <- function(burst_data, bin_size = 0.1) {
  # Validate inputs
  if (!is.data.frame(burst_data)) {
    rlang::abort("burst_data must be a data.frame")
  }

  if (nrow(burst_data) == 0) {
    return(list(synchrony_index = 0, burst_coincidence = data.frame(),
                mean_coincident_bursts = 0))
  }

  required_cols <- c("start_time")
  if (!all(required_cols %in% names(burst_data))) {
    rlang::abort("burst_data must contain 'start_time' column")
  }

  # Handle channel column naming
  if ("burst_chidx" %in% names(burst_data)) {
    burst_data$channel_id <- burst_data$burst_chidx
  } else if (!"channel_id" %in% names(burst_data)) {
    rlang::abort("burst_data must contain 'channel_id' or 'burst_chidx' column")
  }

  if (!is.numeric(bin_size) || bin_size <= 0) {
    rlang::abort("bin_size must be positive")
  }

  # Create time bins
  time_range <- range(burst_data$start_time)
  bins <- seq(time_range[1], time_range[2], by = bin_size)

  # Assign bursts to bins
  burst_data$time_bin <- cut(burst_data$start_time,
                              breaks = bins,
                              labels = FALSE,
                              include.lowest = TRUE)

  # Count coincident bursts per bin
  coincidence <- burst_data %>%
    dplyr::group_by(time_bin) %>%
    dplyr::summarise(
      n_coincident_bursts = dplyr::n(),
      n_channels = dplyr::n_distinct(channel_id),
      .groups = "drop"
    )

  # Calculate synchrony index
  # Higher values indicate more coordinated bursting
  n_channels <- length(unique(burst_data$channel_id))
  synchrony_index <- mean(coincidence$n_channels) / n_channels

  result <- list(
    synchrony_index = synchrony_index,
    burst_coincidence = as.data.frame(coincidence),
    mean_coincident_bursts = mean(coincidence$n_coincident_bursts)
  )

  return(result)
}
