#' Plot spike raster
#'
#' Creates a raster plot showing spike times across channels. Each spike
#' is represented as a vertical tick, with channels on the y-axis and time on the x-axis.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}} with columns:
#'   spike_times, spike_chid
#' @param time_range Numeric vector of length 2. Time range to plot [start, end] in seconds.
#'   If NULL, plots full range (default: NULL)
#' @param channels Integer vector. Specific channels to plot. If NULL, plots all (default: NULL)
#' @param color Character. Color for spike marks (default: "black")
#' @param title Character. Plot title (default: "Spike Raster Plot")
#' @return ggplot2 object
#' @import ggplot2 dplyr rlang
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Plot raster for first 10 seconds
#' plotSpikeRaster(spikes, time_range = c(0, 10))
#'
#' h5$close()
#' }
plotSpikeRaster <- function(spike_data, time_range = NULL, channels = NULL,
                            color = "black", title = "Spike Raster Plot") {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  if (!all(c("spike_times", "spike_chid") %in% names(spike_data))) {
    rlang::abort("spike_data must contain 'spike_times' and 'spike_chid' columns")
  }

  # Filter by time range
  if (!is.null(time_range)) {
    if (!is.numeric(time_range) || length(time_range) != 2) {
      rlang::abort("time_range must be numeric vector of length 2")
    }
    spike_data <- spike_data[spike_data$spike_times >= time_range[1] &
                               spike_data$spike_times <= time_range[2], ]
  }

  # Filter by channels
  if (!is.null(channels)) {
    spike_data <- spike_data[spike_data$spike_chid %in% channels, ]
  }

  if (nrow(spike_data) == 0) {
    rlang::abort("No spikes in specified range")
  }

  # Create plot
  p <- ggplot2::ggplot(spike_data, ggplot2::aes(x = spike_times, y = spike_chid)) +
    ggplot2::geom_point(shape = "|", size = 2, color = color, alpha = 0.7) +
    ggplot2::labs(
      title = title,
      x = "Time (s)",
      y = "Channel ID"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  return(p)
}


#' Plot firing rate heatmap
#'
#' Creates a 2D heatmap showing firing rates across the MEA grid.
#' Requires spatial coordinates (X, Y) in the spike data.
#'
#' @param spike_data data.frame. Spike data with columns: spike_times, spike_chid, X, Y
#' @param bin_size Numeric. Time bin size for rate calculation in seconds. If NULL,
#'   calculates overall firing rate (default: NULL)
#' @param time_range Numeric vector of length 2. Time range for analysis (default: NULL)
#' @param title Character. Plot title (default: "Firing Rate Heatmap")
#' @param color_scale Character. Color scale: "viridis", "plasma", "inferno" (default: "viridis")
#' @return ggplot2 object
#' @import ggplot2 dplyr rlang
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Plot overall firing rate heatmap
#' plotFiringRateHeatmap(spikes)
#'
#' h5$close()
#' }
plotFiringRateHeatmap <- function(spike_data, bin_size = NULL, time_range = NULL,
                                  title = "Firing Rate Heatmap",
                                  color_scale = c("viridis", "plasma", "inferno")) {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  required_cols <- c("spike_times", "spike_chid", "X", "Y")
  if (!all(required_cols %in% names(spike_data))) {
    rlang::abort("spike_data must contain 'spike_times', 'spike_chid', 'X', and 'Y' columns")
  }

  color_scale <- match.arg(color_scale)

  # Remove NA coordinates
  spike_data <- spike_data[!is.na(spike_data$X) & !is.na(spike_data$Y), ]

  if (nrow(spike_data) == 0) {
    rlang::abort("No spikes with valid spatial coordinates")
  }

  # Filter by time range
  if (!is.null(time_range)) {
    spike_data <- spike_data[spike_data$spike_times >= time_range[1] &
                               spike_data$spike_times <= time_range[2], ]
    duration <- time_range[2] - time_range[1]
  } else {
    duration <- max(spike_data$spike_times) - min(spike_data$spike_times)
  }

  # Calculate firing rate per channel
  if (is.null(bin_size)) {
    # Overall firing rate
    rate_data <- spike_data %>%
      dplyr::group_by(spike_chid, X, Y) %>%
      dplyr::summarise(
        firing_rate = dplyr::n() / duration,
        .groups = "drop"
      )
  } else {
    # Binned firing rate (use mean across bins)
    rates <- calculateFiringRate(spike_data, bin_size = bin_size)
    rate_data <- rates %>%
      dplyr::group_by(channel_id) %>%
      dplyr::summarise(
        firing_rate = mean(firing_rate),
        .groups = "drop"
      ) %>%
      dplyr::left_join(
        spike_data %>% dplyr::select(spike_chid, X, Y) %>% dplyr::distinct(),
        by = c("channel_id" = "spike_chid")
      )
  }

  # Create heatmap
  p <- ggplot2::ggplot(rate_data, ggplot2::aes(x = X, y = Y, fill = firing_rate)) +
    ggplot2::geom_tile() +
    ggplot2::labs(
      title = title,
      x = "X Position",
      y = "Y Position",
      fill = "Firing Rate (Hz)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      aspect.ratio = 1
    ) +
    ggplot2::coord_equal()

  # Apply color scale
  if (color_scale == "viridis") {
    p <- p + ggplot2::scale_fill_viridis_c(option = "viridis")
  } else if (color_scale == "plasma") {
    p <- p + ggplot2::scale_fill_viridis_c(option = "plasma")
  } else {
    p <- p + ggplot2::scale_fill_viridis_c(option = "inferno")
  }

  return(p)
}


#' Plot average spike waveforms
#'
#' Visualizes average waveform shapes for specified channels with optional
#' individual waveform overlays.
#'
#' @param waveform_data List. Waveform data from \code{\link{getWaveform}}
#' @param channels Integer vector. Channels to plot (default: plots up to 9 channels)
#' @param show_individual Logical. Whether to show individual waveforms (default: FALSE)
#' @param n_individual Integer. Number of individual waveforms to show if show_individual=TRUE
#'   (default: 50)
#' @param title Character. Plot title (default: "Average Spike Waveforms")
#' @return ggplot2 object
#' @import ggplot2 dplyr rlang
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' waves <- getWaveform(h5, waveform_indices = 1:100)
#'
#' # Plot average waveforms
#' plotWaveforms(waves, channels = c(1, 5, 10))
#'
#' h5$close()
#' }
plotWaveforms <- function(waveform_data, channels = NULL, show_individual = FALSE,
                          n_individual = 50, title = "Average Spike Waveforms") {
  # Validate inputs
  if (!is.list(waveform_data)) {
    rlang::abort("waveform_data must be a list from getWaveform()")
  }

  if (!"waveforms" %in% names(waveform_data) || !"time_axis" %in% names(waveform_data)) {
    rlang::abort("waveform_data must contain 'waveforms' and 'time_axis' elements")
  }

  # For now, create a simple example plot structure
  # This assumes waveforms is a matrix where each row is a waveform
  waveforms <- waveform_data$waveforms
  time_axis <- waveform_data$time_axis

  # Create plot data
  if (is.matrix(waveforms)) {
    # Calculate mean waveform
    mean_wave <- colMeans(waveforms)

    plot_data <- data.frame(
      time = time_axis,
      voltage = mean_wave
    )

    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time * 1000, y = voltage)) +
      ggplot2::geom_line(linewidth = 1.2, color = "blue") +
      ggplot2::labs(
        title = title,
        x = "Time (ms)",
        y = "Voltage"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      )

    # Add individual waveforms if requested
    if (show_individual && nrow(waveforms) > 0) {
      n_show <- min(n_individual, nrow(waveforms))
      sample_idx <- sample(1:nrow(waveforms), n_show)

      for (idx in sample_idx) {
        wave_df <- data.frame(
          time = time_axis * 1000,
          voltage = waveforms[idx, ]
        )
        p <- p + ggplot2::geom_line(data = wave_df,
                                     ggplot2::aes(x = time, y = voltage),
                                     alpha = 0.1, color = "gray50")
      }
    }
  } else {
    rlang::abort("Waveform format not recognized")
  }

  return(p)
}


#' Plot ISI histogram
#'
#' Creates histogram of inter-spike intervals for assessing spike train regularity.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}}
#' @param channel_id Integer. Specific channel to analyze. If NULL, uses all channels
#'   (default: NULL)
#' @param max_isi Numeric. Maximum ISI to display in seconds (default: 1)
#' @param bins Integer. Number of histogram bins (default: 50)
#' @param title Character. Plot title (default: "ISI Distribution")
#' @return ggplot2 object
#' @import ggplot2 rlang
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Plot ISI distribution for all channels
#' plotISIHistogram(spikes, max_isi = 0.5)
#'
#' h5$close()
#' }
plotISIHistogram <- function(spike_data, channel_id = NULL, max_isi = 1,
                             bins = 50, title = "ISI Distribution") {
  # Calculate ISI
  isi_data <- calculateISI(spike_data, channel_id = channel_id)

  # Filter by max ISI
  if ("isi" %in% names(isi_data)) {
    # Single channel
    isi_values <- isi_data$isi[isi_data$isi <= max_isi]
  } else {
    # Multiple channels
    rlang::abort("For multiple channels, specify channel_id for histogram")
  }

  if (length(isi_values) == 0) {
    rlang::abort("No ISI values in specified range")
  }

  plot_data <- data.frame(isi = isi_values * 1000)  # Convert to ms

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = isi)) +
    ggplot2::geom_histogram(bins = bins, fill = "steelblue", color = "white", alpha = 0.8) +
    ggplot2::labs(
      title = title,
      x = "Inter-Spike Interval (ms)",
      y = "Count"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  return(p)
}


#' Plot connectivity matrix
#'
#' Visualizes functional connectivity between channels as a heatmap.
#'
#' @param connectivity_matrix Matrix. Connectivity matrix from \code{\link{buildConnectivityMatrix}}
#' @param title Character. Plot title (default: "Functional Connectivity Matrix")
#' @param cluster Logical. Whether to cluster channels (default: FALSE)
#' @return ggplot2 object
#' @import ggplot2 dplyr rlang
#' @importFrom stats hclust dist
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#' conn <- buildConnectivityMatrix(spikes)
#'
#' # Plot connectivity matrix
#' plotConnectivityMatrix(conn$connectivity_matrix)
#'
#' h5$close()
#' }
plotConnectivityMatrix <- function(connectivity_matrix, title = "Functional Connectivity Matrix",
                                   cluster = FALSE) {
  # Validate input
  if (!is.matrix(connectivity_matrix)) {
    rlang::abort("connectivity_matrix must be a matrix")
  }

  # Reorder if clustering requested
  if (cluster) {
    d <- stats::dist(connectivity_matrix)
    hc <- stats::hclust(d)
    order_idx <- hc$order
    connectivity_matrix <- connectivity_matrix[order_idx, order_idx]
  }

  # Convert to long format for ggplot
  n <- nrow(connectivity_matrix)
  plot_data <- expand.grid(
    Channel1 = 1:n,
    Channel2 = 1:n
  )
  plot_data$Connectivity <- as.vector(connectivity_matrix)

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Channel1, y = Channel2, fill = Connectivity)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, limits = c(-1, 1)
    ) +
    ggplot2::labs(
      title = title,
      x = "Channel ID",
      y = "Channel ID",
      fill = "Correlation"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      aspect.ratio = 1
    ) +
    ggplot2::coord_equal()

  return(p)
}


#' Plot network activity over time
#'
#' Shows population-level activity across the network over time using
#' binned spike counts or firing rates.
#'
#' @param spike_data data.frame. Spike data from \code{\link{bxrSpikedata}}
#' @param bin_size Numeric. Time bin size in seconds (default: 0.1)
#' @param metric Character. Metric to plot: "spike_count" or "firing_rate" (default: "spike_count")
#' @param title Character. Plot title (default: "Network Activity Over Time")
#' @return ggplot2 object
#' @import ggplot2 dplyr rlang
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' spikes <- bxrSpikedata(h5)
#'
#' # Plot network activity
#' plotNetworkActivity(spikes, bin_size = 0.1)
#'
#' h5$close()
#' }
plotNetworkActivity <- function(spike_data, bin_size = 0.1,
                                metric = c("spike_count", "firing_rate"),
                                title = "Network Activity Over Time") {
  # Validate inputs
  if (!is.data.frame(spike_data)) {
    rlang::abort("spike_data must be a data.frame")
  }

  if (!all(c("spike_times") %in% names(spike_data))) {
    rlang::abort("spike_data must contain 'spike_times' column")
  }

  metric <- match.arg(metric)

  # Create time bins
  time_range <- range(spike_data$spike_times)
  bins <- seq(time_range[1], time_range[2], by = bin_size)

  # Bin spikes
  spike_counts <- hist(spike_data$spike_times, breaks = bins, plot = FALSE)$counts
  bin_centers <- (bins[-1] + bins[-length(bins)]) / 2

  if (metric == "firing_rate") {
    values <- spike_counts / bin_size
    y_label <- "Firing Rate (Hz)"
  } else {
    values <- spike_counts
    y_label <- "Spike Count"
  }

  plot_data <- data.frame(
    time = bin_centers,
    value = values
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = value)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 0.8) +
    ggplot2::geom_area(fill = "steelblue", alpha = 0.3) +
    ggplot2::labs(
      title = title,
      x = "Time (s)",
      y = y_label
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )

  return(p)
}
