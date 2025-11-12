#' Generate example MEA spike data
#'
#' Creates realistic synthetic spike data for testing and examples.
#' Simulates a multi-channel MEA recording with diverse firing patterns
#' including active, low-activity, and silent channels, with bursting behavior.
#'
#' @param n_channels Integer. Number of channels (default: 100)
#' @param duration Numeric. Recording duration in seconds (default: 120)
#' @param active_fraction Numeric. Fraction of active channels (default: 0.7)
#' @param base_rate_range Numeric vector. Range of firing rates in Hz (default: c(0.5, 3))
#' @param burst_probability Numeric. Probability of burst events (default: 0.2)
#' @param seed Integer. Random seed for reproducibility (default: NULL)
#' @return data.frame with spike data in BXR format, containing columns:
#'   \describe{
#'     \item{spike_times}{Spike timestamps in seconds}
#'     \item{spike_chid}{Channel index (1-based)}
#'     \item{spike_units}{Unit classification (0 = unsorted)}
#'     \item{peak_amplitude}{Peak spike amplitude}
#'     \item{min_amplitude}{Minimum spike amplitude}
#'     \item{waveform_index}{Index for waveform lookup}
#'     \item{X}{Grid X coordinate}
#'     \item{Y}{Grid Y coordinate}
#'   }
#' @export
#' @examples
#' # Generate default example data
#' spikes <- generateExampleSpikes()
#'
#' # Quick analysis
#' firing_rates <- calculateFiringRate(spikes)
#' print(head(firing_rates))
#'
#' # Generate smaller dataset
#' small_spikes <- generateExampleSpikes(n_channels = 20, duration = 30)
#'
#' # Generate with specific seed for reproducibility
#' spikes1 <- generateExampleSpikes(seed = 123)
#' spikes2 <- generateExampleSpikes(seed = 123)
#' identical(spikes1, spikes2)  # TRUE
generateExampleSpikes <- function(n_channels = 100,
                                  duration = 120,
                                  active_fraction = 0.7,
                                  base_rate_range = c(0.5, 3),
                                  burst_probability = 0.2,
                                  seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate inputs
  if (!is.numeric(n_channels) || n_channels < 1) {
    rlang::abort("n_channels must be a positive integer")
  }

  if (!is.numeric(duration) || duration <= 0) {
    rlang::abort("duration must be positive")
  }

  spike_data_list <- list()

  # Generate spikes for each channel
  for (ch in 1:n_channels) {
    # Determine channel activity level
    activity_rand <- runif(1)

    if (activity_rand < active_fraction) {
      # Active channels
      base_rate <- runif(1, base_rate_range[1], base_rate_range[2])
    } else if (activity_rand < active_fraction + 0.2) {
      # Low activity channels
      base_rate <- runif(1, 0.05, 0.3)
    } else {
      # Nearly silent channels
      base_rate <- runif(1, 0, 0.05)
    }

    # Expected number of spikes
    n_spikes <- rpois(1, base_rate * duration)

    if (n_spikes > 0) {
      # Mix of isolated spikes and bursts
      n_bursts <- rpois(1, base_rate * duration * burst_probability)

      spike_times <- c()

      # Generate bursts
      if (n_bursts > 0) {
        burst_times <- sort(runif(n_bursts, 0, duration - 1))
        for (bt in burst_times) {
          n_spikes_in_burst <- sample(5:15, 1)
          # Short inter-spike intervals within burst
          isi <- rexp(n_spikes_in_burst, rate = 50)
          burst_spikes <- bt + cumsum(isi)
          spike_times <- c(spike_times, burst_spikes)
        }
      }

      # Add isolated spikes to reach target count
      n_isolated <- max(0, n_spikes - length(spike_times))
      if (n_isolated > 0) {
        isolated_times <- runif(n_isolated, 0, duration)
        spike_times <- c(spike_times, isolated_times)
      }

      # Filter valid times and sort
      spike_times <- spike_times[spike_times >= 0 & spike_times <= duration]
      spike_times <- sort(spike_times)

      if (length(spike_times) > 0) {
        # Generate realistic amplitudes
        ch_mean_amp <- runif(1, -150, -80)
        ch_amp_sd <- runif(1, 10, 25)

        ch_data <- data.frame(
          spike_times = spike_times,
          spike_chid = ch,
          spike_units = 0,  # All unsorted
          peak_amplitude = rnorm(length(spike_times),
                                mean = ch_mean_amp,
                                sd = ch_amp_sd),
          min_amplitude = rnorm(length(spike_times),
                               mean = ch_mean_amp - runif(1, 20, 50),
                               sd = ch_amp_sd * 1.5),
          waveform_index = seq_along(spike_times)
        )

        spike_data_list[[length(spike_data_list) + 1]] <- ch_data
      }
    }
  }

  # Combine all spike data
  if (length(spike_data_list) == 0) {
    # Return empty data frame with correct structure
    return(data.frame(
      spike_times = numeric(),
      spike_chid = integer(),
      spike_units = integer(),
      peak_amplitude = numeric(),
      min_amplitude = numeric(),
      waveform_index = integer(),
      X = integer(),
      Y = integer()
    ))
  }

  spikes <- do.call(rbind, spike_data_list)
  rownames(spikes) <- NULL

  # Add grid coordinates (square grid)
  grid_size <- ceiling(sqrt(n_channels))
  spikes$X <- (spikes$spike_chid - 1) %% grid_size + 1
  spikes$Y <- (spikes$spike_chid - 1) %/% grid_size + 1

  # Sort by time
  spikes <- spikes[order(spikes$spike_times), ]

  return(spikes)
}


#' Generate example burst data
#'
#' Creates synthetic burst data compatible with burst analysis functions.
#'
#' @param n_channels Integer. Number of channels (default: 20)
#' @param n_bursts_per_channel Integer. Average bursts per channel (default: 10)
#' @param duration Numeric. Recording duration in seconds (default: 60)
#' @param seed Integer. Random seed (default: NULL)
#' @return data.frame with burst data
#' @export
#' @examples
#' bursts <- generateExampleBursts()
#' burst_stats <- calculateBurstStatistics(bursts)
generateExampleBursts <- function(n_channels = 20,
                                  n_bursts_per_channel = 10,
                                  duration = 60,
                                  seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  bursts <- data.frame()

  for (ch in 1:n_channels) {
    n_bursts <- rpois(1, n_bursts_per_channel)

    if (n_bursts > 0) {
      burst_starts <- sort(runif(n_bursts, 0, duration - 1))
      burst_durations <- rexp(n_bursts, rate = 5)  # Exponential duration
      burst_durations <- pmin(burst_durations, 1)  # Cap at 1 second

      ch_bursts <- data.frame(
        channel_id = ch,
        burst_id = 1:n_bursts,
        start_time = burst_starts,
        end_time = burst_starts + burst_durations,
        duration = burst_durations,
        n_spikes = rpois(n_bursts, 10),
        X = (ch - 1) %% 5 + 1,
        Y = (ch - 1) %/% 5 + 1
      )

      bursts <- rbind(bursts, ch_bursts)
    }
  }

  return(bursts)
}


#' Generate example voltage trace
#'
#' Creates synthetic raw voltage data with spikes and noise for signal
#' processing examples.
#'
#' @param duration Numeric. Duration in seconds (default: 1)
#' @param sampling_rate Numeric. Sampling rate in Hz (default: 10000)
#' @param n_spikes Integer. Number of spikes to add (default: 20)
#' @param noise_level Numeric. Noise standard deviation (default: 10)
#' @param seed Integer. Random seed (default: NULL)
#' @return data.frame with columns: time, voltage
#' @export
#' @examples
#' # Generate voltage trace
#' voltage <- generateExampleVoltage(duration = 2, n_spikes = 30)
#'
#' # Apply filtering
#' filtered <- applyBandpassFilter(voltage, sampling_rate = 10000,
#'                                   low_freq = 300, high_freq = 3000)
#'
#' # Detect spikes
#' timeseries <- list("1" = voltage)
#' detected <- detectSpikes(timeseries, sampling_rate = 10000)
generateExampleVoltage <- function(duration = 1,
                                   sampling_rate = 10000,
                                   n_spikes = 20,
                                   noise_level = 10,
                                   seed = NULL) {

  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_samples <- ceiling(duration * sampling_rate)
  time <- seq(0, duration, length.out = n_samples)

  # Base noise
  voltage <- rnorm(n_samples, mean = 0, sd = noise_level)

  # Add spikes
  if (n_spikes > 0) {
    spike_times <- sample(100:(n_samples - 100), min(n_spikes, n_samples - 200))

    for (st in spike_times) {
      # Spike waveform (simplified)
      spike_window <- 20
      if (st + spike_window <= n_samples) {
        t <- 1:spike_window
        # Negative spike shape
        amplitude <- runif(1, -200, -100)
        spike_shape <- amplitude * exp(-((t - 5)^2) / 10)
        voltage[st:(st + spike_window - 1)] <- voltage[st:(st + spike_window - 1)] + spike_shape
      }
    }
  }

  data.frame(time = time, voltage = voltage)
}
