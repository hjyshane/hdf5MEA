# Generate example spike data for package examples
# This creates realistic synthetic MEA data

set.seed(123)

# Generate example spike data (realistic MEA recording)
n_channels <- 100
recording_duration <- 120  # 2 minutes
sampling_rate <- 10000

# Create realistic spike trains with different firing patterns
spike_data_list <- list()

for (ch in 1:n_channels) {
  # Determine channel activity level
  if (ch <= 70) {
    # Active channels (70%)
    base_rate <- runif(1, 0.5, 3)  # 0.5-3 Hz base rate
  } else if (ch <= 90) {
    # Low activity channels (20%)
    base_rate <- runif(1, 0.05, 0.3)
  } else {
    # Nearly silent channels (10%)
    base_rate <- runif(1, 0, 0.05)
  }

  # Generate spike times with some bursting
  n_spikes <- rpois(1, base_rate * recording_duration)

  if (n_spikes > 0) {
    # Mix of isolated spikes and bursts
    n_bursts <- rpois(1, base_rate * recording_duration * 0.2)

    spike_times <- c()

    # Add bursts
    if (n_bursts > 0) {
      burst_times <- sort(runif(n_bursts, 0, recording_duration - 1))
      for (bt in burst_times) {
        n_spikes_in_burst <- sample(5:15, 1)
        isi <- rexp(n_spikes_in_burst, rate = 50)  # Short ISI
        burst_spikes <- bt + cumsum(isi)
        spike_times <- c(spike_times, burst_spikes)
      }
    }

    # Add isolated spikes
    n_isolated <- max(0, n_spikes - length(spike_times))
    if (n_isolated > 0) {
      isolated_times <- runif(n_isolated, 0, recording_duration)
      spike_times <- c(spike_times, isolated_times)
    }

    # Filter valid times
    spike_times <- spike_times[spike_times >= 0 & spike_times <= recording_duration]
    spike_times <- sort(spike_times)

    if (length(spike_times) > 0) {
      # Create spike data frame
      ch_data <- data.frame(
        spike_times = spike_times,
        spike_chid = ch,
        spike_units = 0,  # All unsorted
        peak_amplitude = rnorm(length(spike_times),
                              mean = runif(1, -150, -80),
                              sd = runif(1, 10, 25)),
        min_amplitude = rnorm(length(spike_times),
                             mean = runif(1, -200, -120),
                             sd = runif(1, 15, 30)),
        waveform_index = seq_along(spike_times)
      )

      spike_data_list[[length(spike_data_list) + 1]] <- ch_data
    }
  }
}

# Combine all spike data
example_spikes <- do.call(rbind, spike_data_list)
rownames(example_spikes) <- NULL

# Add grid coordinates (10x10 grid)
grid_size <- 10
example_spikes$X <- (example_spikes$spike_chid - 1) %% grid_size + 1
example_spikes$Y <- (example_spikes$spike_chid - 1) %/% grid_size + 1

# Sort by time
example_spikes <- example_spikes[order(example_spikes$spike_times), ]

# Save
usethis::use_data(example_spikes, overwrite = TRUE)

cat("Generated example_spikes with:\n")
cat("  - Channels:", length(unique(example_spikes$spike_chid)), "\n")
cat("  - Total spikes:", nrow(example_spikes), "\n")
cat("  - Duration:", max(example_spikes$spike_times), "seconds\n")
cat("  - Mean firing rate:", nrow(example_spikes) / max(example_spikes$spike_times), "Hz\n")
