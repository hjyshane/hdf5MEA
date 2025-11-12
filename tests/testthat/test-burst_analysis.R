# Test burst analysis functions

# Helper function to create mock burst data
create_mock_bursts <- function(n_channels = 5, n_bursts_per_channel = 10, duration = 60) {
  bursts <- data.frame()

  for (ch in 1:n_channels) {
    burst_starts <- sort(runif(n_bursts_per_channel, 0, duration - 1))
    burst_durations <- runif(n_bursts_per_channel, 0.05, 0.5)

    ch_bursts <- data.frame(
      channel_id = ch,
      burst_id = 1:n_bursts_per_channel,
      start_time = burst_starts,
      end_time = burst_starts + burst_durations,
      duration = burst_durations,
      n_spikes = sample(5:20, n_bursts_per_channel, replace = TRUE),
      X = (ch - 1) %% 5 + 1,
      Y = (ch - 1) %/% 5 + 1
    )
    bursts <- rbind(bursts, ch_bursts)
  }

  return(bursts)
}

# Helper to create spike data with bursts
create_bursting_spikes <- function(n_channels = 5, duration = 60) {
  spikes <- data.frame()

  for (ch in 1:n_channels) {
    # Create burst events
    n_bursts <- 10
    burst_times <- sort(runif(n_bursts, 0, duration - 1))

    for (bt in burst_times) {
      # Spikes within burst (short ISI)
      n_spikes_in_burst <- sample(5:15, 1)
      burst_spike_times <- bt + cumsum(runif(n_spikes_in_burst, 0.001, 0.02))

      burst_spikes <- data.frame(
        spike_times = burst_spike_times,
        spike_chid = ch
      )
      spikes <- rbind(spikes, burst_spikes)
    }

    # Add some isolated spikes
    isolated_times <- runif(20, 0, duration)
    isolated_spikes <- data.frame(
      spike_times = isolated_times,
      spike_chid = ch
    )
    spikes <- rbind(spikes, isolated_spikes)
  }

  return(spikes[order(spikes$spike_times), ])
}


# Test detectBursts
test_that("detectBursts finds bursts in spike data", {
  spikes <- create_bursting_spikes(n_channels = 3, duration = 30)

  bursts <- detectBursts(spikes, isi_threshold = 0.1, min_spikes = 5)

  expect_s3_class(bursts, "data.frame")
  if (nrow(bursts) > 0) {
    expect_true(all(c("channel_id", "start_time", "end_time", "duration",
                     "n_spikes", "intra_burst_rate") %in% names(bursts)))
    expect_true(all(bursts$n_spikes >= 5))
    expect_true(all(bursts$duration > 0))
  }
})

test_that("detectBursts respects min_spikes parameter", {
  spikes <- create_bursting_spikes(n_channels = 2, duration = 20)

  bursts_min5 <- detectBursts(spikes, isi_threshold = 0.1, min_spikes = 5)
  bursts_min10 <- detectBursts(spikes, isi_threshold = 0.1, min_spikes = 10)

  if (nrow(bursts_min5) > 0 && nrow(bursts_min10) > 0) {
    expect_true(all(bursts_min5$n_spikes >= 5))
    expect_true(all(bursts_min10$n_spikes >= 10))
    expect_true(nrow(bursts_min10) <= nrow(bursts_min5))
  }
})


# Test calculateBurstStatistics
test_that("calculateBurstStatistics computes correct metrics", {
  bursts <- create_mock_bursts(n_channels = 5, n_bursts_per_channel = 10)

  stats <- calculateBurstStatistics(bursts)

  expect_s3_class(stats, "data.frame")
  expect_true(all(c("channel_id", "n_bursts", "mean_burst_duration",
                   "mean_ibi", "burst_frequency") %in% names(stats)))
  expect_equal(nrow(stats), 5)
  expect_true(all(stats$n_bursts == 10))
})

test_that("calculateBurstStatistics handles BXR naming convention", {
  bursts <- create_mock_bursts(n_channels = 3, n_bursts_per_channel = 5)
  names(bursts)[names(bursts) == "channel_id"] <- "burst_chidx"

  stats <- calculateBurstStatistics(bursts)

  expect_s3_class(stats, "data.frame")
  expect_true("channel_id" %in% names(stats))
})


# Test analyzeBurstPropagation
test_that("analyzeBurstPropagation detects coordinated bursts", {
  # Create synchronized bursts across channels
  bursts <- data.frame()
  event_time <- 10
  for (ch in 1:5) {
    ch_burst <- data.frame(
      channel_id = ch,
      start_time = event_time + rnorm(1, 0, 0.01),
      end_time = event_time + 0.2,
      duration = 0.2,
      X = ch,
      Y = 1
    )
    bursts <- rbind(bursts, ch_burst)
  }

  propagation <- analyzeBurstPropagation(bursts, time_window = 0.05)

  expect_s3_class(propagation, "data.frame")
  if (nrow(propagation) > 0) {
    expect_true(all(c("propagation_id", "n_channels", "duration") %in% names(propagation)))
  }
})


# Test calculateBurstSynchrony
test_that("calculateBurstSynchrony calculates synchrony index", {
  bursts <- create_mock_bursts(n_channels = 5, n_bursts_per_channel = 20, duration = 60)

  synchrony <- calculateBurstSynchrony(bursts, bin_size = 0.5)

  expect_type(synchrony, "list")
  expect_true(all(c("synchrony_index", "burst_coincidence") %in% names(synchrony)))
  expect_true(synchrony$synchrony_index >= 0 && synchrony$synchrony_index <= 1)
  expect_s3_class(synchrony$burst_coincidence, "data.frame")
})


# Test edge cases
test_that("detectBursts handles data with no bursts", {
  # Create spikes with long ISIs
  spikes <- data.frame(
    spike_times = seq(0, 60, by = 2),
    spike_chid = 1
  )

  bursts <- detectBursts(spikes, isi_threshold = 0.1, min_spikes = 5)

  expect_s3_class(bursts, "data.frame")
  expect_equal(nrow(bursts), 0)
})

test_that("calculateBurstStatistics handles empty input", {
  empty_bursts <- data.frame()

  stats <- calculateBurstStatistics(empty_bursts)

  expect_s3_class(stats, "data.frame")
  expect_equal(nrow(stats), 0)
})
