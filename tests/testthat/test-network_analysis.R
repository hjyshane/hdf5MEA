# Test network analysis functions

# Helper to create correlated spike data
create_correlated_spikes <- function(n_channels = 5, duration = 60) {
  spikes <- data.frame()

  # Create some synchronized events
  n_events <- 50
  event_times <- sort(runif(n_events, 0, duration))

  for (event_time in event_times) {
    # Randomly select subset of channels to fire together
    active_channels <- sample(1:n_channels, size = sample(2:n_channels, 1))

    for (ch in active_channels) {
      spike_time <- event_time + rnorm(1, 0, 0.01)  # Small jitter
      spikes <- rbind(spikes, data.frame(
        spike_times = spike_time,
        spike_chid = ch,
        X = ch,
        Y = 1
      ))
    }
  }

  # Add some background spikes
  for (ch in 1:n_channels) {
    bg_times <- runif(20, 0, duration)
    bg_spikes <- data.frame(
      spike_times = bg_times,
      spike_chid = ch,
      X = ch,
      Y = 1
    )
    spikes <- rbind(spikes, bg_spikes)
  }

  return(spikes[order(spikes$spike_times), ])
}


# Test buildConnectivityMatrix
test_that("buildConnectivityMatrix creates valid matrix", {
  spikes <- create_correlated_spikes(n_channels = 5, duration = 30)

  conn <- buildConnectivityMatrix(spikes, bin_size = 0.05, threshold = 0.2)

  expect_type(conn, "list")
  expect_true(all(c("connectivity_matrix", "adjacency_matrix", "network_density") %in% names(conn)))
  expect_true(is.matrix(conn$connectivity_matrix))
  expect_equal(nrow(conn$connectivity_matrix), 5)
  expect_equal(ncol(conn$connectivity_matrix), 5)
  expect_true(conn$network_density >= 0 && conn$network_density <= 1)
})

test_that("buildConnectivityMatrix diagonal is zero", {
  spikes <- create_correlated_spikes(n_channels = 4, duration = 20)

  conn <- buildConnectivityMatrix(spikes, threshold = 0.3)

  expect_true(all(diag(conn$connectivity_matrix) == 0))
})


# Test calculateNetworkSynchrony
test_that("calculateNetworkSynchrony computes synchrony index", {
  spikes <- create_correlated_spikes(n_channels = 6, duration = 30)

  sync <- calculateNetworkSynchrony(spikes, bin_size = 0.1, method = "coincidence")

  expect_type(sync, "list")
  expect_true(all(c("synchrony_index", "temporal_profile") %in% names(sync)))
  expect_true(sync$synchrony_index >= 0 && sync$synchrony_index <= 1)
  expect_s3_class(sync$temporal_profile, "data.frame")
})

test_that("calculateNetworkSynchrony works with correlation method", {
  spikes <- create_correlated_spikes(n_channels = 5, duration = 20)

  sync <- calculateNetworkSynchrony(spikes, bin_size = 0.1, method = "correlation")

  expect_type(sync, "list")
  expect_true("synchrony_index" %in% names(sync))
})


# Test detectSpatialPatterns
test_that("detectSpatialPatterns identifies hotspots", {
  spikes <- create_correlated_spikes(n_channels = 8, duration = 40)

  patterns <- detectSpatialPatterns(spikes, time_window = 0.2)

  expect_type(patterns, "list")
  expect_true(all(c("spatial_events", "hotspots", "wave_events") %in% names(patterns)))
  expect_s3_class(patterns$hotspots, "data.frame")
  if (nrow(patterns$hotspots) > 0) {
    expect_true(all(c("spike_chid", "firing_rate") %in% names(patterns$hotspots)))
  }
})


# Test calculateCrossCorrelation
test_that("calculateCrossCorrelation works with channel pairs", {
  spikes <- create_correlated_spikes(n_channels = 5, duration = 30)

  pairs <- data.frame(ch1 = c(1, 1), ch2 = c(2, 3))
  corr <- calculateCrossCorrelation(spikes, channel_pairs = pairs, max_lag = 0.05)

  expect_s3_class(corr, "data.frame")
  if (nrow(corr) > 0) {
    expect_true(all(c("channel_1", "channel_2", "peak_correlation", "peak_lag") %in% names(corr)))
    expect_equal(nrow(corr), 2)
  }
})


# Test error handling
test_that("buildConnectivityMatrix handles few channels", {
  spikes <- data.frame(
    spike_times = runif(50, 0, 10),
    spike_chid = 1
  )

  expect_error(buildConnectivityMatrix(spikes), "at least 2")
})

test_that("calculateNetworkSynchrony handles insufficient data", {
  spikes <- data.frame(
    spike_times = c(1, 2),
    spike_chid = c(1, 2)
  )

  expect_error(calculateNetworkSynchrony(spikes), "at least 2")
})

test_that("detectSpatialPatterns requires spatial coordinates", {
  spikes <- data.frame(
    spike_times = runif(100, 0, 10),
    spike_chid = sample(1:5, 100, replace = TRUE)
  )

  expect_error(detectSpatialPatterns(spikes), "X.*Y")
})
