# Test spike train analysis functions

# Helper function to create mock spike data
create_mock_spikes <- function(n_channels = 10, n_spikes_per_channel = 100, duration = 60) {
  spikes <- data.frame()

  for (ch in 1:n_channels) {
    spike_times <- sort(runif(n_spikes_per_channel, 0, duration))
    ch_spikes <- data.frame(
      spike_times = spike_times,
      spike_chid = ch,
      peak_amplitude = rnorm(n_spikes_per_channel, mean = -100, sd = 20),
      X = (ch - 1) %% 5 + 1,
      Y = (ch - 1) %/% 5 + 1
    )
    spikes <- rbind(spikes, ch_spikes)
  }

  return(spikes)
}


# Test calculateFiringRate
test_that("calculateFiringRate works with overall rate", {
  spikes <- create_mock_spikes(n_channels = 5, n_spikes_per_channel = 50, duration = 10)

  rates <- calculateFiringRate(spikes)

  expect_s3_class(rates, "data.frame")
  expect_true(all(c("channel_id", "firing_rate", "spike_count", "duration") %in% names(rates)))
  expect_equal(nrow(rates), 5)
  expect_true(all(rates$firing_rate > 0))
  expect_true(all(rates$spike_count == 50))
})

test_that("calculateFiringRate works with binned rate", {
  spikes <- create_mock_spikes(n_channels = 3, n_spikes_per_channel = 100, duration = 10)

  rates <- calculateFiringRate(spikes, bin_size = 1)

  expect_s3_class(rates, "data.frame")
  expect_true(all(c("channel_id", "time_bin_center", "firing_rate") %in% names(rates)))
  expect_true(nrow(rates) > 0)
})

test_that("calculateFiringRate handles time_range", {
  spikes <- create_mock_spikes(duration = 60)

  rates <- calculateFiringRate(spikes, time_range = c(10, 20))

  expect_s3_class(rates, "data.frame")
  expect_equal(rates$duration[1], 10)
})


# Test calculateISI
test_that("calculateISI works for all channels", {
  spikes <- create_mock_spikes(n_channels = 5, n_spikes_per_channel = 50)

  isi <- calculateISI(spikes)

  expect_s3_class(isi, "data.frame")
  expect_true(all(c("channel_id", "mean_isi", "median_isi", "cv_isi") %in% names(isi)))
  expect_equal(nrow(isi), 5)
  expect_true(all(isi$mean_isi > 0, na.rm = TRUE))
})

test_that("calculateISI works for single channel", {
  spikes <- create_mock_spikes(n_channels = 5, n_spikes_per_channel = 50)

  isi <- calculateISI(spikes, channel_id = 1)

  expect_s3_class(isi, "data.frame")
  expect_true("isi" %in% names(isi))
  expect_true(all(isi$channel_id == 1))
})


# Test getSpikeStatistics
test_that("getSpikeStatistics returns expected metrics", {
  spikes <- create_mock_spikes(n_channels = 10, n_spikes_per_channel = 100, duration = 60)

  stats <- getSpikeStatistics(spikes)

  expect_type(stats, "list")
  expect_true(all(c("total_spikes", "n_active_channels", "recording_duration",
                   "mean_firing_rate", "global_firing_rate") %in% names(stats)))
  expect_equal(stats$total_spikes, 1000)
  expect_equal(stats$n_active_channels, 10)
  expect_true(stats$mean_firing_rate > 0)
})


# Test calculateSpikeAutocorrelation
test_that("calculateSpikeAutocorrelation works", {
  spikes <- create_mock_spikes(n_channels = 5, n_spikes_per_channel = 200, duration = 60)

  autocorr <- calculateSpikeAutocorrelation(spikes, channel_id = 1, max_lag = 0.5)

  expect_s3_class(autocorr, "data.frame")
  expect_true(all(c("lag", "autocorr") %in% names(autocorr)))
  expect_true(nrow(autocorr) > 0)
  expect_true(max(autocorr$lag) <= 0.5)
})


# Test error handling
test_that("calculateFiringRate handles empty data", {
  empty_spikes <- data.frame(spike_times = numeric(), spike_chid = integer())

  expect_error(calculateFiringRate(empty_spikes), "empty")
})

test_that("calculateISI handles insufficient spikes", {
  spikes <- data.frame(
    spike_times = c(1),
    spike_chid = c(1)
  )

  expect_error(calculateISI(spikes, channel_id = 1), "fewer than 2")
})
