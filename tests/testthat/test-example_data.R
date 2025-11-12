# Test example data generation functions

test_that("generateExampleSpikes creates valid spike data", {
  spikes <- generateExampleSpikes(n_channels = 20, duration = 30, seed = 123)

  expect_s3_class(spikes, "data.frame")
  expect_true(all(c("spike_times", "spike_chid", "peak_amplitude", "X", "Y") %in% names(spikes)))
  expect_true(nrow(spikes) > 0)
  expect_true(all(spikes$spike_times >= 0 & spikes$spike_times <= 30))
  expect_true(all(spikes$spike_chid >= 1 & spikes$spike_chid <= 20))
})

test_that("generateExampleSpikes is reproducible with seed", {
  spikes1 <- generateExampleSpikes(n_channels = 10, duration = 10, seed = 42)
  spikes2 <- generateExampleSpikes(n_channels = 10, duration = 10, seed = 42)

  expect_equal(spikes1, spikes2)
})

test_that("generateExampleSpikes validates parameters", {
  expect_error(generateExampleSpikes(n_channels = -1), "positive")
  expect_error(generateExampleSpikes(duration = 0), "positive")
})

test_that("generateExampleSpikes creates diverse activity levels", {
  spikes <- generateExampleSpikes(n_channels = 100, duration = 60, seed = 123)

  # Should have variety in firing rates
  rates <- calculateFiringRate(spikes)
  expect_true(sd(rates$firing_rate) > 0)

  # Should have some low-activity or silent channels
  expect_true(min(rates$spike_count) < max(rates$spike_count) * 0.5)
})

test_that("generateExampleSpikes works with different parameters", {
  # Small dataset
  small <- generateExampleSpikes(n_channels = 5, duration = 10)
  expect_true(nrow(small) > 0)

  # Large dataset
  large <- generateExampleSpikes(n_channels = 200, duration = 300)
  expect_true(nrow(large) > nrow(small))
})


test_that("generateExampleBursts creates valid burst data", {
  bursts <- generateExampleBursts(n_channels = 10, n_bursts_per_channel = 5, seed = 123)

  expect_s3_class(bursts, "data.frame")
  expect_true(all(c("channel_id", "start_time", "end_time", "duration", "n_spikes") %in% names(bursts)))
  expect_true(all(bursts$end_time > bursts$start_time))
  expect_true(all(bursts$duration == bursts$end_time - bursts$start_time))
})

test_that("generateExampleBursts is reproducible", {
  bursts1 <- generateExampleBursts(seed = 100)
  bursts2 <- generateExampleBursts(seed = 100)

  expect_equal(bursts1, bursts2)
})


test_that("generateExampleVoltage creates valid voltage data", {
  voltage <- generateExampleVoltage(duration = 1, sampling_rate = 10000, n_spikes = 10, seed = 123)

  expect_s3_class(voltage, "data.frame")
  expect_true(all(c("time", "voltage") %in% names(voltage)))
  expect_equal(nrow(voltage), 10000)
  expect_true(min(voltage$time) >= 0)
  expect_true(max(voltage$time) <= 1)
})

test_that("generateExampleVoltage is reproducible", {
  voltage1 <- generateExampleVoltage(duration = 0.5, seed = 50)
  voltage2 <- generateExampleVoltage(duration = 0.5, seed = 50)

  expect_equal(voltage1, voltage2)
})

test_that("generateExampleVoltage has spikes", {
  voltage <- generateExampleVoltage(duration = 1, n_spikes = 50, seed = 123)

  # Should have negative deflections (spikes)
  expect_true(min(voltage$voltage) < -50)
})


# Integration tests with analysis functions
test_that("generated spike data works with analysis functions", {
  spikes <- generateExampleSpikes(n_channels = 20, duration = 30, seed = 123)

  # Should work with firing rate calculation
  rates <- calculateFiringRate(spikes)
  expect_s3_class(rates, "data.frame")

  # Should work with ISI calculation
  isi <- calculateISI(spikes)
  expect_s3_class(isi, "data.frame")

  # Should work with burst detection
  bursts <- detectBursts(spikes, isi_threshold = 0.1, min_spikes = 5)
  expect_s3_class(bursts, "data.frame")
})

test_that("generated burst data works with burst analysis", {
  bursts <- generateExampleBursts(n_channels = 10, seed = 123)

  stats <- calculateBurstStatistics(bursts)
  expect_s3_class(stats, "data.frame")
})

test_that("generated voltage data works with signal processing", {
  voltage <- generateExampleVoltage(duration = 1, n_spikes = 20, seed = 123)

  # Should work with filtering
  filtered <- applyBandpassFilter(voltage, sampling_rate = 10000,
                                   low_freq = 300, high_freq = 3000)
  expect_s3_class(filtered, "data.frame")

  # Should work with spike detection
  timeseries <- list("1" = voltage)
  detected <- detectSpikes(timeseries, sampling_rate = 10000)
  expect_s3_class(detected, "data.frame")
})
