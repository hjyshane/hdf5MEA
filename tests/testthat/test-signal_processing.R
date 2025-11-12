# Test signal processing functions

# Helper to create mock voltage data
create_mock_voltage <- function(n_samples = 1000, sampling_rate = 10000, add_spikes = TRUE) {
  time <- seq(0, (n_samples - 1) / sampling_rate, by = 1 / sampling_rate)

  # Base noise
  voltage <- rnorm(n_samples, mean = 0, sd = 10)

  # Add spikes if requested
  if (add_spikes) {
    spike_times <- sample(100:(n_samples-100), 20)
    for (st in spike_times) {
      # Spike shape
      spike_window <- 20
      spike_shape <- -200 * exp(-((1:spike_window) - 10)^2 / 10)
      voltage[st:(st + spike_window - 1)] <- voltage[st:(st + spike_window - 1)] + spike_shape
    }
  }

  data.frame(time = time, voltage = voltage)
}

create_mock_timeseries <- function(n_channels = 3, n_samples = 1000) {
  timeseries <- list()
  for (i in 1:n_channels) {
    timeseries[[as.character(i)]] <- create_mock_voltage(n_samples)
  }
  timeseries
}


# Test detectSpikes
test_that("detectSpikes finds spikes in voltage data", {
  timeseries <- create_mock_timeseries(n_channels = 2, n_samples = 2000)

  spikes <- detectSpikes(timeseries, sampling_rate = 10000, method = "mad")

  expect_s3_class(spikes, "data.frame")
  expect_true(all(c("channel_id", "spike_time", "spike_amplitude") %in% names(spikes)))
})

test_that("detectSpikes respects polarity", {
  voltage_data <- create_mock_voltage(n_samples = 1000, add_spikes = TRUE)
  timeseries <- list("1" = voltage_data)

  neg_spikes <- detectSpikes(timeseries, sampling_rate = 10000, polarity = "negative")

  expect_s3_class(neg_spikes, "data.frame")
  if (nrow(neg_spikes) > 0) {
    expect_true(all(neg_spikes$spike_amplitude < 0))
  }
})


# Test applyBandpassFilter
test_that("applyBandpassFilter filters voltage data", {
  voltage_data <- create_mock_voltage(n_samples = 2000)

  filtered <- applyBandpassFilter(voltage_data, sampling_rate = 10000,
                                   low_freq = 300, high_freq = 3000)

  expect_s3_class(filtered, "data.frame")
  expect_equal(nrow(filtered), nrow(voltage_data))
  expect_true("voltage" %in% names(filtered))
})

test_that("applyBandpassFilter works with vector input", {
  voltage_vec <- rnorm(1000)

  filtered <- applyBandpassFilter(voltage_vec, sampling_rate = 10000,
                                   low_freq = 300, high_freq = 3000)

  expect_type(filtered, "double")
  expect_equal(length(filtered), length(voltage_vec))
})


# Test removeBaselineDrift
test_that("removeBaselineDrift removes drift with median method", {
  voltage_data <- create_mock_voltage(n_samples = 1000)
  # Add linear drift
  voltage_data$voltage <- voltage_data$voltage + seq(0, 100, length.out = 1000)

  corrected <- removeBaselineDrift(voltage_data, method = "median",
                                    window_size = 0.1, sampling_rate = 10000)

  expect_s3_class(corrected, "data.frame")
  expect_equal(nrow(corrected), nrow(voltage_data))
  expect_true(mean(corrected$voltage) < mean(voltage_data$voltage))
})

test_that("removeBaselineDrift removes drift with linear method", {
  voltage_data <- create_mock_voltage(n_samples = 1000)
  # Add linear drift
  voltage_data$voltage <- voltage_data$voltage + seq(0, 50, length.out = 1000)

  corrected <- removeBaselineDrift(voltage_data, method = "linear")

  expect_s3_class(corrected, "data.frame")
  expect_equal(nrow(corrected), nrow(voltage_data))
})


# Test detectArtifacts
test_that("detectArtifacts finds and handles artifacts", {
  voltage_data <- create_mock_voltage(n_samples = 1000)
  # Add artifacts
  voltage_data$voltage[100:110] <- 2000  # Large spike
  voltage_data$voltage[500:505] <- -2000

  result <- detectArtifacts(voltage_data, artifact_threshold = 500,
                            sampling_rate = 10000, action = "interpolate")

  expect_type(result, "list")
  expect_true(all(c("cleaned_data", "artifact_times", "n_artifacts") %in% names(result)))
  expect_s3_class(result$cleaned_data, "data.frame")
  expect_true(result$n_artifacts > 0)
})

test_that("detectArtifacts works with list input", {
  timeseries <- create_mock_timeseries(n_channels = 2, n_samples = 1000)
  # Add artifacts to first channel
  timeseries[["1"]]$voltage[100:105] <- 2000

  result <- detectArtifacts(timeseries, artifact_threshold = 500,
                            sampling_rate = 10000, action = "flag")

  expect_type(result, "list")
  expect_true("cleaned_data" %in% names(result))
})


# Test error handling
test_that("detectSpikes handles empty timeseries", {
  empty_timeseries <- list()

  spikes <- detectSpikes(empty_timeseries, sampling_rate = 10000)

  expect_s3_class(spikes, "data.frame")
  expect_equal(nrow(spikes), 0)
})

test_that("applyBandpassFilter validates frequency range", {
  voltage_data <- create_mock_voltage(n_samples = 100)

  expect_error(applyBandpassFilter(voltage_data, sampling_rate = 10000,
                                    low_freq = 3000, high_freq = 300),
               "low_freq must be less")
})

test_that("applyBandpassFilter validates Nyquist frequency", {
  voltage_data <- create_mock_voltage(n_samples = 100)

  expect_error(applyBandpassFilter(voltage_data, sampling_rate = 1000,
                                    low_freq = 100, high_freq = 600),
               "Nyquist")
})
