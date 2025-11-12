# Test quality control functions

# Helper to create spike data with dead channels
create_test_spikes_with_dead <- function() {
  spikes <- data.frame()

  # Active channels with many spikes
  for (ch in 1:5) {
    ch_spikes <- data.frame(
      spike_times = runif(100, 0, 60),
      spike_chid = ch,
      peak_amplitude = rnorm(100, -100, 20)
    )
    spikes <- rbind(spikes, ch_spikes)
  }

  # Low activity channels
  for (ch in 6:7) {
    ch_spikes <- data.frame(
      spike_times = runif(5, 0, 60),
      spike_chid = ch,
      peak_amplitude = rnorm(5, -100, 20)
    )
    spikes <- rbind(spikes, ch_spikes)
  }

  # Dead channels (8-10) have no spikes

  return(spikes)
}


# Test detectDeadChannels
test_that("detectDeadChannels identifies dead channels", {
  spikes <- create_test_spikes_with_dead()
  all_channels <- 1:10

  qc <- detectDeadChannels(spikes, min_spike_count = 10, min_firing_rate = 0.1,
                           all_channels = all_channels)

  expect_type(qc, "list")
  expect_true(all(c("dead_channels", "low_activity_channels", "channel_summary") %in% names(qc)))
  expect_true(qc$n_dead > 0)
  expect_true(all(qc$dead_channels %in% c(6, 7, 8, 9, 10)))
})

test_that("detectDeadChannels classifies channels correctly", {
  spikes <- create_test_spikes_with_dead()

  qc <- detectDeadChannels(spikes, min_spike_count = 10, min_firing_rate = 0.1)

  expect_s3_class(qc$channel_summary, "data.frame")
  expect_true("status" %in% names(qc$channel_summary))
  expect_true(all(qc$channel_summary$status %in% c("active", "low_activity", "dead")))
})

test_that("detectDeadChannels counts active channels", {
  spikes <- create_test_spikes_with_dead()

  qc <- detectDeadChannels(spikes, min_spike_count = 10, all_channels = 1:10)

  expect_true("n_active" %in% names(qc))
  expect_equal(qc$n_active + qc$n_dead + qc$n_low_activity, 10)
})


# Test calculateSNR
test_that("calculateSNR calculates SNR from amplitudes", {
  spikes <- data.frame(
    spike_times = runif(100, 0, 10),
    spike_chid = rep(1:5, each = 20),
    peak_amplitude = c(
      rnorm(20, -100, 10),  # High SNR
      rnorm(20, -50, 5),    # Medium SNR
      rnorm(20, -30, 15),   # Low SNR
      rnorm(20, -80, 8),
      rnorm(20, -60, 12)
    )
  )

  snr <- calculateSNR(spikes, method = "amplitude")

  expect_s3_class(snr, "data.frame")
  expect_true(all(c("channel_id", "signal_mean", "noise_std", "snr_db", "quality") %in% names(snr)))
  expect_equal(nrow(snr), 5)
  expect_true(all(snr$snr_db > 0 | is.na(snr$snr_db)))
})

test_that("calculateSNR classifies quality levels", {
  spikes <- data.frame(
    spike_times = runif(100, 0, 10),
    spike_chid = rep(1:4, each = 25),
    peak_amplitude = rnorm(100, -100, 20)
  )

  snr <- calculateSNR(spikes, method = "amplitude")

  expect_true("quality" %in% names(snr))
  expect_true(all(snr$quality %in% c("poor", "fair", "good", "excellent")))
})


# Test generateQCSummary
test_that("generateQCSummary produces comprehensive report", {
  spikes <- create_test_spikes_with_dead()

  qc <- generateQCSummary(spikes)

  expect_type(qc, "list")
  expect_true(all(c("dead_channels", "snr", "recommendations", "overall_quality") %in% names(qc)))
  expect_true(is.numeric(qc$overall_quality))
  expect_true(qc$overall_quality >= 0 && qc$overall_quality <= 100)
})

test_that("generateQCSummary provides recommendations", {
  spikes <- create_test_spikes_with_dead()

  qc <- generateQCSummary(spikes)

  expect_type(qc$recommendations, "character")
  # Should have recommendations due to dead channels
  expect_true(length(qc$recommendations) > 0)
})

test_that("generateQCSummary quality score reflects data quality", {
  # Good quality data
  good_spikes <- data.frame(
    spike_times = runif(500, 0, 60),
    spike_chid = sample(1:10, 500, replace = TRUE),
    peak_amplitude = rnorm(500, -100, 10)
  )

  qc_good <- generateQCSummary(good_spikes)

  # Poor quality data
  poor_spikes <- data.frame(
    spike_times = runif(20, 0, 60),
    spike_chid = sample(1:2, 20, replace = TRUE),
    peak_amplitude = rnorm(20, -30, 30)
  )

  qc_poor <- generateQCSummary(poor_spikes)

  expect_true(qc_good$overall_quality > qc_poor$overall_quality)
})


# Test error handling
test_that("detectDeadChannels validates input", {
  expect_error(detectDeadChannels(list()), "data.frame")

  invalid_spikes <- data.frame(time = 1:10, value = 1:10)
  expect_error(detectDeadChannels(invalid_spikes), "spike_times.*spike_chid")
})

test_that("calculateSNR requires amplitude column", {
  spikes <- data.frame(
    spike_times = runif(50, 0, 10),
    spike_chid = sample(1:5, 50, replace = TRUE)
  )

  expect_error(calculateSNR(spikes, method = "amplitude"), "peak_amplitude")
})

test_that("generateQCSummary handles minimal data", {
  minimal_spikes <- data.frame(
    spike_times = c(1, 2),
    spike_chid = c(1, 1),
    peak_amplitude = c(-100, -110)
  )

  qc <- generateQCSummary(minimal_spikes)

  expect_type(qc, "list")
  expect_true(qc$overall_quality < 70)  # Should have low quality score
})
