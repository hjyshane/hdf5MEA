# Test file for data extraction and parsing functions

test_that("brwdataParse handles empty binary data", {
  # Test with empty data
  empty_data <- raw(0)
  result <- brwdataParse(empty_data)
  expect_equal(length(result), 0)
})

test_that("brwdataParse validates binary structure", {
  # Test with insufficient data for header
  short_data <- raw(4)  # Less than 8 bytes needed for ChData header
  expect_error(
    brwdataParse(short_data),
    class = "error"
  )
})

test_that("brwtimeseriesConvert handles different modes", {
  # Create mock parsed data
  mock_parsed <- list(
    list(
      channel_id = 1001,
      start_frame = 0,
      end_frame = 100,
      samples = c(0, 0, 50, 100, 0, -50, 0, 0, 200, 0)
    )
  )
  
  sampling_rate <- 1000
  start_frame <- 0
  
  # Test full mode
  result_full <- brwtimeseriesConvert(mock_parsed, sampling_rate, start_frame, mode = "full")
  expect_equal(nrow(result_full[["1001"]]), 10)
  
  # Test events_only mode
  result_events <- brwtimeseriesConvert(mock_parsed, sampling_rate, start_frame, mode = "events_only")
  expect_true(nrow(result_events[["1001"]]) < 10)  # Should have fewer rows
  expect_true(all(result_events[["1001"]]$voltage != 0))  # No zero values
  
  # Test threshold mode
  result_threshold <- brwtimeseriesConvert(mock_parsed, sampling_rate, start_frame, 
                                        mode = "threshold", threshold = 75)
  non_zero_above_threshold <- sum(abs(mock_parsed[[1]]$samples) > 75)
  expect_equal(nrow(result_threshold[["1001"]]), non_zero_above_threshold)
})

test_that("brwtimeseriesConvert calculates time correctly", {
  mock_parsed <- list(
    list(
      channel_id = 2001,
      start_frame = 0,
      end_frame = 5,
      samples = c(1, 2, 3, 4, 5)
    )
  )
  
  sampling_rate <- 1000  # 1000 Hz
  start_frame <- 0
  
  result <- brwtimeseriesConvert(mock_parsed, sampling_rate, start_frame, mode = "full")
  
  # Check time axis
  expected_times <- seq(0, 4) / 1000  # 0, 0.001, 0.002, 0.003, 0.004
  expect_equal(result[["2001"]]$time, expected_times)
  expect_equal(result[["2001"]]$voltage, c(1, 2, 3, 4, 5))
})

test_that("brwtimeseriesConvert handles multiple channels", {
  mock_parsed <- list(
    list(
      channel_id = 3001,
      start_frame = 0,
      end_frame = 3,
      samples = c(10, 20, 30)
    ),
    list(
      channel_id = 3002,
      start_frame = 0,
      end_frame = 2,
      samples = c(100, 200)
    )
  )
  
  result <- brwtimeseriesConvert(mock_parsed, 1000, 0, mode = "full")
  
  expect_equal(length(result), 2)
  expect_true("3001" %in% names(result))
  expect_true("3002" %in% names(result))
  expect_equal(nrow(result[["3001"]]), 3)
  expect_equal(nrow(result[["3002"]]), 2)
})

test_that("eventBased validates target_chunks parameter", {
  skip("Requires mock H5File structure")
  
  # This would test:
  # - Empty target_chunks
  # - Invalid chunk indices
  # - Warning message generation
})

test_that("brwtimeseriesConvert validates input parameters", {
  # Test invalid mode
  mock_parsed <- list(
    list(channel_id = 1, start_frame = 0, end_frame = 1, samples = c(1))
  )
  
  expect_error(
    brwtimeseriesConvert(mock_parsed, 1000, 0, mode = "invalid_mode"),
    "should be one of"
  )
  
  # Test negative threshold
  result <- brwtimeseriesConvert(mock_parsed, 1000, 0, mode = "threshold", threshold = -10)
  expect_s3_class(result[["1"]], "data.frame")  # Should not error
})