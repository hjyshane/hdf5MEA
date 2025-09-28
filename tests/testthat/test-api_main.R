# Test file for main API functions

test_that("get_brw_data validates input parameters", {
  # Test negative start time
  expect_error(
    get_brw_data("dummy.brw", start = -1, duration = 5),
    "Time parameter is wrong"
  )
  
  # Test zero duration
  expect_error(
    get_brw_data("dummy.brw", start = 0, duration = 0),
    "Time parameter is wrong"
  )
  
  # Test negative duration
  expect_error(
    get_brw_data("dummy.brw", start = 0, duration = -1),
    "Time parameter is wrong"
  )
})

test_that("get_brw_data warns about memory usage", {
  # Create a temporary valid file for testing warning
  temp_file <- tempfile(fileext = ".brw")
  file.create(temp_file)
  
  # Test large duration warning (will fail at file opening, but warning should come first)
  # We need to test the warning logic separately
  expect_warning(
    tryCatch(get_brw_data(temp_file, start = 0, duration = 10), error = function(e) NULL),
    "may run out of memory"
  )
  
  file.remove(temp_file)
})

test_that("get_brw_data handles file operations correctly", {
  # Test non-existent file
  expect_error(
    get_brw_data("nonexistent_file.brw", start = 0, duration = 1),
    "File path is invalid"
  )
})

test_that("get_brw_data returns correct structure", {
  skip("Requires test BRW file")
  
  # This test would verify that get_brw_data returns a list with:
  # - binary_chunk
  # - stored_channels  
  # - data_range
  # - sampling_rate
  # - start_frame
  # - start_frames
  # - num_frames
  # - end_frame
})

# Integration tests
test_that("full pipeline works with mock data", {
  skip("Integration test - requires comprehensive mock setup")
  
  # This would test the complete pipeline:
  # get_brw_data() -> dataParse() -> timeseriesConvert()
  # with realistic mock data
})

test_that("sparse modes produce consistent results", {
  # Create test data with known patterns
  test_samples <- c(0, 0, 100, 0, 200, 0, 0, 300, 0, 0)
  mock_parsed <- list(
    list(
      channel_id = 9999,
      start_frame = 0,
      end_frame = length(test_samples),
      samples = test_samples
    )
  )
  
  # Test all modes
  full_result <- brwtimeseriesConvert(mock_parsed, 1000, 0, mode = "full")
  events_result <- brwtimeseriesConvert(mock_parsed, 1000, 0, mode = "events_only")
  threshold_result <- brwtimeseriesConvert(mock_parsed, 1000, 0, mode = "threshold", threshold = 150)
  
  # Verify relationships
  expect_equal(nrow(full_result[["9999"]]), 10)
  expect_equal(nrow(events_result[["9999"]]), 3)  # 100, 200, 300
  expect_equal(nrow(threshold_result[["9999"]]), 2)  # 200, 300 (above 150)
  
  # Check that filtered data maintains correct values
  expect_true(all(events_result[["9999"]]$voltage %in% c(100, 200, 300)))
  expect_true(all(threshold_result[["9999"]]$voltage %in% c(200, 300)))
})

test_that("error handling propagates correctly", {
  # Test that errors from lower-level functions bubble up appropriately
  expect_error(
    get_brw_data("", start = 0, duration = 1),
    class = "error"
  )
})

test_that("default parameters work", {
  # Test default well_id, start, and duration
  expect_error(
    get_brw_data("nonexistent.brw"),  # Uses defaults
    "File path is invalid"  # Should fail on file, not parameters
  )
})