# Test file for BRW v4 I/O functions

test_that("openBRW validates file existence", {
  # Test non-existent file
  expect_error(
    openBRW("nonexistent_file.brw"),
    "File path is invalid"
  )
})

test_that("openBRW handles valid file paths", {
  # Skip if no test data available
  skip_if_not(file.exists("tests/testthat/data/sample.brw"), 
              "No test BRW file available")
  
  # Test with valid file
  h5 <- openBRW("tests/testthat/data/sample.brw")
  expect_s3_class(h5, "H5File")
  h5$close()
})

test_that("timeCheck validates input parameters", {
  # Mock H5File object for testing
  mock_h5 <- list()
  class(mock_h5) <- "H5File"
  
  # Test negative start time
  expect_error(
    timeCheck(mock_h5, start = -1, duration = 5),
    class = "error"
  )
})

test_that("timeCheck calculates frame conversion correctly", {
  skip("Requires mock HDF5 data structure")
  
  # This test would require setting up a mock H5File
  # with proper attributes and TOC structure
})

test_that("selectChunk identifies overlapping chunks", {
  # Test data
  start_frames <- c(0, 1000, 2000, 3000)
  end_frames <- c(1000, 2000, 3000, 4000)
  
  # Test case 1: Request overlaps with single chunk
  chunks <- selectChunk(
    start_frame = 500, 
    num_frames = 200, 
    frame_starts = start_frames, 
    frame_ends = end_frames
  )
  expect_equal(chunks, 1)
  
  # Test case 2: Request overlaps with multiple chunks
  chunks <- selectChunk(
    start_frame = 900, 
    num_frames = 200, 
    frame_starts = start_frames, 
    frame_ends = end_frames
  )
  expect_equal(chunks, c(1, 2))
  
  # Test case 3: No overlapping chunks
  chunks <- selectChunk(
    start_frame = 5000, 
    num_frames = 100, 
    frame_starts = start_frames, 
    frame_ends = end_frames
  )
  expect_equal(length(chunks), 0)
})

test_that("selectChunk handles edge cases", {
  start_frames <- c(0, 1000, 2000)
  end_frames <- c(1000, 2000, 3000)
  
  # Test exact boundary
  chunks <- selectChunk(
    start_frame = 1000, 
    num_frames = 0, 
    frame_starts = start_frames, 
    frame_ends = end_frames
  )
  expect_equal(length(chunks), 0)
  
  # Test single frame request
  chunks <- selectChunk(
    start_frame = 1500, 
    num_frames = 1, 
    frame_starts = start_frames, 
    frame_ends = end_frames
  )
  expect_equal(chunks, 2)
})