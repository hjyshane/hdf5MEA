# Test export functions

# Helper to create temporary directory
create_temp_dir <- function() {
  temp_dir <- tempfile()
  dir.create(temp_dir)
  temp_dir
}

# Helper to create mock spike data
create_test_spikes <- function() {
  data.frame(
    spike_times = runif(100, 0, 10),
    spike_chid = sample(1:5, 100, replace = TRUE),
    peak_amplitude = rnorm(100, -100, 20),
    X = sample(1:5, 100, replace = TRUE),
    Y = sample(1:5, 100, replace = TRUE)
  )
}


# Test exportToCSV
test_that("exportToCSV writes valid CSV file", {
  skip_if_not_installed("readr")

  spikes <- create_test_spikes()
  temp_dir <- create_temp_dir()
  output_file <- file.path(temp_dir, "test_spikes.csv")

  result <- exportToCSV(spikes, output_file)

  expect_true(file.exists(output_file))
  expect_equal(result, output_file)

  # Read back and verify
  read_back <- readr::read_csv(output_file, show_col_types = FALSE)
  expect_equal(nrow(read_back), nrow(spikes))
  expect_equal(ncol(read_back), ncol(spikes))

  unlink(temp_dir, recursive = TRUE)
})

test_that("exportToCSV validates file extension", {
  spikes <- create_test_spikes()

  expect_error(exportToCSV(spikes, "test.txt"), ".csv")
})


# Test exportToParquet
test_that("exportToParquet writes valid Parquet file", {
  skip_if_not_installed("arrow")

  spikes <- create_test_spikes()
  temp_dir <- create_temp_dir()
  output_file <- file.path(temp_dir, "test_spikes.parquet")

  result <- exportToParquet(spikes, output_file)

  expect_true(file.exists(output_file))
  expect_equal(result, output_file)

  # Read back and verify
  read_back <- arrow::read_parquet(output_file)
  expect_equal(nrow(read_back), nrow(spikes))

  unlink(temp_dir, recursive = TRUE)
})

test_that("exportToParquet validates file extension", {
  spikes <- create_test_spikes()

  expect_error(exportToParquet(spikes, "test.csv"), ".parquet")
})


# Test exportSummaryReport
test_that("exportSummaryReport generates report", {
  spikes <- create_test_spikes()
  temp_dir <- create_temp_dir()
  output_file <- file.path(temp_dir, "report.txt")

  report <- exportSummaryReport(spikes, output_file = output_file, include_plots = FALSE)

  expect_true(file.exists(output_file))
  expect_type(report, "character")
  expect_true(length(report) > 0)

  # Check report content
  report_text <- paste(report, collapse = "\n")
  expect_true(grepl("SPIKE TRAIN STATISTICS", report_text))
  expect_true(grepl("Total spikes", report_text))

  unlink(temp_dir, recursive = TRUE)
})

test_that("exportSummaryReport handles burst data", {
  spikes <- create_test_spikes()
  bursts <- data.frame(
    channel_id = c(1, 1, 2),
    start_time = c(1, 5, 3),
    end_time = c(1.5, 5.5, 3.3),
    duration = c(0.5, 0.5, 0.3)
  )

  temp_dir <- create_temp_dir()
  output_file <- file.path(temp_dir, "report_with_bursts.txt")

  report <- exportSummaryReport(spikes, bursts, output_file = output_file)

  expect_true(file.exists(output_file))
  report_text <- paste(report, collapse = "\n")
  expect_true(grepl("BURST STATISTICS", report_text))

  unlink(temp_dir, recursive = TRUE)
})


# Test batchExport
test_that("batchExport exports multiple datasets", {
  skip_if_not_installed("readr")

  spikes <- create_test_spikes()
  bursts <- data.frame(
    channel_id = c(1, 2),
    start_time = c(1, 3),
    end_time = c(1.5, 3.5),
    duration = c(0.5, 0.5)
  )

  data_list <- list(spikes = spikes, bursts = bursts)
  temp_dir <- create_temp_dir()

  files <- batchExport(data_list, output_dir = temp_dir, format = "csv")

  expect_type(files, "character")
  expect_equal(length(files), 2)
  expect_true(all(file.exists(files)))

  unlink(temp_dir, recursive = TRUE)
})

test_that("batchExport creates output directory", {
  skip_if_not_installed("readr")

  spikes <- create_test_spikes()
  data_list <- list(spikes = spikes)
  temp_dir <- file.path(tempdir(), "new_export_dir")

  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)
  }

  files <- batchExport(data_list, output_dir = temp_dir, format = "csv")

  expect_true(dir.exists(temp_dir))
  expect_true(file.exists(files["spikes"]))

  unlink(temp_dir, recursive = TRUE)
})


# Test error handling
test_that("exportToCSV validates input", {
  expect_error(exportToCSV(list(), "test.csv"), "data.frame")
  expect_error(exportToCSV(data.frame(), c("a.csv", "b.csv")), "single character")
})

test_that("exportToParquet validates input", {
  expect_error(exportToParquet(list(), "test.parquet"), "data.frame")
})

test_that("batchExport validates input", {
  expect_error(batchExport(data.frame(), "out"), "named list")
  expect_error(batchExport(list(a = 1, b = 2), "out"), "data.frame")
})
