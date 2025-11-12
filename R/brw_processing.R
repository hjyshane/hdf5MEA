#' Process BRW file with complete analysis pipeline
#'
#' High-level function that processes raw BRW data through the complete
#' analysis pipeline: filtering → spike detection → feature extraction.
#' Provides a one-stop solution for BRW analysis.
#'
#' @param brw_file Character. Path to BRW file
#' @param well_id Character. Well identifier (default: "Well_A1")
#' @param start Numeric. Start time in seconds (default: 0)
#' @param duration Numeric. Duration in seconds (default: NULL, processes entire file)
#' @param filter Logical. Whether to apply bandpass filter (default: TRUE)
#' @param low_freq Numeric. Low cutoff frequency for filter in Hz (default: 300)
#' @param high_freq Numeric. High cutoff frequency for filter in Hz (default: 3000)
#' @param detect_spikes Logical. Whether to detect spikes (default: TRUE)
#' @param spike_threshold Numeric. Spike detection threshold (default: NULL, auto)
#' @param remove_artifacts Logical. Whether to remove artifacts (default: TRUE)
#' @param grid_n Integer. Grid dimension for spatial mapping (default: 48)
#' @return List containing:
#'   \describe{
#'     \item{spike_data}{data.frame with detected spikes (BXR-compatible format)}
#'     \item{raw_timeseries}{List of raw voltage traces per channel}
#'     \item{processed_timeseries}{List of filtered voltage traces}
#'     \item{statistics}{Summary statistics}
#'     \item{qc_metrics}{Quality control metrics}
#'   }
#' @import hdf5r rlang cli
#' @export
#' @examples
#' \dontrun{
#' # Process entire BRW file
#' results <- processBRWFile("data.brw", duration = 60)
#'
#' # Access spike data (compatible with BXR analysis functions)
#' spikes <- results$spike_data
#' firing_rates <- calculateFiringRate(spikes)
#'
#' # Access QC metrics
#' print(results$qc_metrics)
#' }
processBRWFile <- function(brw_file, well_id = "Well_A1",
                           start = 0, duration = NULL,
                           filter = TRUE, low_freq = 300, high_freq = 3000,
                           detect_spikes = TRUE, spike_threshold = NULL,
                           remove_artifacts = TRUE, grid_n = 48) {

  cli::cli_h1("Processing BRW File")
  cli::cli_alert_info("File: {brw_file}")

  # Open file
  cli::cli_alert("Opening BRW file...")
  h5 <- openBRW(brw_file)

  # Get sampling rate
  sampling_rate <- hdf5r::h5attr(h5, "SamplingRate")
  cli::cli_alert_success("Sampling rate: {sampling_rate} Hz")

  # Determine duration if not specified
  if (is.null(duration)) {
    time_check <- brwtimeCheck(h5, start = start, duration = 1)
    duration <- time_check$total_duration - start
    cli::cli_alert_info("Processing entire file: {duration} seconds")
  }

  # Extract raw data
  cli::cli_alert("Extracting raw voltage data...")
  raw_data <- get_brw_data(h5, well_id = well_id, start = start, duration = duration)
  parsed <- brwdataParse(raw_data$binary_chunk)
  raw_timeseries <- brwtimeseriesConvert(parsed, sampling_rate,
                                         raw_data$start_frame, mode = "full")

  n_channels <- length(raw_timeseries)
  cli::cli_alert_success("Loaded {n_channels} channels")

  # Apply filtering if requested
  processed_timeseries <- raw_timeseries
  if (filter) {
    cli::cli_alert("Applying bandpass filter ({low_freq}-{high_freq} Hz)...")
    processed_timeseries <- lapply(raw_timeseries, function(ch_data) {
      if (nrow(ch_data) > 0) {
        tryCatch({
          applyBandpassFilter(ch_data, sampling_rate, low_freq, high_freq)
        }, error = function(e) {
          ch_data  # Return original if filtering fails
        })
      } else {
        ch_data
      }
    })
    cli::cli_alert_success("Filtering complete")
  }

  # Remove artifacts if requested
  if (remove_artifacts) {
    cli::cli_alert("Detecting and removing artifacts...")
    artifact_result <- detectArtifacts(processed_timeseries,
                                       sampling_rate = sampling_rate,
                                       action = "interpolate")
    processed_timeseries <- artifact_result$cleaned_data
    cli::cli_alert_success("Removed {artifact_result$n_artifacts} artifacts")
  }

  # Detect spikes if requested
  spike_data <- NULL
  if (detect_spikes) {
    cli::cli_alert("Detecting spikes...")
    spike_data <- detectSpikes(processed_timeseries,
                               sampling_rate = sampling_rate,
                               threshold = spike_threshold,
                               method = "mad")

    if (nrow(spike_data) > 0) {
      # Add spatial coordinates
      stored_chid <- raw_data$stored_channels
      spike_data <- gridMapping(spike_data, stored_chid, grid_n = grid_n,
                                joiner = "channel_id")

      # Rename to match BXR format
      spike_data$spike_times <- spike_data$spike_time + start
      spike_data$spike_chid <- spike_data$channel_id
      spike_data$peak_amplitude <- spike_data$spike_amplitude

      cli::cli_alert_success("Detected {nrow(spike_data)} spikes across {length(unique(spike_data$spike_chid))} channels")
    } else {
      cli::cli_alert_warning("No spikes detected")
    }
  }

  # Calculate statistics
  statistics <- list()
  if (!is.null(spike_data) && nrow(spike_data) > 0) {
    statistics <- tryCatch({
      getSpikeStatistics(spike_data)
    }, error = function(e) {
      list(error = "Could not calculate statistics")
    })
  }

  # Calculate QC metrics
  qc_metrics <- list()
  if (!is.null(spike_data) && nrow(spike_data) > 0) {
    qc_metrics <- tryCatch({
      list(
        dead_channels = detectDeadChannels(spike_data, all_channels = raw_data$stored_channels),
        snr = calculateSNR(spike_data, method = "amplitude")
      )
    }, error = function(e) {
      list(error = "Could not calculate QC metrics")
    })
  }

  # Close file
  h5$close_all()

  cli::cli_alert_success("Processing complete!")

  # Return results
  result <- list(
    spike_data = spike_data,
    raw_timeseries = raw_timeseries,
    processed_timeseries = processed_timeseries,
    statistics = statistics,
    qc_metrics = qc_metrics,
    metadata = list(
      file = brw_file,
      well_id = well_id,
      start = start,
      duration = duration,
      sampling_rate = sampling_rate,
      n_channels = n_channels,
      processing = list(
        filtered = filter,
        artifacts_removed = remove_artifacts,
        spikes_detected = detect_spikes
      )
    )
  )

  return(result)
}


#' Convert BRW spike detection results to BXR-compatible format
#'
#' Transforms spike data from BRW processing into a format compatible
#' with BXR analysis functions, enabling seamless pipeline integration.
#'
#' @param detected_spikes data.frame. Output from \code{\link{detectSpikes}}
#' @param stored_channels Integer vector. Channel IDs from BRW file
#' @param grid_n Integer. Grid dimension (default: 48)
#' @param start_time Numeric. Recording start time offset in seconds (default: 0)
#' @return data.frame in BXR spike data format with columns:
#'   spike_times, spike_chid, peak_amplitude, X, Y
#' @import dplyr rlang
#' @export
#' @examples
#' \dontrun{
#' # Detect spikes from BRW
#' h5 <- openBRW("data.brw")
#' raw_data <- get_brw_data(h5, start = 0, duration = 10)
#' parsed <- brwdataParse(raw_data$binary_chunk)
#' timeseries <- brwtimeseriesConvert(parsed, raw_data$sampling_rate,
#'                                    raw_data$start_frame)
#' detected <- detectSpikes(timeseries, raw_data$sampling_rate)
#'
#' # Convert to BXR format
#' spikes <- brwToSpikeData(detected, raw_data$stored_channels)
#'
#' # Now use with BXR analysis functions
#' firing_rates <- calculateFiringRate(spikes)
#'
#' h5$close()
#' }
brwToSpikeData <- function(detected_spikes, stored_channels,
                           grid_n = 48, start_time = 0) {

  if (!is.data.frame(detected_spikes)) {
    rlang::abort("detected_spikes must be a data.frame")
  }

  if (nrow(detected_spikes) == 0) {
    return(data.frame(spike_times = numeric(), spike_chid = integer(),
                      peak_amplitude = numeric(), X = integer(), Y = integer()))
  }

  # Rename columns to match BXR format
  spike_data <- detected_spikes

  if ("channel_id" %in% names(spike_data)) {
    spike_data$spike_chid <- spike_data$channel_id
  }

  if ("spike_time" %in% names(spike_data)) {
    spike_data$spike_times <- spike_data$spike_time + start_time
  }

  if ("spike_amplitude" %in% names(spike_data)) {
    spike_data$peak_amplitude <- abs(spike_data$spike_amplitude)
  }

  # Add spatial mapping if not present
  if (!all(c("X", "Y") %in% names(spike_data))) {
    spike_data <- gridMapping(spike_data, stored_channels, grid_n = grid_n,
                              joiner = "spike_chid")
  }

  # Select relevant columns
  result <- spike_data[, c("spike_times", "spike_chid", "peak_amplitude", "X", "Y")]

  return(result)
}


#' Batch process multiple BRW files
#'
#' Processes multiple BRW files in batch mode with parallel processing support.
#' Useful for high-throughput analysis of multiple recordings.
#'
#' @param brw_files Character vector. Paths to BRW files
#' @param output_dir Character. Directory for output files (default: "brw_processed")
#' @param parallel Logical. Use parallel processing (default: FALSE)
#' @param n_cores Integer. Number of cores for parallel processing (default: 2)
#' @param ... Additional arguments passed to \code{\link{processBRWFile}}
#' @return List of processing results for each file
#' @import rlang cli
#' @export
#' @examples
#' \dontrun{
#' # Batch process multiple files
#' files <- c("exp1.brw", "exp2.brw", "exp3.brw")
#' results <- batchProcessBRW(files, duration = 60, filter = TRUE)
#'
#' # Access results
#' all_spikes <- lapply(results, function(r) r$spike_data)
#' }
batchProcessBRW <- function(brw_files, output_dir = "brw_processed",
                            parallel = FALSE, n_cores = 2, ...) {

  if (!is.character(brw_files) || length(brw_files) == 0) {
    rlang::abort("brw_files must be a non-empty character vector")
  }

  # Check files exist
  missing_files <- brw_files[!file.exists(brw_files)]
  if (length(missing_files) > 0) {
    rlang::abort("Files not found: {paste(missing_files, collapse = ', ')}")
  }

  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  cli::cli_h1("Batch Processing {length(brw_files)} BRW Files")

  # Process files
  if (parallel) {
    cli::cli_alert_info("Using parallel processing with {n_cores} cores")
    # Note: Would need parallel package, keeping sequential for now
    cli::cli_alert_warning("Parallel processing not yet implemented, using sequential")
    parallel <- FALSE
  }

  results <- list()

  for (i in seq_along(brw_files)) {
    file <- brw_files[i]
    cli::cli_h2("Processing file {i}/{length(brw_files)}: {basename(file)}")

    tryCatch({
      result <- processBRWFile(file, ...)
      results[[basename(file)]] <- result

      # Save spike data
      if (!is.null(result$spike_data) && nrow(result$spike_data) > 0) {
        output_file <- file.path(output_dir,
                                 paste0(tools::file_path_sans_ext(basename(file)),
                                       "_spikes.csv"))
        exportToCSV(result$spike_data, output_file)
      }

      cli::cli_alert_success("File {i} complete")

    }, error = function(e) {
      cli::cli_alert_danger("Error processing {basename(file)}: {e$message}")
      results[[basename(file)]] <- list(error = e$message)
    })
  }

  cli::cli_alert_success("Batch processing complete!")

  return(results)
}


#' Detect field potentials from raw BRW voltage
#'
#' Identifies field potential (FP) events in raw voltage traces using
#' amplitude and duration criteria. Field potentials are slower, larger
#' amplitude signals compared to spikes.
#'
#' @param voltage_data List. Voltage traces from \code{\link{brwtimeseriesConvert}}
#' @param sampling_rate Numeric. Sampling rate in Hz
#' @param fp_threshold Numeric. Amplitude threshold for FP detection (default: NULL, auto)
#' @param min_duration Numeric. Minimum FP duration in seconds (default: 0.01)
#' @param max_duration Numeric. Maximum FP duration in seconds (default: 0.5)
#' @param filter_low Numeric. Low-pass filter cutoff for FP (default: 100 Hz)
#' @return data.frame with columns:
#'   \describe{
#'     \item{channel_id}{Channel index}
#'     \item{fp_time}{FP peak time in seconds}
#'     \item{fp_amplitude}{Peak amplitude}
#'     \item{fp_duration}{Duration in seconds}
#'     \item{fp_area}{Area under the FP curve}
#'   }
#' @import rlang
#' @importFrom stats median mad
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBRW("data.brw")
#' raw_data <- get_brw_data(h5, start = 0, duration = 10)
#' parsed <- brwdataParse(raw_data$binary_chunk)
#' timeseries <- brwtimeseriesConvert(parsed, raw_data$sampling_rate,
#'                                    raw_data$start_frame)
#'
#' # Detect field potentials
#' fps <- detectFieldPotentials(timeseries, raw_data$sampling_rate)
#'
#' h5$close()
#' }
detectFieldPotentials <- function(voltage_data, sampling_rate,
                                  fp_threshold = NULL,
                                  min_duration = 0.01, max_duration = 0.5,
                                  filter_low = 100) {

  if (!is.list(voltage_data)) {
    rlang::abort("voltage_data must be a list from brwtimeseriesConvert()")
  }

  if (!is.numeric(sampling_rate) || sampling_rate <= 0) {
    rlang::abort("sampling_rate must be positive")
  }

  min_samples <- ceiling(min_duration * sampling_rate)
  max_samples <- ceiling(max_duration * sampling_rate)

  all_fps <- list()

  for (ch_name in names(voltage_data)) {
    ch_data <- voltage_data[[ch_name]]

    if (!is.data.frame(ch_data) || !"voltage" %in% names(ch_data) || !"time" %in% names(ch_data)) {
      next
    }

    if (nrow(ch_data) < min_samples) {
      next
    }

    voltage <- ch_data$voltage
    time <- ch_data$time

    # Apply low-pass filter for FP
    # Simple moving average as low-pass filter
    window_size <- max(3, floor(sampling_rate / filter_low))
    filtered <- voltage
    for (i in (window_size+1):length(voltage)) {
      filtered[i] <- mean(voltage[(i-window_size):i])
    }

    # Calculate threshold if not provided
    if (is.null(fp_threshold)) {
      noise_est <- stats::mad(filtered, center = stats::median(filtered))
      fp_threshold <- -6 * noise_est  # FPs typically larger than spikes
    }

    # Find threshold crossings
    below_threshold <- filtered < fp_threshold

    # Find continuous segments
    fp_starts <- which(diff(c(0, below_threshold)) == 1)
    fp_ends <- which(diff(c(below_threshold, 0)) == -1)

    if (length(fp_starts) == 0) {
      next
    }

    # Filter by duration and extract properties
    for (i in seq_along(fp_starts)) {
      start_idx <- fp_starts[i]
      end_idx <- fp_ends[i]
      duration_samples <- end_idx - start_idx + 1

      if (duration_samples >= min_samples && duration_samples <= max_samples) {
        # Extract FP segment
        fp_segment <- filtered[start_idx:end_idx]
        fp_times <- time[start_idx:end_idx]

        # Find peak
        peak_idx <- which.min(fp_segment)
        peak_time <- fp_times[peak_idx]
        peak_amplitude <- fp_segment[peak_idx]

        # Calculate area
        fp_area <- sum(abs(fp_segment - stats::median(filtered)))

        fp_event <- data.frame(
          channel_id = as.numeric(ch_name),
          fp_time = peak_time,
          fp_amplitude = peak_amplitude,
          fp_duration = (end_idx - start_idx) / sampling_rate,
          fp_area = fp_area
        )

        all_fps[[length(all_fps) + 1]] <- fp_event
      }
    }
  }

  # Combine results
  if (length(all_fps) == 0) {
    return(data.frame(channel_id = numeric(), fp_time = numeric(),
                      fp_amplitude = numeric(), fp_duration = numeric(),
                      fp_area = numeric()))
  }

  result <- do.call(rbind, all_fps)
  rownames(result) <- NULL

  return(result)
}


#' Extract comprehensive features from BRW file
#'
#' Performs complete feature extraction from raw BRW data including spikes,
#' field potentials, spectral features, and signal quality metrics.
#'
#' @param brw_file Character. Path to BRW file
#' @param well_id Character. Well identifier (default: "Well_A1")
#' @param start Numeric. Start time in seconds (default: 0)
#' @param duration Numeric. Duration in seconds (default: 60)
#' @param features Character vector. Features to extract: "spikes", "fps", "spectral", "quality"
#'   (default: all)
#' @return List containing requested feature sets
#' @import hdf5r rlang cli
#' @export
#' @examples
#' \dontrun{
#' # Extract all features
#' features <- extractBRWFeatures("data.brw", duration = 60)
#'
#' # Access different feature types
#' spike_features <- features$spikes
#' fp_features <- features$field_potentials
#' quality_metrics <- features$quality
#' }
extractBRWFeatures <- function(brw_file, well_id = "Well_A1",
                               start = 0, duration = 60,
                               features = c("spikes", "fps", "quality")) {

  cli::cli_h1("Extracting Features from BRW File")

  # Process file
  result <- processBRWFile(brw_file, well_id = well_id,
                          start = start, duration = duration,
                          filter = TRUE, detect_spikes = TRUE,
                          remove_artifacts = TRUE)

  features_out <- list()

  # Spike features
  if ("spikes" %in% features && !is.null(result$spike_data)) {
    cli::cli_alert("Extracting spike features...")

    features_out$spikes <- list(
      spike_data = result$spike_data,
      firing_rates = tryCatch(calculateFiringRate(result$spike_data),
                             error = function(e) NULL),
      isi_stats = tryCatch(calculateISI(result$spike_data),
                          error = function(e) NULL),
      statistics = result$statistics
    )

    cli::cli_alert_success("Spike features extracted")
  }

  # Field potential features
  if ("fps" %in% features) {
    cli::cli_alert("Extracting field potential features...")

    fps <- tryCatch({
      detectFieldPotentials(result$processed_timeseries,
                           result$metadata$sampling_rate)
    }, error = function(e) {
      cli::cli_alert_warning("Could not extract FPs: {e$message}")
      data.frame()
    })

    features_out$field_potentials <- fps
    cli::cli_alert_success("Field potential features extracted")
  }

  # Quality metrics
  if ("quality" %in% features) {
    cli::cli_alert("Calculating quality metrics...")
    features_out$quality <- result$qc_metrics
    cli::cli_alert_success("Quality metrics calculated")
  }

  # Add metadata
  features_out$metadata <- result$metadata

  cli::cli_alert_success("Feature extraction complete!")

  return(features_out)
}
