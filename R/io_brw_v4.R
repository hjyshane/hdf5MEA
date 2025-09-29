#' Open and validate BRW file
#' 
#' Opens a BRW v4 file and validates it's a proper HDF5 file structure
#' according to 3Brain specifications. Supports Raw, EventsBasedSparseRaw,
#' and WaveletBasedEncodedRaw data formats.
#' 
#' @param file Character. Path to BRW file (.brw extension)
#' @return H5File object from hdf5r package for reading BRW data
#' @import hdf5r rlang cli
#' @export
#' @examples
#' \dontrun{
#' # Open BRW file
#' h5 <- openBRW("experiment.brw")
#' 
#' # Check sampling rate
#' sampling_rate <- h5_file$attr_open("SamplingRate")$read()
#' 
#' # Always close when done
#' h5$close()
#' }
openBRW <- function(file) {
  # Validate input parameters
  if (!is.character(file) || length(file) != 1) {
    rlang::abort("file must be a single character string")
  }
  
  if (!file.exists(file)) {
    rlang::abort("File does not exist: {file}")
  }

  # Open H5 file
  h5 <- hdf5r::H5File$new(file, mode = "r")
  
  # Then, validate file is correct file object.
  if (!inherits(h5, "H5File")) {
    rlang::abort("Input file must be a H5File object.")
  }
  
  return(h5)
}

#' Parse EventsBasedSparseRaw binary data
#' 
#' Parses EventsBasedSparseRaw binary data according to 3Brain specification.
#' Each ChData contains channel ID, data size, frame ranges, and sample values.
#' Implements the binary format: ChData header (8 bytes) + Range header (16 bytes) + samples.
#' 
#' @param binary_chunk Raw vector. Binary data from EventsBasedSparseRaw dataset
#' @return List of lists, each containing:
#'   \describe{
#'     \item{channel_id}{Integer. Linear channel index}
#'     \item{start_frame}{Integer. Range start frame (inclusive)}
#'     \item{end_frame}{Integer. Range end frame (exclusive)}
#'     \item{samples}{Integer vector. 16-bit digital sample values}
#'   }
#' @details The function parses the binary format specified in 3Brain documentation:
#'   \itemize{
#'     \item ChData header: 4 bytes channel ID + 4 bytes data size
#'     \item Range header: 8 bytes start frame + 8 bytes end frame  
#'     \item Sample data: 2 bytes per sample (16-bit integers)
#'   }
#' @import rlang cli
#' @export
#' @examples
#' \dontrun{
#' # Get raw binary data
#' raw_data <- get_brw_data("data.brw", start = 10, duration = 2)
#' 
#' # Parse binary data
#' parsed <- brwdataParse(raw_data$binary_chunk)
#' 
#' # Access first channel's data
#' ch1_data <- parsed[[1]]
#' print(ch1_data$channel_id)
#' print(ch1_data$samples[1:10])
#' }
brwdataParse <- function(binary_chunk) {
  # Create empty list for parsed data
  parsed_data <- list()
  
  # Handle empty input
  if (length(binary_chunk) == 0) {
    return(parsed_data)
  }
  
  # Transform into raw vector
  raw <- as.raw(binary_chunk)
  
  # Initialize position counter
  i <- 1
  
  # Loop for parsing ChData structures
  while (i <= length(raw)) {
    # Check if enough data for ChData header (8 bytes)
    if (i + 7 > length(raw)) {
      rlang::abort("Insufficient data for ChData header")
    }
    
    # Channel ID parse (4 bytes)
    ch_id <- readBin(raw[i:(i+3)], "integer", size = 4, endian = "little")
    
    # Data size parse (4 bytes)
    data_size <- readBin(raw[(i+4):(i+7)], "integer", size = 4, endian = "little")
    
    # Check if enough data for Range header (16 bytes)
    if (i + 23 > length(raw)) {
      rlang::abort("Insufficient data for Range header")
    }
    
    # Range header: startFrame(8 bytes) + endFrame(8 bytes)
    start_frame <- readBin(raw[(i+8):(i+15)], "integer", size = 8, endian = "little")
    end_frame <- readBin(raw[(i+16):(i+23)], "integer", size = 8, endian = "little")
    
    # Sample data parsing
    sample_count <- end_frame - start_frame
    sample_start <- i + 24
    
    # Check if enough data for samples
    if (sample_start + sample_count*2 - 1 > length(raw)) {
      rlang::abort("Insufficient data for samples")
    }
    
    # Each sample is 2 bytes
    # Read as 16-bit integer
    samples <- readBin(raw[sample_start:(sample_start + sample_count*2 - 1)], 
                       "integer", size = 2, n = sample_count, endian = "little")
    
    # Save parsed data
    parsed_data[[length(parsed_data) + 1]] <- list(
      channel_id = ch_id,
      start_frame = start_frame,
      end_frame = end_frame,
      samples = samples
    )
    
    # Move to next ChData (skip by data_size)
    i <- i + 8 + data_size  # ChData header(8) + data
  }
  
  return(parsed_data)
}

#' Select data chunks for time range
#' 
#' Finds data chunks that overlap with the requested time frame using TOC information.
#' Uses efficient binary search for chunk selection in large datasets.
#' 
#' @param start_frame Integer. Starting frame index (0-based)
#' @param num_frames Integer. Number of frames requested
#' @param frame_starts Numeric vector. Start frames of all chunks from TOC
#' @param frame_ends Numeric vector. End frames of all chunks from TOC
#' @return Integer vector of chunk indices that overlap with requested range
#' @details Chunks are selected if they have any overlap with the requested range:
#'   \code{(chunk_start < request_end) & (chunk_end > request_start)}
#' @import rlang cli
#' @export
#' @examples
#' \dontrun{
#' # Example TOC data
#' starts <- c(0, 1000, 2000, 3000)
#' ends <- c(1000, 2000, 3000, 4000)
#' 
#' # Find chunks overlapping frames 1500-2500
#' chunks <- brwselectChunk(1500, 1000, starts, ends)
#' print(chunks)  # Should return c(2, 3)
#' }
brwselectChunk <- function(start_frame, num_frames, frame_starts, frame_ends) {
  # Validate input parameters
  if (!is.numeric(start_frame) || length(start_frame) != 1 || start_frame < 0) {
    rlang::abort("start_frame must be a non-negative integer")
  }
  
  if (!is.numeric(num_frames) || length(num_frames) != 1 || num_frames <= 0) {
    rlang::abort("num_frames must be a positive integer")
  }
  
  if (!is.numeric(frame_starts) || !is.numeric(frame_ends)) {
    rlang::abort("frame_starts and frame_ends must be numeric vectors")
  }
  
  if (length(frame_starts) != length(frame_ends)) {
    rlang::abort("frame_starts and frame_ends must have the same length")
  }
  
  if (length(frame_starts) == 0) {
    rlang::abort("No chunks available (empty frame arrays)")
  }
  
  end_frame <- start_frame + num_frames
  
  # Find overlapping chunks
  overlap_condition <- (frame_starts < end_frame) & (frame_ends > start_frame)
  target_chunks <- which(overlap_condition)
  
  return(target_chunks)
}

#' Validate time range and convert to frames
#' 
#' Validates requested time range against available data in BRW file
#' and converts time values to frame indices using sampling rate.
#' Extracts TOC information for chunk management.
#' 
#' @param h5 H5File object from \code{\link{openBRW}}
#' @param start Numeric. Start time in seconds
#' @param duration Numeric. Duration in seconds  
#' @return List containing:
#'   \describe{
#'     \item{sampling_rate}{Numeric. Sampling rate in Hz from file}
#'     \item{start_frame}{Integer. Requested start frame (0-based)}
#'     \item{num_frames}{Integer. Number of frames for requested duration}
#'     \item{end_frame}{Integer. Requested end frame (exclusive, 0-based)}
#'     \item{start_frames}{Integer vector. Start frames of all data chunks}
#'     \item{end_frames}{Integer vector. End frames of all data chunks}
#'   }
#' @import hdf5r rlang cli
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBRW("data.brw")
#' 
#' # Check time range and get frame info
#' time_info <- brwtimeCheck(h5, start = 10, duration = 5)
#' 
#' print(paste("Sampling rate:", time_info$sr))
#' print(paste("Start frame:", time_info$start_frame))
#' 
#' # Use for subsequent data extraction
#' chunks <- brwselectChunk(time_info$start_frame, time_info$num_frames,
#'                         time_info$start_frames, time_info$end_frames)
#' 
#' h5$close()
#' }
brwtimeCheck <- function(h5, start, duration) {
  # Validate input parameters
  if (!inherits(h5, "H5File")) {
    rlang::abort("h5 must be an H5File object from openBRW()")
  }
  
  if (!is.numeric(start) || length(start) != 1 || start < 0) {
    rlang::abort("start must be a non-negative numeric value")
  }
  
  if (!is.numeric(duration) || length(duration) != 1 || duration <= 0) {
    rlang::abort("duration must be a positive numeric value")
  }
  
  # Check for required datasets
  if (!h5$exists("TOC")) {
    rlang::abort("TOC dataset not found in BRW file")
  }
  

  # Get Sampling rate
  sampling_rate <- getAttributes(h5, attr = "SamplingRate")
  
  # Extract whole chunk information
  root_toc <- h5[["TOC"]][, ] # read chunk data from h5 file
  start_frames <- as.numeric(root_toc[1, ]) # all start frames of each chunk
  end_frames <- as.numeric(root_toc[2, ]) # all end frames of each chunk
  
  # Transform time range to frame
  start_frame <- as.integer(start * sampling_rate) # requested starting point
  num_frames <- as.integer(duration * sampling_rate) # number of frames for requested duration
  end_frame <- start_frame + num_frames # request ending point
  
  # check if time range is valid.
  if (start_frame < min(start_frames) || end_frame > max(end_frames)) {
    rlang::abort("Requested time range exceeds available data")
  }  
  
  return(list(
    sampling_rate = sampling_rate,
    start_frame = start_frame,
    num_frames = num_frames,
    end_frame = end_frame,
    start_frames = start_frames,
    end_frames = end_frames
  ))
}

#' Convert parsed data to time series
#' 
#' Converts parsed ChData structures into time series data frames with 
#' multiple filtering options for sparse data handling. Supports different
#' modes for memory efficiency and analysis needs.
#' 
#' @param parsed_data List. Output from \code{\link{brwdataParse}}
#' @param sampling_rate Numeric. Sampling frequency in Hz from file attributes
#' @param start_frame Integer. Reference start frame for time calculation
#' @param mode Character. Data filtering mode:
#'   \describe{
#'     \item{"full"}{Keep all samples (default)}
#'     \item{"events_only"}{Keep only non-zero samples}
#'     \item{"threshold"}{Keep samples above threshold}
#'   }
#' @param threshold Numeric. Threshold value for "threshold" mode (default: 50)
#' @return Named list of data frames, one per channel. Each data frame contains:
#'   \describe{
#'     \item{time}{Numeric. Time in seconds relative to data start}
#'     \item{voltage}{Integer. Digital sample values}
#'   }
#' @details Channel names in returned list correspond to channel IDs from BRW file.
#'   Time axis is calculated relative to the start of the extracted data.
#' @import hdf5r rlang cli dplyr
#' @export
#' @examples
#' \dontrun{
#' # Get and parse data
#' raw_data <- get_brw_data("data.brw", start = 10, duration = 2)
#' parsed <- brwdataParse(raw_data$binary_chunk)
#' 
#' # Convert to time series (all data)
#' ts_full <- brwtimeseriesConvert(parsed, 
#'                                 raw_data$start_frame, mode = "full")
#' 
#' # Convert to time series (events only)
#' ts_events <- brwtimeseriesConvert(parsed,
#'                                   raw_data$start_frame, mode = "events_only")
#' 
#' # Access channel data
#' ch1_data <- ts_full[["1"]]  # Channel ID 1
#' plot(ch1_data$time, ch1_data$voltage, type = "l")
#' }
brwtimeseriesConvert <- function(parsed_data,
                                 sampling_rate,
                                 start_frame,
                                 mode = c("full", "events_only", "threshold"),
                                 threshold = 50) {
  
  # Validate input parameters
  if (!is.list(parsed_data) || length(parsed_data) == 0) {
    rlang::abort("parsed_data must be a non-empty list from brwdataParse()")
  }
  
  if (!is.numeric(sampling_rate) || length(sampling_rate) != 1 || sampling_rate <= 0) {
    rlang::abort("sampling_rate must be a positive numeric value")
  }
  
  if (!is.numeric(start_frame) || length(start_frame) != 1 || start_frame < 0) {
    rlang::abort("start_frame must be a non-negative integer")
  }
  
  start_frame <- as.integer(start_frame)
  
  # Validate and match mode argument
  mode <- match.arg(mode)
  
  if (!is.numeric(threshold) || length(threshold) != 1) {
    rlang::abort("threshold must be a numeric value")
  }
  
  # Initialize list for each channel
  channel_data <- list()
  
  # Loop for conversion
  for (ch_data in parsed_data) {
    
    # Get channel id
    ch_id <- ch_data$channel_id
    
    # Create time axis (relative to start of data)
    time_seconds <- seq(0, length(ch_data$samples) - 1) / sampling_rate
    
    # Create data frame
    df <- data.frame(time = time_seconds, voltage = ch_data$samples)
    
    # Apply filtering based on mode
    if (mode == "events_only") {
      df <- df[df$voltage != 0, ]
    } else if (mode == "threshold") {
      df <- df[abs(df$voltage) > threshold, ]
    } 
    # If mode == "full", keep all data
    
    # Store or append to channel data
    if (is.null(channel_data[[as.character(ch_id)]])) {
      channel_data[[as.character(ch_id)]] <- df
    } else {
      channel_data[[as.character(ch_id)]] <- rbind(channel_data[[as.character(ch_id)]], df)
    }
  }
  
  return(channel_data)
}

#' Extract MEA signal data by time range
#' 
#' Main high-level function for extracting MEA signal data from BRW v4 files.
#' Handles complete workflow: file opening, time validation, chunk selection,
#' and binary data extraction from EventsBasedSparseRaw format.
#' 
#' @param h5 H5File object. Open BXR file handle from \code{\link{openBXR}}
#' @param well_id Character. Well identifier (default: "Well_A1")
#' @param start Numeric. Start time in seconds (default: 0)
#' @param duration Numeric. Duration in seconds (default: 5)
#' @return List containing complete extraction results:
#'   \describe{
#'     \item{binary_chunk}{Raw vector. Binary data from EventsBasedSparseRaw}
#'     \item{stored_channels}{Integer vector. Channel IDs available in well}
#'     \item{data_range}{Integer vector. Start and end positions in raw data}
#'     \item{sampling_rate}{Numeric. Sampling frequency in Hz}
#'     \item{start_frame}{Integer. Requested start frame}
#'     \item{start_frames}{Numeric vector. All chunk start frames}
#'     \item{num_frames}{Integer. Number of requested frames}
#'     \item{end_frame}{Integer. Requested end frame}
#'   }
#' @details This is the main user-facing function for BRW data extraction.
#'   It automatically handles resource management (file closing) and provides
#'   comprehensive output for downstream analysis.
#' @import hdf5r rlang cli
#' @export
#' @seealso \code{\link{brwdataParse}}, \code{\link{brwtimeseriesConvert}}
#' @examples
#' \dontrun{
#' # Extract 2 seconds starting from 10 seconds
#' data <- get_brw_data("experiment.brw", well_id = "Well_A1", 
#'                      start = 10, duration = 2)
#' 
#' # Parse and convert to time series
#' parsed <- brwdataParse(data$binary_chunk)
#' timeseries <- brwtimeseriesConvert(parsed, data$sampling_rate, 
#'                                    data$start_frame)
#' 
#' # Access metadata
#' print(paste("Sampling rate:", data$sampling_rate, "Hz"))
#' print(paste("Channels available:", length(data$stored_channels)))
#' }
get_brw_data <- function(h5,
                         well_id = "Well_A1", 
                         start = 0, 
                         duration = 5) {
  
  # Validate input parameters
  if (!inherits(h5, "H5File")) {
    rlang::abort("h5 must be an H5File object from openBXR()")
  }
  
  if (!is.character(well_id) || length(well_id) != 1) {
    rlang::abort("well_id must be a single character string")
  }
  
  if (!grepl("^Well_[A-Z][0-9]+$", well_id)) {
    cli::cli_warn("well_id '{well_id}' does not follow 3Brain naming convention")
  }
  
  if (!is.numeric(start) || length(start) != 1 || start < 0) {
    rlang::abort("start must be a non-negative numeric value")
  }
  
  if (!is.numeric(duration) || length(duration) != 1 || duration <= 0) {
    rlang::abort("duration must be a positive numeric value")
  }
  
  # Time setup and validation
  time_info <- tryCatch({
    brwtimeCheck(h5, start, duration)
  }, error = function(e) {
    rlang::abort("Time validation failed: {e$message}")
  })
  
  # Select data chunks that overlap with requested time range
  target_chunks <- tryCatch({
    brwselectChunk(time_info$start_frame, 
                   time_info$num_frames, 
                   time_info$start_frames, 
                   time_info$end_frames)
  }, error = function(e) {
    rlang::abort("Chunk selection failed: {e$message}")
  })
  
  # Validate chunk selection
  if (is.null(target_chunks) || length(target_chunks) == 0) {
    rlang::abort("No data chunks found for requested time range")
  }
  
  # Extract binary data from selected chunks
  raw_data <- tryCatch({
    eventBased(h5, well_id, target_chunks)
  }, error = function(e) {
    rlang::abort("Data extraction failed: {e$message}")
  })
  
  # Calculate end frame for completeness
  end_frame <- time_info$start_frame + time_info$num_frames
  
  # Return comprehensive data structure
  result <- list(
    binary_chunk = raw_data$binary_chunk,
    stored_channels = raw_data$stored_channels,
    data_range = raw_data$data_range,
    sampling_rate = time_info$sr,
    start_frame = time_info$start_frame,
    num_frames = time_info$num_frames,
    end_frame = end_frame,
    start_frames = time_info$start_frames,
    end_frames = time_info$end_frames
  )
  
  return(result)
}

#' Extract EventsBasedSparseRaw binary data
#' 
#' Low-level function to extract binary data from EventsBasedSparseRaw dataset
#' for specified chunks. Handles TOC-based data positioning and provides
#' warnings about sparse data characteristics.
#' 
#' @param h5 H5File object from \code{\link{openBRW}}
#' @param well_id Character. Well identifier (e.g., "Well_A1")
#' @param target_chunks Integer vector. Chunk indices to extract (1-based indexing)
#'   from Table of Contents
#' @return List containing extracted binary data and metadata:
#'   \describe{
#'     \item{binary_chunk}{Raw bytes. EventsBasedSparseRaw binary data for parsing}
#'     \item{stored_channels}{Integer vector. Available channel indices (0-based)}
#'     \item{data_range}{Integer vector. Start and end positions in dataset}
#'   }
#' @details EventsBasedSparseRaw format stores only periods around detected events,
#'   not continuous data. The actual time coverage may be less than requested.
#'   Use TOC (Table of Contents) information to understand data organization.
#' @import hdf5r rlang cli
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBRW("data.brw")
#' 
#' # Get chunk information first
#' time_info <- brwtimeCheck(h5, start = 10, duration = 2)
#' chunks <- brwselectChunk(time_info$start_frame, time_info$num_frames,
#'                          time_info$start_frames, time_info$end_frames)
#' 
#' # Extract binary data
#' raw_data <- eventBased(h5, "Well_A1", chunks)
#' 
#' # Parse binary data
#' parsed <- brwdataParse(raw_data$binary_chunk)
#' 
#' h5$close()
#' }
eventBased <- function(h5, well_id, target_chunks) {
  # Validate input parameters
  if (!inherits(h5, "H5File")) {
    rlang::abort("h5 must be an H5File object from openBRW()")
  }
  
  if (!is.character(well_id) || length(well_id) != 1) {
    rlang::abort("well_id must be a single character string")
  }
  
  if (!grepl("^Well_[A-Z][0-9]+$", well_id)) {
    cli::cli_warn("well_id '{well_id}' does not follow 3Brain naming convention")
  }
  
  if (!is.numeric(target_chunks) || any(target_chunks < 1)) {
    rlang::abort("target_chunks must be positive integers (1-based indexing)")
  }
  
  # Read EventsBasedSparseRaw data and TOC
  raw_data <- h5[[paste0(well_id, "/EventsBasedSparseRaw")]][] # saved as chunk size
  events_toc <- h5[[paste0(well_id, "/EventsBasedSparseRawTOC")]][]
  
  # Read channel information
  stored_channels <- h5[[paste0(well_id, "/StoredChIdxs")]][]
  n_channels <- length(stored_channels)
  
  # Calculate data range to extract
  if (length(target_chunks) > 0) {
    start_pos <- events_toc[target_chunks[1]]
    end_pos <- if (max(target_chunks) == length(events_toc)) {
      length(raw_data)
    } else {
      events_toc[max(target_chunks) + 1] # to get the end of the last target_chunks
    }
  }
  
  # Extract binary data
  if (length(target_chunks) > 0) {
    binary_chunk <- raw_data[start_pos:end_pos]
    
    # Warning about EventsBasedSparseRaw
    cli::cli_warn("EventsBasedSparseRaw contains only event periods. Actual time range may differ from request.")
    
  } else {
    rlang::abort("No chunks found for the requested time range")
  }
  
  return(list(
    binary_chunk = binary_chunk,
    stored_channels = stored_channels,
    data_range = c(start_pos, end_pos)
  ))
}