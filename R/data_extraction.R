#' @title Data Extraction and Parsing Functions
#' @description Functions for extracting and parsing EventsBasedSparseRaw data
#' @author hdf5MEA package

#' Extract EventsBasedSparseRaw binary data
#' 
#' Extracts binary data from EventsBasedSparseRaw dataset for specified chunks.
#' 
#' @param h5 H5File object
#' @param well_id Character. Well identifier (e.g., "Well_A1")
#' @param target_chunks Integer vector. Chunk indices to extract
#' @return List containing binary_chunk, stored_channels, and data_range
#' @export
eventBased <- function(h5, well_id, target_chunks) {
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
    warning("EventsBasedSparseRaw contains only event periods. Actual time range may differ from request.")
    
  } else {
    stop("No chunks found for the requested time range")
  }
  
  return(list(
    binary_chunk = binary_chunk,
    stored_channels = stored_channels,
    data_range = c(start_pos, end_pos)
  ))
}

#' Parse binary data into ChData structures
#' 
#' Parses EventsBasedSparseRaw binary data according to 3Brain specification.
#' Each ChData contains channel ID, frame range, and sample values.
#' 
#' @param binary_chunk Raw vector. Binary data from EventsBasedSparseRaw
#' @return List of lists, each containing channel_id, start_frame, end_frame, and samples
#' @export
dataParse <- function(binary_chunk) {
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
      stop("Insufficient data for ChData header")
    }
    
    # Channel ID parse (4 bytes)
    ch_id <- readBin(raw[i:(i+3)], "integer", size = 4, endian = "little")
    
    # Data size parse (4 bytes)
    data_size <- readBin(raw[(i+4):(i+7)], "integer", size = 4, endian = "little")
    
    # Check if enough data for Range header (16 bytes)
    if (i + 23 > length(raw)) {
      stop("Insufficient data for Range header")
    }
    
    # Range header: startFrame(8 bytes) + endFrame(8 bytes)
    start_frame <- readBin(raw[(i+8):(i+15)], "integer", size = 8, endian = "little")
    end_frame <- readBin(raw[(i+16):(i+23)], "integer", size = 8, endian = "little")
    
    # Sample data parsing
    sample_count <- end_frame - start_frame
    sample_start <- i + 24
    
    # Check if enough data for samples
    if (sample_start + sample_count*2 - 1 > length(raw)) {
      stop("Insufficient data for samples")
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

#' Convert parsed data to time series format
#' 
#' Converts parsed ChData structures into time series data frames with 
#' options for sparse data handling.
#' 
#' @param parsed_data List. Output from dataParse()
#' @param sampling_rate Numeric. Sampling rate in Hz
#' @param start_frame Integer. Reference start frame for time calculation
#' @param mode Character. Data filtering mode: "full", "events_only", or "threshold"
#' @param threshold Numeric. Threshold value for "threshold" mode
#' @return List of data frames, one per channel, with time and voltage columns
#' @export
timeseriesConvert <- function(parsed_data, 
                              sampling_rate, 
                              start_frame,
                              mode = c("full", "events_only", "threshold"),
                              threshold = 50) {
  
  sparse_mode <- match.arg(mode)
  
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
    if (sparse_mode == "events_only") {
      df <- df[df$voltage != 0, ]
    } else if (sparse_mode == "threshold") {
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