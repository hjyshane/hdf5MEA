#' @title Main User API Functions
#' @description High-level functions for extracting MEA data by time ranges
#' @author hdf5MEA package

#' Extract signal data by time ranges
#' 
#' Main function for extracting MEA signal data from BRW v4 files.
#' Handles file opening, time validation, chunk selection, and binary data extraction.
#' 
#' @param file Character. Path to BRW file
#' @param well_id Character. Well identifier (default: "Well_A1")
#' @param start Numeric. Start time in seconds (default: 0)
#' @param duration Numeric. Duration in seconds (default: 5)
#' @return List containing:
#'   \itemize{
#'     \item binary_chunk: Raw binary data
#'     \item stored_channels: Vector of channel IDs
#'     \item data_range: Start and end positions in raw data
#'     \item sampling_rate: Sampling frequency in Hz
#'     \item start_frame: Requested start frame
#'     \item start_frames: All chunk start frames
#'     \item num_frames: Number of requested frames
#'     \item end_frames: Requested end frame
#'   }
#' @export
#' @examples
#' \dontrun{
#' # Extract 2 seconds starting from 10 seconds
#' data <- get_brw_data("experiment.brw", well_id = "Well_A1", start = 10, duration = 2)
#' 
#' # Parse and convert to time series
#' parsed <- dataParse(data$binary_chunk)
#' timeseries <- timeseriesConvert(parsed, data$sampling_rate, data$start_frame)
#' }
get_brw_data <- function(file,
                     well_id = "Well_A1", 
                     start = 0, 
                     duration = 5) {
  
  # Validate input parameters
  if (start < 0 | duration <= 0) {
    stop("Time parameter is wrong.")
  } 
  if (duration > 5) {
    warning("If duration is more than 5 seconds, You may run out of memory")
  }
  
  # File operations
  h5 <- openBRW(file)
  
  # Time setup and validation
  time_info <- brwtimeCheck(h5, start, duration)
  
  # Get requested chunks
  target_chunks <- brwselectChunk(time_info$start_frame, 
                               time_info$num_frames, 
                               time_info$start_frames, 
                               time_info$end_frames)
  
  # Data extraction
  raw_data <- eventBased(h5, well_id, target_chunks)
  
  # Clean up resources
  h5$close()
  
  # Return comprehensive data structure
  return(list(
    binary_chunk = raw_data$binary_chunk,
    stored_channels = raw_data$stored_channels,
    data_range = raw_data$data_range,
    sampling_rate = time_info$sr,
    start_frame = time_info$start_frame, 
    start_frames = time_info$start_frames,
    num_frames = time_info$num_frames,
    end_frame = time_info$start_frame + time_info$num_frames
  ))
}