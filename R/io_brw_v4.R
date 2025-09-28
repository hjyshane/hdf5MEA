#' @title BRW v4 File I/O Functions
#' @description Functions for opening BRW v4 files and processing time/chunk information
#' @author hdf5MEA package

#' Open and validate BRW file
#' 
#' Opens a BRW file and validates it's a proper HDF5 file structure.
#' 
#' @param file Character. Path to BRW file
#' @return H5File object from hdf5r package
#' @import hdf5r
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBRW("data.brw")
#' h5$close()
#' }
openBRW <- function(file) {
  # Check file path is correct.
  if (!file.exists(file)) {
    stop("File path is invalid. Check path or file exist.")
  } 
  
  # Open H5 file
  h5 <- hdf5r::H5File$new(file, mode = "r")
  
  # Then, validate file is correct file object.
  if (!inherits(h5, "H5File")) {
    stop("Input file must be a H5File object.")
  }
  
  return(h5)
}

#' Check time range and convert to frames
#' 
#' Validates requested time range against available data and converts
#' time values to frame indices.
#' 
#' @param h5 H5File object
#' @param start Numeric. Start time in seconds
#' @param duration Numeric. Duration in seconds
#' @return List containing sampling rate, frame information, and TOC data
#' @import hdf5r
#' @export
timeCheck <- function(h5, start, duration) {
  # Get Sampling rate and transform time-frame
  sr <- hdf5r::h5attr(h5, "SamplingRate")
  
  # Extract whole chunk information
  root_toc <- h5[["TOC"]][, ] # read chunk data from h5 file
  start_frames <- as.numeric(root_toc[1, ]) # all start frames of each chunk
  end_frames <- as.numeric(root_toc[2, ]) # all end frames of each chunk
  
  # Transform time range to frame
  start_frame <- as.integer(start * sr) # requested starting point
  num_frames <- as.integer(duration * sr) # number of frames for requested duration
  end_frame <- start_frame + num_frames # request ending point
  
  # check if time range is valid.
  if (start_frame < min(start_frames) || end_frame > max(end_frames)) {
    stop("Requested time range exceeds available data")
  }  
  
  return(list(
    sr = sr,
    start_frame = start_frame,
    num_frames = num_frames,
    start_frames = start_frames,
    end_frames = end_frames
  ))
}

#' Select chunks for time frame of interest
#' 
#' Finds chunks that overlap with the requested time frame using TOC data.
#' 
#' @param start_frame Integer. Starting frame index
#' @param num_frames Integer. Number of frames requested
#' @param frame_starts Numeric vector. Start frames of all chunks
#' @param frame_ends Numeric vector. End frames of all chunks
#' @return Integer vector of chunk indices that overlap with requested range
#' @export
selectChunk <- function(start_frame, num_frames, frame_starts, frame_ends) {
  end_frame <- start_frame + num_frames
  
  # Find overlapping chunks
  target_chunks <- which(
    (frame_starts < end_frame) & (frame_ends > start_frame)
  )
  
  return(target_chunks)
}