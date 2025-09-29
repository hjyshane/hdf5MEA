#' Extract HDF5 file attributes
#' 
#' Retrieves attributes from HDF5 file root or specific attribute value.
#' Supports both single attribute extraction and complete attribute listing
#' for metadata inspection and file validation.
#' 
#' @param h5 H5File object from \code{\link{openBRW}} or \code{\link{openBXR}}
#' @param attr Character. Specific attribute name to extract. If NULL, 
#'   returns all attributes (default: NULL)
#' @return If attr specified: attribute value (type depends on attribute).
#'   If attr is NULL: named list of all root attributes
#' @return Attribute value(s):
#'   \describe{
#'     \item{Single attribute}{Atomic value. Type depends on HDF5 attribute (integer, numeric, character)}
#'     \item{All attributes}{Named list. All available attributes with their values}
#'   }
#' @import hdf5r cli rlang
#' @export
#' @examples
#' \dontrun{
#' # Open file and get specific attributes
#' h5_file <- openBRW("recording.brw")
#' 
#' # Get sampling rate for time conversion
#' sampling_rate <- getAttributes(h5_file, attr = "SamplingRate")
#' 
#' # Get voltage conversion factors
#' max_analog <- getAttributes(h5_file, attr = "MaxAnalogValue")
#' min_analog <- getAttributes(h5_file, attr = "MinAnalogValue")
#' 
#' # Get all attributes for file inspection
#' all_attrs <- getAttributes(h5_file)
#' print(names(all_attrs))
#' 
#' h5_file$close_all()
#' }
getAttributes <- function(h5, attr = NULL) {
  # Validate input parameters
  if (!inherits(h5, "H5File")) {
    rlang::abort("h5 must be an H5File object from openBRW() or openBXR()")
  }
  
  if (!is.null(attr)) {
    if (!is.character(attr) || length(attr) != 1) {
      rlang::abort("attr must be a single character string or NULL")
    }
    
    
    # Check if specific attribute exists
    if (!is.null(attr)) {
      if (!is.character(attr) || length(attr) != 1) {
        rlang::abort("attr must be a single character string or NULL")
      }

      # Extract specific attribute
      attr_value <- tryCatch({
        hdf5r::h5attr(h5, attr)
      }, error = function(e) {
        rlang::abort("Failed to read attribute '{attr}': {e$message}")
      })
      
      cli::cli_inform("Retrieved attribute '{attr}': {class(attr_value)} value")
      return(attr_value)
      
    } else {
      # Extract all attributes
      all_attributes <- tryCatch({
        hdf5r::h5attributes(h5)
      }, error = function(e) {
        rlang::abort("Failed to read file attributes: {e$message}")
      })
      
      return(all_attributes)
    }
  }
}

#' Parse experiment settings from BRW/BXR file
#' 
#' Extracts and parses the ExperimentSettings JSON object from BrainWave files.
#' Handles character encoding conversion and JSON parsing to provide structured
#' access to experimental metadata and analysis parameters.
#' 
#' @param h5 H5File object from \code{\link{openBRW}} or \code{\link{openBXR}}
#' @return List containing parsed experimental settings with nested structure.
#' @details The ExperimentSettings object contains comprehensive metadata about:
#'   \itemize{
#'     \item Hardware configuration (plate type, chip layout)
#'     \item Recording parameters (filters, gain, sampling)
#'     \item Analysis settings (thresholds, algorithms)
#'     \item Experimental conditions and protocols
#'   }
#'   
#'   Character encoding is handled automatically, converting from latin1 to UTF-8
#'   and removing non-printable characters for robust JSON parsing.
#' @import hdf5r cli rlang
#' @importFrom jsonlite fromJSON
#' @export
#' @examples
#' \dontrun{
#' h5 <- openBRW("experiment.brw")
#' 
#' # Get experiment settings
#' settings <- getExperimentSettings(h5)
#' 
#' # Access specific information
#' plate_type <- settings$PlateInfo$PlateType
#' sampling_rate <- settings$RecordingInfo$SamplingRate
#' 
#' # Check analysis parameters (if BXR file)
#' if ("AnalysisInfo" %in% names(settings)) {
#'   spike_threshold <- settings$AnalysisInfo$SpikeDetection$Threshold
#' }
#' 
#' h5$close()
#' }
getExperimentSettings <- function(h5) {
  # Validate input parameters
  if (!inherits(h5, "H5File")) {
    rlang::abort("h5 must be an H5File object from openBRW() or openBXR()")
  }
  
  # Get data
  exp <- h5[["ExperimentSettings"]][]
  
  # Convert to JSON format
  json <- iconv(exp, from = "latin1", to = "UTF-8", sub = "")
  
  # Clean non-printable characters
  json_clean <- gsub("[^[:print:]]", "", json)
  
  # Validate JSON string is not empty after cleaning
  if (nchar(json_clean) == 0) {
    rlang::abort("No valid characters remain after cleaning JSON string")
  }
  
  # Parse JSON to a list
  tryCatch({
    exp <- jsonlite::fromJSON(json_clean)
  }, error = function(e) {
    rlang::abort(paste("Failed to parse ExperimentSettings JSON:", e$message))
  })
  
  # return
  return(exp)
}

#' Create grid coordinate mapping for MEA channels
#' 
#' Maps linear channel indices to 2D grid coordinates for spatial analysis
#' and visualization. Supports different chip configurations and flexible
#' data frame integration through customizable join columns.
#' 
#' @param data data.frame. Input data containing channel information
#' @param stored_chid Integer vector. Available channel indices from MEA chip
#' @param grid_n Integer. Grid dimension for square array layout:
#'   \describe{
#'     \item{48}{2304-channel chips (48x48 grid)}
#'     \item{32}{1024-channel chips (32x32 grid)}
#'   }
#' @param joiner Character. Column name in data to join with channel mapping.
#'   Common values: "spike_chid", "fp_chid", "burst_chidx" (default: "spike_chid")
#' @return data.frame with original data plus spatial coordinates:
#'   \describe{
#'     \item{X}{Integer. Grid X coordinate (1-based, left to right)}
#'     \item{Y}{Integer. Grid Y coordinate (1-based, top to bottom)}
#'   }
#' @details Channel indexing follows 3Brain convention:
#'   \itemize{
#'     \item Linear indices are 0-based in file, converted to 1-based grid
#'     \item Grid coordinates use mathematical convention (1,1) = top-left
#'     \item X increases left to right, Y increases top to bottom
#'     \item Formula: X = (index \%\% grid_n) + 1, Y = (index \%/\% grid_n) + 1
#'   }
#' @import dplyr tibble
#' @export
#' @examples
#' \dontrun{
#' # Example spike data
#' spike_data <- data.frame(
#'   spike_times = c(1.2, 1.5, 2.1),
#'   spike_chid = c(0, 47, 2303)
#' )
#' 
#' # Channel indices from file
#' stored_channels <- 0:2303
#' 
#' # Add grid coordinates for 2304-channel chip
#' mapped_data <- gridMapping(spike_data, stored_channels, 
#'                            grid_n = 48, joiner = "spike_chid")
#' 
#' # Now data includes X, Y coordinates
#' print(mapped_data[c("spike_chid", "X", "Y")])
#' #   spike_chid  X  Y
#' #           0   1  1  (top-left)
#' #          47  48  1  (top-right)  
#' #        2303  48 48  (bottom-right)
#' }
gridMapping <- function(data, stored_chid, grid_n = 48, joiner = "spike_chid") {
  # Validate inputs
  if (!is.data.frame(data)) {
    rlang::abort("data must be a data.frame")
  }
  
  if (!joiner %in% names(data)) {
    rlang::abort(paste("Column", joiner, "not found in data"))
  }
  
  if (!is.numeric(grid_n) || length(grid_n) != 1 || grid_n <= 0) {
    rlang::abort("grid_n must be a positive integer")
  }
  
  # idx0 is zero-based index
  idx0 <- as.integer(stored_chid)
  
  # Validate channel indices
  if (any(idx0 < 0)) {
    rlang::abort("Channel indices must be non-negative")
  }
  
  if (any(idx0 >= grid_n^2)) {
    cli::cli_warn("Some channel indices exceed expected grid size ({grid_n}x{grid_n})")
  }
  
  # Create mapping table
  map <- tibble::tibble(
    OriginalChipid = idx0,
    X = (idx0 %% grid_n) + 1,  # Convert 0-based to 1-based
    Y = (idx0 %/% grid_n) + 1  # Convert 0-based to 1-based
  ) %>% 
    dplyr::distinct()
  
  # Validate mapping uniqueness
  if (nrow(map) != length(unique(idx0))) {
    cli::cli_warn("Duplicate channel indices detected in stored_chid")
  }
  
  # Join with input data
  mapped <- data %>%
    dplyr::left_join(map, by = setNames("OriginalChipid", joiner))
  
  # Check for unmapped channels
  unmapped <- sum(is.na(mapped$X))
  if (unmapped > 0) {
    cli::cli_warn("{unmapped} rows could not be mapped to grid coordinates")
  }
  
  return(mapped)
}