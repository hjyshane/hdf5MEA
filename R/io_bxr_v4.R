#' Open and validate BXR file
#' 
#' Opens a BXR v3/v4 file and validates it's a proper HDF5 file structure
#' according to 3Brain specifications.
#' 
#' @param file Character. Path to BXR file (.bxr extension)
#' @return H5File object from hdf5r package for reading BXR data
#' @import hdf5r
#' @export
#' @seealso \code{\link{bxrSpikedata}}, \code{\link{bxrFpdata}}
#' @examples
#' \dontrun{
#' # Open BXR file
#' h5 <- openBXR("experiment.bxr")
#' 
#' # Always close when done
#' h5$close()
#' }
openBXR <- function(file) {
  # Check file path is correct.
  if (!file.exists(file)) {
    rlang::abort("File path is invalid. Check path or file exist.")
  } 
  
  # Open H5 file
  h5 <- hdf5r::H5File$new(file, mode = "r")
  
  # Then, validate file is correct file object.
  if (!inherits(h5, "H5File")) {
    rlang::abort("Input file must be a H5File object.")
  }
  
  return(h5)
}

#' Extract spike event data from BXR file
#' 
#' Reads spike times, channel indices, units, and amplitude information
#' from a BXR file with optional time filtering and grid mapping.
#' 
#' @param h5 H5File object from \code{\link{openBXR}}
#' @param well_id Character. Well identifier (default: "Well_A1")
#' @param grid_n Integer. Grid size for channel mapping (default: 48)
#' @param start Numeric. Start time in seconds for filtering (optional)
#' @param duration Numeric. Duration in seconds for filtering (optional)
#' @return data.frame with columns:
#'   \describe{
#'     \item{spike_times}{Spike timestamps in seconds}
#'     \item{spike_chid}{Channel linear indices}
#'     \item{spike_units}{Unit classification from spike sorting}
#'     \item{peak_amplitude}{Maximum absolute amplitude}
#'     \item{min_amplitude}{Minimum amplitude}
#'     \item{waveform_index}{Index for waveform lookup}
#'     \item{grid_x, grid_y}{Grid coordinates (if mapping applied)}
#'   }
#' @import hdf5r
#' @export
#' @seealso \code{\link{getWaveform}}, \code{\link{bxrSpikeBursts}}
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' 
#' # Get all spikes
#' spikes <- bxrSpikedata(h5)
#' 
#' # Get spikes in specific time window
#' spikes_subset <- bxrSpikedata(h5, start = 10, duration = 30)
#' 
#' h5$close()
#' }
bxrSpikedata <- function(h5, 
                         well_id = "Well_A1",
                         grid_n = 48,
                         start = NULL,
                         duration = NULL) {
  # First get Sampling rate
  sr <- getAttributes(h5, attr = "SamplingRate")
  
  # data frame with SpikeTimes, SpikeChIdxs, SpikeUnits, StoredChIdxs
  spike_data <- data.frame(
    spike_times = h5[[paste0(well_id, "/SpikeTimes")]][] / sr,
    spike_chid = h5[[paste0(well_id, "/SpikeChIdxs")]][],
    spike_units = h5[[paste0(well_id, "/SpikeUnits")]][]
  )
  
  # Extract amplitude data
  wavelength <- hdf5r::h5attr(h5[[paste0(well_id, "/SpikeForms")]], "Wavelength")
  spike_forms <- h5[[paste0(well_id, "/SpikeForms")]][]
  n_spikes <- nrow(spike_data)
  
  # create empty numeric vector to store values same length as spike_data.
  peak_amplitudes <- numeric(n_spikes)
  min_amplitudes <- numeric(n_spikes)
  
  for (i in 1:n_spikes) {
    start_idx <- (i - 1) * wavelength + 1
    end_idx <- i * wavelength
    waveform <- spike_forms[start_idx:end_idx]
    
    peak_amplitudes[i] <- max(abs(waveform))
    min_amplitudes[i] <- min(waveform)
  }
  
  # Add amplitude columns
  spike_data$peak_amplitude <- peak_amplitudes
  spike_data$min_amplitude <- min_amplitudes
  
  # Add waveform index for later lookup
  spike_data$waveform_index <- 1:n_spikes
  
  # Time filtering if specified
  if (!is.null(start) && !is.null(duration)) {
    whole <- range(min(spike_data$spike_times), max(spike_data$spike_times))
    end_time <- start + duration
    
    if (start >= whole[1] & end_time <= whole[2]) {
      spike_data <- spike_data[spike_data$spike_times >= start & 
                                 spike_data$spike_times <= end_time, ]
    } else {
      cli::cli_warn("Specified time range does not match with spike_data time range.")
    } 
  } else {
    cli::cli_inform("Time is not specified, return whole spike data.")
  }
  
  # Chip ID
  stored_chid <- h5[[paste0(well_id, "/StoredChIdxs")]][]
  
  # Create mapping for grid
  spike_data <- gridMapping(spike_data, stored_chid, grid_n = grid_n)
  
  return(spike_data)
}

#' Extract spike burst data from BXR file
#' 
#' Reads spike burst events including start/end times and durations
#' with optional time filtering and grid mapping.
#' 
#' @param h5 H5File object from \code{\link{openBXR}}
#' @param well_id Character. Well identifier (default: "Well_A1")
#' @param grid_n Integer. Grid size for channel mapping (default: 48)
#' @param start Numeric. Start time in seconds for filtering (optional)
#' @param duration Numeric. Duration in seconds for filtering (optional)
#' @return data.frame with columns:
#'   \describe{
#'     \item{burst_id}{Sequential burst identifier}
#'     \item{burst_chidx}{Channel linear index}
#'     \item{start_time}{Burst start time in seconds}
#'     \item{end_time}{Burst end time in seconds}
#'     \item{duration}{Burst duration in seconds}
#'     \item{grid_x, grid_y}{Grid coordinates (if mapping applied)}
#'   }
#' @import hdf5r
#' @export
#' @seealso \code{\link{bxrSpikedata}}, \code{\link{bxrNetworkBursts}}
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' 
#' # Get all spike bursts
#' bursts <- bxrSpikeBursts(h5)
#' 
#' # Get bursts in specific time window
#' bursts_subset <- bxrSpikeBursts(h5, start = 10, duration = 30)
#' 
#' h5$close()
#' }
bxrSpikeBursts <- function(h5,
                           well_id = "Well_A1",
                           grid_n = 48,
                           start = NULL,
                           duration = NULL) {
  # Extract spike burst data
  burst_times <- h5[[paste0(well_id, "/SpikeBurstTimes")]][,]  
  burst_chidx <- h5[[paste0(well_id, "/SpikeBurstChIdxs")]][]
  sampling_rate <- getAttributes(h5, attr = "SamplingRate")
  
  if (length(burst_chidx) == 0) {
    cli::cli_inform("There are no bursts.")
    return(data.frame()) 
  }
  
  # Create burst data frame
  burst_data <- data.frame(
    burst_id = 1:length(burst_chidx),
    burst_chidx = burst_chidx,
    start_time = burst_times[1, ] / sampling_rate,  
    end_time = burst_times[2, ] / sampling_rate,   
    duration = (burst_times[2, ] - burst_times[1, ]) / sampling_rate  
  )
  
  # Time filtering if specified
  if (!is.null(start) && !is.null(duration)) {
    end_time <- start + duration
    
    burst_data <- burst_data[burst_data$start_time < end_time & 
                               burst_data$end_time > start, ]
  } else {
    cli::cli_inform("Time is not specified, return whole burst_data.")
  }
  
  # Chip ID
  stored_chid <- h5[[paste0(well_id, "/StoredChIdxs")]][]
  
  # Create mapping for grid
  burst_data <- gridMapping(burst_data, stored_chid, grid_n = grid_n, joiner = "burst_chidx")
  
  return(burst_data)
}

#' Extract network burst data from BXR file
#' 
#' Reads network-wide burst events that span multiple channels
#' with optional time filtering and grid mapping.
#' 
#' @param h5 H5File object from \code{\link{openBXR}}
#' @param well_id Character. Well identifier (default: "Well_A1")
#' @param grid_n Integer. Grid size for channel mapping (default: 48)
#' @param start Numeric. Start time in seconds for filtering (optional)
#' @param duration Numeric. Duration in seconds for filtering (optional)
#' @return data.frame with columns:
#'   \describe{
#'     \item{burst_id}{Sequential network burst identifier}
#'     \item{network_burst_chidx}{Channel linear index}
#'     \item{start_time}{Network burst start time in seconds}
#'     \item{end_time}{Network burst end time in seconds}
#'     \item{duration}{Network burst duration in seconds}
#'     \item{grid_x, grid_y}{Grid coordinates (if mapping applied)}
#'   }
#' @import hdf5r
#' @export
#' @seealso \code{\link{bxrSpikeBursts}}, \code{\link{bxrSpikedata}}
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' 
#' # Get all network bursts
#' net_bursts <- bxrNetworkBursts(h5)
#' 
#' h5$close()
#' }
bxrNetworkBursts <- function(h5, 
                             well_id = "Well_A1",
                             grid_n = 48,
                             start = NULL,
                             duration = NULL) {
  network_burst_times <- h5[[paste0(well_id, "/SpikeNetworkBurstTimes")]][,]
  network_burst_chidx <- h5[[paste0(well_id, "/SpikeBurstChIdxs")]][]
  sampling_rate <- getAttributes(h5, attr = "SamplingRate")
  
  if (ncol(network_burst_times) == 0) {
    cli::cli_inform("There are no network bursts.")
    return(data.frame())  
  }
  
  # Create network burst data frame
  network_burst_data <- data.frame(
    burst_id = 1:length(network_burst_chidx),
    network_burst_chidx = network_burst_chidx,
    start_time = network_burst_times[1, ] / sampling_rate,  
    end_time = network_burst_times[2, ] / sampling_rate,   
    duration = (network_burst_times[2, ] - network_burst_times[1, ]) / sampling_rate  
  )
  
  # Time filtering if specified
  if (!is.null(start) && !is.null(duration)) {
    end_time <- start + duration
    
    network_burst_data <- network_burst_data[network_burst_data$start_time < end_time & 
                                               network_burst_data$end_time > start, ]
  } else {
    cli::cli_inform("Time is not specified, return whole network_burst_data.")
  }
  
  # Chip ID
  stored_chid <- h5[[paste0(well_id, "/StoredChIdxs")]][]
  
  # Create mapping for grid
  network_burst_data <- gridMapping(network_burst_data, stored_chid, grid_n = grid_n, joiner = "network_burst_chidx")
  
  return(network_burst_data)
}

#' Extract field potential data from BXR file
#' 
#' Reads field potential events including timestamps and amplitude
#' information with optional time filtering and grid mapping.
#' 
#' @param h5 H5File object from \code{\link{openBXR}}
#' @param well_id Character. Well identifier (default: "Well_A1")
#' @param grid_n Integer. Grid size for channel mapping (default: 48)
#' @param start Numeric. Start time in seconds for filtering (optional)
#' @param duration Numeric. Duration in seconds for filtering (optional)
#' @return data.frame with columns:
#'   \describe{
#'     \item{fp_id}{Sequential field potential identifier}
#'     \item{fp_times}{Field potential timestamps in seconds}
#'     \item{fp_chid}{Channel linear index}
#'     \item{peak_amplitude}{Maximum absolute amplitude}
#'     \item{min_amplitude}{Minimum amplitude}
#'     \item{grid_x, grid_y}{Grid coordinates (if mapping applied)}
#'   }
#' @import hdf5r
#' @export
#' @seealso \code{\link{bxrSpikedata}}, \code{\link{getWaveform}}
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' 
#' # Get all field potentials
#' fps <- bxrFpdata(h5)
#' 
#' # Get field potentials in specific time window
#' fps_subset <- bxrFpdata(h5, start = 10, duration = 30)
#' 
#' h5$close()
#' }
bxrFpdata <- function(h5, 
                      well_id = "Well_A1",
                      grid_n = 48,
                      start = NULL,
                      duration = NULL) {
  # Extract field potential data
  fp_times <- h5[[paste0(well_id, "/FpTimes")]][]
  fp_chidx <- h5[[paste0(well_id, "/FpChIdxs")]][]
  sampling_rate <- getAttributes(h5, attr = "SamplingRate")
  
  if (length(fp_times) == 0) {
    cli::cli_inform("There are no field potentials.")
    return(data.frame())
  }
  
  # Create field potential data frame
  fp_data <- data.frame(
    fp_id = 1:length(fp_times),
    fp_times = fp_times / sampling_rate,
    fp_chid = fp_chidx
  )
  
  # Extract amplitude data from FpForms
  wavelength <- hdf5r::h5attr(h5[[paste0(well_id, "/FpForms")]], "Wavelength")
  fp_forms <- h5[[paste0(well_id, "/FpForms")]][]
  n_fps <- nrow(fp_data)
  
  # create empty numeric vector to store values same length as fp_data
  peak_amplitudes <- numeric(n_fps)
  min_amplitudes <- numeric(n_fps)
  
  for (i in 1:n_fps) {
    start_idx <- (i - 1) * wavelength + 1
    end_idx <- i * wavelength
    waveform <- fp_forms[start_idx:end_idx]
    
    peak_amplitudes[i] <- max(abs(waveform))
    min_amplitudes[i] <- min(waveform)
  }
  
  # Add amplitude columns
  fp_data$peak_amplitude <- peak_amplitudes
  fp_data$min_amplitude <- min_amplitudes
  
  # Time filtering if specified
  if (!is.null(start) && !is.null(duration)) {
    end_time <- start + duration
    fp_data <- fp_data[fp_data$fp_times >= start & 
                         fp_data$fp_times <= end_time, ]
  }
  
  # Grid mapping
  stored_chid <- h5[[paste0(well_id, "/StoredChIdxs")]][]
  fp_data <- gridMapping(fp_data, stored_chid, grid_n = grid_n, joiner = "fp_chid")
  
  return(fp_data)
}

#' Extract waveform data from BXR file
#' 
#' Retrieves spike waveforms either for specific indices or all waveforms,
#' with timing information for analysis and visualization.
#' 
#' @param h5 H5File object from \code{\link{openBXR}}
#' @param well_id Character. Well identifier (default: "Well_A1")
#' @param waveform_indices Numeric vector. Specific waveform indices to extract.
#'   If NULL, returns all waveforms (default: NULL)
#' @return List containing:
#'   \describe{
#'     \item{waveforms}{List of waveforms (if indices specified) or raw array (if all)}
#'     \item{wavelength}{Integer. Number of samples per waveform}
#'     \item{sampling_rate}{Numeric. Sampling rate in Hz}
#'     \item{time_axis}{Numeric vector. Time axis for waveforms in seconds}
#'   }
#' @import hdf5r
#' @export
#' @seealso \code{\link{bxrSpikedata}}
#' @examples
#' \dontrun{
#' h5 <- openBXR("data.bxr")
#' sr <- getAttributes(h5, "SamplingRate")
#' 
#' # Get specific waveforms
#' waves <- getWaveform(h5, waveform_indices = c(1, 5, 10))
#' 
#' # Get all waveforms
#' all_waves <- getWaveform(h5)
#' 
#' h5$close()
#' }
getWaveform <- function(h5, well_id = "Well_A1", waveform_indices = NULL) {
  # Basic info
  wavelength <- hdf5r::h5attr(h5[[paste0(well_id, "/SpikeForms")]], "Wavelength")
  spike_forms <- h5[[paste0(well_id, "/SpikeForms")]][]
  sampling_rate <- getAttributes(h5, attr = "SamplingRate")
  
  # Start list
  waveforms <- list()
  
  # retrieve waveforms
  if (!is.null(waveform_indices)){
    for (idx in waveform_indices) {
      start_idx <- (idx - 1) * wavelength + 1
      end_idx <- idx * wavelength
      waveforms[[paste0("index_", idx)]] <- spike_forms[start_idx:end_idx]
    }
    cli::cli_inform("Specified waveforms are returned")
    return(list(
      waveforms = waveforms,
      wavelength = wavelength,
      sampling_rate = sampling_rate,
      time_axis = (0:(wavelength-1)) / sampling_rate
    ))
  } else {
    cli::cli_inform("Waveforms are not specified. All waveform data is returned.")
    return(list(
      waveforms = spike_forms,
      wavelength = wavelength,
      sampling_rate = sampling_rate,
      time_axis = (0:(wavelength-1)) / sampling_rate
    ))
  }
}
