# hdf5MEA

> R I/O toolkit for multi-well HDF5 MEA files (BRW/BXR v4).


## Overview

hdf5MEA is an R package for reading and processing 3Brain BRW/BXR v4 HDF5 multi-well MEA (Multi-Electrode Array) data files.
It provides reliable loading, partial data extraction, and export capabilities with integration for existing analysis algorithms.

## Features

- Reliable I/O: Robust reading of BRW v4 (raw data) and BXR v3 (processed results) files
- Memory Efficient: Streaming and partial loading for large datasets (>1GB)
- Export Support: CSV, Parquet, and HDF5 export formats
- Data Validation: File integrity and format validation
- Integration Ready: Compatible with existing R neuroscience workflows


## Quick Start

### Installation
```r
# Install from GitHub
devtools::install_github("hjyshane/hdf5MEA")

# Load package
library(hdf5MEA)
```

### Reading BXR Files (Raw Results)
```r
library(hdf5MEA)

# Open BRW file
file_info <- openBRW("experiment.brw")

# Extract data
sampling_rate <- getAttributes(data, attr = "SamplingRate")
data <- get_brw_data("experiment.brw", 
                     start = 0,      # seconds
                     duration = 1)  # seconds

# Parse binary data
parsed_data <- dataParse(data$binary_chunk)

# Convert to time series
timeseries <- brwtimeseriesConvert(
  parsed_data,  
  data$start_frame,
  mode = "full"  # or "events_only", "threshold"
)

# Access channel data
channel <- timeseries[[100]]
plot(channel$time, channel$voltage, type = "l")

```

### Sparse Data Modes for BRW Files

```r
# Full data (default)
full_data <- brwtimeseriesConvert(parsed_data, 
                                  data$start_frame,
                                  mode = "full")

# Events only (non-zero values)
events_data <- brwtimeseriesConvert(parsed_data, 
                                    data$start_frame,
                                    mode = "events_only")

# Threshold-based filtering
filtered_data <- brwtimeseriesConvert(parsed_data, 
                                      data$start_frame,
                                      mode = "threshold", 
                                      threshold = 100
```

### Reading BXR Files (Processed Results)
```r
# Open BXR file  
data <- openBXR(bxr_test)
sampling_rate <- getAttributes(data, attr = "SamplingRate")

# Extract spike data
spikes <- bxrSpikedata(data)

# Extract spike waveforms
waveforms <- getWaveform(data)

# Extract burst data (if available)
bursts <- bxrSpikeBursts(data,
                         sampling_rate)

```

## Version

Current version: 0.1.0
