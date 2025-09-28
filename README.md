# hdf5MEA

> R I/O toolkit for multi-well HDF5 MEA files (BRW/BXR v4).


## Overview

hdf5MEA is an R package for reading and processing 3Brain BRW/BXR v4 HDF5 multi-well MEA (Multi-Electrode Array) data files.
It provides reliable loading, partial data extraction, and export capabilities with integration for existing analysis algorithms.

## Features

- **BRW v4 Support**: Read EventsBasedSparseRaw data format
- **Time-based Partial Loading**: Extract specific time ranges without loading entire files
- **Memory Efficient**: Streaming iterators and chunk-based processing
- **Sparse Data Handling**: Multiple modes for handling sparse neural data

## Installation

```r
# Install from GitHub
devtools::install_github("hjyshane/hdf5MEA")

# Load package
library(hdf5MEA)
```

## Requirements

- R >= 4.3.0
- hdf5r >= 1.3.0
- data.table >= 1.14.0

## Quick Start

### Basic Usage

```r
library(hdf5MEA)

# Extract 2 seconds of data starting from 10 seconds
raw_data <- get_data(
  file = "experiment.brw",
  well_id = "Well_A1", 
  start = 10,
  duration = 2
)

# Parse binary data
parsed_data <- dataParse(raw_data$binary_chunk)

# Convert to time series
timeseries <- timeseriesConvert(
  parsed_data, 
  raw_data$sampling_rate, 
  raw_data$start_frame,
  mode = "full"  # or "events_only", "threshold"
)

# Access channel data
channel_133375 <- timeseries[["133375"]]
plot(channel_133375$time, channel_133375$voltage, type = "l")
```

### Sparse Data Modes

```r
# Full data (default)
full_data <- timeseriesConvert(parsed_data, sr, start_frame, mode = "full")

# Events only (non-zero values)
events_data <- timeseriesConvert(parsed_data, sr, start_frame, mode = "events_only")

# Threshold-based filtering
filtered_data <- timeseriesConvert(parsed_data, sr, start_frame, 
                                  mode = "threshold", threshold = 100)
```

## File Format Support

### Currently Supported
- **BRW v4**: EventsBasedSparseRaw format
- **Wells**: Multi-well plates (Well_A1, Well_A2, etc.)
- **Time Ranges**: Partial loading by time intervals

### Planned Support
- **BXR v4**: Spike detection results
- **Export Formats**: CSV, Parquet, HDF5
- **Validation**: File integrity checks

## API Reference

### Main Functions

- `get_data()`: Extract time-based data from BRW files
- `dataParse()`: Parse EventsBasedSparseRaw binary data  
- `timeseriesConvert()`: Convert to time series format

### Core Functions

- `openBRW()`: Open and validate BRW files
- `timeCheck()`: Validate time ranges and convert to frames
- `selectChunk()`: Select data chunks for time ranges
- `eventBased()`: Extract EventsBasedSparseRaw data

## Memory Considerations

- Durations > 5 seconds will trigger memory warnings
- Use sparse modes for large datasets
- EventsBasedSparseRaw contains only event periods, not continuous data

## Legal Notice

This package is an independent implementation based on 3Brain's published file format specifications. It does not contain any proprietary 3Brain code or circumvent access controls.

- **3Brain** is a trademark of 3Brain GmbH
- File format specifications used under fair use for interoperability
- No reverse engineering of proprietary algorithms


## License

MIT License - see LICENSE file for details

## Version

Current version: 0.1.0

### Changelog

#### v0.1.0 (Initial Release)
- BRW v4 EventsBasedSparseRaw support
- Time-based partial loading
- Sparse data handling modes
- Memory-efficient chunk processing