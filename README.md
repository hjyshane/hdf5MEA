# hdf5MEA

> R I/O toolkit for multi-well HDF5 MEA files (BRW/BXR v4).


## Overview

hdf5MEA is an R package for reading and processing 3Brain BRW/BXR v4 HDF5 multi-well MEA (Multi-Electrode Array) data files.
It provides reliable loading, partial data extraction, and export capabilities with integration for existing analysis algorithms.

## Features

### Data I/O
- Reliable I/O: Robust reading of BRW v4 (raw data) and BXR v3 (processed results) files
- Memory Efficient: Streaming and partial loading for large datasets (>1GB)
- Export Support: CSV, Parquet, and report export formats
- Data Validation: File integrity and format validation

### Spike Train Analysis
- Firing rate calculation (overall and time-binned)
- Inter-spike interval (ISI) analysis and distributions
- Spike train statistics and autocorrelation
- Comprehensive spike train metrics

### Burst Analysis
- Threshold-based burst detection from spike trains
- Burst statistics (duration, frequency, IBI)
- Burst propagation and synchrony analysis
- Spatial burst pattern detection

### Network Analysis
- Functional connectivity matrices
- Cross-correlation between channel pairs
- Network-wide synchrony metrics
- Spatial pattern and wave detection
- Hub channel identification

### Signal Processing
- Spike detection from raw voltage traces
- Bandpass filtering for noise reduction
- Baseline drift removal
- Artifact detection and removal

### Visualization
- Spike raster plots
- Firing rate heatmaps (2D MEA grid)
- ISI histograms
- Connectivity matrix plots
- Network activity over time
- Waveform visualization

### Quality Control
- Dead/inactive channel detection
- Signal-to-noise ratio (SNR) calculation
- Data completeness checks
- Automated QC summary reports


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

## Analysis Examples

### Spike Train Analysis
```r
# Calculate firing rates
firing_rates <- calculateFiringRate(spikes)

# Time-binned firing rates
binned_rates <- calculateFiringRate(spikes, bin_size = 1)  # 1-second bins

# ISI analysis
isi_stats <- calculateISI(spikes)
isi_ch1 <- calculateISI(spikes, channel_id = 100)  # Specific channel

# Comprehensive statistics
stats <- getSpikeStatistics(spikes)
print(stats$mean_firing_rate)
```

### Burst Detection and Analysis
```r
# Detect bursts from spike trains
bursts <- detectBursts(spikes, isi_threshold = 0.1, min_spikes = 5)

# Calculate burst statistics
burst_stats <- calculateBurstStatistics(bursts)

# Analyze burst propagation
propagation <- analyzeBurstPropagation(bursts, time_window = 0.05)

# Burst synchrony
synchrony <- calculateBurstSynchrony(bursts, bin_size = 0.1)
```

### Network Analysis
```r
# Build functional connectivity matrix
connectivity <- buildConnectivityMatrix(spikes, threshold = 0.3)
print(connectivity$network_density)
print(connectivity$hub_channels)

# Calculate network synchrony
network_sync <- calculateNetworkSynchrony(spikes, bin_size = 0.01)
print(network_sync$synchrony_index)

# Detect spatial patterns
patterns <- detectSpatialPatterns(spikes, time_window = 0.1)
print(patterns$hotspots)
print(patterns$wave_events)
```

### Signal Processing
```r
# Detect spikes from raw voltage
h5 <- openBRW("data.brw")
raw_data <- get_brw_data(h5, start = 0, duration = 10)
parsed <- brwdataParse(raw_data$binary_chunk)
timeseries <- brwtimeseriesConvert(parsed, raw_data$sampling_rate,
                                   raw_data$start_frame)

detected_spikes <- detectSpikes(timeseries, sampling_rate = raw_data$sampling_rate)

# Apply bandpass filter
filtered <- applyBandpassFilter(timeseries[[1]]$voltage,
                                sampling_rate = raw_data$sampling_rate,
                                low_freq = 300, high_freq = 3000)

# Remove artifacts
cleaned <- detectArtifacts(timeseries, artifact_threshold = 1000,
                           sampling_rate = raw_data$sampling_rate,
                           action = "interpolate")
```

### Visualization
```r
library(ggplot2)

# Spike raster plot
p1 <- plotSpikeRaster(spikes, time_range = c(0, 60))
print(p1)

# Firing rate heatmap
p2 <- plotFiringRateHeatmap(spikes)
print(p2)

# ISI histogram
p3 <- plotISIHistogram(spikes, channel_id = 100, max_isi = 0.5)
print(p3)

# Network activity over time
p4 <- plotNetworkActivity(spikes, bin_size = 0.1)
print(p4)

# Connectivity matrix
connectivity <- buildConnectivityMatrix(spikes)
p5 <- plotConnectivityMatrix(connectivity$connectivity_matrix)
print(p5)
```

### Export and Reporting
```r
# Export to CSV
exportToCSV(spikes, "spike_data.csv")

# Export to Parquet (better compression)
exportToParquet(spikes, "spike_data.parquet", compression = "snappy")

# Generate comprehensive analysis report
exportSummaryReport(spikes, bursts,
                    output_file = "analysis_report.txt",
                    include_plots = TRUE,
                    plot_dir = "plots")

# Batch export multiple datasets
data_list <- list(spikes = spikes, bursts = bursts)
batchExport(data_list, output_dir = "exported_data", format = "parquet")
```

### Quality Control
```r
# Detect dead channels
qc_dead <- detectDeadChannels(spikes, min_spike_count = 10,
                               min_firing_rate = 0.1)
print(qc_dead$n_dead)
print(qc_dead$dead_channels)

# Calculate SNR
snr <- calculateSNR(spikes, method = "amplitude")
poor_snr_channels <- snr$channel_id[snr$quality == "poor"]

# Comprehensive QC summary
qc_summary <- generateQCSummary(spikes, h5 = h5)
print(qc_summary$overall_quality)
print(qc_summary$recommendations)
```

## Version

Current version: 0.1.0
