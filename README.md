# hdf5MEA

R I/O toolkit for multi-well HDF5 MEA files (BRW/BXR v4).

## Status
MVP: BRW/BXR v4 metadata + partial reads + export.

## Scope (v0.1.0)
- Inspect and load .BRW / .BXR files
- Partial load for channel or frame
- Exports to CSV/Parquet/HDF5

## Install
```r
# devtools::install_github("hjyshane/hdf5MEA")