# MicrobiomeBenchmarkData

A resource of datasets with biological ground truth for benchmarking
differential abundance methods.

These datasets are also available through Zenodo at https://sandbox.zenodo.org/deposit/952241
(provisional link at sandbox zenodo).

# Installation

```
if (!"BiocManager" %in% installed.packages()[,"Package"])
    install.packages("BiocManager")

## Development version
BiocManager::install("waldronlab/MicrobiomeBenchmarkData", build_vignettes = TRUE)

```


