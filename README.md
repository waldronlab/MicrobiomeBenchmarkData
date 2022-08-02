# MicrobiomeBenchmarkData

A resource of datasets with biological ground truth for benchmarking
differential abundance methods.

These datasets are also available through Zenodo at https://sandbox.zenodo.org/deposit/952241
(provisional link at sandbox zenodo).

# Installation

```
## Install BioConductor if not installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

## Release version (not yet in Bioc)
BiocManager::install("MicrobiomeBenchmarkData")

## Development version
BiocManager::install("waldronlab/MicrobiomeBenchmarkData", build_vignettes = TRUE)
```

# Links

+ Data on Zenodo: https://zenodo.org/record/6911027
+ R package on Bioconductor: <<<Insert link here when created>>>
+ R package source code: https://github.com/waldronlab/MicrobiomeBenchmarkData
+ R package issues: https://github.com/waldronlab/MicrobiomeBenchmarkData/issues
+ Scripts for preparing the datasets: https://github.com/waldronlab/MicrobiomeBenchmarkDataPrep
+ Code for reproducibility of the analyses: https://waldronlab.io/MicrobiomeBenchmarkDataAnalyses
