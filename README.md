# MicrobiomeBenchmarkData

Access to a resource of datasets with biological ground truth for benchmarking
differential abundance methods. The datasets are available from Zenodo: 
https://doi.org/10.5281/zenodo.6911026

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

+ Data on Zenodo: https://doi.org/10.5281/zenodo.6911026
+ R package on Bioconductor: <<<Insert link here when created>>>
+ R package source code: https://github.com/waldronlab/MicrobiomeBenchmarkData
+ R package issues: https://github.com/waldronlab/MicrobiomeBenchmarkData/issues
+ Scripts for preparing the datasets: https://github.com/waldronlab/MicrobiomeBenchmarkDataPrep
+ Code for reproducibility of the analyses: https://waldronlab.io/MicrobiomeBenchmarkDataAnalyses
