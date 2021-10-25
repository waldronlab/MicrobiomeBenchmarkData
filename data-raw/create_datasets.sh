#! /bin/bash

## Create independent files
Rscript --vanilla Beghini_2019_16S_smoking.R
Rscript --vanilla HMP_2012_16S_gingival.R
Rscript --vanilla HMP_2012_16S_gingival_subset.R
Rscript --vanilla HMP_2012_WMS_gingival.R
Rscript --vanilla Stammler_2016_16S_spikein.R

## Merge and format the sampleMetadata file
Rscript --vanilla sampleMetadata.R
