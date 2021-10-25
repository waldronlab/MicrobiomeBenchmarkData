
# Scripts to generate the data

This directory contains the R scripts to generate the data to upload to
Zenodo, and how to create the sampleMetadata.RData and sysdata.rda files.

The scripts must be run in the current directory with (set with `setwd`) or
from the terminal with the `create_datasets.sh` script
(commenting and uncommenting lines as necessary).

## Summary of the main steps

1. Create and run an R script (*.R) with code to generate the following files
per dataset:
    1. \*_sample_metadata.tsv
    2. \*_count_matrix.tsv
    3. \*_taxonomy_table.tsv
    4. \*_taxonomy_tree.newick
2. Create a combined sampleMetadata.tsv file with the sampleMetadata.R script.
The file must be created to the current directory and also saved with 
`usethis::use_data` so it is also saved to the data directory of
the package.
3. Upload all of the .tsv (including sampleMetadata.tsv) and .newick files to 
Zenodo.org (or [sandbox.zendodo.org](https://sandbox.zenodo.org/record/946131)). 
4. Once the files have been uploaded to Zenodo, add the resource names and
the download url (from Zenodo) to the
[extdata/metadata.csv](https://github.com/waldronlab/MicrobiomeBenchmarkData/blob/main/inst/extdata/metadata.csv)
file.
5. Update R/sysdata.rda with the sysdata.R script.

Note: 
All of the scripts can be run with the create_datasets.sh script by
(un)commenting the necessary lines.
