## code to prepare `sysdata` dataset goes here

# Load packages ------------------------------------------------------------------------------

# library(SummarizedExperiment)
# library(mia)
# library(magrittr)
# library(dplyr)

# Data from Calgaro, 2020 --------------------------------------------------------------------

# The data from Calgaro, 2020 (https://doi.org/10.1186/s13059-020-02104-1) includes six dataset of
# paried samples, well-suited for benchmarking of methos to do differential abundance (DA)
# analyses. The data comes originally from two sources: the HMP16SData and the
# curatedMetagenomicData packages.

# 16S rRNA data
# load(url("https://github.com/mcalgaro93/sc2meta/blob/master/data/16Sdatasets_for_replicability_filtered.RData?raw=true"))
# tse_list_16s <- lapply(ps_list_16S, makeTreeSummarizedExperimentFromphyloseq)
# names_16s <- lapply(tse_list_16s, function(x) list(features = rownames(x), samples = colnames(x))) %>%
#   magrittr::set_names(paste0("16S_", names(.)))

# WMS data
# load(url("https://github.com/mcalgaro93/sc2meta/blob/master/data/WMSdatasets_for_replicability_filtered.RData?raw=true"))
# tse_list_wms <- lapply(ps_list_WMS, makeTreeSummarizedExperimentFromphyloseq)
# names_wms <- lapply(tse_list_wms, function(x) list(features = rownames(x), samples = colnames(x))) %>%
#   magrittr::set_names(paste0("WMS_", names(.)))
#
# calgaro2020_names <- c(names_16s, names_wms)

urls <- list(
  # For the .calgaro2020Datasets function
  calgaro2020_16S = "https://github.com/mcalgaro93/sc2meta/raw/master/data/16Sdatasets_for_replicability_filtered.RData",
  calgaro2020_WMS = "https://github.com/mcalgaro93/sc2meta/raw/master/data/WMSdatasets_for_replicability_filtered.RData",
  # For the .beghini2019Nychanesmicrobiome function
  nychanes_metadata = "https://github.com/waldronlab/nychanesmicrobiome/raw/master/inst/extdata/public_v2_010518.sas7bdat",
  nychanes_otu_table = "https://github.com/waldronlab/nychanesmicrobiome/raw/master/inst/extdata/otu_table_mc10_w_tax.biom",
  nychanes_taxa_tree = "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/rep_set.tre",
  nychanes_rep_set = "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/rep_set.fna",
  nychanes_original_map = "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/original_map.tsv",
  nychanes_smokingsampleselection = "https://raw.githubusercontent.com/waldronlab/nychanesmicrobiome/master/inst/extdata/smokingsampleselection.tsv"
)

usethis::use_data(
  urls,
  internal = TRUE, overwrite = TRUE
  )
