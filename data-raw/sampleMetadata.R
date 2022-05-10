library(dplyr)

# Import metadata ---------------------------------------------------------
## Use yout custom path here
# setwd('/home/samuel/Projects/CUNY/MicrobiomeBenchmarkData/MicrobiomeBenchmarkData/data-raw')


format_sample_metadata <- function(df) {

    ## One of the main formats is that everything is changed to lower case,
    ## althoug this might not be what is required needed for all cases.

    isSampleNameAndCharacter <- function(x) {
        purrr::map_lgl(x, is.character) & !grepl("(SAMPLE_NAME)|(sample_name)", colnames(x))
    }

    vct_lgl <- isSampleNameAndCharacter(df)

    df %>%
        magrittr::set_colnames(., tolower(colnames(.))) %>%
        magrittr::set_colnames(., gsub(" ", "_", colnames(.))) %>%
        purrr::map_if(.x = ., .p = vct_lgl, .f = tolower) %>%
        purrr::map_if(.x = ., .p = is.character, .f = ~gsub(" ", "_", .x)) %>%
        purrr::map_at(.at = "sample_name", .f = ~ as.character(.)) %>%
        tibble::as_tibble()
}

## HMP_2012_16S_gingival (for both v35 and v13) ####
HMP_2012_16S_gingival <- readr::read_tsv(
    file = "data-raw/HMP_2012_16S_gingival_sample_metadata.tsv"
) %>%
    format_sample_metadata()

## HMP_2012_16S_gingival_subset ####
HMP_2012_16S_gingival_V35_subset <- readr::read_tsv(
    file = "data-raw/HMP_2012_16S_gingival_V35_subset_sample_metadata.tsv"
) %>%
    format_sample_metadata()

## HMP_2012_WMS_gingival ####
HMP_2012_WMS_gingival <- readr::read_tsv(
    file = "data-raw/HMP_2012_WMS_gingival_sample_metadata.tsv"
) %>%
    format_sample_metadata()

## Beghini_2019_16S_smoking ####
Beghini_2019_16S_smoking <- readr::read_tsv(
    file = "data-raw/Beghini_2019_16S_smoking_sample_metadata.tsv"
) %>%
    format_sample_metadata()

Stammler_2016_16S_spikein <- readr::read_tsv(
    file = "data-raw/Stammler_2016_16S_spikein_sample_metadata.tsv"
) %>%
    format_sample_metadata()

Ravel_2011_16S_BV <- read.table(
    file = "data-raw/Ravel_2011_16S_BV_sample_metadata.tsv", sep = '\t', header = TRUE
) %>%
    format_sample_metadata() %>%
    rename(sample_name = sample_id)
# Export merged metadata --------------------------------------------------

datasets <- list(
  HMP_2012_16S_gingival = HMP_2012_16S_gingival,
  HMP_2012_16S_gingival_V35_subset = HMP_2012_16S_gingival_V35_subset,
  HMP_2012_WMS_gingival = HMP_2012_WMS_gingival,
  Beghini_2019_16S_smoking = Beghini_2019_16S_smoking,
  Stammler_2016_16S_spikein = Stammler_2016_16S_spikein,
  Ravel_2011_16S_BV = Ravel_2011_16S_BV
)

sampleMetadata <- bind_rows(datasets, .id = "dataset")

## Export file as tsv (this for upload to Zenodo)
readr::write_tsv(x = sampleMetadata, file = "sampleMetadata.tsv")

## Save data
usethis::use_data(sampleMetadata, overwrite = TRUE)

