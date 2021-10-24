library(dplyr)

# Import metadata ---------------------------------------------------------

format_sample_metadata <- function(df) {

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
    file = "HMP_2012_16S_gingival_sample_metadata.tsv"
) %>%
    format_sample_metadata()

## HMP_2012_16S_gingival_subset ####
HMP_2012_16S_gingival_V35_subset <- readr::read_tsv(
    file = "HMP_2012_16S_gingival_V35_subset_sample_metadata.tsv"
) %>%
    format_sample_metadata()

## HMP_2012_WMS_gingival ####
HMP_2012_WMS_gingival <- readr::read_tsv(
    file = "HMP_2012_WMS_gingival_sample_metadata.tsv"
) %>%
    format_sample_metadata()

## Beghini_2019_16S_smoking ####
Beghini_2019_16S_smoking <- readr::read_tsv(
    file = "Beghini_2019_16S_smoking_sample_metadata.tsv"
) %>%
    format_sample_metadata()

Stammler_2016_16S_spikein <- readr::read_tsv(
    file = "Stammler_2016_16S_spikein_sample_metadata.tsv"
) %>%
    format_sample_metadata()

# Export merged metadata --------------------------------------------------

datasets <- list(
  HMP_2012_16S_gingival = HMP_2012_16S_gingival,
  HMP_2012_16S_gingival_V35_subset = HMP_2012_16S_gingival_V35_subset,
  HMP_2012_WMS_gingival = HMP_2012_WMS_gingival,
  Beghini_2019_16S_smoking = Beghini_2019_16S_smoking,
  Stammler_2016_16S_spikein = Stammler_2016_16S_spikein
)

sampleMetadata <- bind_rows(datasets, .id = "dataset")

## Export file as tsv (this for upload to Zenodo)
readr::write_tsv(x = sampleMetadata, file = "sampleMetadata.tsv")

## Save data
usethis::use_data(sampleMetadata, overwrite = TRUE)
