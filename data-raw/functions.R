# Some functions to be used within the scripts

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

import_calgaro_2020 <- function() {
    load(url("https://github.com/mcalgaro93/sc2meta/blob/master/data/16Sdatasets_for_replicability_filtered.RData?raw=true"))
    mia::makeTreeSummarizedExperimentFromPhyloseq(ps_list_16S[["Subgingival_Supragingival"]])
}

#
