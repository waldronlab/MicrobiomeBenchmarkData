
showMockrobiotaDatasets <- function() {
  metadata_tsv <- "https://raw.githubusercontent.com/caporaso-lab/mockrobiota/master/inventory.tsv"
  readr::read_tsv(metadata_tsv, show_col_types = FALSE) %>%
    dplyr::select(study_id = StudyID, study_type = `study-type`,
                  description = `human-readable-description`) %>%
    dplyr::mutate(study_type = dplyr::case_when(
      study_type == "marker-gene" ~ "16S rRNA",
      study_type == "metagenome" ~ "WMS"
    ))
}


getMockrobiotaDatasets <- function(x) {


}
